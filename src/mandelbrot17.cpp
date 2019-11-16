//============================================================================
// Name		: mandelbrot17.cpp
// Author	  : Markus Flad
// Version	 :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <string>
#include <fstream>
#include <array>
#include <vector>
#include <utility>
#include <list>
#include <complex>
#include <algorithm>
#include <map>
#include <mutex>
#include <thread>
#include <iostream>
#include <immintrin.h>

constexpr std::size_t BITS_IN_BYTE = 8;

class PortableBinaryBitmap {
public:
	PortableBinaryBitmap(const std::string& filename, std::size_t width, std::size_t height)
	: _file (filename)
	, _width (width)
	, _height (height)
	, _lineSize (_width / BITS_IN_BYTE)
	, _writingFile (false)
	, _currentY (0) {
		_file << "P4" << '\n';
		_file << _width << ' ' << _height << '\n';
	}
	std::size_t width() const {
		return _width;
	}
	std::size_t height() const {
		return _height;
	}
	std::size_t lineSize() const {
		return _lineSize;
	}
	void setNextLine(std::size_t y, std::vector<char>&& line) {
		std::lock_guard lock (_mtx);
		_lines.emplace_back(HorizontalLine{y, std::move(line)});
		if (!_writingFile) {
			writeFile();
		}
	}
protected:
	struct HorizontalLine {
		bool operator<(const HorizontalLine& other) {
			return y < other.y;
		}
		std::size_t y;
		std::vector<char> line;
	};
	void writeFile() {
		std::vector<std::vector<char>> pendingLines;
		{
			std::lock_guard lock (_mtx);
			_lines.sort();
			std::size_t nextY = _currentY;
			while (!_lines.empty()) {
				const auto& line = _lines.front();
				if (line.y != nextY) {
					break;
				}
				pendingLines.push_back(std::move(line.line));
				_lines.pop_front();
				nextY++;
			}
			_currentY = nextY;
			_writingFile = true;
		}
		for (const auto& line : pendingLines) {
			for (char pixelGroup : line) {
				_file.write(&pixelGroup, sizeof(char));
			}
		}
		_file.flush();
		{
			std::lock_guard lock (_mtx);
			_writingFile = false;
		}
	}
private:
	std::recursive_mutex _mtx;
	std::ofstream _file;
	std::size_t _width;
	std::size_t _height;
	std::size_t _lineSize;
	std::list<HorizontalLine> _lines;
	bool _writingFile;
	std::size_t _currentY;
};

union Simd256DoubleOctet {
	typedef __m256d SimdRegisterType;
	static constexpr std::size_t size() {
		return sizeof(sizeof(val));
	}
	static constexpr std::size_t numberOfDoublesInRegister() {
		return sizeof(__m256d) / sizeof(double);
	}
	static constexpr std::size_t numberOfRegisters() {
		return size() / numberOfDoublesInRegister();
	}
	__m256d reg[2];
	double val[8];
};

template <class SimdDoubleOctet>
class VectorizedComplex {
public:
	static struct imaginary {
	} i;
	struct SquareIntermediateResult {
		double squaredAbs(std::size_t i) const {
			return _squaredAbs.val[i];
		}
		char squaredAbsLessEqualThen(double threshold) {
			char result = 0;
			if (_squaredAbs.val[0] <= threshold) result |= 0b10000000;
			if (_squaredAbs.val[1] <= threshold) result |= 0b01000000;
			if (_squaredAbs.val[2] <= threshold) result |= 0b00100000;
			if (_squaredAbs.val[3] <= threshold) result |= 0b00010000;
			if (_squaredAbs.val[4] <= threshold) result |= 0b00001000;
			if (_squaredAbs.val[5] <= threshold) result |= 0b00000100;
			if (_squaredAbs.val[6] <= threshold) result |= 0b00000010;
			if (_squaredAbs.val[7] <= threshold) result |= 0b00000001;
			return result;
		}
		SimdDoubleOctet _squaredAbs;
	};
	VectorizedComplex() = default;
	VectorizedComplex(double commonRealValue, double commonImagValue) {
		setVectorValues(_reals, commonRealValue);
		setVectorValues(_imags, commonImagValue);
	}
	VectorizedComplex(double commonImagValue, imaginary i) {
		setVectorValues(_imags, commonImagValue);
	}
	VectorizedComplex& setValues(double commonRealValue, double commonImagValue) {
		setVectorValues(_reals, commonRealValue);
		setVectorValues(_imags, commonImagValue);
		return *this;
	}
	VectorizedComplex& setRealValues(double commonRealValue) {
		setVectorValues(_reals, commonRealValue);
		return *this;
	}
	VectorizedComplex& setImagValues(double commonImagValue) {
		setVectorValues(_imags, commonImagValue);
		return *this;
	}
	void real(std::size_t i, double realValue) {
		_reals.val[i] = realValue;
	}
	VectorizedComplex square(SquareIntermediateResult& sir) const {
		VectorizedComplex resultNumbers;
		for (std::size_t i=0; i<SimdDoubleOctet::numberOfRegisters(); i++) {
			auto realSquared = _reals.reg[i] * _reals.reg[i];
			auto imagSquared = _imags.reg[i] * _imags.reg[i];
			auto realTimesImag = _reals.reg[i] * _imags.reg[i];
			resultNumbers._reals.reg[i] = realSquared - imagSquared;
			resultNumbers._imags.reg[i] = realTimesImag + realTimesImag;
			sir._squaredAbs.reg[i] = realSquared + imagSquared;
		}
		return resultNumbers;
	}
	friend VectorizedComplex operator+(const VectorizedComplex& lhs, const VectorizedComplex& rhs) {
		VectorizedComplex resultNumbers;
		for (std::size_t i=0; i<SimdDoubleOctet::numberOfRegisters(); i++) {
			resultNumbers._reals.reg[i] = lhs._reals.reg[i] + rhs._reals.reg[i];
			resultNumbers._imags.reg[i] = lhs._imags.reg[i] + rhs._imags.reg[i];
		}
		return resultNumbers;
	}
protected:
	static void setVectorValues(SimdDoubleOctet& simdDoubleOctet, double v) {
		typedef typename SimdDoubleOctet::SimdRegisterType SimdRegisterType;
		typename SimdDoubleOctet::SimdRegisterType* vValues = simdDoubleOctet.reg;
		for (std::size_t i=0; i<SimdDoubleOctet::size(); i+=SimdDoubleOctet::numberOfDoublesInRegister()) {
			if constexpr (simdDoubleOctet.numberOfDoublesInRegister() == 2) {
				*vValues = SimdRegisterType{v, v};
			} else if constexpr (simdDoubleOctet.numberOfDoublesInRegister() == 4) {
				*vValues = SimdRegisterType{v, v, v, v};
			} else if constexpr (simdDoubleOctet.numberOfDoublesInRegister() == 8) {
				*vValues = SimdRegisterType{v, v, v, v, v, v, v, v};
			}
			vValues++;
		}
	}
private:
	SimdDoubleOctet _reals;
	SimdDoubleOctet _imags;
};

template <class SimdDoubleOctet>
class CalculatorThread {
public:
	typedef VectorizedComplex<SimdDoubleOctet> VComplex;

	CalculatorThread(std::size_t yBegin, std::size_t yRaster, const std::complex<double>& cFirst, const std::complex<double>& cLast,
			std::size_t maxIterations, double pointOfNoReturn, PortableBinaryBitmap& pbm)
	: _yBegin(yBegin)
	, _yRaster(yRaster)
	, _cFirst(cFirst)
	, _cLast(cLast)
	, _maxIterations(maxIterations)
	, _pointOfNoReturn(pointOfNoReturn)
	, _pbm(pbm) {
	}
	void operator()() const {
		double rasterReal = (_cLast.real() - _cFirst.real()) / _pbm.width();
		double rasterImag = (_cLast.imag() - _cFirst.imag()) / _pbm.height();
		double squaredPointOfNoReturn = _pointOfNoReturn * _pointOfNoReturn;
		for (std::size_t y=_yBegin; y<_pbm.height(); y+=_yRaster) {
			double cImagValue = _cFirst.imag() + y*rasterImag;
			std::vector<char> mandelbrotLine(_pbm.lineSize());
			for (std::size_t x=0; x<_pbm.width(); x+=SimdDoubleOctet::size()) {
				VComplex z(0, 0);
				VComplex c(cImagValue, VComplex::i);
				for (std::size_t i=0; i<SimdDoubleOctet::size(); i++) {
					double cRealValue = _cFirst.real() + (x+i)*rasterReal;
					c.real(i, cRealValue);
				}
				char absLessEqualPointOfNoReturn = 0;
				typename VComplex::SquareIntermediateResult sir;
				for (std::size_t i=0; i<_maxIterations; i++) {
					z = z.square(sir) + c;
					absLessEqualPointOfNoReturn = sir.squaredAbsLessEqualThen(squaredPointOfNoReturn);
					if (!absLessEqualPointOfNoReturn) {
						break;
					}
				}
				mandelbrotLine[x/BITS_IN_BYTE] = absLessEqualPointOfNoReturn;
			}
			_pbm.setNextLine(y, std::move(mandelbrotLine));
		}
	}
private:
	std::size_t _yBegin;
	std::size_t _yRaster;
	std::complex<double> _cFirst;
	std::complex<double> _cLast;
	std::size_t _maxIterations;
	double _pointOfNoReturn;
	PortableBinaryBitmap& _pbm;
};

int main() {
	const std::size_t N = 16000;
	const std::complex<double> cFirst (-1.5, -1.0);
	const std::complex<double> cLast (0.5, 1.0);
	const std::size_t maxIterations = 50;
	const double pointOfNoReturn = 4.0;
	PortableBinaryBitmap pbm ("mandelbrot17.pbm", N, N);
	std::size_t numberOfThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads;
	for (std::size_t i=0; i<numberOfThreads; i++) {
		CalculatorThread<Simd256DoubleOctet> calculatorThread(i, numberOfThreads, cFirst, cLast, maxIterations, pointOfNoReturn, pbm);
		threads.push_back(std::thread(calculatorThread));
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
