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
#include <queue>
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
	void setNextLine(std::size_t y, const std::vector<char>& line) {
		std::unique_lock<std::mutex> lock (_mtx);
		if (y != _currentY) {
			_pendingLines[y] = line;
			return;
		}
		setNextLine(line);
		std::vector<int> yOfLinesDone;
		for (auto it = _pendingLines.begin(); it!=_pendingLines.end(); ++it) {
			std::size_t smallestY = it->first;
			if (smallestY == _currentY) {
				setNextLine(it->second);
				yOfLinesDone.push_back(smallestY);
			}
		}
		for (int yOfLineDone : yOfLinesDone) {
			_pendingLines.erase(yOfLineDone);
		}
		_file.flush();
	}
protected:
	void setNextLine(const std::vector<char>& line) {
		for (char pixelGroup : line) {
			_file.write(&pixelGroup, sizeof(char));
		}
		_currentY++;
	}
private:
	std::mutex _mtx;
	std::ofstream _file;
	std::size_t _width;
	std::size_t _height;
	std::size_t _lineSize;
	std::size_t _currentY;
	std::map<int, std::vector<char>> _pendingLines;
};

template <class SimdRegisterType, std::size_t MAX_VECTORIZATION>
class VectorizedComplex {
public:
	static constexpr std::size_t maxVectorization() {
		return MAX_VECTORIZATION;
	}
	static constexpr std::size_t numberOfDoublesInRegister() {
		return sizeof(SimdRegisterType) / sizeof(double);
	}
	static constexpr std::size_t numberOfRegisters() {
		return MAX_VECTORIZATION / numberOfDoublesInRegister();
	}
	static struct imaginary {
	} i;
	struct SquareIntermediateResult {
		double squaredAbs(std::size_t i) const {
			return reinterpret_cast<const double*>(_vSquaredAbs)[i];
		}
		char squaredAbsLessEqualThen(double threshold) {
			static_assert(MAX_VECTORIZATION == 8, "squaredAbsLessEqualThen() "
					"only implemented for MAX_VECTORIZATION == 8");
			double* squaredAbsArray = reinterpret_cast<double*>(_vSquaredAbs);
			char result = 0;
			if (squaredAbsArray[0] <= threshold) result |= 0b10000000;
			if (squaredAbsArray[1] <= threshold) result |= 0b01000000;
			if (squaredAbsArray[2] <= threshold) result |= 0b00100000;
			if (squaredAbsArray[3] <= threshold) result |= 0b00010000;
			if (squaredAbsArray[4] <= threshold) result |= 0b00001000;
			if (squaredAbsArray[5] <= threshold) result |= 0b00000100;
			if (squaredAbsArray[6] <= threshold) result |= 0b00000010;
			if (squaredAbsArray[7] <= threshold) result |= 0b00000001;
			return result;
		}
		SimdRegisterType _vSquaredAbs[numberOfRegisters()];
	};
	VectorizedComplex() = default;
	VectorizedComplex(double commonRealValue, double commonImagValue) {
		setVectorValues(_vReals, commonRealValue);
		setVectorValues(_vImags, commonImagValue);
	}
	VectorizedComplex(double commonImagValue, imaginary i) {
		setVectorValues(_vImags, commonImagValue);
	}
	VectorizedComplex& setValues(double commonRealValue, double commonImagValue) {
		setVectorValues(_vReals, commonRealValue);
		setVectorValues(_vImags, commonImagValue);
		return *this;
	}
	VectorizedComplex& setRealValues(double commonRealValue) {
		setVectorValues(_vReals, commonRealValue);
		return *this;
	}
	VectorizedComplex& setImagValues(double commonImagValue) {
		setVectorValues(_vImags, commonImagValue);
		return *this;
	}
	void real(std::size_t i, double realValue) {
		reinterpret_cast<double*>(_vReals)[i] = realValue;
	}
	VectorizedComplex square(SquareIntermediateResult& sir) const {
		VectorizedComplex resultNumbers;
		for (std::size_t i=0; i<numberOfRegisters(); i++) {
			auto realSquared = _vReals[i] * _vReals[i];
			auto imagSquared = _vImags[i] * _vImags[i];
			auto realTimesImag = _vReals[i] * _vImags[i];
			resultNumbers._vReals[i] = realSquared - imagSquared;
			resultNumbers._vImags[i] = realTimesImag + realTimesImag;
			sir._vSquaredAbs[i] = realSquared + imagSquared;
		}
		return resultNumbers;
	}
	friend VectorizedComplex operator+(const VectorizedComplex& lhs, const VectorizedComplex& rhs) {
		VectorizedComplex resultNumbers;
		for (std::size_t i=0; i<numberOfRegisters(); i++) {
			resultNumbers._vReals[i] = lhs._vReals[i] + rhs._vReals[i];
			resultNumbers._vImags[i] = lhs._vImags[i] + rhs._vImags[i];
		}
		return resultNumbers;
	}
protected:
	static void setVectorValues(SimdRegisterType* vValues, double v) {
		for (std::size_t i=0; i<maxVectorization(); i+=numberOfDoublesInRegister()) {
			if constexpr (numberOfDoublesInRegister() == 2) {
				*vValues = SimdRegisterType{v, v};
			} else if constexpr (numberOfDoublesInRegister() == 4) {
				*vValues = SimdRegisterType{v, v, v, v};
			} else if constexpr (numberOfDoublesInRegister() == 8) {
				*vValues = SimdRegisterType{v, v, v, v, v, v, v, v};
			}
			vValues++;
		}
	}
private:
	SimdRegisterType _vReals[numberOfRegisters()];
	SimdRegisterType _vImags[numberOfRegisters()];
};

template <class SimdRegisterType, std::size_t MAX_VECTORIZATION>
class CalculatorThread {
public:
	typedef VectorizedComplex<SimdRegisterType, MAX_VECTORIZATION> VComplex;

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
			for (std::size_t x=0; x<_pbm.width(); x+=VComplex::maxVectorization()) {
				VComplex z(0, 0);
				VComplex c(cImagValue, VComplex::i);
				for (std::size_t i=0; i<VComplex::maxVectorization(); i++) {
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
			_pbm.setNextLine(y, mandelbrotLine);
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
		CalculatorThread<__m256d, 8> calculatorThread(i, numberOfThreads, cFirst, cLast, maxIterations, pointOfNoReturn, pbm);
		threads.push_back(std::thread(calculatorThread));
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
