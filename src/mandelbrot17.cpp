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

class PortableBinaryBitmap {
public:
	PortableBinaryBitmap(const std::string& filename, std::size_t width, std::size_t height)
	: _file (filename)
	, _width (width)
	, _height (height)
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
	void setNextLine(std::size_t y, const std::vector<bool>& line) {
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
	void setNextLine(const std::vector<bool>& line) {
		int currentBitPos = 0;
		char nextByte = 0;
		for (auto isBlack : line) {
			if (!isBlack) {
				switch(currentBitPos) {
				case 7:
					nextByte |= 0x01;
					break;
				case 6:
					nextByte |= 0x02;
					break;
				case 5:
					nextByte |= 0x04;
					break;
				case 4:
					nextByte |= 0x08;
					break;
				case 3:
					nextByte |= 0x10;
					break;
				case 2:
					nextByte |= 0x20;
					break;
				case 1:
					nextByte |= 0x40;
					break;
				case 0:
					nextByte |= 0x80;
					break;
				}
			}
			currentBitPos++;
			if (currentBitPos >= 8) {
				_file.write(&nextByte, sizeof(char));
				currentBitPos = 0;
				nextByte = 0;
			}
		}
		_currentY++;
	}
private:
	std::mutex _mtx;
	std::ofstream _file;
	std::size_t _width;
	std::size_t _height;
	std::size_t _currentY;
	std::map<int, std::vector<bool>> _pendingLines;
};

template <class SimdRegisterType, std::size_t MAX_VECTORIZATION>
class VectorizedComplex {
public:
	static constexpr std::size_t numberOfDoublesInRegister() {
		return sizeof(SimdRegisterType) / sizeof(double);
	}
	static constexpr std::size_t iterations() {
		return MAX_VECTORIZATION / numberOfDoublesInRegister();
	}
	static constexpr std::size_t maxVectorization() {
		return MAX_VECTORIZATION;
	}
	static struct imaginary {
	} i;
	VectorizedComplex() = default;
	VectorizedComplex(const VectorizedComplex& other) = default;
	VectorizedComplex(double commonRealValue, double commonImagValue) {
		setVectorValues(_reals, commonRealValue);
		setVectorValues(_imags, commonImagValue);
	}
	VectorizedComplex(double commonImagValue, imaginary i) {
		setVectorValues(_imags, commonImagValue);
	}
	VectorizedComplex& operator=(const VectorizedComplex& other) = default;
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
    	_reals[i] = realValue;
    }
    VectorizedComplex square() const {
		VectorizedComplex resultNumbers;
		auto reals = reinterpret_cast<const SimdRegisterType*>(_reals);
		auto imags = reinterpret_cast<const SimdRegisterType*>(_imags);
		auto resultReals = reinterpret_cast<SimdRegisterType*>(resultNumbers._reals);
		auto resultImags = reinterpret_cast<SimdRegisterType*>(resultNumbers._imags);
		auto resultSquaredAbs = reinterpret_cast<SimdRegisterType*>(resultNumbers._squaredAbs);
		for (std::size_t i=0; i<iterations(); i++) {
			auto realSquared = reals[i] * reals[i];
			auto imagSquared = imags[i] * imags[i];
			auto realTimesImag = reals[i] * imags[i];
			resultReals[i] = realSquared - imagSquared;
			resultImags[i] = realTimesImag + realTimesImag;
			resultSquaredAbs[i] = realSquared + imagSquared;
		}
		return resultNumbers;
    }
    double squaredAbs(std::size_t i) const {
    	return _squaredAbs[i];
    }
	friend VectorizedComplex operator+(const VectorizedComplex& lhs, const VectorizedComplex& rhs) {
		VectorizedComplex resultNumbers;
		auto lhsReals = reinterpret_cast<const SimdRegisterType*>(lhs._reals);
		auto lhsImags = reinterpret_cast<const SimdRegisterType*>(lhs._imags);
		auto rhsReals = reinterpret_cast<const SimdRegisterType*>(rhs._reals);
		auto rhsImags = reinterpret_cast<const SimdRegisterType*>(rhs._imags);
		auto resultReals = reinterpret_cast<SimdRegisterType*>(resultNumbers._reals);
		auto resultImags = reinterpret_cast<SimdRegisterType*>(resultNumbers._imags);
		for (std::size_t i=0; i<iterations(); i++) {
			resultReals[i] = lhsReals[i] + rhsReals[i];
			resultImags[i] = lhsImags[i] + rhsImags[i];
		}
		return resultNumbers;
	}
protected:
	static void setVectorValues(double* values, double v) {
		auto vValues = reinterpret_cast<SimdRegisterType*>(values);
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
	double _reals[MAX_VECTORIZATION];
	double _imags[MAX_VECTORIZATION];
	double _squaredAbs[MAX_VECTORIZATION];
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
		double rasterImag = (_cLast.imag() - _cFirst.imag()) / _pbm.width();
		double squaredPointOfNoReturn = _pointOfNoReturn * _pointOfNoReturn;
		for (std::size_t y=_yBegin; y<_pbm.height(); y+=_yRaster) {
			double cImagValue = _cFirst.imag() + y*rasterImag;
			std::vector<bool> mandelbrotLine(_pbm.width());
			std::array<std::size_t, VComplex::maxVectorization()> numberOfIterations;
			std::fill(numberOfIterations.begin(), numberOfIterations.end(), 0);
			for (std::size_t x=0; x<_pbm.width(); x+=VComplex::maxVectorization()) {
				VComplex z(0, 0);
				VComplex c(cImagValue, VComplex::i);
				for (std::size_t i=0; i<VComplex::maxVectorization(); i++) {
					double cRealValue = _cFirst.real() + (x+i)*rasterReal;
					c.real(i, cRealValue);
				}
				std::size_t i=0;
				bool anyZNotExceeded = true;
				while (anyZNotExceeded && i < _maxIterations) {
					z = z.square() + c;
					i++;
					anyZNotExceeded = false;
					for (std::size_t j=0; j<VComplex::maxVectorization(); j++) {
						if (z.squaredAbs(j) < squaredPointOfNoReturn) {
							numberOfIterations[j] = i;
							anyZNotExceeded = true;
						}
					}
				}
				for (std::size_t j=0; j<VComplex::maxVectorization(); j++) {
					mandelbrotLine[x+j] = (numberOfIterations[j] < _maxIterations);
				}
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
