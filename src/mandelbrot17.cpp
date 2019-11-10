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

template <class NumberType, std::size_t N>
class VectorizedNumbers {
public:
	VectorizedNumbers() = default;
	VectorizedNumbers(const VectorizedNumbers& other) = default;
	VectorizedNumbers& operator=(const VectorizedNumbers& other) = default;
	VectorizedNumbers(const NumberType* firstNumber)	{
		const NumberType* currentNumber = firstNumber;
		for (size_t i = 0; i < N; i++) {
			_numbers[i] = *currentNumber;
			currentNumber++;
		}
	}
    const NumberType& operator[](std::size_t i) const {
        return _numbers[i];
    }
    NumberType& operator[](std::size_t i) {
        return _numbers[i];
    }
    const std::array<NumberType, N>& numbers() const {
    	return _numbers;
    }
	friend VectorizedNumbers operator+(const VectorizedNumbers& lhs, const VectorizedNumbers& rhs) {
		VectorizedNumbers resultNumbers;
		for (std::size_t i=0; i<N; i++) {
			const NumberType& lhsNumber = lhs._numbers[i];
			const NumberType& rhsNumber = rhs._numbers[i];
			resultNumbers._numbers[i] = lhsNumber + rhsNumber;
		}
		return resultNumbers;
	}
	friend VectorizedNumbers operator*(const VectorizedNumbers& lhs, const VectorizedNumbers& rhs) {
		VectorizedNumbers resultNumbers;
		for (std::size_t i=0; i<N; i++) {
			const NumberType& lhsNumber = lhs._numbers[i];
			const NumberType& rhsNumber = rhs._numbers[i];
			resultNumbers._numbers[i] = lhsNumber * rhsNumber;
		}
		return resultNumbers;
	}
protected:
	std::array<NumberType, N> _numbers;
};

template <std::size_t N>
class VectorizedComplexWithSquaredAbs : public VectorizedNumbers<std::complex<double>, N> {
public:
	VectorizedComplexWithSquaredAbs() = default;
	VectorizedComplexWithSquaredAbs(const VectorizedComplexWithSquaredAbs& other) = default;
	VectorizedComplexWithSquaredAbs& operator=(const VectorizedComplexWithSquaredAbs& other) = default;
	VectorizedComplexWithSquaredAbs& operator=(const VectorizedNumbers<std::complex<double>, N>& other) {
		VectorizedNumbers<std::complex<double>, N>::_numbers = other.numbers();
		return *this;
	}
	VectorizedComplexWithSquaredAbs(const VectorizedNumbers<std::complex<double>, N>& number,
			const VectorizedNumbers<double, N>& numbersSquaredAbs)
	: VectorizedNumbers<std::complex<double>, N>(number)
	, _numbersSquaredAbs(numbersSquaredAbs) {
	}
	VectorizedComplexWithSquaredAbs square() {
		VectorizedNumbers<double, N> squaredAbs;
		for (std::size_t i=0; i<N; i++) {
			std::complex<double>& number = VectorizedNumbers<std::complex<double>, N>::_numbers[i];
			double r2 = number.real() * number.real();
			double i2 = number.imag() * number.imag();
			double ri = number.real() * number.imag();
			squaredAbs[i] = r2 + i2;
			number.real(r2 - i2);
			number.imag(ri + ri);
		}
		return VectorizedComplexWithSquaredAbs(*this, squaredAbs);
	}
	VectorizedNumbers<double, N> squaredAbs() const {
		return _numbersSquaredAbs;
	}
private:
	VectorizedNumbers<double, N> _numbersSquaredAbs;
};

class CalculatorThread {
public:
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
			constexpr std::size_t vectorizationConcurrency = 8;
			VectorizedComplexWithSquaredAbs<vectorizationConcurrency> z;
			VectorizedNumbers<std::complex<double>, vectorizationConcurrency> c;
			VectorizedNumbers<std::size_t, vectorizationConcurrency> numberOfIterations;
			std::size_t startX = 0;
			for (std::size_t x=0; x<_pbm.width(); x++) {
				double cRealValue = _cFirst.real() + x*rasterReal;
				std::size_t i = x%vectorizationConcurrency;
				z[i] = std::complex<double>(0, 0);
				c[i] = std::complex<double>(cRealValue, cImagValue);
				numberOfIterations[i] = 0;
				if (i==7) {
					std::size_t i=0;
					bool anyZNotExceeded = true;
					while (anyZNotExceeded && i < _maxIterations) {
						auto zSquare = z.square();
						z = zSquare + c;
						i++;
						anyZNotExceeded = false;
						for (std::size_t j=0; j<vectorizationConcurrency; j++) {
							if (zSquare.squaredAbs()[j] < squaredPointOfNoReturn) {
								numberOfIterations[j] = i;
								anyZNotExceeded = true;
							}
						}
					}
					for (std::size_t j=0; j<vectorizationConcurrency; j++) {
						mandelbrotLine[startX+j] = (numberOfIterations[j] < _maxIterations);
					}
					startX += vectorizationConcurrency;
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
		CalculatorThread calculatorThread(i, numberOfThreads, cFirst, cLast, maxIterations, pointOfNoReturn, pbm);
		threads.push_back(std::thread(calculatorThread));
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
