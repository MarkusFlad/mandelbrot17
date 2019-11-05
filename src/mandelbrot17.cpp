//============================================================================
// Name		: mandelbrot17.cpp
// Author	  : Markus Flad
// Version	 :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <string>
#include <fstream>
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
	PortableBinaryBitmap(const std::string& filename, int width, int height)
	: _file (filename)
	, _width (width)
	, _height (height)
	, _currentY (0) {
		_file << "P4" << '\n';
		_file << _width << ' ' << _height << '\n';
	}
	int width() const {
		return _width;
	}
	int height() const {
		return _height;
	}
	void setNextLine(int y, const std::vector<bool>& line) {
		std::unique_lock<std::mutex> lock (_mtx);
		if (y != _currentY) {
			_pendingLines[y] = line;
			return;
		}
		setNextLine(line);
		std::vector<int> yOfLinesDone;
		for (auto it = _pendingLines.begin(); it!=_pendingLines.end(); ++it) {
			int smallestY = it->first;
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
	int _width;
	int _height;
	int _currentY;
	std::map<int, std::vector<bool>> _pendingLines;
};

double squaredAbs(const std::complex<double>& c) {
	double cReal = c.real();
	double cImag = c.imag();
	return cReal * cReal + cImag * cImag;
}

class CalculatorThread {
public:
	CalculatorThread(int yBegin, int yRaster, const std::complex<double>& cFirst, const std::complex<double>& cLast,
			int maxIterations, double pointOfNoReturn, PortableBinaryBitmap& pbm)
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
		for (int y=_yBegin; y<_pbm.height(); y+=_yRaster) {
			double cImagValue = _cFirst.imag() + y*rasterImag;
			std::vector<bool> mandelbrotLine(_pbm.width());
			for (int x=0; x<_pbm.width(); x++) {
				std::complex<double> z(0, 0);
				double cRealValue = _cFirst.real() + x*rasterReal;
				std::complex<double> c(cRealValue, cImagValue);
				int i=0;
				while (squaredAbs(z) <= squaredPointOfNoReturn && i < _maxIterations) {
					z = z*z + c;
					i++;
				}
				mandelbrotLine[x] = (i < _maxIterations);
			}
			_pbm.setNextLine(y, mandelbrotLine);
		}
	}
private:
	int _yBegin;
	int _yRaster;
	std::complex<double> _cFirst;
	std::complex<double> _cLast;
	int _maxIterations;
	double _pointOfNoReturn;
	PortableBinaryBitmap& _pbm;
};

int main() {
	const int N = 16000;
	const std::complex<double> cFirst (-1.5, -1.0);
	const std::complex<double> cLast (0.5, 1.0);
	const int maxIterations = 50;
	const int pointOfNoReturn = 4.0;
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
