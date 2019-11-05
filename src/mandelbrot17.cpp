//============================================================================
// Name        : mandelbrot17.cpp
// Author      : Markus Flad
// Version     :
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

using std::complex;

class PortableBinaryBitmap {
public:
	PortableBinaryBitmap(const std::string& filename, int width, int HEIGHT)
	: _file (filename)
	, _width (width)
	, _height (HEIGHT)
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

double squaredAbs(const complex<double>& c) {
	double cReal = c.real();
	double cImag = c.imag();
	return cReal * cReal + cImag * cImag;
}

class CalculatorThread {
public:
	CalculatorThread(int yBegin, int yRaster, int maxY, double minReal, double rasterReal, double minImag, double rasterImag,
			double pointOfNoReturn, int maxIterations, PortableBinaryBitmap& pbm)
	: _yBegin(yBegin)
	, _yRaster(yRaster)
	, _maxY(maxY)
	, _minReal(minReal)
	, _rasterReal(rasterReal)
	, _minImag(minImag)
	, _rasterImag(rasterImag)
	, _pointOfNoReturn(pointOfNoReturn)
	, _maxIterations(maxIterations)
	, _pbm(pbm) {
	}
    void operator()() const {
    	double squaredPointOfNoReturn = _pointOfNoReturn * _pointOfNoReturn;
    	for (int y=_yBegin; y<_maxY; y+=_yRaster) {
			double cImagValue = _minImag + y*_rasterImag;
			std::vector<bool> mandelbrotLine(_pbm.width());
			for (int x=0; x<_pbm.width(); x++) {
				complex<double> z(0, 0);
				double cRealValue = _minReal + x*_rasterReal;
				complex<double> c(cRealValue, cImagValue);
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
    int _maxY;
    double _minReal;
    double _rasterReal;
    double _minImag;
    double _rasterImag;
    double _pointOfNoReturn;
    int _maxIterations;
    PortableBinaryBitmap& _pbm;
};

int main() {
	constexpr int N = 16000;
	constexpr int WIDTH = N;
	constexpr int HEIGHT = N;
	constexpr double MIN_REAL = -1.5;
	constexpr double MAX_REAL = 0.5;
	constexpr double MIN_IMAG = -1;
	constexpr double MAX_IMAG = 1;
	constexpr double RASTER_REAL = (MAX_REAL - MIN_REAL) / WIDTH;
	constexpr double RASTER_IMAG = (MAX_IMAG - MIN_IMAG) / HEIGHT;
	constexpr int MAXITERATIONS = 50;
	constexpr int POINT_OF_NO_RETURN = 4.0;
	PortableBinaryBitmap pbm ("mandelbrot17.pbm", WIDTH, HEIGHT);
	std::size_t numberOfThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads;
	for (std::size_t i=0; i<numberOfThreads; i++) {
		CalculatorThread calculatorThread(i, numberOfThreads, HEIGHT, MIN_REAL, RASTER_REAL, MIN_IMAG, RASTER_IMAG, POINT_OF_NO_RETURN, MAXITERATIONS, pbm);
		threads.push_back(std::thread(calculatorThread));
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
