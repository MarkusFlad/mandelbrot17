//============================================================================
// Name		   : mandelbrot17.cpp
// Author	   : Markus Flad
// Version	   : 1.0.0
// Description : Calculate mandelbrot in C++17
//============================================================================

#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <algorithm>
#include <thread>
#include <sstream>
#if defined(__AVX512BW__) || defined(__AVX__) || defined(__SSE__)
#include <immintrin.h>
#endif

constexpr std::size_t BITS_IN_BYTE = 8;

std::size_t roundToMultiple (std::size_t number, std::size_t base) {
	return number + ((number % base) ? (base - number % base) : 0);
}

class PortableBinaryBitmap {
public:
	PortableBinaryBitmap(const std::string& filename, std::size_t width, std::size_t height)
	: _file (filename)
	, _width (roundToMultiple(width, BITS_IN_BYTE))
	, _height (roundToMultiple(height, std::thread::hardware_concurrency()))
	, _pixelData ((_width * _height) / BITS_IN_BYTE) {
		_file << "P4" << '\n';
		_file << _width << ' ' << _height << '\n';
	}
	~PortableBinaryBitmap() {
		_file.write(_pixelData.data(), _pixelData.size());
	}
	std::size_t width() const {
		return _width;
	}
	std::size_t height() const {
		return _height;
	}
	struct Line {
		std::size_t y;
		std::size_t width;
		char* data;
	};
	class InterlacedIterator {
	public:
		InterlacedIterator(std::size_t y, std::size_t _width, char* data,
				std::size_t interlaceFactor, std::size_t dataPointerIncrement)
		: _il {y, _width, data}
		, _interlaceFactor (interlaceFactor)
		, _dataPointerIncrement (dataPointerIncrement) {
		}
		Line& operator*() {
			return _il;
		}
		bool operator!=(const InterlacedIterator& other) const {
			return _il.data != other._il.data;
		}
		InterlacedIterator& operator++() {
			_il.y += _interlaceFactor;
			_il.data += _dataPointerIncrement;
			return *this;
		}
	private:
		Line _il;
		std::size_t _interlaceFactor;
		std::size_t _dataPointerIncrement;
	};
	class InterlacedCanvas {
	public:
		InterlacedCanvas(PortableBinaryBitmap& pbm, std::size_t yStart, std::size_t interlacedFactor)
		: _pbm (pbm)
		, _yStart (yStart)
		, _interlacedFactor (interlacedFactor)
		, _widthInBytes (pbm._width / BITS_IN_BYTE)
		, _dataPointerIncrement (interlacedFactor * _widthInBytes) {
		}
		InterlacedIterator begin() {
			return InterlacedIterator(_yStart, _pbm._width,
					_pbm._pixelData.data() + _yStart * _widthInBytes,
					_interlacedFactor, _dataPointerIncrement);
		}
		InterlacedIterator end() {
			return InterlacedIterator(_yStart + _pbm._height, _pbm._width,
					_pbm._pixelData.data() + _pbm._pixelData.size() + _yStart * _widthInBytes,
					_interlacedFactor, _dataPointerIncrement);
		}
		std::size_t yStart() const {
			return _yStart;
		}
	private:
		PortableBinaryBitmap& _pbm;
		std::size_t _yStart;
		std::size_t _interlacedFactor;
		std::size_t _widthInBytes;
		std::size_t _dataPointerIncrement;
	};
	friend class InterlacedCanvas;
	std::vector<InterlacedCanvas> provideInterlacedCanvas(std::size_t interlacedFactor) {
		std::vector<InterlacedCanvas> canvasVector;
		for (std::size_t i=0; i<interlacedFactor; i++) {
			canvasVector.emplace_back(InterlacedCanvas(*this, i, interlacedFactor));
		}
		return canvasVector;
	}
private:
	std::ofstream _file;
	std::size_t _width;
	std::size_t _height;
	std::vector<char> _pixelData;
};

struct NoSimdUnion {
	typedef double NumberType;
	typedef double SimdRegisterType;
	NoSimdUnion()
	: reg(val) {
	}
	NoSimdUnion(const NoSimdUnion& other)
	: reg(val) {
		std::copy(std::begin(other.val), std::end(other.val), std::begin(val));
	}
	NoSimdUnion& operator=(const NoSimdUnion& other) {
		std::copy(std::begin(other.val), std::end(other.val), std::begin(val));
		return *this;
	}
	SimdRegisterType* reg;
	NumberType val[8];
};

#if defined(__AVX512BW__) || defined(__AVX__) || defined(__SSE__)
union Simd128DUnion {
	typedef double NumberType;
	typedef __m128d SimdRegisterType;
	SimdRegisterType reg[4];
	NumberType val[8];
};

union Simd256DUnion {
	typedef double NumberType;
	typedef __m256d SimdRegisterType;
	SimdRegisterType reg[2];
	NumberType val[8];
};

union Simd512DUnion {
	typedef double NumberType;
	typedef __m512d SimdRegisterType;
	SimdRegisterType reg[1];
	NumberType val[8];
};
#endif // defined(__AVX512BW__) || defined(__AVX__) || defined(__SSE__)

template<class SimdUnion>
static constexpr std::size_t size() {
	return sizeof(sizeof(SimdUnion::val));
}
template<class SimdUnion>
static constexpr std::size_t numberOfNumbersInRegister() {
	return sizeof(typename SimdUnion::SimdRegisterType) / sizeof(typename SimdUnion::NumberType);
}
template<class SimdUnion>
static constexpr std::size_t numberOfRegisters() {
	return size<SimdUnion>() / numberOfNumbersInRegister<SimdUnion>();
}

template <class SimdUnion>
class VectorizedComplex {
public:
	typedef typename SimdUnion::NumberType NumberType;
	static struct Imaginary {
	} i;
	class SquareIntermediateResult {
	public:
		SquareIntermediateResult() = default;
		void simdReg(std::size_t i, const typename SimdUnion::SimdRegisterType& reg) {
			_squaredAbs.reg[i] = reg;
		}
		bool squaredAbsLessEqualThen(NumberType threshold) const {
			static_assert(size<SimdUnion>() == 8, "squaredAbsLessEqualThen() "
					"is only implemented for SIMD with size of 8.");
			return std::any_of(std::begin(_squaredAbs.val), std::end(_squaredAbs.val),
					[&threshold](auto v) { return v<=threshold; });
		}
		char squaredAbsLessEqualToPixels(NumberType threshold) const {
			static_assert(size<SimdUnion>() == 8, "squaredAbsLessEqualToPixels() "
					"is only implemented for SIMD with size of 8.");
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
	private:
		SimdUnion _squaredAbs;
	};
	VectorizedComplex() = default;
	VectorizedComplex(NumberType commonRealValue, NumberType commonImagValue) {
		setVectorValues(_reals, commonRealValue);
		setVectorValues(_imags, commonImagValue);
	}
	VectorizedComplex(NumberType commonImagValue, Imaginary i) {
		setVectorValues(_imags, commonImagValue);
	}
	void real(std::size_t i, NumberType realValue) {
		_reals.val[i] = realValue;
	}
	VectorizedComplex square(SquareIntermediateResult& sir) const {
		VectorizedComplex resultNumbers;
		for (std::size_t i=0; i<numberOfRegisters<SimdUnion>(); i++) {
			auto realSquared = _reals.reg[i] * _reals.reg[i];
			auto imagSquared = _imags.reg[i] * _imags.reg[i];
			auto realTimesImag = _reals.reg[i] * _imags.reg[i];
			resultNumbers._reals.reg[i] = realSquared - imagSquared;
			resultNumbers._imags.reg[i] = realTimesImag + realTimesImag;
			sir.simdReg(i, realSquared + imagSquared);
		}
		return resultNumbers;
	}
	friend VectorizedComplex operator+(const VectorizedComplex& lhs, const VectorizedComplex& rhs) {
		VectorizedComplex resultNumbers;
		for (std::size_t i=0; i<numberOfRegisters<SimdUnion>(); i++) {
			resultNumbers._reals.reg[i] = lhs._reals.reg[i] + rhs._reals.reg[i];
			resultNumbers._imags.reg[i] = lhs._imags.reg[i] + rhs._imags.reg[i];
		}
		return resultNumbers;
	}
protected:
	static void setVectorValues(SimdUnion& simdUnion, NumberType v) {
		typedef typename SimdUnion::SimdRegisterType SimdRegisterType;
		SimdRegisterType* vValues = simdUnion.reg;
		constexpr auto numbersInReg = numberOfNumbersInRegister<SimdUnion>();
		for (std::size_t i=0; i<size<SimdUnion>(); i+=numbersInReg) {
			if constexpr (numbersInReg == 1) {
				*vValues = v;
			} else if constexpr (numbersInReg == 2) {
				*vValues = SimdRegisterType{v, v};
			} else if constexpr (numbersInReg == 4) {
				*vValues = SimdRegisterType{v, v, v, v};
			} else if constexpr (numbersInReg == 8) {
				*vValues = SimdRegisterType{v, v, v, v, v, v, v, v};
			}
			vValues++;
		}
	}
private:
	SimdUnion _reals;
	SimdUnion _imags;
};

template <class SimdUnion>
class CalculatorThread {
public:
	typedef VectorizedComplex<SimdUnion> VComplex;
	typedef typename SimdUnion::NumberType NumberType;

	CalculatorThread(const std::complex<NumberType>& cFirst, NumberType rasterReal, NumberType rasterImag,
			std::size_t maxIterations, std::size_t iterationsWithoutCheck, NumberType pointOfNoReturn,
			PortableBinaryBitmap::InterlacedCanvas& canvas)
	: _cFirst(cFirst)
	, _rasterReal(rasterReal)
	, _rasterImag(rasterImag)
	, _maxOuterIterations(maxIterations/iterationsWithoutCheck)
	, _iterationsWithoutCheck(iterationsWithoutCheck)
	, _pointOfNoReturn(pointOfNoReturn)
	, _canvas(canvas) {
	}
	void operator()() {
		NumberType squaredPointOfNoReturn = _pointOfNoReturn * _pointOfNoReturn;
		for (PortableBinaryBitmap::Line& line : _canvas) {
			char* nextPixelGroup = line.data;
			NumberType cImagValue = _cFirst.imag() + line.y*_rasterImag;
			for (std::size_t x=0; x<line.width; x+=size<SimdUnion>()) {
				VComplex z(0, 0);
				VComplex c(cImagValue, VComplex::i);
				for (std::size_t i=0; i<size<SimdUnion>(); i++) {
					c.real(i, _cFirst.real() + (x+i)*_rasterReal);
				}
				typename VComplex::SquareIntermediateResult sir;
				for (std::size_t i=0; i<_maxOuterIterations; i++) {
					for (std::size_t j=0; j<_iterationsWithoutCheck; j++) {
						z = z.square(sir) + c;
					}
					if (!sir.squaredAbsLessEqualThen(squaredPointOfNoReturn)) {
						break;
					}
				}
				*nextPixelGroup = sir.squaredAbsLessEqualToPixels(squaredPointOfNoReturn);
				nextPixelGroup++;
			}
		}
	}
private:
	std::complex<NumberType> _cFirst;
	NumberType _rasterReal;
	NumberType _rasterImag;
	std::size_t _maxOuterIterations;
	std::size_t _iterationsWithoutCheck;
	NumberType _pointOfNoReturn;
	PortableBinaryBitmap::InterlacedCanvas _canvas;
};

#if defined(__AVX512BW__)
typedef Simd512DUnion SystemSimdUnion;
#elif defined __AVX__
typedef Simd256DUnion SystemSimdUnion;
#elif defined __SSE__
typedef Simd128DUnion SystemSimdUnion;
#else
typedef NoSimdUnion SystemSimdUnion;
#endif

int main(int argc, char** argv) {
	std::size_t n = 16000;
	if (argc>=2) {
		std::stringstream nss (argv[1]);
		nss >> n;
	}
	typedef std::complex<SystemSimdUnion::NumberType> ComplexNumber;
	const ComplexNumber cFirst (-1.5, -1.0);
	const ComplexNumber cLast (0.5, 1.0);
	const std::size_t maxIterations = 50;
	const std::size_t iterationsWithoutCheck = 5;
	const SystemSimdUnion::NumberType pointOfNoReturn = 2.0;
	PortableBinaryBitmap pbm ("mandelbrot17.pbm", n, n);
	std::size_t numberOfThreads = std::thread::hardware_concurrency();
	auto canvasVector = pbm.provideInterlacedCanvas(numberOfThreads);
	std::vector<std::thread> threads;
	SystemSimdUnion::NumberType rasterReal = (cLast.real() - cFirst.real()) / pbm.width();
	SystemSimdUnion::NumberType rasterImag = (cLast.imag() - cFirst.imag()) / pbm.height();
	for (auto& canvas : canvasVector) {
		CalculatorThread<SystemSimdUnion> calculatorThread(cFirst, rasterReal, rasterImag,
				maxIterations, iterationsWithoutCheck, pointOfNoReturn, canvas);
		threads.push_back(std::thread(calculatorThread));
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
