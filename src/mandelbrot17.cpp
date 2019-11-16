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
#include <complex>
#include <algorithm>
#include <thread>
#include <immintrin.h>

constexpr std::size_t BITS_IN_BYTE = 8;

class PortableBinaryBitmap {
public:
	PortableBinaryBitmap(const std::string& filename, std::size_t width, std::size_t height)
	: _file (filename)
	, _width (width)
	, _height (height) {
		_file << "P4" << '\n';
		_file << _width << ' ' << _height << '\n';
	}
	~PortableBinaryBitmap() {
		for (const auto& canvas : _canvasVector) {
			_file.write(canvas.data().data(), canvas.data().size());
		}
	}
	std::size_t width() const {
		return _width;
	}
	std::size_t height() const {
		return _height;
	}
	class Canvas {
	public:
		Canvas(std::size_t width, std::size_t height)
		: _width (width)
		, _height (height)
		, _data((width * height) / BITS_IN_BYTE, 0)
		, _nextPixelGroupIndex(0) {
		}
		std::size_t width() const {
			return _width;
		}
		std::size_t height() const {
			return _height;
		}
		const std::vector<char>& data() const {
			return _data;
		}
		void writePixelGroup(char pixelGroup) {
			_data[_nextPixelGroupIndex] = pixelGroup;
			_nextPixelGroupIndex++;
		}
	private:
		std::size_t _width;
		std::size_t _height;
		std::vector<char> _data;
		std::size_t _nextPixelGroupIndex;
	};
	std::vector<Canvas>& provideParallelCanvas(std::size_t number) {
		std::size_t canvasHeight = _height/number;
		for (std::size_t i=0; i<number; i++) {
			_canvasVector.emplace_back(Canvas(_width, canvasHeight));
		}
		return _canvasVector;
	}
private:
	std::ofstream _file;
	std::size_t _width;
	std::size_t _height;
	std::vector<Canvas> _canvasVector;
};

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
	struct SquareIntermediateResult {
		bool squaredAbsLessEqualThen(NumberType threshold) const {
			static_assert(size<SimdUnion>() == 8, "squaredAbsLessEqualThen() "
					"is only implemented for SIMD with size of 8.");
			return (_squaredAbs.val[0] <= threshold ||
					_squaredAbs.val[1] <= threshold ||
					_squaredAbs.val[2] <= threshold ||
					_squaredAbs.val[3] <= threshold ||
					_squaredAbs.val[4] <= threshold ||
					_squaredAbs.val[5] <= threshold ||
					_squaredAbs.val[6] <= threshold ||
					_squaredAbs.val[7] <= threshold);
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
			sir._squaredAbs.reg[i] = realSquared + imagSquared;
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
			if constexpr (numbersInReg == 2) {
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
			std::size_t maxIterations, std::size_t iterationsWithoutCheck,
			NumberType pointOfNoReturn, PortableBinaryBitmap::Canvas& canvas)
	: _cFirst(cFirst)
	, _rasterReal(rasterReal)
	, _rasterImag(rasterImag)
	, _maxOuterIterations(maxIterations/iterationsWithoutCheck)
	, _iterationsWithoutCheck(iterationsWithoutCheck)
	, _pointOfNoReturn(pointOfNoReturn)
	, _canvas(canvas) {
	}
	void operator()() const {
		NumberType squaredPointOfNoReturn = _pointOfNoReturn * _pointOfNoReturn;
		std::vector<NumberType> cRealValues;
		for (std::size_t x=0; x<_canvas.width(); x++) {
			cRealValues.push_back (_cFirst.real() + x*_rasterReal);
		}
		for (std::size_t y=0; y<_canvas.height(); y++) {
			NumberType cImagValue = _cFirst.imag() + y*_rasterImag;
			for (std::size_t x=0; x<_canvas.width(); x+=size<SimdUnion>()) {
				VComplex z(0, 0);
				VComplex c(cImagValue, VComplex::i);
				for (std::size_t i=0; i<size<SimdUnion>(); i++) {
					c.real(i, cRealValues[x+i]);
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
				_canvas.writePixelGroup(sir.squaredAbsLessEqualToPixels(squaredPointOfNoReturn));
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
	PortableBinaryBitmap::Canvas& _canvas;
};

#if defined(__AVX512BW__)
typedef Simd512DUnion SystemSimdUnion;
#elif defined __AVX__
typedef Simd256DUnion SystemSimdUnion;
#else
typedef Simd128DUnion SystemSimdUnion;
#endif

int main() {
	const std::size_t N = 16000;
	typedef std::complex<SystemSimdUnion::NumberType> ComplexNumber;
	const ComplexNumber cFirst (-1.5, -1.0);
	const ComplexNumber cLast (0.5, 1.0);
	const std::size_t maxIterations = 50;
	const std::size_t iterationsWithoutCheck = 5;
	const SystemSimdUnion::NumberType pointOfNoReturn = 2.0;
	PortableBinaryBitmap pbm ("mandelbrot17.pbm", N, N);
	std::size_t numberOfThreads = std::thread::hardware_concurrency();
	auto& canvasVector = pbm.provideParallelCanvas(numberOfThreads);
	std::vector<std::thread> threads;
	SystemSimdUnion::NumberType rasterReal = (cLast.real() - cFirst.real()) / pbm.width();
	SystemSimdUnion::NumberType rasterImag = (cLast.imag() - cFirst.imag()) / pbm.height();
	SystemSimdUnion::NumberType nextImag = cFirst.imag();
	for (auto& canvas : canvasVector) {
		ComplexNumber cNextFirst(cFirst.real(), nextImag);
		CalculatorThread<SystemSimdUnion> calculatorThread(cNextFirst, rasterReal, rasterImag,
				maxIterations, iterationsWithoutCheck, pointOfNoReturn, canvas);
		threads.push_back(std::thread(calculatorThread));
		nextImag += rasterImag * canvas.height();
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
