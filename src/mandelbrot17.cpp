//============================================================================
// Name		   : mandelbrot17.cpp
// Author	   : Markus Flad
// Description : Calculate mandelbrot in C++17
//============================================================================

#include <string>
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <thread>
#include <climits>
#if defined(__AVX512BW__) || defined(__AVX__) || defined(__SSE__)
#include <immintrin.h>
#endif

const auto numberOfCpuCores = std::thread::hardware_concurrency();

class PortableBinaryBitmap {
public:
	typedef std::size_t Size;
	PortableBinaryBitmap(std::ostream& ostr, Size width, Size height)
	: _ostr (ostr)
	, _width (roundToMultiple(width, CHAR_BIT))
	, _height (roundToMultiple(height, numberOfCpuCores))
	, _data ((_width * _height) / CHAR_BIT) {
		_ostr << "P4" << '\n';
		_ostr << _width << ' ' << _height << '\n';
	}
	~PortableBinaryBitmap() {
		_ostr.write(_data.data(), _data.size());
	}
	Size width() const {
		return _width;
	}
	Size height() const {
		return _height;
	}
	Size widthInBytes() const {
		return _width / CHAR_BIT;
	}
	struct Line {
		Size y;
		Size width;
		char* data;
	};
	class InterlacedCanvas {
	public:
		class Iterator {
		public:
			Iterator(Size y, Size _width, char* data,
					Size interlaceIncrement, Size dataPointerIncrement)
			: _il {y, _width, data}
			, _interlaceIncrement (interlaceIncrement)
			, _dataPointerIncrement (dataPointerIncrement) {
			}
			Line& operator*() {
				return _il;
			}
			bool operator!=(const Iterator& other) const {
				return _il.data != other._il.data;
			}
			Iterator& operator++() {
				_il.y += _interlaceIncrement;
				_il.data += _dataPointerIncrement;
				return *this;
			}
		private:
			Line _il;
			Size _interlaceIncrement;
			Size _dataPointerIncrement;
		};
		InterlacedCanvas(PortableBinaryBitmap& pbm, Size yStart, Size increment)
		: _pbm (pbm)
		, _yStart (yStart)
		, _increment (increment)
		, _dataStart (yStart * pbm.widthInBytes())
		, _dataPointerIncrement (increment * pbm.widthInBytes()) {
		}
		Size width() const {
			return _pbm.width();
		}
		Size height() const {
			return _pbm.height();
		}
		Iterator begin() {
			return Iterator(_yStart, _pbm.width(),
					_pbm._data.data() + _dataStart,
					_increment, _dataPointerIncrement);
		}
		Iterator end() {
			return Iterator(_yStart + _pbm.height(), _pbm.width(),
					_pbm._data.data() + _pbm._data.size() + _dataStart,
					_increment, _dataPointerIncrement);
		}
	private:
		PortableBinaryBitmap& _pbm;
		Size _yStart;
		Size _increment;
		Size _dataStart;
		Size _dataPointerIncrement;
	};
	std::vector<InterlacedCanvas> provideInterlacedCanvas(Size increment) {
		std::vector<InterlacedCanvas> result;
		for (Size yStart=0; yStart<increment; yStart++) {
			result.emplace_back(InterlacedCanvas(*this, yStart, increment));
		}
		return result;
	}
	static Size roundToMultiple (Size number, Size base) {
		return number + ((number % base) ? (base - number % base) : 0);
	}
private:
	std::ostream& _ostr;
	Size _width;
	Size _height;
	std::vector<char> _data;
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
constexpr std::size_t size() {
	return sizeof(sizeof(SimdUnion::val));
}
template<class SimdUnion>
constexpr std::size_t numberOfNumbersInRegister() {
	return sizeof(typename SimdUnion::SimdRegisterType) /
			sizeof(typename SimdUnion::NumberType);
}
template<class SimdUnion>
constexpr std::size_t numberOfRegisters() {
	return size<SimdUnion>() / numberOfNumbersInRegister<SimdUnion>();
}

template<class SimdUnion>
void setValue(SimdUnion& simdUnion, typename SimdUnion::NumberType v) {
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

template <class SimdUnion>
class VectorizedComplexBase {
public:
	typedef typename SimdUnion::NumberType NumberType;
	typedef typename SimdUnion::SimdRegisterType SimdRegisterType;
	typedef std::size_t Size;
	class SquaredAbs {
	public:
		void simdReg(Size i, const SimdRegisterType& reg) {
			_squaredAbs.reg[i] = reg;
		}
		bool operator>(NumberType threshold) const {
			const auto& sqrdAbsVals = _squaredAbs.val;
			return std::all_of(std::begin(sqrdAbsVals), std::end(sqrdAbsVals),
					[&threshold](auto v) { return v>threshold; });
		}
		char lteToPixels(NumberType threshold) const {
			static_assert(size<SimdUnion>() == 8, "lteToPixels() "
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
};

template <class SimdUnion>
class VectorizedComplexCommonImag : public VectorizedComplexBase<SimdUnion> {
public:
	typedef VectorizedComplexBase<SimdUnion> BT;
	using typename BT::NumberType;
	using typename BT::SimdRegisterType;
	using typename BT::SquaredAbs;
	using typename BT::Size;
	VectorizedComplexCommonImag(const SimdUnion& reals,
			SimdRegisterType& commonImagReg)
	: _reals(reals)
	, _commonImagReg(commonImagReg) {
	}
	const SimdUnion& reals() const {
		return _reals;
	}
	const SimdRegisterType& imagReg() const {
		return _commonImagReg;
	}
private:
	SimdUnion _reals;
	SimdRegisterType& _commonImagReg;
};

template <class SimdUnion>
class VectorizedComplex : public VectorizedComplexBase<SimdUnion> {
public:
	typedef VectorizedComplexBase<SimdUnion> BT;
	using typename BT::NumberType;
	using typename BT::SquaredAbs;
	using typename BT::Size;
	VectorizedComplex() = default;
	VectorizedComplex(NumberType commonRealValue,
			NumberType commonImagValue) {
		setValue(_reals, commonRealValue);
		setValue(_imags, commonImagValue);
	}
	VectorizedComplex& square(SquaredAbs& squaredAbs) {
		for (Size i=0; i<numberOfRegisters<SimdUnion>(); i++) {
			auto realSquared = _reals.reg[i] * _reals.reg[i];
			auto imagSquared = _imags.reg[i] * _imags.reg[i];
			auto realTimesImag = _reals.reg[i] * _imags.reg[i];
			_reals.reg[i] = realSquared - imagSquared;
			_imags.reg[i] = realTimesImag + realTimesImag;
			squaredAbs.simdReg(i, realSquared + imagSquared);
		}
		return *this;
	}
	friend VectorizedComplex operator+(const VectorizedComplex& lhs,
			const VectorizedComplexCommonImag<SimdUnion>& rhs) {
		VectorizedComplex resultNumbers;
		for (Size i=0; i<numberOfRegisters<SimdUnion>(); i++) {
			resultNumbers._reals.reg[i] = lhs._reals.reg[i] + rhs.reals().reg[i];
			resultNumbers._imags.reg[i] = lhs._imags.reg[i] + rhs.imagReg();
		}
		return resultNumbers;
	}
private:
	SimdUnion _reals;
	SimdUnion _imags;
};

template <class SimdUnion>
class MandelbrotCalculator {
public:
	typedef VectorizedComplex<SimdUnion> VComplex;
	typedef VectorizedComplexCommonImag<SimdUnion> VComplexCI;
	typedef typename SimdUnion::NumberType NumberType;
	typedef std::size_t Size;

	MandelbrotCalculator(const std::complex<NumberType>& cFirst,
			const std::complex<NumberType>& cLast,
			Size maxIterations, PortableBinaryBitmap::InterlacedCanvas& canvas,
			NumberType pointOfNoReturn = 2.0, Size iterationsWithoutCheck = 5)
	: _cFirst(cFirst)
	, _cLast(cLast)
	, _maxIterations(maxIterations)
	, _canvas(canvas)
	, _pointOfNoReturn(pointOfNoReturn)
	, _iterationsWithoutCheck(iterationsWithoutCheck) {
	}
	void operator()() {
		const NumberType realRange = _cLast.real() - _cFirst.real();
		const NumberType imagRange = _cLast.imag() - _cFirst.imag();
		const NumberType rasterReal = realRange / _canvas.width();
		const NumberType rasterImag = imagRange / _canvas.height();
		const NumberType squaredPointOfNoReturn =
				_pointOfNoReturn * _pointOfNoReturn;
		const Size maxOuterIterations =
				_maxIterations / _iterationsWithoutCheck;
		std::vector<SimdUnion> cRealValues;
		cRealValues.reserve(_canvas.width() / size<SimdUnion>());
		for (Size x=0; x<_canvas.width(); x+=size<SimdUnion>()) {
			SimdUnion cReals;
			for (Size i=0; i<size<SimdUnion>(); i++) {
				cReals.val[i] = _cFirst.real() + (x+i)*rasterReal;
			}
			cRealValues.push_back(cReals);
		}
		for (PortableBinaryBitmap::Line& line : _canvas) {
			char* nextPixels = line.data;
			SimdUnion cImags;
			setValue(cImags, _cFirst.imag() + line.y*rasterImag);
			for (const SimdUnion& cReals : cRealValues) {
				VComplex z(0, 0);
				VComplexCI c(cReals, cImags.reg[0]);
				typename VComplex::SquaredAbs squaredAbs;
				for (Size i=0; i<maxOuterIterations; i++) {
					for (Size j=0; j<_iterationsWithoutCheck; j++) {
						z = z.square(squaredAbs) + c;
					}
					if (squaredAbs > squaredPointOfNoReturn) {
						break;
					}
				}
				*nextPixels = squaredAbs.lteToPixels(squaredPointOfNoReturn);
				nextPixels++;
			}
		}
	}
private:
	std::complex<NumberType> _cFirst;
	std::complex<NumberType> _cLast;
	Size _maxIterations;
	PortableBinaryBitmap::InterlacedCanvas _canvas;
	NumberType _pointOfNoReturn;
	Size _iterationsWithoutCheck;
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
	typedef SystemSimdUnion::NumberType NumberType;
	typedef std::complex<NumberType> ComplexNumber;
	std::size_t n = 16000;
	if (argc>=2) {
		std::stringstream nss (argv[1]);
		nss >> n;
	}
	const std::size_t maxIterations = 50;
	PortableBinaryBitmap pbm(std::cout, n, n);
	auto canvasVector = pbm.provideInterlacedCanvas(numberOfCpuCores);
	std::vector<std::thread> threads;
	for (auto& canvas : canvasVector) {
		threads.push_back(std::thread(MandelbrotCalculator<SystemSimdUnion> (
				ComplexNumber(-1.5, -1.0), ComplexNumber(0.5, 1.0),
				maxIterations, canvas)));
	}
	for (auto& t : threads) {
		t.join();
	}
	return 0;
}
