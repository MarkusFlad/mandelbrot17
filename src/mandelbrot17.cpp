// The Computer Language Benchmarks Game
// https://salsa.debian.org/benchmarksgame-team/benchmarksgame/
//
// Contributed by Markus Flad
//
// compile with following g++ flags
//  -std=c++17 -O3 -Wall -fomit-frame-pointer -march=native -mfpmath=sse -msse2 -mno-fma

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

// Put everything in a namespace forces inlining
namespace {

const auto numberOfCpuCores = std::thread::hardware_concurrency();

// The PortableBinaryBitmap manages access to the pbm output file and provides
// interlaced canvases that allow threads to write to the bitmap in parallel.
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
        constexpr static Size pixelsPerWrite() {
            return sizeof(data);
        }
        Size y;
        Size width;
        char* data;
    };
    // The InterlacedCanvas provides interlaced access to the bitmap data. Each
    // thread must use its own InterlacedCanvas to write to the bitmap.
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
        std::vector<InterlacedCanvas> interlacedCanvasVector;
        for (Size yStart=0; yStart<increment; yStart++) {
            interlacedCanvasVector.emplace_back(*this, yStart, increment);
        }
        return interlacedCanvasVector;
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

// If the system does not support SIMD, NoSimdUnion can be used.
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
constexpr std::size_t numberOfNumbers() {
    return sizeof(SimdUnion::val) / sizeof(typename SimdUnion::NumberType);
}
template<class SimdUnion>
constexpr std::size_t numberOfNumbersInRegister() {
    return sizeof(typename SimdUnion::SimdRegisterType) /
            sizeof(typename SimdUnion::NumberType);
}
template<class SimdUnion>
constexpr std::size_t numberOfRegisters() {
    return numberOfNumbers<SimdUnion>() /
            numberOfNumbersInRegister<SimdUnion>();
}
template<class SimdUnion>
void setValue(SimdUnion& simdUnion, typename SimdUnion::NumberType v) {
    typedef typename SimdUnion::SimdRegisterType SimdRegisterType;
    SimdRegisterType* vValues = simdUnion.reg;
    constexpr auto numbersInReg = numberOfNumbersInRegister<SimdUnion>();
    for (std::size_t i=0; i<numberOfNumbers<SimdUnion>(); i+=numbersInReg) {
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

// VectorizedComplex provides a convenient interface to deal with complex
// numbers and uses the power of SIMD for high execution speed.
template <class SimdUnion>
class VectorizedComplex {
public:
    typedef typename SimdUnion::NumberType NumberType;
    typedef typename SimdUnion::SimdRegisterType SimdRegisterType;
    typedef std::size_t Size;

    // SquaredAbs is passed to special VectorizedComplex methods that calculate
    // the squared absolute value of the complex number as an intermediate.
    // This means that the calculation does not have to be done twice.
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
            static_assert(numberOfNumbers<SimdUnion>() == 8, "lteToPixels() "
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
    VectorizedComplex(const VectorizedComplex&) = default;
    VectorizedComplex(const SimdUnion& reals, NumberType commonImagValue)
    : _reals(reals) {
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
            const VectorizedComplex& rhs) {
        VectorizedComplex resultNumbers;
        for (Size i=0; i<numberOfRegisters<SimdUnion>(); i++) {
            resultNumbers._reals.reg[i] = lhs._reals.reg[i] + rhs._reals.reg[i];
            resultNumbers._imags.reg[i] = lhs._imags.reg[i] + rhs._imags.reg[i];
        }
        return resultNumbers;
    }
private:
    SimdUnion _reals;
    SimdUnion _imags;
};

// The ComplexPlaneCalculator performs function f(c), with c as a
// VectorizedComplex and a byte as the return value. Due to its eightfold
// vectorization, each returned bit can return a Boolean value from the
// calculation f(c). The full byte is then written to the canvas. This is done
// until the whole bitmap is filled.
template <class SimdUnion, class Functor>
class ComplexPlaneCalculator {
public:
    typedef VectorizedComplex<SimdUnion> VComplex;
    typedef typename SimdUnion::NumberType NumberType;
    typedef typename PortableBinaryBitmap::Line Line;
    typedef std::size_t Size;

    ComplexPlaneCalculator(const std::complex<NumberType>& cFirst,
            const std::complex<NumberType>& cLast,
            PortableBinaryBitmap::InterlacedCanvas& canvas, Functor f)
    : _cFirst(cFirst)
    , _cLast(cLast)
    , _canvas(canvas)
    , _f(f) {
        static_assert(numberOfNumbers<SimdUnion>() == Line::pixelsPerWrite());
    }
    void operator()() {
        const NumberType realRange = _cLast.real() - _cFirst.real();
        const NumberType imagRange = _cLast.imag() - _cFirst.imag();
        const NumberType rasterReal = realRange / _canvas.width();
        const NumberType rasterImag = imagRange / _canvas.height();
        std::vector<SimdUnion> cRealValues;
        cRealValues.reserve(_canvas.width() / Line::pixelsPerWrite());
        for (Size x=0; x<_canvas.width(); x+=Line::pixelsPerWrite()) {
            SimdUnion cReals;
            for (Size i=0; i<Line::pixelsPerWrite(); i++) {
                cReals.val[i] = _cFirst.real() + (x+i)*rasterReal;
            }
            cRealValues.push_back(cReals);
        }
        for (Line& line : _canvas) {
            char* nextPixels = line.data;
            char lastPixels = *nextPixels;
            const NumberType cImagValue = _cFirst.imag() + line.y*rasterImag;
            for (const SimdUnion& cReals : cRealValues) {
                const VComplex c(cReals, cImagValue);
                *nextPixels = _f(c, lastPixels);
                lastPixels = *nextPixels;
                nextPixels++;
            }
        }
    }
private:
    std::complex<NumberType> _cFirst;
    std::complex<NumberType> _cLast;
    PortableBinaryBitmap::InterlacedCanvas _canvas;
    Functor _f;
};

// Functor calculating a Mandelbrot iteration for a VectorizedComplex. This
// means that for eight complex numbers the Mandelbrot calculation is
// (potentially) executed in parallel. The result is a byte that contains a 1
// for each bit if the corresponding complex number is in the Mandelbrot set,
// and a 0 if it is not.
template <class SimdUnion>
class MandelbrotFunction {
public:
    typedef VectorizedComplex<SimdUnion> VComplex;
    typedef typename SimdUnion::NumberType NumberType;
    typedef std::size_t Size;
    constexpr static Size ITERATIONS_WITHOUT_CHECK = 5;
    constexpr static char NO_PIXEL_IN_MANDELBROT_SET = 0x0;

    MandelbrotFunction(Size maxIterations, NumberType pointOfNoReturn = 2.0)
    : _maxOuterIterations(maxIterations / ITERATIONS_WITHOUT_CHECK)
    , _squaredPointOfNoReturn(pointOfNoReturn * pointOfNoReturn) {
    }
    static void mandelbrotIterationsWithoutCheck(VComplex& z, const VComplex& c,
    		typename VComplex::SquaredAbs& squaredAbs) {
		for (Size j=0; j<ITERATIONS_WITHOUT_CHECK; j++) {
			z = z.square(squaredAbs) + c;
		}
    }
    char operator()(const VComplex& c, char lastPixels) const {
        VComplex z = c;
        typename VComplex::SquaredAbs squaredAbs;
        if (lastPixels == NO_PIXEL_IN_MANDELBROT_SET) {
			for (Size i=0; i<_maxOuterIterations; i++) {
				mandelbrotIterationsWithoutCheck(z, c, squaredAbs);
				if (squaredAbs > _squaredPointOfNoReturn) {
					return 0;
				}
			}
        } else {
			for (Size i=0; i<_maxOuterIterations; i++) {
				mandelbrotIterationsWithoutCheck(z, c, squaredAbs);
			}
        }
        return squaredAbs.lteToPixels(_squaredPointOfNoReturn);
    }
private:
    Size _maxOuterIterations;
    NumberType _squaredPointOfNoReturn;
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

} // end namespace

int main(int argc, char** argv) {
    typedef SystemSimdUnion::NumberType NumberType;
    typedef std::complex<NumberType> ComplexNumber;
    typedef ComplexPlaneCalculator<SystemSimdUnion,
            MandelbrotFunction<SystemSimdUnion>> MandelbrotCalculator;
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
        threads.emplace_back(MandelbrotCalculator (ComplexNumber(-1.5, -1.0),
                ComplexNumber(0.5, 1.0), canvas,
                MandelbrotFunction<SystemSimdUnion> (maxIterations)));
    }
    for (auto& t : threads) {
        t.join();
    }
    return 0;
}
