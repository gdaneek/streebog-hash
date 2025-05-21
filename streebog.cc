/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.2-beta
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#include "constants.hh"
#include "streebog.hh"
#include <type_traits> // for metaprog templates
#include <utility>     // for index sequences
#include <string.h>    // for memset memcpy


template<typename T>    // compile-time mmul prcalc
consteval T mmul_A_precalc(T out = T{}, T permutted = T{}) {
    for(int i{};i < 8;i++)for(int j{};j < 256;j++)for(int t{};t < 8;t++) // precalc mmul
    out.v[i][j] ^= ((!(pi[j]&(1<<t))-1)&m_A[63-t-(i<< 3)]);       // flattened mmul if with pi S-box
    return out;
}

inline static constexpr struct{alignas(64) uint64_t v[8][256];} // i dont want use std::array only for this
mmul_lut = mmul_A_precalc<std::decay_t<decltype(mmul_lut)>>();


template<uint64_t... I> using is = std::index_sequence<I...>;
template<uint64_t I> using make_is = std::make_index_sequence<I>;


void Streebog::reset() {
    memcpy((void*)h, (void*)(IV + (mode == Streebog::Mode::H512? 0 : 8)), 64);
    memset(n, 0, sizeof(n)); memset(sum, 0, sizeof(sum));
}


Streebog::Streebog(const Mode _mode) : mode{_mode} {
    this->reset();
}

#if defined(__AVX2__) && !defined(DISABLE_MANUAL_VECTORIZATION) // Manual AVX2 using

void i512_sum(__m256i* _a, __m256i* _b, __m256i* _dst) {    // idk how convert it to AVX2.
                                                            // This, by the way, generates a vector assembler anyway
    uint64_t* a = (uint64_t*)_a, *b = (uint64_t*)_b, *dst = (uint64_t*)_dst;
    for(int i = 0, carry = 0; i < 8; ++i) {
        uint64_t tmp = a[i] + b[i] + carry;
        carry = (tmp < a[i]) || (carry && (tmp == a[i] + b[i]));
        dst[i] = tmp;
    }
}

void LPS(__m256i* in) {
     __m256i bl{}, bh{};

    for(int i = 0;i < 8;i++) {
        __m128i v = _mm_loadl_epi64((const __m128i*)((uint64_t * const)in + i));
        __m128i idx_lo = _mm_cvtepu8_epi32(v);
        __m128i idx_hi = _mm_cvtepu8_epi32(_mm_srli_si128(v, 4));
        __m256i lo_vals = _mm256_i32gather_epi64((const long long*)&mmul_lut.v[i], idx_lo, 8);
        __m256i hi_vals = _mm256_i32gather_epi64((const long long*)&mmul_lut.v[i], idx_hi, 8);
        bl = _mm256_xor_si256(bl, lo_vals);
        bh = _mm256_xor_si256(bh, hi_vals);
     }

   _mm256_store_si256(in, bl);
   _mm256_store_si256(in + 1, bh);
}


void X(__m256i* l, __m256i* r, __m256i* o) {
    _mm256_store_si256(o, _mm256_xor_si256(*l, *r));
    _mm256_store_si256(o+1, _mm256_xor_si256(*++l, *++r));
}


void Streebog::G(__m256i* m, bool is_zero) {
    __m256i K[2], tmp[2];

    memcpy(K, h, 64);
    if(!is_zero) X(K, n, K);
    LPS(K);

    X(K, m, tmp); LPS(tmp);
    X(K, (__m256i*)C, K); LPS(K);

    [&]<uint64_t... I>(is<I...>) {
        ((X(K, tmp, tmp), LPS(tmp),
          X(K, (__m256i*)(C+((I + 1) << 3)), K), LPS(K)), ...);
    } (make_is<11>());

    X(K, tmp, tmp);  // last iteration X only
    X(tmp, h, tmp); X(tmp, m, h);
}

#else   // Metaprogramming cross-platform impl

void i512_sum(const uint64_t * const a, uint64_t const * const b, uint64_t * const dst, int carry = 0) {
    [&]<uint64_t... Is>(is<Is...>) {
        (([&](const uint64_t I) {
            uint64_t tmp = a[I] + b[I] + carry;
            carry = (tmp < a[I]) || (carry && (tmp == a[I] + b[I]));
            dst[I] = tmp;
        } (Is)), ...);
    } (make_is<8>());
}


void LPS(uint64_t * const in) {
    uint64_t b[8]{};
    [&]<uint64_t... I>(is<I...>) __attribute__((always_inline)) {
        (([&]<uint64_t... J>(const uint64_t idx, is<J...>) {
            ((b[J] ^= mmul_lut.v[idx][((uint8_t*)in)[(idx << 3) + J]]), ...);
        } (I, make_is<8>())), ...);
     } (make_is<8>());

     memcpy(in, b, 64);
}


inline void X(uint64_t const * const l, uint64_t const * const r,uint64_t * const o) {
    [&]<uint64_t... I>(is<I...>) {
        ((o[I] = l[I] ^ r[I]), ...);
    } (make_is<8>());
}

void Streebog::G(uint64_t const * const m, bool is_zero) {
    alignas(32) uint64_t K[8], tmp[8];
    memcpy(K, h, 64);
    if(!is_zero) X(K, n, K);
    LPS(K);
    X(K, m, tmp); LPS(tmp);
    X(K, C, K); LPS(K);

    [&]<uint64_t... I>(is<I...>) {
        ((X(K, tmp, tmp), LPS(tmp),
          X(K, C+((I + 1) << 3), K), LPS(K)), ...);
    } (make_is<11>());

    X(K, tmp, tmp);  // last iteration X only
    X(tmp, h, tmp); X(tmp, m, h);
}


#endif

void Streebog::update(void* m, const uint64_t size) {
    #if defined(__AVX2__) && !defined(DISABLE_MANUAL_VECTORIZATION)
    auto m256 = reinterpret_cast<__m256i*>(m);
    for(uint64_t i{};i < (size >> 6);i++, *(uint64_t*)n += 0x200)
        G(m256 + (i << 1)), i512_sum(sum, m256 + (i << 1), sum);
    #else
    auto m64 = reinterpret_cast<uint64_t const * const>(m);
    for(uint64_t i{};i < (size >> 6);i++, *n += 0x200)
        G(m64 + (i << 3)), i512_sum(sum, m64 + (i << 3), sum);
    #endif

}


uint64_t const * const Streebog::finalize(void* m, const uint64_t size) {
     #if defined(__AVX2__) && !defined(DISABLE_MANUAL_VECTORIZATION)
    __m256i buff[2]{};
    #else
    alignas(32) uint64_t buff[8]{};
    #endif
    uint64_t _d = size & ~0x3FULL;
    auto rem = size - _d;
    this->update(m, _d);                                            // process whole chunks
    memcpy((void*)buff, (char*)m+_d, rem);                          // pad input
    ((uint8_t*)buff)[rem] = 0x01;

    G(buff), *(uint64_t*)n += (rem << 3), i512_sum(sum, buff, sum); // last step
    G(n, true), G(sum, true);

    return (uint64_t const * const)(this->h);

}

uint64_t const * const Streebog::operator()(void* m, const uint64_t size, void* out) {
    auto ret = this->finalize(m, size);
    if(out != nullptr)
        memcpy(out, ret + (mode == Mode::H512? 0 : 4), (mode == Mode::H512? 8 : 4) << 4);

    return ret;
}
