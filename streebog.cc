/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.2
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#include "constants.hh"
#include "streebog.hh"
#include <type_traits> // for metaprog templates
#include <utility>     // for index sequences
#include <string.h>    // for memset memcpy

template<typename T>                                                        // compile-time mmul prcalc
consteval T mmul_A_precalc(T out = T{}, T permutted = T{}) {
    for(int i{};i < 8;i++)for(int j{};j < 256;j++)for(int t{};t < 8;t++)    // precalc mmul
    out.v[i][j] ^= ((!(pi[j]&(1<<t))-1)&m_A[63-t-(i<< 3)]);                 // flattened mmul if with pi S-box
    return out;
}

inline static constexpr struct{alignas(64) uint64_t v[8][256];} // i dont want use std::array only for this
mmul_lut = mmul_A_precalc<std::decay_t<decltype(mmul_lut)>>();

template<uint64_t... I> using is = std::index_sequence<I...>;
template<uint64_t I> using make_is = std::make_index_sequence<I>;
using ui64 = uint64_t;
using llp = const long long * __restrict;

void Streebog::reset() {
    memcpy((void*)h, (void*)(IV + (mode == Streebog::Mode::H512? 0 : 8)), 64);
    memset(n, 0, sizeof(n)); memset(sum, 0, sizeof(sum));
}

Streebog::Streebog(const Mode _mode) : mode{_mode} {
    this->reset();
}


inline void i512_sum(void* _a, void* _b, void* _dst, int carry = 0) {   // idk how convert it to AVX2. fully vectorized assembler is generated anyway.
    ui64 *a = (ui64*)_a, *b = (ui64*)_b, *dst = (ui64*)_dst;
    [&]<ui64... Is>(is<Is...>) {
        (([&](const ui64 I) {
            ui64 tmp = a[I] + b[I] + carry;
            carry = (tmp < a[I]) || (carry && (tmp == a[I] + b[I]));
            dst[I] = tmp;
        } (Is)), ...);
    } (make_is<8>());
}

#if defined(__AVX2__) && !defined(DISABLE_MANUAL_AVX) // Manual AVX2 using

__attribute__((always_inline))
inline uint64_t reduce(__m256i v) {
    __m128i x = _mm_xor_si128(_mm256_castsi256_si128(v),
                              _mm256_extracti128_si256(v, 1));
    __m128i s = _mm_unpackhi_epi64(x, x);

    return _mm_cvtsi128_si64(_mm_xor_si128(x, s));
}


void LPSX(__m256i const * __restrict l, __m256i const * __restrict r, __m256i* __restrict out) {
    const __m256i mask = _mm256_set1_epi64x(0xFF), perm = _mm256_set_epi32(0,0,0,0,6,4,2,0);
    const auto lo = _mm_set_epi32(768, 512, 256, 0), ro =  _mm_set_epi32(1792, 1536, 1280, 1024);

    __m256i in0 = _mm256_xor_si256(l[0], r[0]), in1 = _mm256_xor_si256(l[1], r[1]);

    for(int i = 0;i < 8;i++) {
        __m256i bytes0 = _mm256_and_si256(in0, mask), bytes1 = _mm256_and_si256(in1, mask);

        __m128i idx0 = _mm_add_epi32(_mm256_castsi256_si128(_mm256_permutevar8x32_epi32(bytes0, perm)), lo),
                idx1 = _mm_add_epi32(_mm256_castsi256_si128(_mm256_permutevar8x32_epi32(bytes1, perm)), ro);

        ((uint64_t* __restrict)out)[i] = reduce(_mm256_xor_si256(_mm256_i32gather_epi64((llp)mmul_lut.v, idx0, 8),
                                                                 _mm256_i32gather_epi64((llp)mmul_lut.v, idx1, 8)));

        in0 = _mm256_srli_epi64(in0, 8), in1 = _mm256_srli_epi64(in1, 8);
    }
}

__attribute__((always_inline))
inline void X(__m256i const * __restrict l, __m256i const * __restrict r, __m256i* __restrict out) {
    _mm256_store_si256(out, _mm256_xor_si256(l[0], r[0]));
    _mm256_store_si256(out+1, _mm256_xor_si256(l[1], r[1]));
}


void Streebog::G(__m256i * __restrict m, bool is_zero) {
    __m256i K[2], tmp[2];
    static __m256i zeros[2]{};
    memcpy(K, h, 64);

    is_zero? LPSX(K, zeros, K) : LPSX(K, n, K);
    LPSX(K, m, tmp);
    LPSX(K, (__m256i*)C, K);

    for(int i = 0;i < 11;i++) {
        LPSX(K, tmp, tmp);
        LPSX(K, (__m256i*)(C+((i + 1) << 3)), K);
    }

    X(K, tmp, tmp);  // last iteration X only
    X(tmp, h, tmp);
    X(tmp, m, h);
}

#else   // Metaprogramming cross-platform impl


inline void LPSX(ui64 const* __restrict lhs, ui64 const* __restrict rhs,  ui64* __restrict out) {
    ui64 r[8];
    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((r[I] = lhs[I] ^ rhs[I],
          out[I] = 0), ...);
    } (make_is<8>());

    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((out[I >> 3] ^= mmul_lut.v[I % 8][r[I & 7] & 0xff],
          r[I & 7] >>= 8), ...);
    } (make_is<64>());
}

void Streebog::G(ui64 const * __restrict m, bool is_zero) {
    alignas(32) ui64 K[8], tmp[8], zeros[8]{};
    memcpy(K, h, 64);

    is_zero? LPSX(K, zeros, K) : LPSX(K, n, K);
    LPSX(K, m, tmp);
    LPSX(K, C, K);

    [&]<ui64... I>(is<I...>) {
        ((LPSX(K, tmp, tmp),
          LPSX(K, C+((I + 1) << 3), K)), ...);
    } (make_is<11>());

    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((h[I] ^= (tmp[I] ^ K[I] ^ m[I])), ...);
    } (make_is<8>());
}

#endif

void Streebog::update(void* m, const ui64 size) {
    for(ui64 i{};i < (size >> 6);i++, *(uint64_t*)n += 0x200)
        #if defined(__AVX2__) && !defined(DISABLE_MANUAL_AVX)
        G((__m256i*)m + (i << 1)), i512_sum(sum, (__m256i*)m + (i << 1), sum);
        #else
        G((ui64*)m + (i << 3)), i512_sum(sum, (ui64*)m + (i << 3), sum);
        #endif
}

ui64 const * const Streebog::finalize(void* m, const ui64 size) {
    #if defined(__AVX2__) && !defined(DISABLE_MANUAL_AVX)
    __m256i buff[2]{};
    #else
    alignas(32) uint64_t buff[8]{};
    #endif
    ui64 _d = size & ~0x3FULL;
    auto rem = size - _d;
    this->update(m, _d);                                            // process whole chunks
    memcpy((void*)buff, (char*)m+_d, rem);                          // pad input
    ((uint8_t*)buff)[rem] = 0x01;
    G(buff), *(ui64*)n += (rem << 3), i512_sum(sum, buff, sum);     // last step
    G(n, true), G(sum, true);

    return (ui64 const * const)(this->h);
}

ui64 const * const Streebog::operator()(void* m, const ui64 size, void* out) {
    auto ret = this->finalize(m, size);
    if(out != nullptr)
        memcpy(out, ret + (mode == Mode::H512? 0 : 4), (mode == Mode::H512? 8 : 4) << 4);

    return ret;
}
