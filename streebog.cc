/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.3
 * @see https://github.com/gdaneek/streebog-hash
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

void Streebog::reset() {
    memcpy((void*)h, (void*)(IV + (mode == Streebog::Mode::H512? 0 : 8)), 64);
    memset(n, 0, sizeof(n)); memset(sum, 0, sizeof(sum));
}

Streebog::Streebog(const Mode _mode) : mode{_mode} {
    this->reset();
}

inline void vadd512(void* _a, void* _b, void* __restrict _dst) {   // idk how convert it to AVX2. fully vectorized assembler is generated anyway.
    ui64 *a = (ui64*)_a, *b = (ui64* __restrict)_b, *dst = (ui64*)_dst;
    bool carry{};
    [&]<ui64... Is>(is<Is...>) __attribute__((always_inline)) {
        (([&](const ui64 I) __attribute__((always_inline)) {
            ui64 tmp = a[I] + b[I] + carry;
            carry = (tmp < a[I]) || (carry && (tmp == a[I] + b[I]));
            dst[I] = tmp;
        } (Is)), ...);
    } (make_is<8>());
}

#if defined(USE_MANUAL_AVX)
#include "experimental/avx_streebog.cc"
#else                                   // Metaprogramming cross-platform impl

inline void LPSX(ui64 const* __restrict lhs, ui64 const* __restrict rhs,  ui64* __restrict out) {
    alignas(32) ui64 r[8];
    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((r[I] = lhs[I] ^ rhs[I],
          out[I] = 0), ...);
    } (make_is<8>());

    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((out[I >> 3] ^= mmul_lut.v[I % 8][(uint8_t)r[I & 7]],
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
        ((tmp[I] ^= K[I]), ...);
    } (make_is<8>());

    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((tmp[I] ^= m[I]), ...);
    } (make_is<8>());

    [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
        ((h[I] ^= tmp[I]), ...);
    } (make_is<8>());
}

#endif

void Streebog::update(void * __restrict m, const ui64 size) {
    for(ui64 i{};i < (size >> 6);i++, *(uint64_t*)n += 0x200)
        #if defined(USE_MANUAL_AVX)
        G((__m256i*)m + (i << 1)), vadd512(sum, (__m256i*)m + (i << 1), sum);
        #else
        G((ui64 *)m + (i << 3)), vadd512(sum, (ui64*)m + (i << 3), sum);
        #endif
}

ui64 const * const Streebog::finalize(void * __restrict m, const ui64 size) {
    #if defined(USE_MANUAL_AVX)
    __m256i buff[2]{};
    #else
    alignas(32) uint64_t buff[8]{};
    #endif
    const ui64 _d = size & ~0x3FULL;
    auto rem = size - _d;
    this->update(m, _d);                                            // process whole chunks
    memcpy((void*)buff, (char*)m+_d, rem);                          // pad input
    ((uint8_t*)buff)[rem] = 0x01;
    G(buff), *(ui64*)n += (rem << 3), vadd512(sum, buff, sum);     // last step
    G(n, true), G(sum, true);

    return (ui64 const * const)(this->h);
}

ui64 const * const Streebog::operator()(void* m, const ui64 size, void* out) {
    auto ret = this->finalize(m, size);
    if(out != nullptr)
        memcpy(out, ret + (mode == Mode::H512? 0 : 4), (mode == Mode::H512? 8 : 4) << 3);

    return ret;
}
