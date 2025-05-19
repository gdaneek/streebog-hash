/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 1.1
 *
 * @details
 *  In this code, I will often use constructions that, in their meaning, are a loop unrolling at the compilation stage.
 *  This is necessary in order to generate a lot of SIMD transformations
 *  and apply AVX or SSE optimization to them with a significant increase in performance.
 *  I know that this is possible using compiler extensions (i.e. openMP `loop unroll` pragmas),
 *  but I want to have full control over the unrolling without compiler extensions using standard C++ mechanisms only.
 *
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#include "constants.hh"
#include "streebog.hh"
#include <type_traits>
#include <utility>

#include <string.h>


template<typename T>
consteval T b() {
        T out{};
        for(auto i = 0;i < 8;i++) {
            uint8_t byte = 0;
            for(auto j = 0;j < 256;j++, byte++) {
                for(auto t = 0;t < 8;t++)
                    out.data[i][j] ^= ((!(byte&(1<<t))-1)&m_A[63-t-(i<< 3)]);
            }
        }

        return out;
}
constexpr struct {uint64_t data[8][256]; } v = b<std::remove_cv_t<decltype(v)>>();


constexpr inline auto get_precalc_mmul(const uint8_t byte_num, const uint8_t byte_val) {
    return v.data[byte_num][byte_val];
}


inline constexpr auto get_IV(const MODE mode) {
    return IV +( mode == MODE::H512? 0 : 8);
}


inline constexpr auto get_IV(const MODE mode, uint64_t * const out) {
    auto IV = get_IV(mode);
    for(auto i = 0;i < 8;i++)
        out[i] = IV[i];
}


void Streebog::i512_sum(const uint64_t * const a, uint64_t const * const b, uint64_t * const dst) {
    uint8_t carry = 0;
    for(int i = 0; i < 8; ++i) {
        uint64_t tmp = a[i] + b[i] + carry;
        carry = (tmp < a[i]) || (carry && (tmp == a[i] + b[i]));
        dst[i] = tmp;
    }
}


 Streebog::Streebog(MODE mode) : n{}, sum{}, K{}, tmp2{}, tmp3{} {
    get_IV(mode, h); // compile-time filling
}


void Streebog::LPS(uint64_t* in) {

    for(auto i = 0;i < 64;i++)
        ((uint8_t*)in)[i] = pi[((uint8_t*)in)[i]];

    for(uint8_t i{};i < 8; ++i)
        for(uint8_t j{};j < 8;++j)
            ((uint8_t*)tmp3)[i*8 +j] = ((uint8_t*)in)[(j << 3) + i];

    for(auto i = 0;i < 8;i++) {
        uint64_t bf{};
        auto i8 = reinterpret_cast<uint8_t*>(tmp3+i);
        for(auto j = 0;j < 8;j++)
            bf ^= get_precalc_mmul(j, i8[j]);
        in[i] = bf;
    }
}


void Streebog::X(uint64_t const * const l, uint64_t const * const r,uint64_t * const o) {
    for(auto i = 0 ;i < 8;i++)
        o[i] = l[i] ^ r[i];
}

void Streebog::G(uint64_t const * const m, bool is_zero) {
    if (is_zero) for(auto i = 0;i < 8;K[i++] = h[i]);
    else X(h, n, K);

    LPS(K); // K
    X(K, m, tmp2); // buff
    LPS(tmp2);
    X(K, C, K);   // curr_K?
    LPS(K);

    for(auto i = 1;i < 12;++i) {
        X(K, tmp2, tmp2);
        LPS(tmp2);
        X(K, C+(i << 3), K);
        LPS(K);
    }

    X(K, tmp2, tmp2);  // last iteration X only
    X(tmp2, h, tmp2);
    X(tmp2, m, h);   // answer there
}


void  Streebog::update(void* m, const uint64_t size) {   // make void*
    auto m64 = reinterpret_cast<uint64_t const * const>(m);
    const auto full_blocks_n = (size >> 6);
    for(uint64_t i{};i < full_blocks_n;i++, m64 += 8) {      // size / block == size >> 6
        G(m64);
        i512_sum(sum, m64, sum);
        n[0] += 0x200;
    }
}

uint64_t const * const Streebog::finalize(uint8_t const * const m, const uint64_t size) {
    auto blocks_n = size / block;

    uint64_t buff[8]{};

    uint64_t _d = blocks_n * block;
    this->update((void*)m, _d);

    auto rem = size - _d;
    auto _b = (uint8_t*)buff;

    for(auto i = 0;i < rem; ++i)
        _b[i] = m[i + _d];
    _b[rem] = 0x01;

    G(buff);

    n[0] += (rem << 3);

    i512_sum(sum, buff, sum);

    G(n, true);
    G(sum, true);

    return this->h;
}
