/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.0
 *
 *
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#include "constants.hh"
#include "streebog.hh"
#include <type_traits>
#include <utility>
#include <string.h>


template<typename T>    // compile-time mmul prcalc
consteval T mmul_A_precalc(T out = T{}) {
        for(int i = 0;i < 8;i++)
        for(int j = 0;j < 256;j++)
        for(int t = 0;t < 8;t++)
        out.v[i][j] ^= ((!(j&(1<<t))-1)&m_A[63-t-(i<< 3)]);

        return out;
}

constexpr struct {uint64_t v[8][256];} mmul_A_LUT = mmul_A_precalc<std::remove_cv_t<decltype(mmul_A_LUT)>>();


inline constexpr auto get_IV(const Streebog::Mode mode, uint64_t * const out) {
     memcpy((void*)out, (void*)(IV +( mode ==  Streebog::Mode::H512? 0 : 8)), 64);
}


void Streebog::i512_sum(const uint64_t * const a, uint64_t const * const b, uint64_t * const dst) {
    uint8_t carry = 0;
    for(int i = 0; i < 8; ++i) {
        uint64_t tmp = a[i] + b[i] + carry;
        carry = (tmp < a[i]) || (carry && (tmp == a[i] + b[i]));
        dst[i] = tmp;
    }
}


Streebog::Streebog(const Mode mode){ this->reset(mode); }

void Streebog::reset(const Mode mode) {
    get_IV(mode, h);
    memset(n, 0, sizeof(n));
    memset(sum, 0, sizeof(sum));
}

void Streebog::LPS(uint64_t* in) {
    uint64_t bb[8]{}, w;
    for(int i = 0;i < 8;i++) {
        w = in[i];
        const uint64_t* p = mmul_A_LUT.v[i];
        for(int j = 0;j < 8;j++) {
            bb[j] ^= p[pi [ ((unsigned char*)&w)[j]] ];
        }
    }

    memcpy(in, bb, 64);
}


void Streebog::X(uint64_t const * const l, uint64_t const * const r,uint64_t * const o) {
    for(auto i = 0;i < 8;i++)
        o[i] = l[i] ^ r[i];
}


void Streebog::G(uint64_t const * const m, bool is_zero) {
    uint64_t K[8], tmp[8];

    is_zero? (void)memcpy(K, h, 64) : X(h, n, K);
    LPS(K);
    X(K, m, tmp);
    LPS(tmp);
    X(K, C, K);
    LPS(K);

    for(auto i = 1;i < 12;++i) {
        X(K, tmp, tmp);
        LPS(tmp);
        X(K, C+(i << 3), K);
        LPS(K);
    }

    X(K, tmp, tmp);  // last iteration X only
    X(tmp, h, tmp);
    X(tmp, m, h);
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
    auto blocks_n = size >> 6;
    uint64_t _d = blocks_n << 6;

    uint64_t buff[8]{};

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
