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
    // precalc mmul
    for(int i = 0;i < 8;i++)
    for(int j = 0;j < 256;j++)
    for(int t = 0;t < 8;t++)
    out.v[i][j] ^= ((!(j&(1<<t))-1)&m_A[63-t-(i<< 3)]);

    // make pi
    T permutted{};

    for(int i = 0;i < 8;i++)
    for(int j = 0;j < 256;j++)
    permutted.v[i][j] = out.v[i][pi[j]];

    return permutted;
}

constexpr struct alignas(32) {  alignas(32) uint64_t v[8][256];} mmul_A_LUT = mmul_A_precalc<std::remove_cv_t<decltype(mmul_A_LUT)>>();


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


void Streebog::reset() {
    get_IV(this->mode, h);
    memset(n, 0, sizeof(n));
    memset(sum, 0, sizeof(sum));
}

#include <immintrin.h>

void Streebog::LPS(uint64_t* in) {

    __m256i bl{}, bh{};
     for(int i = 0;i < 8;i++) {
        uint64_t w = in[i];

        const uint64_t* p = mmul_A_LUT.v[i];

        uint8_t b0 =  w        & 0xFF;
        uint8_t b1 = (w >> 8)  & 0xFF;
        uint8_t b2 = (w >> 16) & 0xFF;
        uint8_t b3 = (w >> 24) & 0xFF;
        uint8_t b4 = (w >> 32) & 0xFF;
        uint8_t b5 = (w >> 40) & 0xFF;
        uint8_t b6 = (w >> 48) & 0xFF;
        uint8_t b7 = (w >> 56) & 0xFF;

        __m256i pi_vl = _mm256_set_epi64x(p[b3], p[b2], p[b1], p[b0]);
        __m256i pi_vh = _mm256_set_epi64x(p[b7], p[b6], p[b5], p[b4]);

        bl = _mm256_xor_si256(bl, pi_vl);
        bh = _mm256_xor_si256(bh, pi_vh);

     }

   _mm256_store_si256((__m256i*)(in + 0), bl);
   _mm256_store_si256((__m256i*)(in + 4), bh);
}


void Streebog::X(uint64_t const * const l, uint64_t const * const r,uint64_t * const o) {
    [&]<uint64_t... I>(std::index_sequence<I...>) {
        ((o[I] = l[I] ^ r[I]), ...);
    } (std::make_index_sequence<8>());
}

void Streebog::G(uint64_t const * const m, bool is_zero) {
    alignas(32) uint64_t K[8], tmp[8];

    is_zero? (void)memcpy(K, h, 64) : X(h, n, K);
    LPS(K);
    X(K, m, tmp);
    LPS(tmp);
    X(K, C, K);
    LPS(K);

    [&]<uint64_t... I>(std::index_sequence<I...>) {
        ((X(K, tmp, tmp),
        LPS(tmp),
        X(K, C+((I+1) << 3), K),
        LPS(K)), ...);
    } (std::make_index_sequence<11>());


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

uint64_t const * const Streebog::finalize(void* m, const uint64_t size) {

    auto mb = (uint8_t const * const)m;
    auto blocks_n = size >> 6;
    uint64_t _d = blocks_n << 6;

    alignas(32) uint64_t buff[8]{};

    this->update((void*)mb, _d);

    auto rem = size - _d;
    auto _b = (uint8_t*)buff;

    for(auto i = 0;i < rem; ++i)
        _b[i] = mb[i + _d];
    _b[rem] = 0x01;

    G(buff);

    n[0] += (rem << 3);

    i512_sum(sum, buff, sum);

    G(n, true);
    G(sum, true);

    return this->h;
}

uint64_t const * const Streebog::operator()(void* m, const uint64_t size, void* out) {
    auto ret = this->finalize(m, size);
    if(out == nullptr)
        return ret;

    for(auto i = 0;i < (mode == Mode::H512? 8 : 4);i++)
        ((uint64_t*)out)[i] = ret[i + (mode == Mode::H512? 0 : 4)];

    return ret;
}
