/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.1
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#include "constants.hh"
#include "streebog.hh"
#include <type_traits> // for metaprog templates
#include <utility>     // for index sequences
#include <string.h>    // for memset memcpy


template<typename T>    // compile-time mmul prcalc
consteval T mmul_A_precalc(T out = T{}, T permutted = T{}) {
    // precalc mmul
    for(int i = 0;i < 8;i++)
    for(int j = 0;j < 256;j++)
    for(int t = 0;t < 8;t++)
    out.v[i][j] ^= ((!(j&(1<<t))-1)&m_A[63-t-(i<< 3)]); // flattened if
    // S-box:
    for(int i = 0;i < 8;i++)
    for(int j = 0;j < 256;j++)
    permutted.v[i][j] = out.v[i][pi[j]];

    return permutted;
}

constexpr struct {alignas(32) uint64_t v[8][256];}
mmul_lut = mmul_A_precalc<std::decay_t<decltype(mmul_lut)>>();


void Streebog::reset() {
    memcpy((void*)h, (void*)(IV + (mode == Streebog::Mode::H512? 0 : 8)), 64);
    memset(n, 0, sizeof(n));
    memset(sum, 0, sizeof(sum));
}


void Streebog::i512_sum(const uint64_t * const a, uint64_t const * const b, uint64_t * const dst) {
    for(int i = 0, carry = 0; i < 8; ++i) {
        uint64_t tmp = a[i] + b[i] + carry;
        carry = (tmp < a[i]) || (carry && (tmp == a[i] + b[i]));
        dst[i] = tmp;
    }
}

#if defined(__AVX2__)

void Streebog::LPS(uint64_t* in) {
    __m256i bl{}, bh{};
     for(int i = 0;i < 8;i++) {
        uint64_t w = in[i];
        const uint64_t* p = mmul_lut.v[i];
        uint8_t b0 =  w & 0xFF, b1 = (w >> 8)  & 0xFF, b2 = (w >> 16) & 0xFF, b3 = (w >> 24) & 0xFF,
        b4 = (w >> 32) & 0xFF, b5 = (w >> 40) & 0xFF, b6 = (w >> 48) & 0xFF, b7 = (w >> 56) & 0xFF;

        bl = _mm256_xor_si256(bl, _mm256_set_epi64x(p[b3], p[b2], p[b1], p[b0]));
        bh = _mm256_xor_si256(bh, _mm256_set_epi64x(p[b7], p[b6], p[b5], p[b4]));
     }

   _mm256_store_si256((__m256i*)(in + 0), bl);
   _mm256_store_si256((__m256i*)(in + 4), bh);
}

#else

void Streebog::LPS(uint64_t * const in) {
    uint64_t b[8]{};

    for(int i = 0;i < 8;i++) {
        uint8_t* wb = (uint8_t*)(in + i);
        auto p = (const long long*)mmul_lut.v[i];
        for(int j = 0;j < 8;j++)
            b[j] ^= p[wb[j]];
     }

     memcpy(in, b, 64);
}

#endif


void Streebog::X(uint64_t const * const l, uint64_t const * const r,uint64_t * const o) {
    [&]<uint64_t... I>(std::index_sequence<I...>) {
        ((o[I] = l[I] ^ r[I]), ...);
    } (std::make_index_sequence<8>());
}


void Streebog::G(uint64_t const * const m, bool is_zero) {
    alignas(32) uint64_t K[8], tmp[8];

    is_zero? (void)memcpy(K, h, 64) : X(h, n, K); LPS(K);
    X(K, m, tmp); LPS(tmp);
    X(K, C, K); LPS(K);

    [&]<uint64_t... I>(std::index_sequence<I...>) {
        ((X(K, tmp, tmp), LPS(tmp),
          X(K, C+((I + 1) << 3), K), LPS(K)), ...);
    } (std::make_index_sequence<11>());

    X(K, tmp, tmp);  // last iteration X only
    X(tmp, h, tmp); X(tmp, m, h);
}


void Streebog::update(void* m, const uint64_t size) {
    auto m64 = reinterpret_cast<uint64_t const * const>(m);
    for(uint64_t i{};i < (size >> 6);i++, m64 += 8, *n += 0x200)
        G(m64), i512_sum(sum, m64, sum);
}


uint64_t const * const Streebog::finalize(void* m, const uint64_t size) {
    alignas(32) uint64_t buff[8]{};
    uint64_t _d = size & ~0x3FULL;
    auto rem = size - _d;
    this->update(m, _d);                                    // process whole chunks
    memcpy(buff, (char*)m+_d, rem);                         // pad input
    ((uint8_t*)buff)[rem] = 0x01;

    G(buff), n[0] += (rem << 3), i512_sum(sum, buff, sum);  // last step
    G(n, true), G(sum, true);

    return this->h;
}


uint64_t const * const Streebog::operator()(void* m, const uint64_t size, void* out) {
    auto ret = this->finalize(m, size);
    if(out != nullptr)
        memcpy(out, ret + (mode == Mode::H512? 0 : 4), (mode == Mode::H512? 8 : 4) << 4);

    return ret;
}
