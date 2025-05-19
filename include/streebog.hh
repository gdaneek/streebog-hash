/**
 * @file    streebog.hh
 * @brief   GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 1.1
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#pragma once
#include <stdint.h>

/**
 * @brief list of available operating modes
 * @note `Hxxx`means a xxx-bit hash algorithm
 */
enum class MODE {
    H512,
    H256,
    __COUNT__
};


class Streebog {
    uint64_t n[8];
    uint64_t sum[8];
    uint64_t h[8];
    uint64_t K[8], tmp2[8];
    uint64_t tmp3[8];

    void i512_sum(const uint64_t * const l, uint64_t const * const r, uint64_t * const o);

    void LPS(uint64_t* in);
    void X(uint64_t const * const l, uint64_t const * const r, uint64_t * const o);
    void G(uint64_t const * const m, bool is_zero = false);

public:

    inline static constexpr auto block = 64; // 64 bytes

     explicit Streebog(MODE mode);
     void update(void* m, const uint64_t size);
     uint64_t const * const finalize(uint8_t const * const m, const uint64_t size);
    auto operator()(auto&&... args);
};
