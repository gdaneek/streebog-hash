/**
 * @file    streebog.hh
 * @brief   GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.0
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#pragma once
#include <stdint.h>

class Streebog {
    uint64_t n[8];
    uint64_t sum[8];
    uint64_t h[8];

    void i512_sum(const uint64_t * const l, uint64_t const * const r, uint64_t * const o);

    void LPS(uint64_t* in);
    void X(uint64_t const * const l, uint64_t const * const r, uint64_t * const o);
    void G(uint64_t const * const m, bool is_zero = false);

public:

    enum class Mode {
        H512,
        H256,
        __COUNT__
    };

    void reset(const Mode mode);

    explicit Streebog(const Mode mode);
    void update(void* m, const uint64_t size);
    uint64_t const * const finalize(uint8_t const * const m, const uint64_t size);
};
