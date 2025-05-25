/**
 * @file    streebog.hh
 * @brief   GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.3
 * @see https://github.com/gdaneek/streebog-hash
 */

#pragma once
#include <stdint.h>

#if defined(USE_MANUAL_AVX)
#include <immintrin.h>
#include <xmmintrin.h>
#endif

/**
 * @brief GOST 34.11-2018 (34.11-2012 - `Streebog`) implementation with block-by-block hash calculation support
 * @details
 * The standard does not explicitly declare endianness,
 * so this implementation follows the little-endian approach as the most optimal in terms of programming and performance.
 * Instantiate an object and calculate big data hashes using a combination of update (updates the state of h, n, sum)
 * and finalize (implements final processing)
 * @warning by default, the resulting hash is written in little endian (i.e., back to how it is presented in the control examples)
 */
class Streebog {

    #if defined(USE_MANUAL_AVX)
    __m256i n[2];
    __m256i sum[2];
    __m256i h[2];
    void G(__m256i* m, bool is_zero = false);
    #else
    alignas(32) uint64_t n[8];   ///< N variable (number of bits)
    alignas(32) uint64_t sum[8]; ///< Σ variable (sum of all data blocks)
    alignas(32) uint64_t h[8];   ///< h variable (output hash)
    void G(uint64_t const * const m, bool is_zero = false); ///< implementation of G transformation
    #endif

public:

    /**
     * @brief operating modes of the hash function
     * @note Hxxx means xxx-bit hash function mode
     */
    enum class Mode {
        H512,
        H256,
        __COUNT__
    };

    const Mode mode; ///< current algo mode (512-bit | 256-bit)

    /**
     * @brief forcibly resets the state of the class
     * @note h variable takes the IV (init vector) value corresponding to the operating mode
     */
    void reset();
    explicit Streebog(const Mode _mode);


    /**
     * @brief calculates the partial hash of a chunk of data
     * @param m input data
     * @param size data size in bytes
     * @warning do not use this method if you can immediately provide all the data that you need to calculate the hash from.
     * Instead, use operator() or finalize().
     */
    void update(const void* m, const uint64_t size);


    /**
     * @brief processes the last chunk of data, calculating the resulting hash
     * @param m input data
     * @param size data size in bytes
     * @note use operator() for a more convenient call if you can provide all the data at once
     */
    uint64_t const * const finalize(const void* m, const uint64_t size);

    /**
     * @brief alias for finalize() method
     * @param m input data
     * @param size data size in bytes
     * @param out array for writing output; it may not be provided,
     * and then the function will work exactly the same as finalize()
     */
    uint64_t const * const operator()(const void* m, const uint64_t size, void* out = nullptr);

};

#ifdef STREEBOG_ENABLE_WRAPPERS

#include <array>

inline auto streebog512{[](const void* in, const uint64_t in_sz) {
    std::array<uint64_t, 8> out;
    Streebog{Streebog::Mode::H512}(in, in_sz, out.data());

    return out;
}};


inline auto streebog256{[](const void* in, const uint64_t in_sz) {
    std::array<uint64_t, 4> out;
    Streebog{Streebog::Mode::H256}(in, in_sz, out.data());

    return out;
}};

#endif
