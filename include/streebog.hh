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


/**
 * @brief implementation of 256-bit hash algorithm
 * @param in pointer to input data
 * @param bytes_n input data size in bytes
 * @param out array to write out the result (256-bit hash)
 */
void streebog256(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out);


/**
 * @brief implementation of 512-bit hash algorithm
 * @param in pointer to input data
 * @param bytes_n input data size in bytes
 * @param out array to write out the result (256-bit hash)
 */
void streebog512(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out);


/**
 * @brief implementation of  both 256-bit and 512-bit hash algorithms
 * @param in pointer to input data
 * @param bytes_n input data size in bytes
 * @param out array to write out the result (256-bit or 512-bit hash)
 * @note use MODE enumeration to select hash algorithm
 */
void streebog(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out, const MODE mode = MODE::H256);





#ifdef ENABLE_STREEBOG_WRAPPERS

#include<array>
#include <utility>

template<typename T>
inline auto streebog256(T const * const in, const uint64_t bytes_n) noexcept {
    std::array<uint8_t, 32> out;
    streebog256((uint8_t const * const)in, sizeof(T) * bytes_n, out.data());
    return out;
}

template<typename T>
inline auto streebog512(T const * const in, const uint64_t bytes_n) noexcept {
    std::array<uint8_t, 64> out;
    streebog512((uint8_t const * const)in, sizeof(T) * bytes_n, out.data());
    return out;
}


inline auto streebog256(auto&& cc) {    // cc - continuous (in-memory) container
    return streebog256(cc.data(), cc.size());
}


inline auto streebog512(auto&& cc) {    // cc - continuous (in-memory) container
    return streebog512(cc.data(), cc.size());
}


inline auto streebog(const MODE mode, auto&&... args) {
    return (mode == MODE::H512? streebog512(std::forward<decltype(args)>(args)...)
                              : streebog256(std::forward<decltype(args)>(args)...));
}

#endif
