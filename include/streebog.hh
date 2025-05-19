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


/**
 * @brief GOST 34.11-2018 (34.11-2012 - `Streebog`) implementation with block-by-block hash calculation support
 * @details
 * The standard does not explicitly declare endianness,
 * so this implementation follows the little-endian approach as the most optimal in terms of programming and performance.
 * Instantiate an object and calculate big data hashes using a combination of update (updates the state of h, n, sum)
 * and finalize (implements final processing). For simple cases (when all the data has already been placed in the buffer),
 * it is rational to use the class as a functor called immediately after creation, or the more convenient functional wrapper streebog below.
 *
 * @warning by default, the resulting hash is written in little endian (i.e., back to how it is presented in the control examples)
 */
class Streebog {
    uint64_t n[8]; ///< N variable (number of bits)
    uint64_t sum[8]; ///< Σ variable (sum of all data blocks)
    uint64_t h[8]; ///< h variable (output hash)

    void i512_sum(const uint64_t * const l, uint64_t const * const r, uint64_t * const o);

    void LPS(uint64_t* in);
    void X(uint64_t const * const l, uint64_t const * const r, uint64_t * const o);
    void G(uint64_t const * const m, bool is_zero = false);

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

    const Mode mode;

    /**
     * @brief forcibly resets the state of the class
     * @note h variable takes the IV (init vector) value corresponding to the operating mode
     */
    void reset();
    explicit Streebog(const Mode _mode) : mode{_mode} { this->reset(); }


    /**
     * @brief calculates the partial hash of a chunk of data
     * @param m input data
     * @param size data size in bytes
     * @warning do not use this method if you can immediately provide all the data that you need to calculate the hash from.
     * Instead, use operator() or finalize().
     */
    void update(void* m, const uint64_t size);


    /**
     * @brief processes the last chunk of data, calculating the resulting hash
     * @param m input data
     * @param size data size in bytes
     * @note use operator() for a more convenient call if you can provide all the data at once
     */
    uint64_t const * const finalize(void* m, const uint64_t size);

    /**
     * @brief alias for finalize() method
     * @param m input data
     * @param size data size in bytes
     * @param out array for writing output; it may not be provided,
     * and then the function will work exactly the same as finalize()
     */
    uint64_t const * const operator()(void* m, const uint64_t size, void* out = nullptr);

};


/**
 * @brief A wrapper function for more convenient class invocation in hash calculations based on data provided immediately
 */
inline uint64_t const * const streebog(const Streebog::Mode mode, void* in, const uint64_t in_sz, void* out) {
    return Streebog{mode}(in, in_sz, out);
}


#ifdef STREEBOG_ENABLE_WRAPPERS

inline uint64_t const * const streebog256(void* in, const uint64_t in_sz, void* out = nullptr) {
    return Streebog{Streebog::Mode::H256}(in, in_sz, out);
}


inline uint64_t const * const streebog512(void* in, const uint64_t in_sz, void* out = nullptr) {
    return Streebog{Streebog::Mode::H512}(in, in_sz, out);
}


inline auto streebog(void* in, const uint64_t in_sz, void* out = nullptr, const Streebog::Mode mode = Streebog::Mode::H256) {
    return (mode == Streebog::Mode::H256? streebog256(in, in_sz, out) : streebog512(in, in_sz, out));
}


#include <array>

inline auto streebog512(void* in, const uint64_t in_sz) {
    std::array<uint64_t, 8> out;
    streebog512(in, in_sz, out.data());

    return out;
}


inline auto streebog256(void* in, const uint64_t in_sz) {
    std::array<uint64_t, 4> out;
    streebog256(in, in_sz, out.data());

    return out;
}



#endif
