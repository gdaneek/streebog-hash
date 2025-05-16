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

/**
 * @brief
 *
 */
template<uint64_t F, typename funcT, typename... Args>
inline constexpr auto unroll(funcT&& func, Args&&... args) {
    [&]<uint64_t... I>(std::integer_sequence<uint64_t, I...>) __attribute__((always_inline)) {
        ((func.template operator()<I>(std::forward<Args>(args)...)), ... );
    } (std::make_integer_sequence<uint64_t, F>());
}

template<uint64_t F, typename PWFuncT, typename ReduceFuncT, typename... Args>
inline constexpr auto unroll(PWFuncT&& func, ReduceFuncT&& reduce, Args&&... args) {
    return [&]<uint64_t... I>(std::integer_sequence<uint64_t, I...>) __attribute__((always_inline)) {
        return reduce(func.template operator()<I>(std::forward<Args>(args)...)...);
    } (std::make_integer_sequence<uint64_t, F>());
}


/**
 * @brief Implementation of the S transformation (pi substitution)
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @warning If you are not using a UNROLL, only a single-byte data type is allowed as input.
 * use reinterpret_cast to uint8_t* before function call to avoid mistakes
 */
template<typename T>
constexpr inline auto S(T * const in) {

#ifdef METAPROG_UNROLL
    unroll<0x40 / sizeof(T)>([in]<uint64_t i>() {
        in[i] = unroll<sizeof(T)>(
                []<uint64_t j>(uint8_t* v) { return (T)pi[v[j]] << (j << 3); },
                [](auto&&... args) { return (args | ...); },
                (uint8_t*)(in + i));
    });
#else
    for(uint8_t i{};i < 0x40; ++i)
        in[i] = pi[in[i]];              // substitution by pi definition
#endif

    return in;
}


/**
 * @brief Implementation of the P transformation (r permutation)
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @warning If you are not using a UNROLL, only a single-byte data type is allowed as input.
 * use reinterpret_cast to uint8_t* before function call to avoid mistakes
 */
template<typename T>
constexpr inline auto P(T * const in) {
    constexpr auto SZ = 64 / sizeof(T);
    T b[SZ]{};

#ifdef METAPROG_UNROLL
    unroll<SZ>(
    [&]<uint64_t i>() {
        b[i] = unroll<sizeof(T)>(
            []<uint64_t j>(T* in) { return (in[j]>>(56-(i<<3))&0xFF)<<(56-(j<<3)); },
            [](auto&&... args) { return (args | ...); },
            in);
    });

    unroll<SZ>([&]<uint64_t i>(){ in[i] = b[i]; });
#else
    for(uint8_t i{};i < 8; ++i)
        for(uint8_t j{};j < 8;++j)
            b[63-i] = in[(j << 3) + i];

    for(uint8_t i{};i < SZ; ++i)
        in[i] = b[i];
#endif

    return in;
}


/**
 * @brief Implementation of the L transformation (GF2 mmul)
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @warning If you are not using a UNROLL, only a eight-byte data type is allowed as input.
 * use reinterpret_cast to uint64_t* before function call to avoid mistakes
 */
template<typename T>
constexpr inline auto L(T * const in) {

#ifdef METAPROG_UNROLL
    unroll<8>([in]<uint64_t i>() {
        in[i] = unroll<64>(
            []<uint64_t j>(T v) { return ((T)!(v&1ull<<j)-1)&m_A[63-j]; },
            [](auto&&... args){return (args ^ ...);},
            in[i]);
    });
#else
    for(T i = 0, r = 0;i < 8;in[i++] = r, r = 0)
        for(auto j = 0;j < 64;j++)
            r ^= ((T)!(in[i] & (1ull << j))-1) & m_A[63-j];   // flattened if
#endif

    return in;
}


/**
 * @brief Wrapper for calls to frequently used LPS transformations
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 */
template<typename T> constexpr inline auto LPS(T* in) { return L(P(S(in))); }


/**
 * @brief Implementation of the X transformation
 * @param l sizeof(T) * Len bytes array on which the transformation must be performed
 * @param r sizeof(T) * Len bytes array on which the transformation must be performed
 * @param o sizeof(T) * Len bytes array to write out the result
 * @return pointer to the output array
 * @note   it's no different from xor operation, so it is often called instead of it
 */
template<typename T, uint64_t Len = 8>
constexpr inline auto X(T const * const l, T const * const r, T * const o) {

#ifdef METAPROG_UNROLL
    unroll<Len>([&]<uint64_t i>(){ o[i] = l[i] ^ r[i]; });
#else
    for(uint8_t i{};i < Len;++i)
        o[i] = l[i] ^ r[i];
#endif

    return o;
}


/**
 * @brief Wrapper for calls to frequently used LPSX transformations
 * @param l sizeof(T) * Len bytes array on which the transformation must be performed
 * @param r sizeof(T) * Len bytes array on which the transformation must be performed
 * @param o sizeof(T) * Len bytes array to write out the result
 */
template<typename T, uint64_t Len = 8>
constexpr inline auto LPSX(T const * const l, T const * const r, T * const o) {
    return LPS<T>(X<T, Len>(l, r, o));
}


/**
 * @brief a function that provides a pointer to the desired initialization vector (IV)
 * @param mode hashing mode of operation (H256 -- 256bit or H512 -- 512bit)
 * @return pointer to 512 bit IV
 */
inline constexpr auto get_IV(const MODE mode) {
    return IV +( mode == MODE::H512? 0 : 8);
}


inline constexpr auto get_IV(const MODE mode, uint64_t * const out) {
    auto IV = get_IV(mode);
#ifdef METAPROG_UNROLL
    unroll<8>([&]<uint64_t i>(){ out[i] = IV[i]; });
#else
    for(auto i = 0;i < 8;out[i++] = IV[i]);
#endif

}


/**
 * @brief Implementation of the E transformation
 * @param K pointer to initial K value (to 512 bit array)
 * @param m operand passed with ()
 * @param out array to write out the result
 * @return pointer to the output array
 */
template<typename T, uint64_t Len = 8>
constexpr inline auto E(const T * const K, T const * const m, T * const out) {
    T b[64 / sizeof(T)]{}, curr_K[64 / sizeof(T)]{};

    LPSX<T, Len>(K, m, b);
    LPSX<T, Len>(K, C, curr_K);

#ifdef METAPROG_UNROLL
    unroll<11>([&]<uint64_t i>() {
        LPSX<T, Len>(curr_K, b, b);
        LPSX<T, Len>(curr_K, C+((i+1) << 3), curr_K);
    });
#else
    for(auto i = 1;i < 12;++i)
        LPSX<T, Len>(curr_K, b, b),
        LPSX<T, Len>(curr_K, C+(i << 3), curr_K);
#endif

    return X<T, Len>(curr_K, b, out);  // last iteration X only
}


/**
 * @brief Implementation of the G transformation
 * @param N pointer to current N value (to 512 bit array)
 * @param h pointer to current h value (to 512 bit array)
 * @param out array to write out the result
 * @return pointer to the output array
 */
template<typename T, uint64_t Len = 8>
constexpr inline auto G(const T * const N, T const * const h, T const * const m, T * const out) {
    T b[64 / sizeof(T)]{};
    LPSX<T, Len>(h, N, b);
    E<T, Len>(b, m, b);
    return X<T, Len>(X<T, Len>(b, h, b), m, out);
}


/**
 * @brief Implementation of 512-bit addition
 * @param l lhs 512-bit operand
 * @param r rhs 512-bit operand
 * @param o array to write out the result
 * @return pointer to the output array
 */
constexpr inline auto i512_sum(const uint64_t * const l, uint64_t const * const r, uint64_t * const o) {
    uint64_t res[8]{};
    for (int i = 7; i > -1; --i) {
        uint64_t ai = l[i], bi = r[i];
        uint64_t sum = ai + bi;
        auto carry = (sum < ai) || (sum < bi) || (sum == (uint64_t)-1);
        res[i] += sum;
        if(i != 0)
            res[i-1] += carry;

    }

#ifdef METAPROG_UNROLL
    unroll<8>([&]<uint64_t i>(){ o[i] = res[i]; });
#else
    for(auto i = 0;i < 8;o[i++] = res[i]);
#endif

}

constexpr inline bool is_little_endian() { return ((uint32_t)0x01020304 & 0xFF) == 0x04; }


constexpr inline auto endianess_swap(uint64_t * const dst, uint64_t const * const src) {
#ifdef METAPROG_UNROLL
    unroll<8>([&]<uint64_t i>(){ dst[i] = __builtin_bswap64(src[i]); });
#else
    for(auto i = 0;i < 8;dst[i++] = __builtin_bswap64(src[i]));
#endif

}

// TODO: rewrite code to remove swaps. idk why devs of the algorithm chose such an inconvenient endianess
//  in addition, it would be more logical to take the blocks for processing not from the end, but from the beginning of the message
template<MODE mode>
inline void streebog(uint8_t const * const m, const uint64_t size, uint8_t * const out) {
    constexpr auto block = 64; // 64 bytes
    const auto blocks_n = size / block;
    uint64_t n[8]{}, sum[8]{}, h[8]{}, buff[8]{}, it{}, *out64 = (uint64_t * const)out;
    constexpr static uint64_t zero[8]{};

    get_IV(mode, h); // compile-time filling

    for(;it + block - 1 < size;it += block) {
        auto m64 = (uint64_t const * const)m + (size - block - it >> 3);
        if constexpr (is_little_endian()) { // hash algo created for big endian
            endianess_swap(buff, m64);
            G(n, h, buff, h);
            i512_sum(sum, buff, sum);
        } else {
            G(n, h, m64, h);
            i512_sum(sum, m64, sum);
        }
        n[7] += 0x200; // ATTENTION: the index (7) depends on endianess?
    }

    const auto rem = size - it;
    auto _b = (uint8_t*)buff;

#ifdef METAPROG_UNROLL
    unroll<8>([&]<uint64_t i>(){buff[i] = 0;});
#else
    for(auto i = 0;i < 8;buff[i++] = 0);
#endif

    for(auto i = 0;i < rem; ++i)
        _b[i+64-rem] = m[i];

    _b[64 - rem - 1] = 0x01;

    endianess_swap(buff, buff);
    G(n, h, buff, h);

    n[7] += (rem << 3);

    i512_sum(sum, buff, sum);
    G(zero, h, n, h);

    if constexpr (mode == MODE::H512)
        G(zero, h, sum, out64);

    else if constexpr (mode == MODE::H256) {
        G(zero, h, sum, h);
        out64[0] = h[0], out64[1] = h[1], out64[2] = h[2], out64[3] = h[3];
    }
}


void streebog256(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out) {
    streebog<MODE::H256>(in, bytes_n, out);
}


void streebog512(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out) {
    streebog<MODE::H512>(in, bytes_n, out);
}


void streebog(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out, const MODE mode) {
    return mode == MODE::H256? streebog256(in, bytes_n, out) : streebog512(in, bytes_n, out);
}

