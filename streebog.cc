/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 1.0
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

// TODO: all with unrolled without make_is and is

/**
 * @brief an abbreviated form from index sequence
 */
template<size_t... I> using is = std::index_sequence<I...>;

/**
 * @brief an abbreviated form from make index sequence
 */
template<size_t N> using make_is = std::make_index_sequence<N>; //<

/**
 * @brief
 *
 */
template<size_t F, typename funcT, typename... Args>
inline constexpr auto unrolled(funcT&& func, Args&&... args) {
    [&]<size_t... I>(std::index_sequence<I...>) __attribute__((always_inline)) {
        ((func(I, std::forward<Args>(args)...)), ... );
    } (std::make_index_sequence<F>());
}

/**
 * @brief Implementation of the S transformation (pi substitution)
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @note Use UNROLL = true parameter to get a unroll of all cycles of the compilation stage
 * @warning If you are not using a UNROLL, only a single-byte data type is allowed as input.
 * use reinterpret_cast to uint8_t* before function call to avoid mistakes
 */
template<typename T, bool UNROLL = true> requires (UNROLL || (sizeof(T) == 1))
constexpr inline auto S(T * const in) {
    if constexpr (UNROLL) {
        [in]<size_t... Is>(is<Is...>) {
            ((in[Is] = []<size_t... I>(/* T* */uint8_t* v, is<I...>)  {  // __attribute__((always_inline))
                //return (((T)pi[(v >> (I << 3)) & 0xFF] << (I << 3)) | ...);
                return (((T)pi[v[I]] << (I << 3)) | ...);
            } ((uint8_t*)(in+Is), make_is<sizeof(T)>())), ...);
        } (make_is<0x40 / sizeof(T)>());

        return in;
    }

    for(uint8_t i{};i < 0x40; ++i)
        in[i] = pi[in[i]];              // substitution by pi definition

    return in;
}

/**
 * @brief Implementation of the P transformation (r permutation)
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @note Use UNROLL = true parameter to get a unroll of all cycles of the compilation stage
 * @warning If you are not using a UNROLL, only a single-byte data type is allowed as input.
 * use reinterpret_cast to uint8_t* before function call to avoid mistakes
 */
template<typename T, bool UNROLL = true> requires (UNROLL || (sizeof(T) == 1))
constexpr inline auto P(T * const in) {
    constexpr auto SZ = 64 / sizeof(T);
    T b[SZ]{};
    if constexpr (UNROLL) {
        [in]<size_t... I>(T * const out, is<I...>) {
            ([]<size_t... B>(T* in, T* out, const size_t i, is<B...>) {
                out[i] = ((((in[B] >> (56 - 8*i)) & 0xFF) << (56 - 8*B)) | ...);
                //out[i] = ((((in[B] >> (8*i)) & 0xFF) << (56 - 8*B)) | ...);
            } (in, out, I, make_is<sizeof(T)>()), ...);
        } (b, make_is<SZ>());

        unrolled<SZ>([&](const size_t i){ in[i] = b[i]; });

        return in;
    }

    for(uint8_t i{};i < 8; ++i)
        for(uint8_t j{};j < 8;++j)
            b[63-i] = in[(j << 3) + i];

    for(uint8_t i{};i < SZ; ++i)
        in[i] = b[i];

    return in;
}


/**
 * @brief Implementation of the L transformation (GF2 mmul)
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @note Use UNROLL = true parameter to get a unroll of all cycles of the compilation stage
 * @warning If you are not using a UNROLL, only a eight-byte data type is allowed as input.
 * use reinterpret_cast to uint64_t* before function call to avoid mistakes
 */
template<typename T, bool UNROLL = true> requires (UNROLL || (sizeof(T) == 1))
constexpr inline auto L(T * const in) {
    if constexpr (UNROLL) {
        [in]<size_t... Is>(is<Is...>) {
            ((in[Is] = []<size_t... I>(T v, is<I...>) {
                return ((((T)!(v & (1ull << I))-1) & m_A[63-I]) ^ ...);
            } (in[Is], make_is<64>())), ...);
        } (make_is<8>());

        return in;
    }

    for(auto i = 0;i < 8;i++) {
        T res{};
        for(auto j = 0;j < 64;j++)
            res ^= ((T)!(in[i] & (1ull << j))-1) & m_A[63-j];
        in[i] = res;
    }

    return in;
}

/**
 * @brief Wrapper for calls to frequently used LPS transformations
 * @param in 512 bit array on which the transformation must be performed
 * @return pointer to the converted array
 * @note Use UNROLL = true parameter to get a unroll of all cycles of the compilation stage
 */
template<typename T, bool UNROLL = true>
constexpr inline auto LPS(T* in) {
    return L<T, UNROLL>(P<T, UNROLL>(S<T, UNROLL>(in)));
}


/**
 * @brief Implementation of the X transformation
 * @param l sizeof(T) * Len bytes array on which the transformation must be performed
 * @param r sizeof(T) * Len bytes array on which the transformation must be performed
 * @param o sizeof(T) * Len bytes array to write out the result
 * @return pointer to the output array
 * @note Use UNROLL = true parameter to get a unroll of all cycles of the compilation stage
 * @note in fact, it's no different from xor operation, so it is often called instead of it
 */
template<typename T, bool UNROLL = true, size_t Len = 8>
constexpr inline auto X(T const * const l, T const * const r, T * const o) {
    if constexpr (UNROLL) {
        unrolled<Len>([&](const size_t i){ o[i] = l[i] ^ r[i]; });
        return o;
    }

    for(uint8_t i{};i < Len;++i)
        o[i] = l[i] ^ r[i];

    return o;
}


/**
 * @brief Wrapper for calls to frequently used LPSX transformations
 * @param l sizeof(T) * Len bytes array on which the transformation must be performed
 * @param r sizeof(T) * Len bytes array on which the transformation must be performed
 * @param o sizeof(T) * Len bytes array to write out the result
 * @note Use UNROLL = true parameter to get a unroll of all cycles of the compilation stage
 */
template<typename T, bool UNROLL = true, size_t Len = 8>
constexpr inline auto LPSX(T const * const l, T const * const r, T * const o) {
    return LPS<T, UNROLL>(X<T, UNROLL, Len>(l, r, o));
}

/**
 * @brief a function that provides a pointer to the desired initialization vector (IV)
 * @param mode hashing mode of operation (H256 -- 256bit or H512 -- 512bit)
 * @return pointer to 512 bit IV
 */
inline constexpr auto get_IV(const MODE mode) {
    switch(mode) {
        case MODE::H512: return IV;
        case MODE::H256: return IV + 8;
    };

    return IV;
}


inline constexpr auto get_IV(const MODE mode, uint64_t * const out) {
    auto IV = get_IV(mode);
    unrolled<8>([&](const size_t i){ out[i] = IV[i]; });
}


/**
 * @brief Implementation of the E transformation
 * @param K pointer to initial K value (to 512 bit array)
 * @param m operand passed with ()
 * @param out array to write out the result
 * @return pointer to the output array
 */
template<typename T, bool UNROLL = true, size_t Len = 8>
constexpr inline auto E(const T * const K, T const * const m, T * const out) {
    T b[64 / sizeof(T)]{}, curr_K[64 / sizeof(T)]{};

    LPSX<T, UNROLL, Len>(K, m, b);
    LPSX<T, UNROLL, Len>(K, C, curr_K);

    if constexpr (UNROLL) {
        unrolled<11>([&](const size_t i) {
            LPSX<T, UNROLL, Len>(curr_K, b, b);
            LPSX<T, UNROLL, Len>(curr_K, C+((i+1) << 3), curr_K);
        });
    } else {
        for(auto i = 1;i < 12;++i) {
            LPSX<T, UNROLL, Len>(curr_K, b, b);
            LPSX<T, UNROLL, Len>(curr_K, C+(i << 3), curr_K);
        }
    }

    return X<T, UNROLL, Len>(curr_K, b, out);  // last iteration X only
}


/**
 * @brief Implementation of the G transformation
 * @param N pointer to current N value (to 512 bit array)
 * @param h pointer to current h value (to 512 bit array)
 * @param out array to write out the result
 * @return pointer to the output array
 */
template<typename T, bool UNROLL = true, size_t Len = 8>
constexpr inline auto G(const T * const N, T const * const h, T const * const m, T * const out) {
    T b[64 / sizeof(T)]{};
    LPSX<T, UNROLL, Len>(h, N, b);
    E<T, UNROLL, Len>(b, m, b);
    return X<T, UNROLL, Len>(X<T, UNROLL, Len>(b, h, b), m, out);
}


/**
 * @brief Implementation of 512-bit addition
 * @param l lhs 512-bit operand
 * @param r rhs 512-bit operand
 * @param o array to write out the result
 * @return pointer to the output array
 */
constexpr inline auto i512_sum(const uint64_t * const l, uint64_t const * const r, uint64_t * const o) {    // TODO: ATTENTION: remove  all size_t and it lib
    uint64_t res[8]{};
    for (int i = 7; i > -1; --i) {
        uint64_t ai = l[i], bi = r[i];
        uint64_t sum = ai + bi;
        auto carry = (sum < ai) || (sum < bi) || (sum == (uint64_t)-1);
        res[i] += sum;
        if(i != 0)
            res[i-1] += carry;

    }

    unrolled<8>([&](const size_t i){ o[i] = res[i]; });
}

constexpr inline bool is_little_endian() { return (((uint32_t)0x01020304 & 0xFF) == 0x04); }


// TODO: rewrite code to remove swaps. idk why devs of the algorithm chose such an inconvenient endianess
//  in addition, it would be more logical to take the blocks for processing not from the end, but from the beginning of the message
template<MODE mode>
inline void streebog(uint8_t const * const m, const uint64_t size, uint8_t * const out) noexcept {
    constexpr auto block_sz = 64; // 64 bytes
    const auto blocks_n = size / block_sz;
    uint64_t n[8]{}, sum[8]{}, h[8]{}, buff[8]{}, m_it{};
    constexpr static uint64_t zero[8]{};

    get_IV(mode, h); // compile-time filling

    for(;m_it + (block_sz-1) < size;m_it += block_sz) {
        // __builtin_prefetch(&array[i + 4]);
        auto m64 = reinterpret_cast<uint64_t const * const>(m) + ((size - block_sz - m_it) >> 3);
        if constexpr (is_little_endian()) {
            unrolled<8>([&](const size_t i){ buff[i] = __builtin_bswap64(m64[i]); });
            G(n, h, buff, h);
            i512_sum(sum, buff, sum);
        } else {
            G(n, h, m64, h);
            i512_sum(sum, m64, sum);
        }

        n[7] += 0x200;
    }

    const auto rem = size - m_it;
    auto _buff = reinterpret_cast<uint8_t*>(buff);

    unrolled<8>([&](const size_t i){ buff[i] = 0; });
    for(auto i = 0;i < rem;_buff[i + 64-rem] = m[i], ++i);

    _buff[64 - rem - 1] = 0x01;

    unrolled<8>([&](char i){ buff[i] = __builtin_bswap64(buff[i]); });
    G(n, h, buff, h);

    n[7] += (rem << 3);

    i512_sum(sum, buff, sum);
    G(zero, h, n, h);

    if constexpr (mode == MODE::H512)
        G(zero, h, sum, reinterpret_cast<uint64_t * const>(out));
    else if constexpr (mode == MODE::H256) {
        G(zero, h, sum, h);
        unrolled<4>([&](char i){ ((uint64_t * const)(out))[i] = h[i]; });
    }
}


void streebog256(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out) noexcept {
    streebog<MODE::H256>(in, bytes_n, out);
}


void streebog512(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out) noexcept {
    streebog<MODE::H512>(in, bytes_n, out);
}


void streebog(uint8_t const * const in, const uint64_t bytes_n, uint8_t * const out, const MODE mode) noexcept {
    if(mode == MODE::H256)
        return streebog256(in, bytes_n, out);

   return streebog512(in, bytes_n, out);
}

