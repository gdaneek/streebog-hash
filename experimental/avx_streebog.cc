 
using llp = const long long * __restrict;

__attribute__((always_inline))
inline uint64_t reduce(__m256i v) {
    __m128i x = _mm_xor_si128(_mm256_castsi256_si128(v),
                              _mm256_extracti128_si256(v, 1));
    __m128i s = _mm_unpackhi_epi64(x, x);

    return _mm_cvtsi128_si64(_mm_xor_si128(x, s));
}


void LPSX(__m256i const * __restrict l, __m256i const * __restrict r, __m256i* __restrict out) {
    const __m256i mask = _mm256_set1_epi64x(0xFF), perm = _mm256_set_epi32(0,0,0,0,6,4,2,0);
    const auto lo = _mm_set_epi32(768, 512, 256, 0), ro =  _mm_set_epi32(1792, 1536, 1280, 1024);

    __m256i in0 = _mm256_xor_si256(l[0], r[0]), in1 = _mm256_xor_si256(l[1], r[1]);

    for(int i = 0;i < 8;i++) {
        __m256i bytes0 = _mm256_and_si256(in0, mask), bytes1 = _mm256_and_si256(in1, mask);

        __m128i idx0 = _mm_add_epi32(_mm256_castsi256_si128(_mm256_permutevar8x32_epi32(bytes0, perm)), lo),
                idx1 = _mm_add_epi32(_mm256_castsi256_si128(_mm256_permutevar8x32_epi32(bytes1, perm)), ro);

        ((uint64_t* __restrict)out)[i] = reduce(_mm256_xor_si256(_mm256_i32gather_epi64((llp)mmul_lut.v, idx0, 8),
                                                                 _mm256_i32gather_epi64((llp)mmul_lut.v, idx1, 8)));

        in0 = _mm256_srli_epi64(in0, 8), in1 = _mm256_srli_epi64(in1, 8);
    }
}

__attribute__((always_inline))
inline void X(__m256i const * __restrict l, __m256i const * __restrict r, __m256i* __restrict out) {
    _mm256_store_si256(out, _mm256_xor_si256(l[0], r[0]));
    _mm256_store_si256(out+1, _mm256_xor_si256(l[1], r[1]));
}


void Streebog::G(__m256i * __restrict m, bool is_zero) {
    __m256i K[2], tmp[2];
    static __m256i zeros[2]{};
    memcpy(K, h, 64);

    is_zero? LPSX(K, zeros, K) : LPSX(K, n, K);
    LPSX(K, m, tmp);
    LPSX(K, (__m256i*)C, K);

    for(int i = 0;i < 11;i++) {
        LPSX(K, tmp, tmp);
        LPSX(K, (__m256i*)(C+((i + 1) << 3)), K);
    }

    X(K, tmp, tmp);  // last iteration X only
    X(tmp, h, tmp);
    X(tmp, m, h);
}
