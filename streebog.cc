/**
 * @file    streebog.cc
 * @brief   Implementation of GOST 34.11-2018 hash functions 256 and 512 bits
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.3.1
 * @see https://github.com/gdaneek/streebog-hash
 */

#include "streebog.hh"
#include "constants.hh"

#include <string.h>  // for memset memcpy

#include <array>
#include <type_traits>  // for metaprog templates
#include <utility>      // for index sequences


using ui64 = uint64_t;

template <uint64_t... I>
using is = std::index_sequence<I...>;

template <uint64_t I>
using make_is = std::make_index_sequence<I>;

consteval auto mmul_A_precalc() {
  std::array<std::array<ui64, 256>, 8> out{};
  for (int i{}; i < 8; i++) {
    for (int j{}; j < 256; j++) {
      for (int t{}; t < 8; t++) {
        out[i][j] ^= ((!(pi[j] & (1 << t)) - 1) & m_A[63 - t - (i << 3)]);
      }
    }
  }

  return out;
}

inline constexpr auto mmul_lut = mmul_A_precalc();

void Streebog::reset() {
  memcpy((void*)h, (void*)(IV + (mode == Streebog::Mode::H512 ? 0 : 8)), 64);
  memset(n, 0, sizeof(n));
  memset(sum, 0, sizeof(sum));
}

Streebog::Streebog(const Mode _mode) : mode{_mode} { this->reset(); }

inline void vadd512(void* _a, void* _b, void* __restrict _dst) {
  ui64 *a = (ui64*)_a, *b = (ui64*)_b, *dst = (ui64*__restrict)_dst;
  bool carry{};
  [&]<ui64... Is>(is<Is...>) {
    (([&](const ui64 I) {
       ui64 tmp = a[I] + b[I] + carry;
       carry = (tmp < a[I]) | ((tmp == a[I]) & carry);
       dst[I] = tmp;
     } (Is)), ...);
  } (make_is<8>());
}

inline void LPSX(ui64 const* __restrict lhs, ui64 const* __restrict rhs, ui64* __restrict out) {
  alignas(32) ui64 r[8];
  [&]<ui64... I>(is<I...>) {
    ((r[I] = lhs[I] ^ rhs[I], out[I] = 0), ...);
  } (make_is<8>());

  [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
    ((out[I >> 3] ^= mmul_lut[I & 7][(uint8_t)r[I & 7]],
      r[I & 7] >>= 8), ...);
  } (make_is<64>());
}

void Streebog::G(ui64 const* __restrict m, bool is_zero) {
  alignas(32) ui64 K[8], tmp[8], zeros[8]{};
  memcpy(K, h, 64);

  is_zero ? LPSX(K, zeros, K) : LPSX(K, n, K);
  LPSX(K, m, tmp);
  LPSX(K, C, K);

  [&]<ui64... I>(is<I...>) {
    ((LPSX(K, tmp, tmp),
      LPSX(K, C + ((I + 1) << 3), K)), ...);
  } (make_is<11>());

  [&]<ui64... I>(is<I...>) __attribute__((always_inline)) {
    ((tmp[I] ^= K[I]), ...);
    ((tmp[I] ^= m[I]), ...);
    ((h[I] ^= tmp[I]), ...);
  } (make_is<8>());
}

void Streebog::update(void* __restrict m, const ui64 size) {
  for (ui64 i{}; i < (size >> 6); i++) {
    G((ui64*)m + (i << 3));
    vadd512(sum, (ui64*)m + (i << 3), sum);
    *(uint64_t*)n += 0x200;
  }
}

ui64 const* const Streebog::finalize(void* __restrict m, const ui64 size) {
  alignas(32) uint64_t buff[8]{};
  const ui64 _d = size & ~0x3FULL;
  auto rem = size - _d;
  this->update(m, _d);                      // process whole chunks
  memcpy((void*)buff, (char*)m + _d, rem);  // pad input
  ((uint8_t*)buff)[rem] = 0x01;
  G(buff);
  *(ui64*)n += (rem << 3);
  vadd512(sum, buff, sum);  // last step
  G(n, true), G(sum, true);

  return (ui64 const* const)(this->h);
}

ui64 const* const Streebog::operator()(void* m, const ui64 size, void* out) {
  auto ret = this->finalize(m, size);
  if (out != nullptr) {
    auto ret_offset = (mode == Mode::H512 ? 0 : 4);
    auto bytes_n = (mode == Mode::H512 ? 8 : 4);
    memcpy(out, ret + ret_offset, bytes_n << 3);
  }

  return ret;
}
