/**
 * @file    streebog_test.cc
 * @brief   GOST 34.11-2018 hash functions 256 and 512 bits tests based on control examples from the document application
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 2.0
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "streebog.hh"

const uint8_t small_m[] = { // input msg - control examples 1, 3 from the document application
       0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
       0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
       0x36, 0x37, 0x38, 0x39, 0x30, 0x31, 0x32, 0x33,
       0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x30, 0x31,
       0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39,
       0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
       0x38, 0x39, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
       0x36, 0x37, 0x38, 0x39, 0x30, 0x31, 0x32
};


const uint8_t big_m[] = { // input msg - control examples 2, 4 from the document application
  0xd1, 0xe5, 0x20, 0xe2, 0xe5, 0xf2, 0xf0, 0xe8,
  0x2c, 0x20, 0xd1, 0xf2, 0xf0, 0xe8, 0xe1, 0xee,
  0xe6, 0xe8, 0x20, 0xe2, 0xed, 0xf3, 0xf6, 0xe8,
  0x2c, 0x20, 0xe2, 0xe5, 0xfe, 0xf2, 0xfa, 0x20,
  0xf1, 0x20, 0xec, 0xee, 0xf0, 0xff, 0x20, 0xf1,
  0xf2, 0xf0, 0xe5, 0xeb, 0xe0, 0xec, 0xe8, 0x20,
  0xed, 0xe0, 0x20, 0xf5, 0xf0, 0xe0, 0xe1, 0xf0,
  0xfb, 0xff, 0x20, 0xef, 0xeb, 0xfa, 0xea, 0xfb,
  0x20, 0xc8, 0xe3, 0xee, 0xf0, 0xe5, 0xe2, 0xfb
};


template<typename T>
bool equal(T const * const lhs, const uint64_t lhs_sz, T const * const rhs) {   // i dont want use <algorithm> for only this func
    bool flag = true;
    for(auto i = 0ull;i < lhs_sz;++i)
        flag &= (lhs[i] == rhs[i]);

    return flag;
}


TEST_SUITE("hash 512-bit") {

    TEST_CASE("small message") {

        uint64_t expected[] = {
            0xd5b9f54a1ad0541b,
            0x6254288dd6863dcc,
            0x352f227524bc9ab1,
            0xfa1fbae42b1285c0,
            0x823a7b76f830ad00,
            0x11c324f074654c38,
            0x7fef082b3381a4e2,
            0x486f64c191787941
        };

        Streebog stbg{Streebog::Mode::H512};
        auto res = stbg.finalize(small_m, sizeof(small_m));

        REQUIRE(equal(expected, 8, res));
    }

    TEST_CASE("big message") {

        uint64_t expected[] = {
            0x6fcabf2622e6881e,
            0xe06915d5f2f19499,
            0x1ae60f3b5a47f8da,
            0x7613966de4ee0053,
            0xb8a2ad4935e85f03,
            0xb3e56c497ccd0f62,
            0x60642bdcddb90c3f,
            0x28fbc9bada033b14
        };

        Streebog stbg{Streebog::Mode::H512};
        auto res = stbg.finalize(big_m, sizeof(big_m));

        REQUIRE(equal(expected, 8, res));
    }
}


TEST_SUITE("hash 256-bit") {

    TEST_CASE("small message") {

        uint64_t expected[] = {
            0x890b59d8ef1e159d,
            0x27f94ab76cbaa6da,
            0xa449b16b0251d05d,
            0x00557be5e584fd52
        };

        Streebog stbg{Streebog::Mode::H256};
        auto res = stbg.finalize(small_m, sizeof(small_m));

        REQUIRE(equal(expected, 4, res+4));
    }

    TEST_CASE("big message") {

        uint64_t expected[] = {
            0x5d9e40904efed29d,
            0xb005746d97537fa8,
            0x749a66fc28c6cac0,
            0x508f7e553c06501d
        };

        Streebog stbg{Streebog::Mode::H256};
        auto res = stbg.finalize(big_m, sizeof(big_m));

        REQUIRE(equal(expected, 4, res+4));
    }
}


