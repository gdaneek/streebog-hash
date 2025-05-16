/**
 * @file    streebog_test.cc
 * @brief   GOST 34.11-2018 hash functions 256 and 512 bits tests based on control examples from the document application
 * @author  https://github.com/gdaneek
 * @date    30.05.2025
 * @version 1.0
 * @see https://github.com/gdaneek/GOST-34.11-2018
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "streebog.hh"

const uint8_t small_m[] = { // input msg - control examples 1, 3 from the document application
         0x32,0x31,0x30,0x39,0x38,0x37,0x36,
    0x35,0x34,0x33,0x32,0x31,0x30,0x39,0x38,
    0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,
    0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,
    0x31,0x30,0x39,0x38,0x37,0x36,0x35,0x34,
    0x33,0x32,0x31,0x30,0x39,0x38,0x37,0x36,
    0x35,0x34,0x33,0x32,0x31,0x30,0x39,0x38,
    0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30
};


const uint8_t big_m[] = { // input msg - control examples 2, 4 from the document application
    0xfb,0xe2,0xe5,0xf0,0xee,0xe3,0xc8,0x20,
    0xfb,0xea,0xfa,0xeb,0xef,0x20,0xff,0xfb,
    0xf0,0xe1,0xe0,0xf0,0xf5,0x20,0xe0,0xed,
    0x20,0xe8,0xec,0xe0,0xeb,0xe5,0xf0,0xf2,
    0xf1,0x20,0xff,0xf0,0xee,0xec,0x20,0xf1,
    0x20,0xfa,0xf2,0xfe,0xe5,0xe2,0x20,0x2c,
    0xe8,0xf6,0xf3,0xed,0xe2,0x20,0xe8,0xe6,
    0xee,0xe1,0xe8,0xf0,0xf2,0xd1,0x20,0x2c,
    0xe8,0xf0,0xf2,0xe5,0xe2,0x20,0xe5,0xd1
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
            0x486f64c191787941,
            0x7fef082b3381a4e2,
            0x11c324f074654c38,
            0x823a7b76f830ad00,
            0xfa1fbae42b1285c0,
            0x352f227524bc9ab1,
            0x6254288dd6863dcc,
            0xd5b9f54a1ad0541b
        }, out[8]{};

        streebog512(small_m, sizeof(small_m), (uint8_t*)out);
        REQUIRE(equal(expected, 8, out));
    }

    TEST_CASE("big message") {

        uint64_t expected[] = {
            0x28fbc9bada033b14,
            0x60642bdcddb90c3f,
            0xb3e56c497ccd0f62,
            0xb8a2ad4935e85f03,
            0x7613966de4ee0053,
            0x1ae60f3b5a47f8da,
            0xe06915d5f2f19499,
            0x6fcabf2622e6881e
        }, out[8]{};

        streebog512((const uint8_t*)big_m, sizeof(big_m), (uint8_t*)out);
        REQUIRE(equal(expected, 8, out));
    }
}


TEST_SUITE("hash 256-bit") {

    TEST_CASE("small message") {

        uint64_t expected[] = {
            0x00557be5e584fd52,
            0xa449b16b0251d05d,
            0x27f94ab76cbaa6da,
            0x890b59d8ef1e159d
        }, out[4]{};

        streebog256((const uint8_t*)small_m, sizeof(small_m), (uint8_t*)out);
        REQUIRE(equal(expected, 4, out));
    }

    TEST_CASE("big message") {

        uint64_t expected[] = {
            0x508f7e553c06501d,
            0x749a66fc28c6cac0,
            0xb005746d97537fa8,
            0x5d9e40904efed29d
        }, out[4]{};

        streebog256(big_m, sizeof(big_m), (uint8_t*)out);
        REQUIRE(equal(expected, 4, out));
    }
}




