cmake_minimum_required(VERSION 3.20)

project(streebog)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(stbg streebog.cc example/canonical.cc)
target_include_directories(stbg PUBLIC include/)


add_executable(stbg512 streebog.cc example/canonical.cc)
target_include_directories(stbg512 PUBLIC include/)
target_compile_options(stbg512 PRIVATE)


add_executable(stbg256 streebog.cc example/canonical.cc)
target_include_directories(stbg256 PUBLIC include/)
target_compile_options(stbg256 PRIVATE)


add_library(streebog STATIC streebog.cc)
target_include_directories(streebog PUBLIC include/)
target_compile_options(streebog PRIVATE -DSTREEBOG_ENABLE_WRAPPERS)


set(TARGETS stbg stbg512 stbg256)

foreach(target IN LISTS TARGETS streebog)
    target_compile_options(${target} PRIVATE
        -std=c++20
        -O3 -Ofast
        -march=native
        -flto
        -fno-exceptions
        -fno-rtti
        -nostdlib
    )
endforeach()

enable_testing()

add_executable(streebog_test streebog.cc test/streebog_test.cc )
target_include_directories(streebog_test PUBLIC include/)
add_test(NAME streebog_tests COMMAND streebog_test)
target_compile_options(streebog_test PRIVATE)


