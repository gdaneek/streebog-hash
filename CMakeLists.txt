cmake_minimum_required(VERSION 3.20)

project(streebog)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(stbg streebog.cc examples/example.cc)
target_include_directories(stbg PUBLIC include/)


add_executable(stbg512 streebog.cc examples/example.cc)
target_include_directories(stbg512 PUBLIC include/)
target_compile_options(stbg512 PRIVATE)


add_executable(stbg256 streebog.cc examples/example.cc)
target_include_directories(stbg256 PUBLIC include/)
target_compile_options(stbg256 PRIVATE)


add_library(streebog STATIC streebog.cc)
target_include_directories(streebog PUBLIC include/)
target_compile_options(streebog PRIVATE -DSTREEBOG_ENABLE_WRAPPERS)


set(TARGETS stbg stbg512 stbg256)

foreach(target IN LISTS TARGETS streebog)
    target_compile_options(${target} PRIVATE
        -std=c++20
        -O3
        -march=native
        -flto
        -fno-exceptions
        -fno-rtti
        -DDISABLE_MANUAL_AVX # remove this for get full-avx implementation
    )
endforeach()

enable_testing()

add_executable(tests streebog.cc tests/streebog_test.cc )
target_include_directories(tests PUBLIC include/)
add_test(NAME streebog_tests COMMAND tests)
target_compile_options(tests PRIVATE)


