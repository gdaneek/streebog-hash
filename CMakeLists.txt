cmake_minimum_required(VERSION 3.20)

project(streebog)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGETS stbg streebog512 streebog256 stbg512-le)

foreach(target IN LISTS TARGETS)
add_executable(${target} streebog.cc examples/example.cc)
target_include_directories(${target} PUBLIC include/)
endforeach()


add_library(streebog STATIC streebog.cc)
target_include_directories(streebog PUBLIC include/)

target_compile_options(streebog256 PRIVATE -DH256_ONLY)
target_compile_options(streebog512 PRIVATE -DH512_ONLY)
target_compile_options(stbg512-le PRIVATE -DH512_ONLY -DSTANDARD_AS_LITTLE_ENDIAN)

add_executable(tests streebog.cc tests/streebog_test.cc )
target_include_directories(tests PUBLIC include/)

target_compile_options(tests PRIVATE -DMETAPROG_UNROLL)

enable_testing()
add_test(NAME streebog_tests COMMAND tests)

foreach(target IN LISTS TARGETS streebog)
    target_compile_options(${target} PRIVATE
        -std=c++20
        -O3
        -Ofast
        # -Os
        -march=native
      #  -funroll-loops
        -flto
        -fno-exceptions
        -fno-rtti
        -fno-semantic-interposition
        -fno-asynchronous-unwind-tables
        -fno-unwind-tables
        -nostdlib           # для streebog lib надо оставить. и USE_WRAPPERS разрешить
        -Wl,--gc-sections,--strip-all
        -DMETAPROG_UNROLL
    )
endforeach()



