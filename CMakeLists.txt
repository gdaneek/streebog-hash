cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TARGETS streebog streebog512 streebog256)

foreach(target IN LISTS TARGETS)
add_executable(${target} streebog.cc examples/example.cc)
target_include_directories(${target} PUBLIC include/)
endforeach()


add_library(gost341118 STATIC streebog.cc)
target_include_directories(gost341118 PUBLIC include/)

target_compile_options(streebog256 PRIVATE -DH256_ONLY)
target_compile_options(streebog512 PRIVATE -DH512_ONLY)

add_executable(tests streebog.cc tests/streebog_test.cc )
target_include_directories(tests PUBLIC include/)


enable_testing()
add_test(NAME streebog_tests COMMAND tests)

foreach(target IN LISTS TARGETS gost341118)
    target_compile_options(${target} PRIVATE
        -std=c++20
        -O3
        -Ofast
        -s
        -march=native
        -flto
        -fno-exceptions
        -fno-rtti
        -fno-semantic-interposition
    )
endforeach()



