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
#[[

add_library(streebog STATIC streebog.cc)
target_include_directories(streebog PUBLIC include/)
]]

set(TARGETS stbg stbg512 stbg256)

foreach(target IN LISTS TARGETS streebog)
    target_compile_options(${target} PRIVATE
        -std=c++20
        -O3
        #-Ofast
        #-g
        -march=native
        -mavx2
        -flto
        -fno-exceptions
        -fno-rtti

        #-fno-semantic-interposition
        #-fno-asynchronous-unwind-tables
        #-fno-unwind-tables
        #-nostdlib           # для streebog lib надо оставить. и USE_WRAPPERS разрешить
        #-Wl,--gc-sections,--strip-all
    )
endforeach()

#[[
enable_testing()

add_executable(tests streebog.cc tests/streebog_test.cc )
target_include_directories(tests PUBLIC include/)
add_test(NAME streebog_tests COMMAND tests)
target_compile_options(tests PRIVATE)]]


