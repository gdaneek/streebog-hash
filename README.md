# Streebog hash

High-performance implementation of GOST 34.11-2012 standards (also GOST 34.11-2018).

### Installation

Clone the repository and use cmake to build the project:

```
git clone https://github.com/gdaneek/streebog-hash.git
cd streebog-hash
mkdir build
cd build
cmake ..
cmake --build .
```

After the build, you will receive
- executable *stbg512* containing only a 512-bit hash, 
- executable *stbg256* containing only a 256-bit hash, 
- executable *stbg* with both modes of operation
- static library streeboglib.a, which you can connect to your project. 
- Tests based on control examples are compiled into the *tests* executable file.

### Developer documentation

See `docs/code` subfolder or generate it yourself using Doxygen:
```
cd docs
doxygen
```

### Benchmarks

See [benchmarks](docs/benchmarks.md)

### Advanced settings

The header file contains some wrappers for the class that simplify working with the hash as an STL container. Use the `STREEBOG_ENABLE_WRAPPERS` macro to allow them.

Also use the macro `USE_MANUAL_AVX` to get a fully AVX/AVX2 vectorized version of the programs (experimental  in the implementations of versions v2.2+)

