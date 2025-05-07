 
#include "streebog.hh"

#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

int main(int argc, char* argv[]) {

    int fd = open(argv[1], O_RDONLY);
    if (fd == -1)
        return 1;

    uint64_t file_size = lseek(fd, 0, SEEK_END);
    lseek(fd, 0, SEEK_SET);

    void* mapped = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        return 1;
    }

    uint64_t out[8];
    #ifdef H512_ONLY
    streebog512((const uint8_t*)mapped, file_size, (uint8_t*)out);
    constexpr int real_sz = 8;
    #elif defined(H256_ONLY)
    streebog256((const uint8_t*)mapped, file_size, (uint8_t*)out);
    constexpr int real_sz = 4;
    #else
    auto mode = (MODE)atoi(argv[2]);
    streebog((const uint8_t*)mapped, file_size, (uint8_t*)out, mode);
    int real_sz = (mode == MODE::H256)? 4 : 8;
    #endif

    for(int i = 0; i < real_sz; i++)
        printf("%016" PRIx64, out[i]);
    printf("\n");

    return 0;
}
