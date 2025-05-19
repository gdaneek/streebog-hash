 
#include "streebog.hh"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <cinttypes>

#define CHUNK_SIZE 1048576

int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "Использование: %s <файл>\n", argv[0]);
        return 1;
    }

    int fd = open(argv[1], O_RDONLY | O_DIRECT);
    if (fd == -1) {
        fprintf(stderr, "Ошибка открытия файла: %s\n", strerror(errno));
        return 1;
    }

    char* buf;
    if (posix_memalign((void**)&buf, 32, CHUNK_SIZE)) {
        perror("posix_memalign");
        close(fd);
        return 1;
    }

    Streebog stbg(Streebog::Mode::H512);

    ssize_t bytes_read;

    uint64_t hash[8];

    while ((bytes_read = read(fd, buf, CHUNK_SIZE)) > 0) {
        if(bytes_read < CHUNK_SIZE) {
            auto ret = stbg.finalize((const uint8_t*)buf, bytes_read);
            for(auto i = 0;i < 8;hash[i++] = ret[i]);
        } else stbg.update(buf, bytes_read);


    }

    if (bytes_read == -1) {
        perror("read");
        free(buf);
        close(fd);
        return 1;
    }

    for(int i = 7; i > -1; i--)
        printf("%016" PRIx64, hash[i]);
    printf("\n");

    free(buf);
    close(fd);
    return 0;
}
