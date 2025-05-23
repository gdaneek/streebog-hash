
#include "streebog.hh"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cinttypes>

int main(int argc, char** argv) {

    int fd = open(argv[1], O_RDONLY);
    if (fd == -1) {
        perror("Ошибка при открытии файла");
        exit(EXIT_FAILURE);
    }

    off_t file_size = lseek(fd, 0, SEEK_END);
    if (file_size == -1) {
        perror("Ошибка при получении размера файла");
        close(fd);
        exit(EXIT_FAILURE);
    }

    if (lseek(fd, 0, SEEK_SET) == -1) {
        close(fd);
        exit(EXIT_FAILURE);
    }

    void *mapped_data = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);

    if (mapped_data == MAP_FAILED) {
        close(fd);
        exit(EXIT_FAILURE);
    }

    close(fd);

    alignas(32) uint64_t hash[8];
    Streebog{Streebog::Mode::H512}(mapped_data, file_size, hash);

    for(int i = 7; i > -1; i--)
        printf("%016" PRIx64, hash[i]);
    printf("\n");

    if (munmap(mapped_data, file_size) == -1) {
        perror("Ошибка при освобождении памяти");
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
