// DEV: This previously allowed swapping `TFILE` etc with their disk-based equivalents but we're avoiding those explicitly
#include "MemoryFileSystem.h"

#define TFILE MemoryFileSystem::FILE
#define Tfopen MemoryFileSystem::fopen
#define Tfclose MemoryFileSystem::fclose
#define Tfread MemoryFileSystem::fread
#define Tfwrite MemoryFileSystem::fwrite
#define Tfseek MemoryFileSystem::fseek
#define Tftell MemoryFileSystem::ftell
#define Tfprintf MemoryFileSystem::fprintf
#define Tfgets MemoryFileSystem::fgets
#define Tmktemp MemoryFileSystem::_mktemp

#define FILE fprintf(stderr, "Reached FILE unexpectedly"); exit(1); FILE
#define fopen fprintf(stderr, "Reached fopen unexpectedly"); exit(1); fopen
#define fclose fprintf(stderr, "Reached fclose unexpectedly"); exit(1); fclose
#define fread fprintf(stderr, "Reached fread unexpectedly"); exit(1); fread
#define fwrite fprintf(stderr, "Reached fwrite unexpectedly"); exit(1); fwrite
#define fseek fprintf(stderr, "Reached fseek unexpectedly"); exit(1); fseek
#define ftell fprintf(stderr, "Reached ftell unexpectedly"); exit(1); ftell
#define fprintf fprintf(stderr, "Reached fprintf unexpectedly"); exit(1); fprintf
#define fgets fprintf(stderr, "Reached fgets unexpectedly"); exit(1); fgets
#define mktemp fprintf(stderr, "Reached mktemp unexpectedly"); exit(1); mktemp
