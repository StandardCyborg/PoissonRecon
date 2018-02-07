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
