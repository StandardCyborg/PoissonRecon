#ifndef FILES_INCLUDED
#define FILES_INCLUDED

#define USE_MEMORY_FILE_SYSTEM

#ifdef USE_MEMORY_FILE_SYSTEM

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

#else

#include <stdio.h>

#define TFILE FILE
#define Tfopen fopen
#define Tfclose fclose
#define Tfread fread
#define Tfwrite fwrite
#define Tfseek fseek
#define Tftell ftell
#define Tfprintf fprintf
#define Tfgets fgets
#define Tmktemp _mktemp

#endif // USE_MEMORY_FILE_SYSTEM

#endif // !FILES_INCLUDED
