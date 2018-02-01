#ifndef FILES_INCLUDED
#define FILES_INCLUDED

#ifdef USE_MEMORY_FILE_SYSTEM

#include "MemoryFileSystem.h"

#define TFILE MemoryFileSystem::FILE
#define Tfopen MemoryFileSystem::fopen
#define Tfclose MemoryFileSystem::fclose
#define Tfread MemoryFileSystem::fread
#define Tfwrite MemoryFileSystem::fwrite
#define Tfseek MemoryFileSystem::fseek
#define Tfprintf MemoryFileSystem::fprintf
#define Tfgets MemoryFileSystem::fgets

#else

#include <stdio.h>

#define TFILE FILE
#define Tfopen fopen
#define Tfclose fclose
#define Tfread fread
#define Tfwrite fwrite
#define Tfseek fseek
#define Tfprintf fprintf
#define Tfgets fgets

#endif // USE_MEMORY_FILE_SYSTEM

#endif // !FILES_INCLUDED
