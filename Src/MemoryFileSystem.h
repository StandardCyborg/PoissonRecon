#ifndef MEMORY_FILE_SYSTEM_HEADER
#define MEMORY_FILE_SYSTEM_HEADER

#include <stdlib.h>
#include <string>
#include <map>

namespace MemoryFileSystem
{
	class FILE
	{
		public:
			FILE(std::string name);
    
			std::string GetName();

		private:
			std::string Name;
	};

	class InternalFILE
	{
		public:
			InternalFILE();
			~InternalFILE();

			unsigned char *buffer;
			size_t iBufferSize;
			size_t iLocation;
			size_t iFileSize;

			bool Reading;
			bool Writing;
			bool Binary;

			FILE *fp; // NULL if not open
	};

	class MemoryFileDrive
	{
		public:
			MemoryFileDrive();
			~MemoryFileDrive();

			// Creates a file that will reside in memory
			FILE *CreateVirtualFile(const char *filename, const char *mode);

			// will allocate memory for a virtual file or realloc if there is already a memory block
			unsigned char *ReAllocBuffer(unsigned char *, size_t size, size_t new_size);

			bool HasFile(const char *filename);

			int CloseFile(FILE *fp);

			// this is used to keep track of the files and to make sure we don't have any memory leaks
			std::map<std::string, InternalFILE> m_memoryFileMap;

			static const int defaultbufferSize = 1024 * 1024;

		private:
			void CreateNewInternalFile(std::string FileName);
	};

	FILE *fopen(const char *filename, const char *mode);

	int fclose(FILE *);

	size_t fwrite(const void * ptr, size_t size, size_t count, FILE * stream);

	size_t fread(void * ptr, size_t size, size_t count, FILE * stream);

	int fseek(FILE * stream, long int offset, int origin);

	int ftell(FILE * stream);

	int fprintf(FILE * stream, const char * format, ...);

	char * fgets(char * str, int num, FILE * stream);
	
	char *_mktemp(char *);

	void WriteFileInMemoryToDisc(const char *filename);

	void WriteFileInMemoryToStdout(const char *filename);
};

#endif