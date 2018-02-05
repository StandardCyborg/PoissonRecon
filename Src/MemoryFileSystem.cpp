#include <map>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>

#if !defined(WIN32)
#include <unistd.h>
#else
#include <Fcntl.h>
#include <io.h>
#endif

#include "MemoryFileSystem.h"

using namespace std;

namespace MemoryFileSystem
{
	int GetSizeNeeded(const char * format, va_list args);

	int tempCount = 0;
};

// This is global so it will be created when the program starts and auto-deleted when it exits freeing memory
MemoryFileSystem::MemoryFileDrive GlobalDrive;

// -------------------------------------------------------------------------

MemoryFileSystem::FILE::FILE(string name)
	:Name(name)
{
}

string MemoryFileSystem::FILE::GetName()
{
	return Name;
}

// ------------------------------------------------------------------
MemoryFileSystem::InternalFILE::InternalFILE()
{
	buffer = NULL;
	iBufferSize = 0;
	iLocation = 0;
	iFileSize = 0;
	Reading = false;
	Writing = false;
	fp = NULL;
}

MemoryFileSystem::InternalFILE::~InternalFILE()
{
	if (NULL != buffer) delete[] buffer;
}

// ------------------------------------------------------------------


MemoryFileSystem::MemoryFileDrive::MemoryFileDrive()
{
}

MemoryFileSystem::MemoryFileDrive::~MemoryFileDrive()
{
}

bool MemoryFileSystem::MemoryFileDrive::HasFile(const char *filename)
{
	std::string FileName(filename);

	return m_memoryFileMap.count(FileName) > 0;
}

MemoryFileSystem::FILE *MemoryFileSystem::MemoryFileDrive::CreateVirtualFile(const char *filename, const char *mode)
{
	std::string FileName(filename);
	std::string Mode(mode);

	// make sure the file is not currently open
	if (m_memoryFileMap.count(FileName) > 0)
	{
		if (m_memoryFileMap[FileName].Reading || m_memoryFileMap[FileName].Writing)
		{
			// the file is already open
			return NULL;
		}
	}

	// reading so file must exist
	if (Mode.find('r') != string::npos)
	{
		if (m_memoryFileMap.count(FileName) > 0)
		{
			// file is open for reading
			m_memoryFileMap[FileName].iLocation = 0;
			m_memoryFileMap[FileName].Reading = true;
		}
		else
		{
			// see if the file exist on the disk and open it if so
			::FILE *fp = ::fopen(filename, "rb");

			if (fp != NULL)
			{
				m_memoryFileMap.insert(std::make_pair(FileName, InternalFILE()));

				::fseek(fp, 0, SEEK_END);
				size_t fileSize = ftell(fp);
				::fseek(fp, 0, SEEK_SET);

				unsigned char *buffer = new unsigned char[fileSize];
				
				for (size_t i = 0; i < fileSize; ++i)
				{
					buffer[i] = ::fgetc(fp);
				}

				m_memoryFileMap[FileName].iBufferSize = fileSize;
				m_memoryFileMap[FileName].buffer = buffer;
				m_memoryFileMap[FileName].iLocation = 0;
				m_memoryFileMap[FileName].iFileSize = fileSize;
				m_memoryFileMap[FileName].Reading = true;
				m_memoryFileMap[FileName].iLocation = 0;
			}
			else
			{
				// file not in memory or on the disk
				return NULL;
			}
		}
	}
	// writing, overwrite if exist
	else if (Mode.find('w') != string::npos)
	{
		if (m_memoryFileMap.count(FileName) > 0)
		{
			// leave the buffer alone, just change the file size which will 'delete' the data
			m_memoryFileMap[FileName].iFileSize = 0;
			m_memoryFileMap[FileName].iLocation = 0;
		}
		else
		{
			if (Mode.find('x') != string::npos) return NULL; // 'x' means fail if file doesn't exist

			CreateNewInternalFile(FileName);
		}

		m_memoryFileMap[FileName].Writing = true;
	}
	// append, create if it doesn't exist
	else if (Mode.find('a') != string::npos)
	{
		if (m_memoryFileMap.count(FileName) > 0)
		{
			size_t newbufferSize = m_memoryFileMap[FileName].iBufferSize + defaultbufferSize;
			m_memoryFileMap[FileName].iLocation = m_memoryFileMap[FileName].iFileSize;

			// grow the buffer a little
			m_memoryFileMap[FileName].buffer = ReAllocBuffer(m_memoryFileMap[FileName].buffer,
															 m_memoryFileMap[FileName].iBufferSize,
															 newbufferSize);
		}
		else
		{
			CreateNewInternalFile(FileName);
		}

		m_memoryFileMap[FileName].Writing = true;
	}

	// should this be read/write
	if (Mode.find('+') != string::npos)
	{
		m_memoryFileMap[FileName].Writing = true;
		m_memoryFileMap[FileName].Reading = true;
	}

	// binary
	if (Mode.find('b') != string::npos)
	{
		m_memoryFileMap[FileName].Binary = true;
	}

	// at this we have an internal file and need to create a FILE *to point to it
	FILE *fp = new FILE(FileName);

	// store it so we can free memory if the programmer forgot
	m_memoryFileMap[FileName].fp = fp;

	return fp;
}

void MemoryFileSystem::MemoryFileDrive::CreateNewInternalFile(std::string FileName)
{
	// create a new file
	unsigned char *buffer = new unsigned char[defaultbufferSize];

	m_memoryFileMap.insert(std::make_pair(FileName, InternalFILE()));

	m_memoryFileMap[FileName].buffer = buffer;
	m_memoryFileMap[FileName].iBufferSize = defaultbufferSize;
	m_memoryFileMap[FileName].iFileSize = 0;
	m_memoryFileMap[FileName].iLocation = 0;
}

unsigned char *MemoryFileSystem::MemoryFileDrive::ReAllocBuffer(unsigned char *buffer, size_t size, size_t new_size)
{
	unsigned char *new_buffer = new unsigned char[new_size];

	for (size_t i = 0; i < size; ++i)
	{
		new_buffer[i] = buffer[i];
	}

	delete [] buffer;

	return new_buffer;
}

int MemoryFileSystem::MemoryFileDrive::CloseFile(FILE *fp)
{
	if (m_memoryFileMap.count(fp->GetName()) > 0)
	{
		// file is open for reading
		m_memoryFileMap[fp->GetName()].Reading = false;
		m_memoryFileMap[fp->GetName()].Writing = false;
		m_memoryFileMap[fp->GetName()].fp = NULL;

		delete fp;

		return 0;
	}
	else
	{
		return EOF;
	}
}

// ------------------------------------------------------------------

MemoryFileSystem::FILE *MemoryFileSystem::fopen(const char *filename, const char *mode)
{
	return GlobalDrive.CreateVirtualFile(filename, mode);
}

int MemoryFileSystem::fclose(FILE *fp)
{
	return GlobalDrive.CloseFile(fp);
}

size_t MemoryFileSystem::fwrite(const void * ptr, size_t size, size_t count, FILE *stream)
{
	unsigned char *buffer = GlobalDrive.m_memoryFileMap[stream->GetName()].buffer;
	size_t iBufferSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iBufferSize;
	size_t iLocation = GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation;
	size_t iFileSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize;

	size_t sizeTocopy = size * count;

	// is the buffer big enough
	if ((iLocation + sizeTocopy) > iBufferSize)
	{
		int extra_needed = ((iLocation + sizeTocopy) - iBufferSize) + MemoryFileDrive::defaultbufferSize;

		buffer = GlobalDrive.ReAllocBuffer(buffer, iBufferSize, iBufferSize + extra_needed);
		iBufferSize += extra_needed;
	}

	unsigned char *srcptr = (unsigned char *)ptr;

	// copy the data
	for (size_t i = 0; i < count; ++i)
	{
		memcpy(&buffer[iLocation], srcptr, size);
		iLocation += size;
		srcptr += size;
	}

	if (iLocation > iFileSize) iFileSize = iLocation;

	GlobalDrive.m_memoryFileMap[stream->GetName()].buffer = buffer;
	GlobalDrive.m_memoryFileMap[stream->GetName()].iBufferSize = iBufferSize;
	GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation = iLocation;
	GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize = iFileSize;

	return count;
}

size_t MemoryFileSystem::fread(void * ptr, size_t size, size_t count, FILE * stream)
{
	unsigned char *destptr = (unsigned char *)ptr;

	unsigned char *buffer = GlobalDrive.m_memoryFileMap[stream->GetName()].buffer;
	size_t iLocation = GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation;
	size_t iFileSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize;

	for (size_t i = 0; i < count; ++i)
	{
		if ((iLocation + size) > iFileSize) return i;

		memcpy(destptr, &buffer[iLocation], size);
		destptr += size;
		iLocation += size;
	}

	GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation = iLocation;

	return count;
}

int MemoryFileSystem::ftell(FILE *stream)
{
	return GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation;
}

int MemoryFileSystem::fseek(FILE * stream, long int offset, int origin)
{
	size_t iLocation = GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation;
	size_t iFileSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize;

	size_t new_location;

	if (origin == SEEK_END) new_location = iFileSize;
	if (origin == SEEK_SET) new_location = 0;
	new_location += offset;

	if (offset > iFileSize) return -1;

	GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation = new_location;

	return 0;
}

int MemoryFileSystem::fprintf(MemoryFileSystem::FILE * stream, const char * format, ...)
{
	va_list args;

	// get the amount of space needed
	va_start(args, format);
	int space_needed = GetSizeNeeded(format, args);
	va_end(args);

	// write to the file
	unsigned char *buffer = GlobalDrive.m_memoryFileMap[stream->GetName()].buffer;
	size_t iBufferSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iBufferSize;
	size_t iLocation = GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation;
	size_t iFileSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize;

	if (iLocation + space_needed > iBufferSize)
	{
		int extra_needed = ((iLocation + space_needed) - iBufferSize) + MemoryFileDrive::defaultbufferSize;

		buffer = GlobalDrive.ReAllocBuffer(buffer, iBufferSize, iBufferSize + extra_needed);
		iBufferSize += extra_needed;
	}

	va_start(args, format);

	char *startLocation = (char *)&buffer[iLocation];

	int count = vsnprintf(startLocation, iBufferSize - iLocation, format, args);
	iLocation += count;

	if (iFileSize < iLocation) iFileSize = iLocation;
	
	va_end(args);

	GlobalDrive.m_memoryFileMap[stream->GetName()].buffer = buffer;
	GlobalDrive.m_memoryFileMap[stream->GetName()].iBufferSize = iBufferSize;
	GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation = iLocation;
	GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize = iFileSize;

	return count;
}

char * MemoryFileSystem::fgets(char * str, int num, FILE * stream)
{
	unsigned char *buffer = GlobalDrive.m_memoryFileMap[stream->GetName()].buffer;
	size_t iBufferSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iBufferSize;
	size_t iLocation = GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation;
	size_t iFileSize = GlobalDrive.m_memoryFileMap[stream->GetName()].iFileSize;

	int characters_left = iFileSize - iLocation;

	if (characters_left <= 0) return NULL; // end of file

	char *src = (char *)&buffer[iLocation];
	char *dest = str;

	--num; // decrement to account for the '\0'
	int count = num < characters_left ? num : characters_left;
	int numread = 0;

	for (int i = 0; i < count; ++i)
	{
		int c = *src++;

		*dest++ = c;
		++numread;

		if (c == '\n') break;
	}

	str[numread] = '\0';

	iLocation += numread;

	GlobalDrive.m_memoryFileMap[stream->GetName()].iLocation = iLocation;

	return str;
}

char *MemoryFileSystem::_mktemp(char *_template){
	char buffer[32];

	sprintf(buffer, "%i", tempCount++);

	strcat(_template, buffer);

	return _template;
}

// this is temporary, hence the poor implementation
void MemoryFileSystem::WriteFileInMemoryToDisc(const char *filename)
{
	::FILE *fp = ::fopen(filename, "wb");

	if (GlobalDrive.HasFile(filename) && fp != NULL)
	{
		unsigned char *buffer = GlobalDrive.m_memoryFileMap[filename].buffer;
		size_t iBufferSize = GlobalDrive.m_memoryFileMap[filename].iBufferSize;
		size_t iLocation = GlobalDrive.m_memoryFileMap[filename].iLocation;
		size_t iFileSize = GlobalDrive.m_memoryFileMap[filename].iFileSize;

		unsigned char *c = buffer;

		::fwrite(c, 1, iFileSize, fp);
	}

	::fclose(fp);
}

// this is temporary, hence the poor implementation
void MemoryFileSystem::WriteFileInMemoryToStdout(const char *filename, bool includeheader, int padding)
{
#if defined(WIN32)
	setmode(fileno(stdout), O_BINARY);
#endif

	const int block_size = 128;

	if (GlobalDrive.HasFile(filename))
	{
		unsigned char *buffer = GlobalDrive.m_memoryFileMap[filename].buffer;
		size_t iBufferSize = GlobalDrive.m_memoryFileMap[filename].iBufferSize;
		size_t iLocation = GlobalDrive.m_memoryFileMap[filename].iLocation;
		size_t iFileSize = GlobalDrive.m_memoryFileMap[filename].iFileSize;

		unsigned char *c = buffer;

		if (includeheader)
		{
			// write a 32 char string with the length of the file
			char buffer[32];
			sprintf(buffer, "%u", iFileSize);
			::fwrite(buffer, 1, 32, stdout);
		}

		::fwrite(c, 1, iFileSize, stdout);

		// add extra padding if required
		if (includeheader && padding > 0)
		{
			for (int i = 0; i < padding; ++i)
			{
				putchar('\n');
			}
		}

		fflush(stdout);
	}
}


int MemoryFileSystem::GetSizeNeeded(const char * format, va_list args)
{
	char smallbuffer[32]; // use a buffer we know will be too small

	return vsnprintf(smallbuffer, 31, format, args);
}

// ----------------------------------------------------

