// http://www.manticmoo.com/articles/jeff/programming/cpp-stdin.php
// Never mind... stringstream isn't what we want; it's only for double READING IN stdin
// https://stackoverflow.com/a/17959055
// http://cpp.indi.frih.net/blog/2014/09/how-to-read-an-entire-file-into-memory-in-cpp/
// http://pubs.opengroup.org/onlinepubs/9699919799/functions/fmemopen.html
// http://www.cplusplus.com/reference/cstdio/fread/
// http://www.cplusplus.com/reference/cstdio/fwrite/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>

static char buffer[] = "foobar";

int main() {
  // Load in all content into a traversable stream
  // std::stringstream content;
  // content << std::cin.rdbuf();

  FILE *foo = fmemopen(buffer, strlen(buffer), "w+b");
  if (foo == NULL) {
    std::cerr << "Unable to open memory stream";
    exit(1);
  }

  char buffer2[] = { 'x' , 'y' , 'z' };
  fwrite(buffer2, sizeof(char), sizeof(buffer2), foo);

  // Allocate buffer for ourselfs
  // TODO: Can we do sizeof buffer
  char *buffer = (char*) malloc (sizeof(char)*strlen(buffer));
  if (buffer == NULL) {
    std::cerr << "Memory error";
    exit(2);
  }

  fread(buffer, 1, 5, foo);

  std::cout << buffer;
  // // Output our content twice
  // std::string line;
  // std::cout << "rewind" << content.tellp() << std::endl;
  // while(content) {
  //   getline(content, line);
  //   std::cout << line << std::endl;
  // };
  // std::cout << "rewind" << content.tellp() << std::endl;
  // content.seekp(166);
  // std::cout << "rewind" << content.tellp() << std::endl;
  // while(content) {
  //   getline(content, line);
  //   std::cout << line << std::endl;
  // };

  // Exit our program
  return 0;
}