// http://www.manticmoo.com/articles/jeff/programming/cpp-stdin.php
// Never mind... stringstream isn't what we want; it's only for double READING IN stdin
// https://stackoverflow.com/a/17959055
// http://cpp.indi.frih.net/blog/2014/09/how-to-read-an-entire-file-into-memory-in-cpp/
#include <iostream>
#include <sstream>
#include <string>

int main() {
  // Load in all content into a traversable stream
  std::stringstream content;
  content << std::cin.rdbuf();


  // Output our content twice
  std::string line;
  content.close();
  std::cout << "rewind" << content.tellp() << std::endl;
  while(content) {
    getline(content, line);
    std::cout << line << std::endl;
  };
  std::cout << "rewind" << content.tellp() << std::endl;
  content.seekp(166);
  std::cout << "rewind" << content.tellp() << std::endl;
  while(content) {
    getline(content, line);
    std::cout << line << std::endl;
  };

  // Exit our program
  return 0;
}