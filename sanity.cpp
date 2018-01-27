// http://www.manticmoo.com/articles/jeff/programming/cpp-stdin.php
// https://stackoverflow.com/a/17959055
#include <iostream>
#include <sstream>
#include <string>

int main() {
  // Load in all content into a traversable stream
  std::stringstream content;
  content << std::cin.rdbuf();


  // Output our content twice
  std::string line;
  while(content) {
    getline(content, line);
    std::cout << line << std::endl;
  };
  std::cout << "rewind" << std::endl;
  content.seekg(100);
  while(content) {
    getline(content, line);
    std::cout << line << std::endl;
  };

  // Exit our program
  return 0;
}