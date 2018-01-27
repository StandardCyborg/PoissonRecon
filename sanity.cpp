// http://www.manticmoo.com/articles/jeff/programming/cpp-stdin.php
// https://stackoverflow.com/a/17959055
#include <iostream>
#include <sstream>
#include <string>

int main() {
  std::stringstream content;
  content << std::cin.rdbuf();
  std::string input_line;

  while(content) {
    getline(content, input_line);
    std::cout << input_line << std::endl;
  };
  content.seekg(0);
  while(content) {
    getline(content, input_line);
    std::cout << input_line << std::endl;
  };

  return 0;
}