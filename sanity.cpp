#include <iostream>
#include <string>

int main() {
  // https://stackoverflow.com/a/17959055
  std::string input_line;
  while(std::cin) {
    getline(std::cin, input_line);
    std::cout << input_line << std::endl;
  };

  return 0;
}