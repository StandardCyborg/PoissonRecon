#include <iostream>
#include <stdio.h>
#include <string>

int main(int argc, char* argv[]) {
  int foo;
  foo << std::cin;
  std::out << foo;
  return 0;
}