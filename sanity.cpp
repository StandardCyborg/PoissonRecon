#include <iostream>
#include <string>
using namespace std;

int main() {
  string input_line;
  while(cin) {
    getline(cin, input_line);
    cout << input_line << endl;
  };

  return 0;
}