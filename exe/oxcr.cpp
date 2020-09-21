//
// Copyright (C) 2020 Francesco Miniati <francesco.miniati@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include <cstring>
#include <string>
#include "FileParser.h"
#include "CRTransport.h"

// main code
int main(int argc, char* argv[]) {

  // get input file
  std::string input_file="unknown";

  for (int i=1; i<argc; ++i)
  {
      if (strncmp(argv[i], "-f", 3) == 0)
          input_file = argv[i + 1];
  }
  std::cout << " input file name is " << input_file << std::endl;

  // read input file
  fm::FileParser input(input_file);

    // cr transport object
  fm::cr_transport::CRTransport crt;

  std::cout << " setting up cr_transport object " << std::endl;
  crt.setup(input_file);

  std::cout << " solving cr-transport equation " << std::endl;
  crt.solve();

  std::cout << " done " << std::endl;
}
