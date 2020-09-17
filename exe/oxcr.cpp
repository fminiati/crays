//
// integrator for cr-transport in momentum space
// F. Miniati, (ZH) 08.12.2018
//

#include <cstdlib>
#include <cstring>
#include "FileParser.H"


#include "CRTransport.H"

//using namespace cr_transport;

// main code
int main(int argc, char* argv[]) {

  // get input file
  string input_file="unknown";

  for (int i=1; i<argc; ++i) {
    if (strncmp(argv[i],"-f",3) == 0) input_file = argv[i+1];
  }
  std::cout << " input file name is " << input_file << std::endl;

  // read input file
  FileParser input(input_file);

    // cr transport object
  cr_transport::CRTransport crt;

  std::cout << " setting up cr_transport object " << std::endl;
  crt.setup(input_file);

  std::cout << " solving cr-transport equation " << std::endl;
  crt.solve();

  std::cout << " done " << std::endl;
}
