#include <fstream>
#include <iostream>
#include "crypto.h"
#include "util.h"

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    cout << "Usage: GenerateKey out.key" << endl;
    return 1; 
  }

  cout << "Generating a random 128-bit key for AES/GCM ... ";
  unsigned char key[PRF_KEY_BYTES];
  int nbytes = randread(key, PRF_KEY_BYTES);
  if (nbytes < PRF_KEY_BYTES) {
    cout << endl << "Error: Failed to sample from /dev/urandom" << endl;
    return 1;
  }

  ofstream ofs(argv[1], ios::binary);
  if (!ofs.is_open()) {
    cout << endl << "Error: Failure to open file" << endl;
    return 1;
  }

  ofs.write((const char *)&key, PRF_KEY_BYTES);
  ofs.close();

  cout << "done." << endl;

  return 0;
}
