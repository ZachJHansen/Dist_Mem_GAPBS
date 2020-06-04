#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char** argv) {
  if (argc < 3) {
    printf("Proper call: ./harness.txt shmem_output.txt normal_output.txt\n");
    exit(1);
  }
  string shmem_line;
  string norm_line;
  ifstream shmem; 
  ifstream normal;
  shmem.open(argv[1]);
  if (shmem.is_open()) {
    normal.open(argv[2]);
    if (normal.is_open()) {
      int i = 0;
      while (getline(shmem, shmem_line)) {
        if (getline(normal, norm_line)) {
          i++;
          if (shmem_line.compare(norm_line) != 0) {
            cout << "Mismatch! Shmem (" << shmem_line << ") != Normal (" << norm_line << ")" << " at line " << i << endl;
            exit(0);
          }
        } else {
          printf("Mismatch! Shmem file longer than Normal file\n");
          exit(0);
        }
      } 
      if (getline(normal, norm_line)) {
        printf("Mismatch! Normal file longer than Shmem file\n");
        exit(0);
      } else {
        printf("Files match\n");
      }
      normal.close();
    } else {
      printf("Couldnt open file %s\n", argv[2]);
    }
    normal.close();
  } else {
    printf("Couldnt open file %s\n", argv[1]);
  }
  return 0;
}
