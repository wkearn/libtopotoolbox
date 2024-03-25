extern "C" {
#include "topotoolbox.h"
}

#include <iostream>

int main(int argc, char *argv[]) {
  std::cout << "topotoolbox v" << TOPOTOOLBOX_VERSION_MAJOR << "."
            << TOPOTOOLBOX_VERSION_MINOR << "." << TOPOTOOLBOX_VERSION_PATCH
            << std::endl;
}
