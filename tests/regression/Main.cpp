#include "Emplode.hpp"

using namespace emplode;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: Emplode <script>" << std::endl;
        return 1;
    }

    Emplode emplode;
    emplode.Load(argv[1]);

    return 0;
}