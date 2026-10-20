#include <cstring>
#include <iostream>

int index_main(int argc, char* argv[]);
int search_main(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    if (argc == 2 && std::strcmp(argv[1], "--version") == 0) {
        std::cout << "skippy v1\n";
        return 0;
    }
    if (argc == 3 && std::strcmp(argv[2], "--version") == 0 &&
        (std::strcmp(argv[1], "index") == 0 || std::strcmp(argv[1], "search") == 0)) {
        std::cout << "skippy v1\n";
        return 0;
    }

    if (argc > 1) {
        if (std::strcmp(argv[1], "index") == 0)
            return index_main(argc - 1, argv + 1);
        if (std::strcmp(argv[1], "search") == 0)
            return search_main(argc - 1, argv + 1);
        if ((std::strcmp(argv[1], "--help") == 0 || std::strcmp(argv[1], "-help") == 0) &&
            argc == 2) {
            std::cout << "Usage: skippy <index|search> [options]\n"
                      << "  skippy index --help   Index a reference genome\n"
                      << "  skippy search --help  Search reads against an index\n"
                      << "  skippy --version      Show version\n";
            return 0;
        }
        std::cerr << "Error: unknown command or option: " << argv[1] << '\n';
    }
    else {
        std::cerr << "Error: expected index or search command\n";
    }

    std::cerr << "Usage: skippy <index|search> [options]\n"
              << "  skippy index --help   Index a reference genome\n"
              << "  skippy search --help  Search reads against an index\n"
              << "  skippy --version      Show version\n";
    return 2;
}
