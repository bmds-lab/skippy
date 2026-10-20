#include <cstring>
#include <iostream>

int index_main(int argc, char* argv[]);
int search_main(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    if (argc > 1) {
        if (std::strcmp(argv[1], "index") == 0)
            return index_main(argc - 1, argv + 1);
        if (std::strcmp(argv[1], "search") == 0)
            return search_main(argc - 1, argv + 1);
        if (std::strcmp(argv[1], "--help") != 0 && std::strcmp(argv[1], "-help") != 0)
            std::cerr << "Unknown command: " << argv[1] << '\n';
    }

    std::cerr << "Usage: skippy <index|search> [options]\n"
              << "  skippy index --help   Index a reference genome\n"
              << "  skippy search --help  Search reads against an index\n";
    return argc > 1 && (std::strcmp(argv[1], "--help") == 0 || std::strcmp(argv[1], "-help") == 0)
               ? 0
               : 1;
}
