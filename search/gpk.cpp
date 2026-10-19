/*  GPK - invoke KISS on multiple directories
 *  takes input directory (must be absolute path because using symbolic links
 *  (see code below) splits it into temporary directories - one per process -
 *  and then runs it as a separate pricess on the machine.
 *  to compile:  g++ -std=c++17 -O2 -march=native -Wall -Wextra -pedantic gpk.cpp -o gpk
 */

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

struct FileInfo
{
    std::string filename;
    uintmax_t size;
};

// ---------------------------------------------------------------
// MAKE_ABSOLUTE  working with absolute paths in this program
// ---------------------------------------------------------------
static std::string make_absolute(const std::string &path)
{
    if (path.empty())
        return path;

    fs::path p = path;
    if (p.is_absolute())
    {
        return path;
    }
    else
    {
        return (fs::current_path() / p).string();
    }
}

// Function to sort files by size in descending order
bool compareFileSize(const FileInfo &a, const FileInfo &b)
{
    return a.size > b.size;
}

void createSymbolicLinkDirsBySize(const std::string &dataDir, int np,
                                  std::vector<std::string> &tempDirs)
{
    std::vector<FileInfo> files;

    // Get the list of files and their sizes
    for (const auto &entry : fs::directory_iterator(dataDir))
    {
        if (fs::is_regular_file(entry))
        {
            files.push_back({entry.path().string(), fs::file_size(entry.path())});
        }
    }

    // Sort files by size in descending order
    std::sort(files.begin(), files.end(), compareFileSize);

    // Initialize partitions
    std::vector<std::vector<std::string>> partitions(np);
    std::vector<uintmax_t> partitionSizes(np, 0);

    // Greedily distribute files based on size
    for (const auto &file : files)
    {
        // Find the partition with the smallest total size
        auto minElement = std::min_element(partitionSizes.begin(), partitionSizes.end());
        int partitionIndex = std::distance(partitionSizes.begin(), minElement);

        // Add file to the selected partition
        partitions[partitionIndex].push_back(file.filename);
        partitionSizes[partitionIndex] += file.size;
    }

    // Create symbolic links for the files in the temp directories
    for (int i = 0; i < np; ++i)
    {
        fs::create_directory(tempDirs[i]); // Create temp directory
        for (const auto &file : partitions[i])
        {
            std::string linkPath = tempDirs[i] + "/" + fs::path(file).filename().string();
            fs::create_symlink(file, linkPath); // Create symlink in the temp directory
        }
    }

    // Print out the partition sizes for debugging
    for (int i = 0; i < np; ++i)
    {
        std::cout << "Partition " << i + 1 << " total size: " << partitionSizes[i] << " bytes"
                  << std::endl;
    }
}

void invokeKISS(const std::string &tempDir, const std::vector<std::string> &kissArgs)
{
    std::vector<char *> args;
    // KISS program path
    const char *KISS = "./filter";

    // Add the KISS executable as the first argument
    args.push_back(const_cast<char *>("KISS"));

    // Add the command-line arguments passed to GPK (except the data directory)
    for (const auto &arg : kissArgs)
    {
        args.push_back(const_cast<char *>(arg.c_str()));
    }

    // Replace the -data argument with the temp directory
    args.push_back(const_cast<char *>("-data"));
    args.push_back(const_cast<char *>(tempDir.c_str()));

    // Add the null pointer to mark the end of the arguments
    args.push_back(nullptr);

    // Fork the process and execute KISS
    pid_t pid = fork();
    if (pid == 0)
    {
        // In the child process, invoke KISS
        execvp(KISS, args.data());
        // If execvp returns, it means there was an error
        std::cerr << "Failed to execute KISS in " << tempDir << std::endl;
        exit(EXIT_FAILURE);
    }
    else if (pid < 0)
    {
        std::cerr << "Failed to fork for tempDir: " << tempDir << std::endl;
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: GPK -np <number_of_processes> -data <dataDir> [KISS_args...]"
                  << std::endl;
        std::cerr << "Example: time ./gpk -np 16 -data $(pwd)/data/human_10_err_0_to_10_fasta "
                     "-reference GRCh38.fasta -min_matches 2 -cut_off 50 -threads 96 -num_seeds 16 "
                     "-out_dir junkGPK -np <number_of_processes> -data <dataDir> [KISS_args...]"
                  << std::endl;
        std::cerr << "Note: using $(pwd)/ to make it a full path - must not be a relative path"
                  << std::endl;

        return EXIT_FAILURE;
    }
    /*
        // Parse the command line arguments
        int np = 0;
        std::string dataDir;
        std::vector<std::string> kissArgs;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-np" && i + 1 < argc) {
                np = std::stoi(argv[++i]);
            } else if (arg == "-data" && i + 1 < argc) {
                dataDir = argv[++i];
            } else {
                kissArgs.push_back(arg);
            }
        }
    */
    // Parse the command line arguments
    int np = 0;
    std::string dataDir;
    std::string referenceFile; // optional: also fix reference path
    std::vector<std::string> kissArgs;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-np" && i + 1 < argc)
        {
            np = std::stoi(argv[++i]);
        }
        else if (arg == "-data" && i + 1 < argc)
        {
            dataDir = make_absolute(argv[++i]); // <<< FIXED
        }
        else if (arg == "-reference" && i + 1 < argc)
        {
            referenceFile = make_absolute(argv[++i]); // <<< also fixed if you want
        }
        else
        {
            kissArgs.push_back(arg);
        }
    }

    if (np <= 0 || dataDir.empty())
    {
        std::cerr << "Error: Invalid number of processes or data directory." << std::endl;
        return EXIT_FAILURE;
    }

    // Create temp directories
    std::vector<std::string> tempDirs(np);
    for (int i = 0; i < np; ++i)
    {
        tempDirs[i] = "./tempDir" + std::to_string(i + 1);
    }

    // Distribute files to temp directories based on file sizes
    createSymbolicLinkDirsBySize(dataDir, np, tempDirs);

    // Invoke KISS for each temp directory in parallel
    std::vector<pid_t> pids;
    for (const auto &tempDir : tempDirs)
    {
        invokeKISS(tempDir, kissArgs);
    }

    // Wait for all KISS processes to complete
    int status;
    while (wait(&status) > 0)
        ;

    // Cleanup: Optionally, delete the temp directories
    for (const auto &tempDir : tempDirs)
    {
        fs::remove_all(tempDir); // Remove the temp directory
    }

    return 0;
}
