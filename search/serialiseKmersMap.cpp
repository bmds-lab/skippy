//
//  serialiseKmersMap.cpp
//  indexReference
//
//  Created by Shlomo Geva on 19/7/2023.
//
#include "serialiseKmersMap.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>

#include "kmerUtilities.hpp"

void serializeMap(const std::vector<std::vector<uint32_t>>& kmersMap,
                  const std::string& innerMapFilename, const std::string& outerMapFilename)
{
    // Set the buffer size to 8192 bytes (for example)
    constexpr std::streamsize bufferSize = 8192;

    // Open the files with custom buffer sizes
    std::ofstream innerMapFile(innerMapFilename, std::ios::binary);
    innerMapFile.rdbuf()->pubsetbuf(nullptr, bufferSize);

    std::ofstream outerMapFile(outerMapFilename, std::ios::binary);
    outerMapFile.rdbuf()->pubsetbuf(nullptr, bufferSize);

    uint32_t offset = 0;
    for (const auto& innerVector : kmersMap) {
        // Write inner vector data
        innerMapFile.write(reinterpret_cast<const char*>(innerVector.data()),
                           innerVector.size() * sizeof(uint32_t));

        // Write the offset for the outer map
        outerMapFile.write(reinterpret_cast<const char*>(&offset), sizeof(uint32_t));

        // Calculate the offset for the next inner vector
        offset += static_cast<uint32_t>(innerVector.size());
    }

    innerMapFile.close();
    outerMapFile.close();
}

void deserializeMap(const std::string& innerMapFilename, const std::string& outerMapFilename,
                    std::span<uint32_t>& innerMapBlob, std::span<uint32_t>& outerMapBlob)
{
    int fd;
    struct stat sb;
    void* addr;
    size_t length;
    fd = open(innerMapFilename.c_str(), O_RDONLY);
    fstat(fd, &sb);
    length = sb.st_size;
    addr = mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0);
    innerMapBlob = std::span<uint32_t>((uint32_t*)addr, length / sizeof(uint32_t));
    close(fd);
    fd = open(outerMapFilename.c_str(), O_RDONLY);
    fstat(fd, &sb);
    length = sb.st_size;
    addr = mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0);
    outerMapBlob = std::span<uint32_t>((uint32_t*)addr, length / sizeof(uint32_t));
    close(fd);
}

// Helper function to access the index
double getInnerVector(std::vector<store_position>& result,
                      const std::vector<uint32_t>& innerMapBlob,
                      const std::vector<uint32_t>& outerMapBlob, uint32_t index, uint32_t maxLen)
{
    if (index < outerMapBlob.size()) {
        uint32_t startOffset = outerMapBlob[index];
        uint32_t endOffset = (index + 1 < outerMapBlob.size()) ? outerMapBlob[index + 1]
                                                               : (uint32_t)innerMapBlob.size();
        double range = endOffset - startOffset;
        if (range > maxLen) {
            return -range; // ignore crazy long inverted lists
        }
        if (endOffset > startOffset) {
            result.push_back(
                store_position{innerMapBlob.data() + startOffset, endOffset - startOffset});
            return range + 1;
        }
        else
            return 0; // not found
    }
    return 0; // beyond size of outer blob
}

bool writeTextBlobToFile(const char* text, std::size_t length, const std::string& filename)
{
    std::ofstream outputFile(filename, std::ios::binary);
    if (!outputFile.is_open()) {
        return false; // Failed to open the file
    }

    outputFile.write(text, length);

    outputFile.close();

    return true; // File write successful
}

std::span<char> readTextBlobFromFile(const std::string& filename)
{
    int fd;
    struct stat sb;
    void* addr;
    size_t length;

    fd = open(filename.c_str(), O_RDONLY);
    fstat(fd, &sb);
    length = sb.st_size;
    addr = mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0);
    std::span<char> textBlob((char*)addr, length);
    close(fd);
    return textBlob;
    /*
        std::ifstream inputFile(filename, std::ios::binary | std::ios::ate);
        if (!inputFile.is_open()) {
            return std::span<char>(); // Failed to open the file
        }

        std::size_t length = inputFile.tellg(); // Get the file size
        inputFile.seekg(0); // Move the file pointer back to the beginning

        char *buffer = (char *)malloc(length); // Create a buffer to hold the
       content

        if (!inputFile.read(buffer, length)) {
            return std::span<char>(); // Failed to read the content
        }

        inputFile.close();

        //buffer[length]='\0'; // terminate with '\0' to make vaild string
        // Allocate memory for the text blob and copy the content from the buffer
        std::span<char> textBlob(buffer, length);
        return textBlob; // Return the textBlob string
    */
}

QueryTable readTableFromFile(const std::string& filename)
{
    QueryTable table;
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return table;
    }
    //    std::cout << "Table from " << filename << std::endl;
    uint32_t key;
    std::string value;
    std::string line; // Used to read the entire line
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        if (iss >> key) {
            value = line.substr(line.find(' ') + 1);
            //    std::cout << key << " : " << value << std::endl;
            table.entries[key] = value;
        }
    }

    inFile.close();

    return table;
}

// int mainTest() {
//     std::map<uint32_t, std::string> referenceIDMap = {
//         {100, "Value100"},
//         {200, "Value200"},
//         {300, "Value300"},
//         // Add more entries as needed
//     };
//
//     std::string filename = "output.txt";
//     writeMapToFile(filename, referenceIDMap);
//
//     // Read the file back into the query table
//     QueryTable queryTable = readTableFromFile(filename);
//
//     // Perform a lookup and get the key
//     uint32_t queryKey = 250;
//     int key = queryTable.query(queryKey);
//
//     if (key != -1) {
//         // Access the value using the key
//         std::string value = queryTable.entries[key];
//
//         // Display the key and value
//         std::cout << "Key: " << key << ", Value: " << value << std::endl;
//     } else {
//         std::cout << "Key not found." << std::endl;
//     }
//
//     return 0;
// }
