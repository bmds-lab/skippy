//
//  serialiseKmersMap.hpp
//  indexReference
//
//  Created by Shlomo Geva on 19/7/2023.
//

#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "kmerUtilities.hpp"

struct QueryTable {
    std::map<uint32_t, std::string> entries;

    uint32_t query(uint32_t queryKey) const
    {
        auto it = entries.lower_bound(queryKey);
        if (it != entries.begin()) {
            // Check if the exact key was found
            if (it == entries.end() || it->first != queryKey) {
                --it; // Move the iterator one step back to the highest entry not
                      // exceeding the queryKey
            }
            return it->first; // Return the key
        }
        else {
            if (!entries.empty()) {
                // If the queryKey is smaller than the first entry, return the key of
                // the first entry
                return entries.begin()->first;
            }
        }
        return static_cast<uint32_t>(-1); // Return -1 to indicate key not found for an empty table
    }

    // Function to print out the keys of the 'entries' map
    void printKeys() const
    {
        std::cout << "Table entries" << std::endl;
        int n = 0;
        for (const auto& entry : entries) {
            std::cout << n++ << ":  " << entry.first << " " << entry.second << std::endl;
        }
    }

    // Function to access the elements of the underlying map
    const std::string& getValue(uint32_t key) const
    {
        // Use the at() method of the map to access the value
        // It throws an exception if the key is not found
        return entries.at(key);
    }
};

// struct QueryTable {
//     std::map<uint32_t, std::string> entries;
//
//     uint32_t query(uint32_t queryKey) const {
//         auto it = entries.lower_bound(queryKey);
//         if (it != entries.begin()) {
//             --it; // Move the iterator one step back to the highest entry not
//             exceeding the queryKey return it->first; // Return the key
//         } else {
//             if (!entries.empty()) {
//                 // If the queryKey is smaller than the first entry, return
//                 the key of the first entry return entries.begin()->first;
//             }
//         }
//         return -1; // Return -1 to indicate key not found for an empty table
//     }
//
//     // Function to print out the keys of the 'entries' map
//     void printKeys() const {
//         std::cout << "Table entries" << std::endl;
//         int n = 0;
//         for (const auto& entry : entries) {
//             std::cout << n++ << ":  " << entry.first << " "
//                       << entry.second << std::endl;
//         }
//     }
//
//     // Function to access the elements of the underlying map
//     const std::string& getValue(uint32_t key) const {
//         // Use the at() method of the map to access the value
//         // It throws an exception if the key is not found
//         return entries.at(key);
//     }
//
// };

QueryTable readTableFromFile(const std::string& filename);

void serializeMap(const std::vector<std::vector<uint32_t>>& kmersMap,
                  const std::string& innerMapFilename, const std::string& outerMapFilename);

void deserializeMap(const std::string& innerMapFilename, const std::string& outerMapFilename,
                    std::span<uint32_t>& innerMapBlob, std::span<uint32_t>& outerMapBlob);

double getInnerVector(std::vector<store_position>& result,
                      const std::vector<uint32_t>& innerMapBlob,
                      const std::vector<uint32_t>& outerMapBlob, uint32_t index, uint32_t maxLen);

bool writeTextBlobToFile(const char* text, std::size_t length, const std::string& filename);

std::span<char> readTextBlobFromFile(const std::string& filename);
