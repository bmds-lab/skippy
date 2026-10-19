//
//  smithWaterman.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 15/7/2023.
//

#include "main.hpp"
#include "smithWaterman.hpp"
#include "ssw.h"
#include "ssw_cpp.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

extern long long DEBUG_THIS;
/*
     COMPUTE SCORE AND LENGTH (from CIGAR)
 */

extern uint32_t KMERSIZE;
// extern uint32_t MIN_MATCHES;   // minimum number of seed matches required
// extern uint32_t NUM_SEEDS; // seed skipping distancee
extern bool MATCH_ALL; // find all mathces flag
extern bool GET_SW;    // compute strict SW (defualt is naive sw)

#include <algorithm>
#include <cctype>
#include <span>
#include <stdexcept>
#include <string_view>

// Define DEBUG to enable logging, undefine to disable
#define DEBUG

int align(const std::span<char> s, const std::string_view t)
{
#ifdef DEBUG2
    std::cout << "align: Starting with s=\n"
              << std::string_view(s.data(), s.size()) << "\nt=\n"
              << t << ", s.size=" << s.size() << ", t.size=" << t.size() << std::endl;
#endif

    size_t n = s.size();
    int score = 0;
    size_t i = 0;
    size_t j = 0;
    const size_t lookahead =
        10; // 10 is not expensive because we do not expect such long mismatches

    while (i < n && j < n)
    {
#ifdef DEBUG2
        std::cout << "align: i=" << i << ", j=" << j << ", score=" << score
                  << ", s[i]=" << (i < n ? std::string(1, s[i]) : "END")
                  << ", t[j]=" << (j < n ? std::string(1, t[j]) : "END") << std::endl;
#endif

        char s_char = s[i];
        char t_char = t[j];
        bool is_match =
            (s_char == t_char) || (std::toupper(s_char) == std::toupper(t_char) &&
                                   (std::toupper(s_char) == 'A' || std::toupper(s_char) == 'C' ||
                                    std::toupper(s_char) == 'G' || std::toupper(s_char) == 'T'));

        if (is_match)
        {
#ifdef DEBUG2
            std::cout << "align: Match found (case-insensitive or exact), incrementing score"
                      << std::endl;
#endif
            score++;
            i++;
            j++;
        }
        else
        {
#ifdef DEBUG2
            std::cout << "align: No direct match, entering lookahead" << std::endl;
#endif
            bool found_match = false;
            size_t best_di = 0;
            size_t best_dj = 0;
            size_t min_skip = n;

            for (size_t di = 0; di <= lookahead && i + di < n; di++)
            {
                for (size_t dj = 0; dj <= lookahead && j + dj < n; dj++)
                {
                    char s_look = s[i + di];
                    char t_look = t[j + dj];
                    if (di + dj > 0 &&
                        ((s_look == t_look) ||
                         (std::toupper(s_look) == std::toupper(t_look) &&
                          (std::toupper(s_look) == 'A' || std::toupper(s_look) == 'C' ||
                           std::toupper(s_look) == 'G' || std::toupper(s_look) == 'T'))))
                    {
                        if (di + dj < min_skip)
                        {
                            min_skip = di + dj;
                            best_di = di;
                            best_dj = dj;
                            found_match = true;
#ifdef DEBUG2
                            std::cout << "align: Lookahead match found at s[" << i + di
                                      << "]=" << s_look << ", t[" << j + dj << "]=" << t_look
                                      << ", di=" << di << ", dj=" << dj << std::endl;
#endif
                        }
                    }
                }
            }

            if (found_match)
            {
#ifdef DEBUG2
                std::cout << "align: Advancing to i=" << i + best_di + 1
                          << ", j=" << j + best_dj + 1 << std::endl;
#endif
                i += best_di + 1;
                j += best_dj + 1;
                score++;
            }
            else
            {
#ifdef DEBUG2
                std::cout << "align: No match in lookahead, advancing both pointers" << std::endl;
#endif
                i++;
                j++;
            }
        }
    }

#ifdef DEBUG2
    std::cout << "align: Exiting with score=" << score << ", final i=" << i << ", final j=" << j
              << std::endl;
#endif
    return score;
}

std::pair<int, int> computeScoreAndLength(const std::string &cigarString)
{
    /*
     The SAM spec offers us this table of CIGAR operations which indicates which ones "consume" the
     query or the reference, complete with explicit instructions on how to calculate sequence length
     from a CIGAR string: Consumes  Consumes Op  BAM Description query  reference M   0   alignment
     match (can be a sequence match or mismatch)   yes   yes I   1   insertion to the reference yes
     no D   2   deletion from the reference                             no    yes N   3   skipped
     region from the reference                       no    yes S   4   soft clipping (clipped
     sequences present in SEQ)        yes   no H   5   hard clipping (clipped sequences NOT present
     in SEQ)    no    no P   6   padding (silent deletion from padded reference)         no    no =
     7   sequence match                                          yes   yes X   8   sequence mismatch
     yes   yes “Consumes query” and “consumes reference” indicate whether the CIGAR operation causes
     the alignment to step along the query sequence and the reference sequence respectively.
     ...
     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
     */
    int score = 0;
    int length = 0;
    std::stringstream ss(cigarString);
    char c;
    int currentNumber = 0;

    while (ss >> currentNumber >> c)
    {
        // this loops around the string stream, extracting a number and a symbol from the CIGAR
        // this could be a simple if statement, but using a switch to allow
        //     different penalty/reward for alignment operations in score calculations
        switch (c)
        {
        case '=':
            score += currentNumber;
            length += currentNumber; // Add length for the matches to the total length
            break;
        case 'M':
        case 'I':
        case 'S':
        case 'X':
        case 'D':
        case 'N':
        case 'H':
        case 'P':
            score -= currentNumber; // penalise all variations due to misalignment
        default:
            break;
        }
    }
    return {score, length};
}

// Function to print alignment in BLAST-like format
void print_blast_like_output(const StripedSmithWaterman::Alignment &alignment,
                             const std::string &ref, const std::string &query)
{
    int ref_start = alignment.ref_begin;
    int ref_end = alignment.ref_end;
    int query_start = alignment.query_begin;
    int query_end = alignment.query_end;
    std::string cigar_string = alignment.cigar_string;

    std::cout << "*********\nRef: (length= " << ref.size() << ")\n" << ref << "\n";
    std::cout << "Query: (length= " << query.size() << ")\n" << query << "\n";

    std::cout << "Alignment:\n";
    std::cout << "Query start: " << query_start << "\n";
    std::cout << "Query end: " << query_end << "\n";
    std::cout << "Reference start: " << ref_start << "\n";
    std::cout << "Reference end: " << ref_end << "\n";
    std::cout << "Score: " << alignment.sw_score << "\n\n";

    // Parse CIGAR string to get the aligned sequences and match line
    std::string ref_aligned, query_aligned, match_line;
    int ref_idx = ref_start;
    int query_idx = query_start;
    int cigar_len = (int)cigar_string.length();
    for (int i = 0; i < cigar_len; ++i)
    {
        char op = cigar_string[i];
        if (isdigit(op))
        {
            int length = op - '0';
            while (i + 1 < cigar_len && isdigit(cigar_string[i + 1]))
            {
                length = length * 10 + (cigar_string[++i] - '0');
            }
            for (int j = 0; j < length; ++j)
            {
                if (cigar_string[i + 1] == '=')
                { // Match/Mismatch
                    ref_aligned += ref[ref_idx++];
                    query_aligned += query[query_idx++];
                    if (ref[ref_idx - 1] == query[query_idx - 1])
                    {
                        match_line += '|';
                    }
                    else
                    {
                        match_line += ' ';
                    }
                }
                else if (cigar_string[i + 1] == 'I')
                { // Insertion to the reference
                    ref_aligned += '-';
                    query_aligned += query[query_idx++];
                    match_line += ' ';
                }
                else if (cigar_string[i + 1] == 'D')
                { // Deletion from the reference
                    ref_aligned += ref[ref_idx++];
                    query_aligned += '-';
                    match_line += ' ';
                }
            }
            ++i; // Move past the operation character
        }
    }

    std::cout << "Ref Aligned:  " << ref_aligned << "\n";
    std::cout << "Match Line:   " << match_line << "\n";
    std::cout << "Query Aligned:" << query_aligned << "\n\n";
}

// Function to convert a string to uppercase
std::string toUpperCase(const std::string &str)
{
    std::string upperStr = str;
    std::transform(upperStr.begin(), upperStr.end(), upperStr.begin(), ::toupper);
    return upperStr;
}
std::string toUpperCase(const std::span<char> str)
{
    // Create a string from the span
    std::string upperStr(str.begin(), str.end());
    // Transform the string to uppercase
    std::transform(upperStr.begin(), upperStr.end(), upperStr.begin(), ::toupper);
    return upperStr;
}

void printAlignedSequences(const std::span<char> ref, const std::string query, std::string &cigar,
                           int &sw_score, uint32_t &refBegin, uint32_t &refEnd)
{
    int32_t maskLen = (int32_t)query.size() / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the query to the ref
    aligner.Align(query.c_str(), ref.data(), (int)ref.size(), filter, &alignment, maskLen);

    sw_score = alignment.sw_score;
    cigar = alignment.cigar_string;
    refBegin = alignment.ref_begin;
    refEnd = alignment.ref_end;

    // Print the BLAST-like output (for debugging)
    std::string upperRef(ref.begin(), ref.end());
    upperRef = toUpperCase(upperRef);
    print_blast_like_output(alignment, upperRef, query);
}

void alignSequences(const std::span<char> ref, const std::string &query, std::string &cigar,
                    int &sw_score, uint32_t &refBegin, uint32_t &refEnd)
{
    int32_t maskLen = (int32_t)query.size() / 2;
    //    	int32_t maskLen = (int32_t) strlen(query.c_str()) / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the query to the ref
    aligner.Align(query.c_str(), ref.data(), (int)ref.size(), filter, &alignment, maskLen);

    sw_score = alignment.sw_score;
    cigar = alignment.cigar_string;
    refBegin = alignment.ref_begin;
    refEnd = alignment.ref_end;
}

/*
     SSW Aligner
 */
/*
void alignSequences(const std::string ref, const std::string query,std::string& cigar, int&
sw_score, uint32_t& refBegin, uint32_t& refEnd) { int32_t maskLen = (int32_t)
strlen(query.c_str())/2; maskLen = maskLen < 15 ? 15 : maskLen;

    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the query to the ref
    aligner.Align(query.c_str(), ref.c_str(), (int) ref.size(), filter, &alignment, maskLen);

    sw_score = alignment.sw_score;
    cigar = alignment.cigar_string;
    refBegin = alignment.ref_begin;
    refEnd = alignment.ref_end;
}
*/

/*
    REVERSE COMPLEMENT (string)
 */
std::string reverseComplement(const std::string &sequence)
{
    std::string reverseSeq;
    for (char nucleotide : sequence)
    {
        switch (nucleotide)
        {
        case 'A':
            reverseSeq += 'T';
            break;
        case 'T':
            reverseSeq += 'A';
            break;
        case 'C':
            reverseSeq += 'G';
            break;
        case 'G':
            reverseSeq += 'C';
            break;
        case 'a':
            reverseSeq += 'T';
            break;
        case 't':
            reverseSeq += 'A';
            break;
        case 'c':
            reverseSeq += 'G';
            break;
        case 'g':
            reverseSeq += 'C';
            break;
        default:
            reverseSeq += nucleotide;
        }
    }
    std::reverse(reverseSeq.begin(), reverseSeq.end());
    return reverseSeq;
}

// Inline function for efficiency in tight loops
static inline int containsN(const std::string &sequence)
{
    //    if (sequence.size() != 32) return false;

    for (int i = 0; i < sequence.size(); i++)
    {
        if (sequence[i] == 'N')
            return i; // index of 'N'
    }
    return -1; // no 'N' found
}

bool hasLongMatch(const std::string &cigar)
{
    size_t i = 0;
    while (i < cigar.size())
    {
        int count = 0;

        // Parse the number
        while (i < cigar.size() && std::isdigit(cigar[i]))
        {
            count = count * 10 + (cigar[i] - '0');
            ++i;
        }

        // Check if it's a '=' and count >= 16
        if (i < cigar.size() && cigar[i] == '=' && count >= 16)
        {
            return true;
        }

        ++i; // Skip the operation character
    }
    return false;
}

/*
     GETSW
 */

// Function to compute the Smith-Waterman scores between pairs in the results
// In this code, the resultsSW vector is defined as std::vector<std::tuple<uint32, uint32, uint32>>
// to store the index pair and the Smith-Waterman score. The assignment of sampleIndex
// and referenceIndex is done directly from the result tuple, and the Smith-Waterman score
// is computed using the sample and reference sequences. Finally, the index pair
// and the SW score are stored in the resultsSW vector using std::make_tuple.

void getSW(std::span<const std::string> chunk, size_t startIndex, const std::span<char> &genome,
           SafeQueue<Result> &queue,
           const std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> &candidate_results,
           int matchScore, int mismatchScore, int gapScore, int cutoff)
{
    //    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string, int>>
    //    localResultsSW;
    struct Result result;
    float score;
    bool skipMulti = false; // used to skip multiple matches per read when MATCH_ALL is false
    size_t skipCount = 0;
    uint32_t currentSample = UINT_MAX; // sentinel
    // Iterate through each result in the results vector
    std::string cigar, cigar2;
    uint32_t sequenceLength;
    size_t SAMPLE_INDEX; // for index of actual read (not index in the chunk)
                         /****************************************/
    uint32_t debug = UINT32_MAX;
    //    debug=376193;
    /****************************************/
    //    for (const auto& result : candidate_results) {
    for (int n = 0; n < candidate_results.size(); n++)
    {
        uint32_t sampleIndex, revcomp, minPos, maxPos, nSeed;
        std::tie(sampleIndex, minPos, nSeed, revcomp) = candidate_results[n];
        if (currentSample == sampleIndex)
        {
            if (skipMulti)
            {
                ++skipCount;
                continue;
            }
        }
        else
        {
            currentSample = sampleIndex;
            skipMulti = false;
        }
        std::string sampleSequence = chunk[sampleIndex];
        SAMPLE_INDEX = sampleIndex + startIndex; // correct to original read index
        if (SAMPLE_INDEX >= debug - 10 && SAMPLE_INDEX < debug + 9)
            std::cout << "sampleIndex " << SAMPLE_INDEX << "\n";
        if (debug == SAMPLE_INDEX)
            std::cout << "debugging " << debug << "\n";

        sequenceLength = (uint32_t)sampleSequence.size();
        //        std::string exact = std::to_string(sampleSequence.size()) + "=";

        uint32_t minPosF, maxPosF;
        //        int expand = (uint32_t) sampleSequence.size()-(maxPos-minPos); // not covered by
        //        seeds

        if (SAMPLE_INDEX == debug)
        {
            std::cout << SAMPLE_INDEX << ", minPos " << minPos << ", nseed " << nSeed;
            if (revcomp)
                std::cout << ", REVERSED\n";
            else
                std::cout << ", FORWARD\n";
        }

        int swScore = 0;
        uint32_t refBegin = 0;
        uint32_t refEnd = 0;
        std::span<char> referenceSequence;
        std::string revComplement = "";
        uint32_t slack = 0; //******* to allow for insertion and gaps in aligment ********

        if (GET_SW)
        {
            ///// USING SMITHWATERMAN TO ALIGN SEQUENCES /////
            minPosF = std::fmax(((double)minPos - slack), 0);
            maxPosF = std::fmin(minPos + sequenceLength + slack, (uint32_t)genome.size());

            referenceSequence = std::span<char>(&genome[minPosF], maxPosF - minPosF);

            // use SSW lib code to compute SW
            if (revcomp == 0)
            {
                // try alignment in the forward direction
                alignSequences(referenceSequence, sampleSequence, cigar, swScore, refBegin, refEnd);
            }
            else
            {
                // try alignment in the reverse direction
                revComplement = reverseComplement(sampleSequence);
                alignSequences(referenceSequence, revComplement, cigar, swScore, refBegin, refEnd);
            }
            refBegin += minPos; // adjust to genome position because SW only reports position
                                // relative to referenceSequence
            refEnd += minPos;   //    - " -
        }
        else
        {
            minPosF = std::fmax(((double)minPos - slack), 0);
            maxPosF = std::fmin(minPos + sequenceLength + slack, (uint32_t)genome.size());

            referenceSequence = std::span<char>(&genome[minPosF], maxPosF - minPosF);

            int naiveScore;
            if (revcomp == 0)
            {
                naiveScore = align(referenceSequence, sampleSequence);
            }
            else
            {
                revComplement = reverseComplement(sampleSequence);
                naiveScore = align(referenceSequence, revComplement);
            }
            // Dynamic alpha: scales inversely with cutoff to reduce overestimation at low cutoffs
            float alpha =
                0.9f - 0.004f * cutoff; // Adjusts from ~0.82 at cutoff=20 to ~0.7 at cutoff=50
            alpha = std::max(0.6f, std::min(1.0f, alpha)); // Ensure 0.6 <= alpha <= 1.0
            int estimatedMatches = std::max(
                0, std::min(static_cast<int>(sequenceLength),
                            static_cast<int>(alpha * (4.0f * naiveScore - sequenceLength) / 3.0f)));
            swScore = 2 * estimatedMatches;
            refBegin = minPos;
            refEnd = minPos;
        }

        // Store the result index tuple and the SW score in the resultsSW vector

        int cutOFF = (cutoff * sequenceLength + 99) / 100;
        score = static_cast<float>(swScore) / 2.0f;
        if (static_cast<int>(std::ceil(score)) >= cutOFF)
        {
            if (!hasLongMatch(cigar))
            {
                std::cout << " OOPS!!!\n";
                std::cout << "     read: " << sampleSequence << "\n";
                std::cout << "rev compl: " << revComplement << "\n";
                std::cout << "reference: "
                          << std::string(referenceSequence.begin(), referenceSequence.end())
                          << "\n";
                std::cout << "cigar: " << cigar << "\n";
                if (revcomp == 0)
                    printAlignedSequences(referenceSequence, sampleSequence, cigar, swScore,
                                          refBegin, refEnd);
                else
                    printAlignedSequences(referenceSequence, reverseComplement(sampleSequence),
                                          cigar, swScore, refBegin, refEnd);
                exit(1);
            }
            result.sampleIndex = (uint32_t)SAMPLE_INDEX;
            result.cigar = cigar;
            result.refBegin = refBegin;
            result.score = score;
            result.revComp = revcomp;
            //            std::cout << "enqueue " << SAMPLE_INDEX << "\n";
            queue.enqueue(result); // write out the result
            //            std::string exact = std::to_string(sampleSequence.size()) + "=";
            //            std::cout << "Seq Length " << sequenceLength << ", #seeds " << nSeed << ",
            //            " << cigar << "\n";
        }
        if (SAMPLE_INDEX == debug)
        {
            if (revcomp == 0)
                printAlignedSequences(referenceSequence, sampleSequence, cigar, swScore, refBegin,
                                      refEnd);
            else
                printAlignedSequences(referenceSequence, reverseComplement(sampleSequence), cigar,
                                      swScore, refBegin, refEnd);
        }
        if (!MATCH_ALL)
        {
            skipMulti = true;
            currentSample = sampleIndex;
        }
    }
}
