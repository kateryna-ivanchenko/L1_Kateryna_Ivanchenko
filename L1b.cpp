#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include "json.hpp"

using json = nlohmann::json;
using namespace std;

struct SequenceItem {
    string sequenceSegment;
    int readDepth;
    double phredScore;

    SequenceItem() : sequenceSegment(""), readDepth(0), phredScore(0.0) {}
    SequenceItem(string seq, int depth, double phred)
        : sequenceSegment(seq), readDepth(depth), phredScore(phred) {}
};

struct ResultItem {
    SequenceItem item;
    int alignmentScore;

    ResultItem() : item(), alignmentScore(0) {}
    ResultItem(const SequenceItem& inputItem, int score)
        : item(inputItem), alignmentScore(score) {}
};

const string BASE_REFERENCE_SEQUENCE =
    "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC";

struct AlgoScores {
    int matchScore = 2;
    int mismatchScore = -1;
    int gapScore = -1;
} algoScores;

int compute_alignment_score(const string& query, const string& reference) {
    int rows = query.size() + 1;
    int cols = reference.size() + 1;
    vector<vector<int>> scoring_matrix(rows, vector<int>(cols, 0));

    int maxScore = 0;

    for (int row = 1; row < rows; ++row) {
        for (int col = 1; col < cols; ++col) {
            int n1 = (query[row - 1] == reference[col - 1])
                     ? scoring_matrix[row - 1][col - 1] + algoScores.matchScore
                     : scoring_matrix[row - 1][col - 1] + algoScores.mismatchScore;

            int n2 = scoring_matrix[row - 1][col] + algoScores.gapScore;
            int n3 = scoring_matrix[row][col - 1] + algoScores.gapScore;

            int score = max({0, n1, n2, n3});
            scoring_matrix[row][col] = score;

            if (score > maxScore) {
                maxScore = score;
            }
        }
    }
    return maxScore;
}

int main(int argc, char* argv[]) {
    string inputFilePath = "data.json";
    string outputFilePath = "Group_LastnameF_L1_rez.txt";
    int numThreads = 10;

    if (argc > 1) inputFilePath = argv[1];
    if (argc > 2) outputFilePath = argv[2];
    if (argc > 3) numThreads = stoi(argv[3]);

    auto start_time = chrono::high_resolution_clock::now();

    ifstream inputFile(inputFilePath);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open input file " << inputFilePath << endl;
        return 1;
    }

    json inputJson;
    inputFile >> inputJson;
    inputFile.close();

    vector<SequenceItem> loadedItems;
    for (const auto& el : inputJson) {
        loadedItems.emplace_back(
            el.at("sequenceSegment").get<string>(),
            el.at("readDepth").get<int>(),
            el.at("phredScore").get<double>()
        );
    }

    const auto total_data = loadedItems.size();
    cout << "Loaded " << total_data << " items from " << inputFilePath << endl;

    omp_lock_t result_lock;
    omp_init_lock(&result_lock);

    vector<ResultItem> sortedResults;
    long long global_int_sum = 0;
    double global_float_sum = 0.0;

    omp_set_num_threads(numThreads);

    #pragma omp parallel default(none) shared(loadedItems, sortedResults, result_lock, total_data, BASE_REFERENCE_SEQUENCE) reduction(+:global_int_sum, global_float_sum)
    {
        const auto num_threads = static_cast<size_t>(omp_get_num_threads());
        const auto per_thread = total_data / num_threads;
        const auto thread_num = static_cast<size_t>(omp_get_thread_num());

        const auto start = thread_num * per_thread;
        const auto end = (thread_num == num_threads - 1) ? total_data : (thread_num + 1) * per_thread;

        vector<ResultItem> local_matches;

        for (size_t i = start; i < end; ++i) {
            const auto& item = loadedItems[i];
            int score = compute_alignment_score(item.sequenceSegment, BASE_REFERENCE_SEQUENCE);

            if (item.phredScore >= 30.0 && score > 50) {
                local_matches.emplace_back(item, score);
                global_int_sum += item.readDepth;
                global_float_sum += item.phredScore;
            }
        }

        omp_set_lock(&result_lock);
        for (const auto& match : local_matches) {
            auto it = upper_bound(sortedResults.begin(), sortedResults.end(), match,
                [](const ResultItem& a, const ResultItem& b) {
                    return a.alignmentScore > b.alignmentScore;
                });
            sortedResults.insert(it, match);
        }
        omp_unset_lock(&result_lock);
    }

    omp_destroy_lock(&result_lock);

    ofstream outputFile(outputFilePath);
    if (!outputFile.is_open()) {
        cerr << "Error: Could not open output file " << outputFilePath << endl;
        return 1;
    }

    const size_t MAX_PREFIX_LEN = 10;
    const int col1 = 18, col2 = 12, col3 = 12, col4 = 16;

    outputFile << "| " << left << setw(col1) << "SequenceSegment"
               << " | " << left << setw(col2) << "ReadDepth"
               << " | " << left << setw(col3) << "PhredScore"
               << " | " << left << setw(col4) << "AlignmentScore" << " |\n";

    outputFile << "|-" << string(col1, '-')
               << "-|- " << string(col2, '-')
               << "-|- " << string(col3, '-')
               << "-|- " << string(col4, '-') << " -|\n";

    for (const auto& res : sortedResults) {
        outputFile << "| " << left << setw(col1) << res.item.sequenceSegment.substr(0, MAX_PREFIX_LEN)
                   << " | " << left << setw(col2) << res.item.readDepth
                   << " | " << left << setw(col3) << fixed << setprecision(1) << res.item.phredScore
                   << " | " << left << setw(col4) << res.alignmentScore
                   << " |\n";
    }

    outputFile << "\nSums of filtered items:\n";
    outputFile << "Total ReadDepth (int sum): " << global_int_sum << "\n";
    outputFile << "Total PhredScore (float sum): " << fixed << setprecision(2) << global_float_sum << "\n";
    outputFile.close();

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "\n==========================================" << endl;
    cout << "OpenMP Execution Completed." << endl;
    cout << "Worker Threads: " << numThreads << endl;
    cout << "Results saved to: " << outputFilePath << endl;
    cout << "Total Execution Time: " << elapsed.count() << " seconds" << endl;
    cout << "==========================================" << endl;

    return 0;
}