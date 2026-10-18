#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iomanip>
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
    "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
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
    string outputFilePath = "results.txt";

    if (argc > 1) inputFilePath = argv[1];
    if (argc > 2) outputFilePath = argv[2];

    auto start_time = chrono::high_resolution_clock::now();

    ifstream inputFile(inputFilePath);
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open input file " << inputFilePath << endl;
        return 1;
    }

    json inputJson;
    inputFile >> inputJson;
    inputFile.close();

    vector<SequenceItem> items;
    for (const auto& el : inputJson) {
        SequenceItem item(
            el.at("sequenceSegment").get<string>(),
            el.at("readDepth").get<int>(),
            el.at("phredScore").get<double>()
        );
        items.push_back(item);
    }

    cout << "Loaded " << items.size() << " sequence items from " << inputFilePath << endl;

    vector<ResultItem> results;
    results.reserve(items.size());

    for (const auto& item : items) {
        int score = compute_alignment_score(item.sequenceSegment, BASE_REFERENCE_SEQUENCE);
        results.emplace_back(item, score);
    }

    ofstream outputFile(outputFilePath);
    if (!outputFile.is_open()) {
        cerr << "Error: Could not open output file " << outputFilePath << endl;
        return 1;
    }

    const size_t MAX_PREFIX_LEN = 10;

    const int col1 = 18;
    const int col2 = 12;
    const int col3 = 12;
    const int col4 = 16;

    outputFile << "| " << left << setw(col1) << "SequenceSegment"
               << " | " << left << setw(col2) << "ReadDepth"
               << " | " << left << setw(col3) << "PhredScore"
               << " | " << left << setw(col4) << "AlignmentScore" << " |\n";

    outputFile << "|-" << string(col1, '-')
               << "-|- " << string(col2, '-')
               << "-|- " << string(col3, '-')
               << "-|- " << string(col4, '-') << " -|\n";

    for (const auto& res : results) {
        bool accepted = (res.item.phredScore >= 30.0 && res.alignmentScore > 50);
        if (accepted) {
            outputFile << "| " << left << setw(col1) << res.item.sequenceSegment.substr(0, MAX_PREFIX_LEN)
                       << " | " << left << setw(col2) << res.item.readDepth
                       << " | " << left << setw(col3) << fixed << setprecision(1) << res.item.phredScore
                       << " | " << left << setw(col4) << res.alignmentScore
                       << " |\n";
        }
    }
    outputFile.close();

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "\n==========================================" << endl;
    cout << "Sequential Execution Completed." << endl;
    cout << "Processed Items: " << results.size() << endl;
    cout << "Results saved to: " << outputFilePath << endl;
    cout << "Total Execution Time: " << elapsed.count() << " seconds" << endl;
    cout << "==========================================" << endl;

    return 0;
}