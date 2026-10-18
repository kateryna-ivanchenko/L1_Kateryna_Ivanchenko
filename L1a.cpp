#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
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

struct ThreadStats {
    int threadId = 0;
    int processedCount = 0;
    int matchedCount = 0;
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

class DataMonitor {
private:
    vector<SequenceItem> buffer;
    size_t capacity;
    size_t head = 0;
    size_t tail = 0;
    size_t count = 0;
    bool finished = false;

    mutex lock;
    condition_variable cv_not_full;
    condition_variable cv_not_empty;

public:
    explicit DataMonitor(size_t maxCapacity) : capacity(maxCapacity) {
        buffer.resize(capacity);
    }

    void addItem(const SequenceItem& item) {
        unique_lock<mutex> guard(lock);
        cv_not_full.wait(guard, [this]() { return count < capacity; });

        buffer[tail] = item;
        tail = (tail + 1) % capacity;
        count++;

        cv_not_empty.notify_one();
    }

    bool removeItem(SequenceItem& item) {
        unique_lock<mutex> guard(lock);
        cv_not_empty.wait(guard, [this]() { return count > 0 || finished; });

        if (count == 0 && finished) {
            return false;
        }

        item = buffer[head];
        head = (head + 1) % capacity;
        count--;

        cv_not_full.notify_one();
        return true;
    }

    void setFinished() {
        unique_lock<mutex> guard(lock);
        finished = true;
        cv_not_empty.notify_all();
    }
};

class SortedResultMonitor {
private:
    vector<ResultItem> results;
    mutable mutex lock;

public:
    void addItemSorted(const ResultItem& item) {
        lock_guard<mutex> guard(lock);
        auto it = upper_bound(results.begin(), results.end(), item,
            [](const ResultItem& a, const ResultItem& b) {
                return a.alignmentScore > b.alignmentScore;
            });
        results.insert(it, item);
    }

    vector<ResultItem> getItems() const {
        lock_guard<mutex> guard(lock);
        return results;
    }
};

void workerFunction(DataMonitor& dataMon, SortedResultMonitor& resultMon, ThreadStats& stats) {
    SequenceItem item;
    while (dataMon.removeItem(item)) {
        stats.processedCount++;

        int score = compute_alignment_score(item.sequenceSegment, BASE_REFERENCE_SEQUENCE);

        if (item.phredScore >= 30.0 && score > 50) {
            resultMon.addItemSorted(ResultItem(item, score));
            stats.matchedCount++;
        }
    }
}

int main(int argc, char* argv[]) {
    string inputFilePath = "data1.json";
    string outputFilePath = "IF4_IvanchenkoK_L1_rez.txt";
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

    size_t totalData = loadedItems.size();
    cout << "Loaded " << totalData << " items from " << inputFilePath << endl;

    size_t monitorCapacity = max<size_t>(1, totalData / 2);
    DataMonitor dataMonitor(monitorCapacity);
    SortedResultMonitor resultMonitor;

    vector<thread> workers;
    vector<ThreadStats> threadStats(numThreads);

    for (int i = 0; i < numThreads; ++i) {
        threadStats[i].threadId = i + 1;
        workers.emplace_back(workerFunction, ref(dataMonitor), ref(resultMonitor), ref(threadStats[i]));
    }

    for (const auto& item : loadedItems) {
        dataMonitor.addItem(item);
    }
    dataMonitor.setFinished();

    for (auto& worker : workers) {
        worker.join(); 
    }

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

    for (const auto& res : resultMonitor.getItems()) {
        outputFile << "| " << left << setw(col1) << res.item.sequenceSegment.substr(0, MAX_PREFIX_LEN)
                   << " | " << left << setw(col2) << res.item.readDepth
                   << " | " << left << setw(col3) << fixed << setprecision(1) << res.item.phredScore
                   << " | " << left << setw(col4) << res.alignmentScore
                   << " |\n";
    }

    outputFile << "\n";
    for (const auto& stats : threadStats) {
        outputFile << "Thread " << stats.threadId << ": "
                   << stats.matchedCount << " items passed, "
                   << stats.processedCount << " total\n";
    }
    outputFile.close();

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "\n==========================================" << endl;
    cout << "Parallel Execution Completed." << endl;
    cout << "Worker Threads: " << numThreads << endl;
    cout << "Results saved to: " << outputFilePath << endl;
    cout << "Total Execution Time: " << elapsed.count() << " seconds" << endl;
    cout << "==========================================" << endl;

    return 0;
}