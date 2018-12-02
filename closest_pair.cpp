#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <chrono>

#define OUTPUT_TO_FILE true

static constexpr double doublePrecisionLossBound = 9007199254740991; // 2^53
static constexpr auto precision = std::numeric_limits<double>::max_digits10;

// Basic point struct
struct Point {
    double x;
    double y;

    // Calculates squared distance
    inline double dist2(const Point& p) const {
        double subX2 = pow(x - p.x, 2), subY2 = pow(y - p.y, 2);

        assert(subX2 <= doublePrecisionLossBound && subY2 <= doublePrecisionLossBound && "Double precision loss asserted");

        return subX2 + subY2;
    }

    friend std::ostream& operator<<(std::ostream& s, const Point& p) {
        s << "x=" << p.x << ",y=" << p.y;
        return s;
    }
};

// Points pair and distance between them pack structure
struct PointPack {
    double dist2 = std::numeric_limits<double>::infinity();
    Point p1, p2;
};

// Helper predicates for sorting
inline bool cmpXAsc(const Point& p1, const Point& p2) {
    return p1.x < p2.x;
}

inline bool cmpYAsc(const Point& p1, const Point& p2) {
    return p1.y < p2.y;
}

// Reads data from the input file
std::vector<Point> readData(const std::string& filename) {
    std::ifstream is(filename);

    std::vector<Point> data;

    if (!is.is_open()) {
        throw std::runtime_error("Invalid  input file!");
    }
    else {
        Point p;

        while (is >> p.x >> p.y) {
            data.push_back(p);
        }
    }

    is.close();

    return data;
}

// Writes results to the output file
void writeResults(const std::string& filename, const PointPack& minPack) {
    std::ofstream os(filename);
    os.precision(precision);

    if (!os.is_open()) {
        throw std::runtime_error("Invalid  output file!");
    }
    else {
        os << minPack.dist2 << std::endl;
        os << minPack.p1.x << " " << minPack.p1.y << std::endl;
        os << minPack.p2.x << " " << minPack.p2.y << std::endl;
    }

    os.close();
}

struct DataWrapper {
private:
    std::vector<Point>& data;
    size_t dataSize;
    std::vector<Point> tmpClosest;

public:
    PointPack minPack;

    DataWrapper(std::vector<Point>& data): data(data) {
        this->dataSize = data.size();

        // 1) Sort data by X coordinate
        std::sort(data.begin(), data.end(), cmpXAsc);

        // 2) Allocate memory for temporary indexes storage
        this->tmpClosest.reserve(dataSize);
    }

    // Check if given pair closet to each other than old one
    inline void updMinPack(const Point& p1, const Point& p2) {
        auto dist2 = p1.dist2(p2);

        if (dist2 < minPack.dist2) {
            minPack.dist2 = dist2;
            minPack.p1 = p1;
            minPack.p2 = p2;
        }
    }

    // Brute forces O(N^2) closest pair with the points in given range
    inline void bruteForce(size_t startInd, size_t endInd) {
        for (size_t i = startInd; i < endInd; ++i) {
            for (size_t j = i + 1; j < endInd; ++j) {
                updMinPack(data[i], data[j]);
            }
        }
    }

    // Brute forces closest pair with the points with given indexes
    inline void checkClosest() {
        auto tmpClosestSize = tmpClosest.size();

        std::sort(tmpClosest.begin(), tmpClosest.end(), cmpYAsc);

        for (size_t i = 0; i < tmpClosestSize; ++i) {
            for (size_t j = i + 1;
                 j < tmpClosestSize && pow((tmpClosest[j].y - tmpClosest[i].y), 2) < minPack.dist2;
                 ++j) {
                updMinPack(tmpClosest[i], tmpClosest[j]);
            }
        }

        tmpClosest.clear();
    }

    // Recursively finds the closest pair in given set
    void findClosestPair(const size_t startInd,
                         const size_t endInd) {

        // 1) If there are few points, just brute force them and
        if (endInd - startInd <= 3) {
            bruteForce(startInd, endInd);
            return;
        }

        size_t mid = (startInd + endInd) / 2;
        auto& midPoint = data[mid];

        findClosestPair(startInd, mid);
        findClosestPair(mid, endInd);

        for (size_t i = 0; i < dataSize; ++i) {
            if (fabs(data[i].x - midPoint.x) < minPack.dist2) {
                tmpClosest.push_back(data[i]);
            }
        }

        checkClosest();
    }
};

int main() {
    // 1) Read data
    //auto data = readData("data1.dat");

    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0,1000000); // guaranteed unbiased

    std::vector<Point> data;

    for (int i = 0; i < 10000; ++i) {
        Point p;
        p.x = uni(rng);
        p.y = uni(rng);
        data.push_back(p);
    }

    auto dataSize = data.size();

    // 2) Validate data and init data wrapper
    if (dataSize < 2) {
        std::cout << "Too few points" << std::endl;
        return 1;
    }

    DataWrapper dataWrapper(data);

    using clock = std::chrono::system_clock;
    using ms = std::chrono::milliseconds;

    auto before = clock::now();
    // 3) Find closest pair
    dataWrapper.findClosestPair(0, dataSize - 1);
    auto duration = std::chrono::duration_cast<ms>(clock::now() - before);

    std::cout << "It took " << duration.count()/1000.0 << "s" << std::endl;
    // 3) Retrieve the result
    auto& minPack = dataWrapper.minPack;
    minPack.dist2 = sqrt(minPack.dist2);

    // 4) Output result
    std::cout.precision(precision);
    std::cout << "d=" << minPack.dist2 << std::endl;
    std::cout << "p1: " << minPack.p1 << std::endl;
    std::cout << "p2: " << minPack.p2 << std::endl;

    DataWrapper dataWrapper1(data);

    before = clock::now();
    // 3) Find closest pair
    dataWrapper1.bruteForce(0, dataSize - 1);
    duration = std::chrono::duration_cast<ms>(clock::now() - before);

    std::cout << "It took " << duration.count()/1000.0 << "s" << std::endl;

    // 3) Retrieve the result
    auto& minPack1 = dataWrapper1.minPack;
    minPack1.dist2 = sqrt(minPack1.dist2);

    // 4) Output result
    std::cout.precision(precision);
    std::cout << "d=" << minPack1.dist2 << std::endl;
    std::cout << "p1: " << minPack1.p1 << std::endl;
    std::cout << "p2: " << minPack1.p2 << std::endl;

    if (OUTPUT_TO_FILE) {
        writeResults("result.dat", minPack);
    }

    return 0;
}
