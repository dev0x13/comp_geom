#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <cmath>

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

    // Finds closest pair in the given band
    inline void checkClosest() {
        auto tmpClosestSize = tmpClosest.size();

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
    std::vector<Point> findClosestPair(const size_t startInd,
                                       const size_t endInd) {

        std::vector<Point> res(endInd - startInd);

        // 1) If there are few points, just brute force them and sort array piece by Y coordinate
        if (endInd - startInd <= 3) {
            bruteForce(startInd, endInd);
            std::sort(data.begin() + startInd, data.begin() + endInd, cmpYAsc);
            res.assign(data.begin() + startInd, data.begin() + endInd);
            return res;
        }

        // 2) Find the middle (by X coordinate) point
        size_t mid = (startInd + endInd) / 2;
        auto midPoint = data[mid];

        // 3) Recursively find the closest pair in two halves
        auto part1 = findClosestPair(startInd, mid);
        auto part2 = findClosestPair(mid, endInd);

        // 4) Merge sorted arrays
        std::merge(part1.begin(),
                   part1.end(),
                   part2.begin(),
                   part2.end(),
                   res.begin(), cmpYAsc);

        // 45 Check points in the band with the minimal distance width
        for (const auto& m: res) {
            if (pow(m.x - midPoint.x, 2) < minPack.dist2) {
                tmpClosest.push_back(m);
            }
        }

        checkClosest();

        return res;
    }
};

int main() {
    // 1) Read data
    auto data = readData("input.dat");

    auto dataSize = data.size();

    // 2) Validate data and init data wrapper
    if (dataSize < 2) {
        std::cout << "Too few points" << std::endl;
        return 1;
    }

    DataWrapper dataWrapper(data);

    // 3) Find closest pair
    dataWrapper.findClosestPair(0, dataSize);

    // 3) Retrieve the result
    auto& minPack = dataWrapper.minPack;
    minPack.dist2 = sqrt(minPack.dist2);

    // 4) Output result
    std::cout.precision(precision);
    std::cout << "d=" << minPack.dist2 << std::endl;
    std::cout << "p1: " << minPack.p1 << std::endl;
    std::cout << "p2: " << minPack.p2 << std::endl;

    if (OUTPUT_TO_FILE) {
        writeResults("result.dat", minPack);
    }

    return 0;
}
