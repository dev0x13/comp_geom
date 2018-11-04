#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>

#define OUTPUT_TO_FILE true

static constexpr double doubleBound = 0.67 * 10E159;
static constexpr auto precision = std::numeric_limits<double>::max_digits10;

// Basic point struct
struct Point {
    double x;
    double y;

    // Calculates squared distance
    inline double dist2(const Point& p) const {
        double subX = x - p.x, subY = y - p.y;
        assert(fabs(subX) <= doubleBound && fabs(subY) <= doubleBound && "Double overflow asserted");

        return pow(subX, 2) + pow(subY, 2);
    }

    // Calculates triangle area
    static inline double area(const Point& p0, const Point& p1, const Point& p2) {
        return fabs((p1.x - p0.x) * (p2.y - p0.y) -
                    (p1.y - p0.y) * (p2.x - p0.x));
    }

    friend std::ostream& operator<<(std::ostream& s, const Point& p) {
        s << "x=" << p.x << ",y=" << p.y;
        return s;
    }
};

// Points pair and distance between them pack structure
struct PointPack {
    double dist2 = 0;
    const Point *p1 = nullptr, *p2 = nullptr;
};

// Data wrapper struct. Provides wise index increment and diameter candidates processing.
struct DataWrapper {
private:
    size_t dataSize;
    const std::vector<Point>& data;

public:
    PointPack maxPack;

    DataWrapper(const std::vector<Point>& data): data(data) {
        dataSize = data.size();
    }

    // Wise looped index increment
    inline const Point& operator[](size_t i) const {
        return data[i % dataSize];
    }

    // Checks diameter points pair candidate
    inline void checkPointPair(size_t ind1, size_t ind2) {
        const Point &p1 = (*this)[ind1],
                &p2 = (*this)[ind2];
        auto d = p1.dist2(p2);

        if (d > maxPack.dist2) {
            maxPack.dist2 = d;
            maxPack.p1 = &p1;
            maxPack.p2 = &p2;
        }
    }
};

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
void writeResults(const std::string& filename, const PointPack& maxPack) {
    std::ofstream os(filename);
    os.precision(precision);

    if (!os.is_open()) {
        throw std::runtime_error("Invalid  output file!");
    }
    else {
        os << maxPack.dist2 << std::endl;
        os << maxPack.p1->x << " " << maxPack.p1->y << std::endl;
        os << maxPack.p2->x << " " << maxPack.p2->y << std::endl;
    }

    os.close();
}

int main() {
    // 1) Read data
    auto data = readData("data.dat");

    auto dataSize = data.size();

    // 2) Validate data and init data wrapper
    if (dataSize < 3) {
        std::cout << "Too few points" << std::endl;
        return 1;
    }

    DataWrapper dw(data);

    // 3) Find diameter

    // 3.1) Init calipers
    size_t j = 1, j0;
    while (Point::area(dw[dataSize - 1], dw[0], dw[j + 1]) >= Point::area(dw[dataSize - 1], dw[0], dw[j])) {
        j++;
    }

    j0 = j;

    // 3.2) Process all the antipodal pairs
    for (size_t i = 0; i <= j0 && j <= dataSize; ++i) {
        dw.checkPointPair(i, j);

        while (j < dataSize && Point::area(dw[i], dw[i + 1], dw[j + 1]) > Point::area(dw[i], dw[i + 1], dw[j])) {
            dw.checkPointPair(i, ++j);
        }

        // handle parallel enges case
        if (j < dataSize && Point::area(dw[i], dw[i + 1], dw[j + 1]) == Point::area(dw[i], dw[i + 1], dw[j])) {
            dw.checkPointPair(i, j + 1);
        }
    }

    // 3) Retrieve the result and validate it
    auto& maxPack = dw.maxPack;
    assert(maxPack.p1 && maxPack.p2);
    maxPack.dist2 = sqrt(maxPack.dist2);

    // 4) Output result
    std::cout.precision(precision);
    std::cout << "d=" << maxPack.dist2 << std::endl;
    std::cout << "p1: " << *maxPack.p1 << std::endl;
    std::cout << "p2: " << *maxPack.p2 << std::endl;

    if (OUTPUT_TO_FILE) {
        writeResults("result.dat", maxPack);
    }

    return 0;
}
