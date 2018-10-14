#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

#define OUTPUT_TO_FILE true

static constexpr double doubleBound = 0.67 * 10E159;

// basic point struct
struct Point {
  double x;
  double y;

  inline double dist2(const Point& p) const {
    double subX = x - p.x, subY = y - p.y;
    assert(abs(subX) <= doubleBound && abs(subY) <= doubleBound && "Double overflow asserted");

    return pow(subX, 2) + pow(subY, 2);
  }

  inline double length2() const {
    return pow(x, 2) + pow(y, 2);
  }

  friend std::ostream& operator<<(std::ostream& s, const Point& p) {
    s << "x=" << p.x << ",y=" << p.y;
    return s;
  }
};

// basic point struct
struct Caliper {
  size_t pivotInd;
  Point vec;

  void rotate(const std::pair<double, double>& angleCosSin) {
    double x = vec.x, y = vec.y;
    vec.x = x * angleCosSin.first - y * angleCosSin.second;
    vec.y = y * angleCosSin.first + x * angleCosSin.second;
  }

  // finds angle between vectors (0 < angle < M_PI)
  inline auto angleCosSin(const Point& p1, const Point& p2) {
    double x1 = p2.x - p1.x, y1 = p2.y - p1.y;
    assert(abs(x1) <= doubleBound && abs(y1) <= doubleBound && "Double overflow asserted");
    
    double length = sqrt(vec.length2() * (pow(x1, 2) + pow(y1, 2)));
    assert(length != 0 && "Division by zero asserted");

    return std::make_pair(
      (x1 * vec.x + y1 * vec.y) / length,
      (y1 * vec.x - x1 * vec.y) / length);
  }
};

struct DataWrapper {
private:
  Caliper caliper1, caliper2;
  size_t dataSize;
  std::vector<Point> data;

  // looped index increment
  inline size_t wiseIndIncr(size_t ind) {
    if (++ind == dataSize) {
      return 0;
    }

    return ind;
  }

public:
  DataWrapper(size_t ind1, size_t ind2, const std::vector<Point>& data) :
              data(data),
              dataSize(data.size()) {
    // init calipers with parallel lines
    caliper1.pivotInd = ind1;
    caliper1.vec = Point{1, 0};
    caliper2.pivotInd = ind2;
    caliper2.vec = Point{-1, 0};
  }

  // yeilds next antipodal pair
  inline std::pair<size_t, size_t> nextPair() {
    auto angle1 = caliper1.angleCosSin(data[caliper1.pivotInd],
      data[wiseIndIncr(caliper1.pivotInd)]);
    auto angle2 = caliper2.angleCosSin(data[caliper2.pivotInd],
      data[wiseIndIncr(caliper2.pivotInd)]);

    std::pair<double, double> angle;

    if (angle1.first > angle2.first) {
      caliper1.pivotInd = wiseIndIncr(caliper1.pivotInd);
      angle = angle1;
    }
    else if (angle1.first < angle2.first) {
      caliper2.pivotInd = wiseIndIncr(caliper2.pivotInd);
      angle = angle2;
    }
    else {
      caliper1.pivotInd = wiseIndIncr(caliper1.pivotInd);
      caliper2.pivotInd = wiseIndIncr(caliper2.pivotInd);
      angle = angle1;
    }
 
    caliper1.rotate(angle);
    caliper2.rotate(angle);

    return std::make_pair(caliper1.pivotInd, caliper2.pivotInd);
  }
};

// reads data from the input file
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

// write results to the output file
void writeResults(const std::string& filename, double dist,
                  const Point& p1, const Point& p2) {
  std::ofstream os(filename);

  if (!os.is_open()) {
    throw std::runtime_error("Invalid  output file!");
  }
  else {
    os << dist << std::endl;
    os << p1.x << " " << p1.y << std::endl;
    os << p2.x << " " << p2.y << std::endl;
  }

  os.close();
}

int main() {
  // 1) Read data
  auto data = readData("data.dat");

  if (data.size() < 3) {
    std::cout << "Too few points" << std::endl;
    return 1;
  }

  // 2) Find highest and lowest points
  Point *lowest = &data[0], *highest = &data[0];

  size_t lowestInd = 0, highestInd = 0;

  for (size_t i = 1; i < data.size(); ++i) {
    // the first condition handles the special case with a square
    if (data[i].y == highest->y) {
      if (data[i].x > highest->x) {
        highest = &data[i];
        highestInd = i;
      }
    }
    else if (data[i].y < lowest->y) {
      lowest = &data[i];
      lowestInd = i;
    }
    else if (data[i].y > highest->y) {
      highest = &data[i];
      highestInd = i;
    }
  }

  // 3) Find pairs and max distance
  double maxDist2 = data[lowestInd].dist2(data[highestInd]), dist2;
  size_t maxInd1 = lowestInd, maxInd2 = highestInd;
  DataWrapper dataWrapper(lowestInd, highestInd, data);

  auto nextPair = std::make_pair(lowestInd, highestInd);

  do {
    dist2 = data[nextPair.first].dist2(data[nextPair.second]);

    if (dist2 > maxDist2) {
      maxDist2 = dist2;
      maxInd1 = nextPair.first;
      maxInd2 = nextPair.second;
    }

    nextPair = dataWrapper.nextPair();
  } while (!(nextPair.first == lowestInd && nextPair.second == highestInd));

  maxDist2 = sqrt(maxDist2);

  // 4) Output results
  std::cout << "d=" << maxDist2 << std::endl;
  std::cout << "p1: " << data[maxInd1] << std::endl;
  std::cout << "p2: " << data[maxInd2] << std::endl;

  if (OUTPUT_TO_FILE) {
    writeResults("result.dat", maxDist2, data[maxInd1], data[maxInd2]);
  }

  return 0;
}