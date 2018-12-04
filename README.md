# Example Usage
```cpp
#define _USE_MATH_DEFINES
#include <iostream>
#include <PSO.h>
#include <vector>
using namespace std;

// benchmark function
double calculateEasom(const vector<double> &x) {
    return -cos(x[0]) * cos(x[1]) * exp(-pow(x[0] - M_PI, 2) - pow(x[1] - M_PI, 2));
}

// utility function to print vector
void printVector(const vector<double> &vec) {
    for (double f : vec) {
        cout << f << ' ';
    }
    cout << endl;
}

int main() {
    // create bounds of search space
    vector<vector<double>> bounds;
    vector<double> d1;
    d1.push_back(-100);
    d1.push_back(100);
    bounds.push_back(d1);
    bounds.push_back(d1);

    // initialise PSO with benchmark, bounds, num of dimensions, num of particles, num of iterations
    PSO pso = PSO(calculateEasom, bounds, 2, 25, 500);

    // run PSO and get result
    vector<double> result = pso.run();

    // print result
    printVector(result);

    return 0;
}
```