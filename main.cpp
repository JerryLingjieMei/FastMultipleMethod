#include <iostream>
#include <random>
#include "FMM.h"

using namespace std;

class BruteForce {
public:
    double *begin, *end;
    double *weightBegin, *weightEnd;
    double h;

    BruteForce(double *begin, double *end, double *weightBegin, double *weightEnd, double h) :
            begin(begin), end(end), weightBegin(weightBegin), weightEnd(weightEnd), h(h) {}

    double evaluate(double x) {
        double result = 0;
        for (double *it = begin, *weightIt = weightBegin; it != end; it++, weightIt++)
            result += *weightIt * sqrt((x - *it) * (x - *it) + h * h);
        return result;
    }
};

int main() {
    int NList[10] = {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000};
    for (int N:NList) {
        cout << "N: " << N << endl;
        double *place = new double[N], *weight = new double[N];
        for (int i = 0; i < N; i++)
            place[i] = (double) i / N;
        double h = 1. / N;
        default_random_engine generator;
        normal_distribution<double> distribution(0, 1. / N);
        for (int i = 0; i < N; i++)
            weight[i] = distribution(generator);
        uniform_real_distribution<double> uniform(0, 1);
        double *test = new double[1000];
        for (int i = 0; i < 1000; i++)
            test[i] = uniform(generator);

        auto time = chrono::high_resolution_clock::now();
        BruteForce MyBruteForce(place, place + N, weight, weight + N, h);
        cout << "Bruteforce build time: " << (chrono::high_resolution_clock::now() - time).count() << endl;
        double *bruteForceResult = new double[1000];
        time = chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000; i++)
            bruteForceResult[i] = MyBruteForce.evaluate(test[i]);
        cout << "Bruteforce run time: " << (chrono::high_resolution_clock::now() - time).count() / 1000
             << endl, time = chrono::high_resolution_clock::now();

        FMM1d MyFMM1(place, place + N, weight, weight + N, h, 0.01);
        cout << "Error <1e-2 build time: " << (chrono::high_resolution_clock::now() - time).count() << endl;
        double *FMMResult1 = new double[1000];
        time = chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000; i++)
            FMMResult1[i] = MyFMM1.evaluate(test[i]);
        cout << "Error <1e-2 run time: " << (chrono::high_resolution_clock::now() - time).count() / 1000
             << endl, time = chrono::high_resolution_clock::now();

        FMM1d MyFMM2(place, place + N, weight, weight + N, h, 0.0001);
        cout << "Error <1e-4 build time: " << (chrono::high_resolution_clock::now() - time).count() << endl;
        double *FMMResult2 = new double[1000];
        time = chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000; i++)
            FMMResult2[i] = MyFMM2.evaluate(test[i]);
        cout << "Error <1e-4 run time: " << (chrono::high_resolution_clock::now() - time).count() / 1000
             << endl, time = chrono::high_resolution_clock::now();

        FMM1d MyFMM3(place, place + N, weight, weight + N, h, 0.000001);
        cout << "Error <1e-6 build time: " << (chrono::high_resolution_clock::now() - time).count() << endl;
        double *FMMResult3 = new double[1000];
        time = chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000; i++)
            FMMResult3[i] = MyFMM3.evaluate(test[i]);
        cout << "Error <1e-6 run time: " << (chrono::high_resolution_clock::now() - time).count() / 1000
             << endl, time = chrono::high_resolution_clock::now();

        FMM1d MyFMM4(place, place + N, weight, weight + N, h, 0.00000001);
        cout << "Error <1e-8 build time: " << (chrono::high_resolution_clock::now() - time).count() << endl;
        double *FMMResult4 = new double[1000];
        time = chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000; i++)
            FMMResult4[i] = MyFMM4.evaluate(test[i]);
        cout << "Error <1e-8 run time: " << (chrono::high_resolution_clock::now() - time).count() / 1000
             << endl, time = chrono::high_resolution_clock::now();

        double error1 = 0, error2 = 0, error3 = 0, error4 = 0;
        for (int i = 0; i < 1000; i++)
            error1 = max(error1, abs(bruteForceResult[i] - FMMResult1[i])),
            error2 = max(error2, abs(bruteForceResult[i] - FMMResult2[i])),
            error3 = max(error3, abs(bruteForceResult[i] - FMMResult3[i])),
            error4 = max(error4, abs(bruteForceResult[i] - FMMResult4[i]));
        cout << "Designed error: 1e-2, Actual error: " << error1 << endl << "Designed error: 1e-4, Actual error: "
             << error2 << endl
             << "Designed error: 1e-6, Actual error: " << error3 << endl << "Designed error: 1e-8, Actual error: "
             << error4 << endl;
    }
    return 0;
}