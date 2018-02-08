//
// Created by Jerry Mei on 2/8/18.
//

#ifndef FINAL_PROJECT_FFM_H
#define FINAL_PROJECT_FFM_H

#endif //FINAL_PROJECT_FFM_H

#include <cmath>
#include <algorithm>

class FMM1dNode {
public:
    double lower, upper;
    double *begin, *end;
    double *weightBegin, *weightEnd;
    double *coeff, h;
    int p;
    FMM1dNode *lChild = nullptr, *rChild = nullptr;

    FMM1dNode(double lower, double upper, double *begin, double *end, double *weightBegin, double *weightEnd,
              double h, int p) {
        this->lower = lower, this->upper = upper, this->begin = begin, this->end = end,
        this->weightBegin = weightBegin, this->weightEnd = weightEnd, this->h = h, this->p = p, coeff = new double[p +
                                                                                                                   2]();
    }

    void buildCoeff() {
        if (lChild == nullptr)
            for (double *it = begin, *weightIt = weightBegin; it != end; it++, weightIt++) {
                coeff[0] += *weightIt;
                coeff[1] -= *weightIt * (*it - (lower + upper) / 2);
                int *binom = new int[p]();
                for (int i = 0; i < p; i++)
                    binom[i] = 1;
                double alpha = *weightIt, t = *it - (lower + upper) / 2;
                for (int i = 1; i <= p; i += 2) {
                    alpha *= (.5 - (i >> 1)) / ((i >> 1) + 1) * h * h;
                    double temp = alpha;
                    for (int j = 0; i + j <= p; j++)
                        coeff[i + j + 1] += binom[j] * temp, temp *= t;
                    for (int j = 0; j < 2; j++)
                        for (int k = 1; k < p; k++)
                            binom[k] += binom[k - 1];
                }
            }
        else {
            double t1 = (lower - upper) / 4, t2 = (upper - lower) / 4;
            coeff[0] = lChild->coeff[0] + rChild->coeff[0];
            coeff[1] = lChild->coeff[1] + rChild->coeff[1] - t1 * lChild->coeff[0] - t2 * rChild->coeff[0];
            int *binom = new int[p]();
            for (int i = 0; i < p; i++)
                binom[i] = 1;
            for (int i = 1; i <= p; i++) {
                double temp1 = lChild->coeff[i + 1], temp2 = rChild->coeff[i + 1];
                for (int j = 0; i + j <= p; j++)
                    coeff[i + j + 1] += binom[j] * (temp1 + temp2), temp1 *= t1, temp2 *= t2;
                for (int k = 1; k < p; k++)
                    binom[k] += binom[k - 1];
            }
        }

    }

    double evaluate(double x) {
        if (upper * 2 - lower <= x || lower * 2 - upper >= x) {
            double middle = (upper + lower) / 2;
            double result = 0, x0 = 1. / (x - middle);
            for (int i = p + 1; i >= 1; i--)
                result *= x0, result += coeff[i];
            result += coeff[0] * (x - middle);
            if (x > middle)
                return result;
            else
                return -result;
        } else if (lChild == nullptr) {
            double result = 0;
            for (double *it = begin, *weightIt = weightBegin; it != end; it++, weightIt++)
                result += *weightIt * sqrt((x - *it) * (x - *it) + h * h);
            return result;
        } else
            return lChild->evaluate(x) + rChild->evaluate(x);
    }

    ~FMM1dNode() {
        delete coeff;
    }
};

class FMM1d {
public:
    int M, p;
    double eps;
    FMM1dNode *nodes;

    FMM1d(double *begin, double *end, double *weightBegin, double *weightEnd, double h, double eps) :
            M((int) round(log2(end - begin))), eps(eps) {
        double weightAll = 0;
        for (auto iter = weightBegin; iter != weightEnd; iter++)
            weightAll += abs(*iter);
        p = int(log2(8 * weightAll / eps) / log2(3 / sqrt(2)));
        nodes = (FMM1dNode *) operator new(sizeof(FMM1dNode) * (1 << (M + 1)));
        for (int i = 1; i < 1 << (M + 1); i++) {
            FMM1dNode *node = &nodes[i];
            if (i == 1)
                new(node) FMM1dNode(0, 1, begin, end, weightBegin, weightEnd, h, p);
            else {
                FMM1dNode *parent = &nodes[i / 2];
                if (i & 1)
                    new(node) FMM1dNode((parent->lower + parent->upper) / 2, parent->upper, (node - 1)->end,
                                        parent->end, (node - 1)->weightEnd, parent->weightEnd,
                                        h, p), parent->rChild = node;
                else {
                    double *mid = std::lower_bound(parent->begin, parent->end, (parent->lower + parent->upper) / 2);
                    new(node) FMM1dNode(parent->lower, (parent->lower + parent->upper) / 2, parent->begin, mid,
                                        parent->weightBegin, mid - parent->begin + parent->weightBegin,
                                        h, p), parent->lChild = node;
                }
            }
        }
        for (int i = (1 << (M + 1)) - 1; i >= 1; i--)
            nodes[i].buildCoeff();
    }

    double evaluate(double x) {
        return nodes[1].evaluate(x);
    }

    ~FMM1d() {
        delete nodes;
    }

};