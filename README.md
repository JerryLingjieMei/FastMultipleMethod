# Fast Multiple Method

[Fast Multiple Method](https://en.wikipedia.org/wiki/Fast_multipole_method) (FMM) is a fast technique to calculate forces in multiple positions of a n-body system. This project is an implementation of FMM in 1D. An brief introduction to this method can be found in [this short course](https://www.math.nyu.edu/faculty/greengar/shortcourse_fmm.pdf). 

## File Descriptions

### FMM.h

[FMM.h](FMM.h) is the header containing necessary classes for implementing FFM in 1D. The ```FMM1d``` gives a class to intialize a n-body system with an array of weights in differnt places, as well as a parameter ```h``` and maximal error ```eps```. 

```cpp
FMM1d(double *begin, double *end, double *weightBegin, double *weightEnd, double h, double eps)
```

This class has equipped with a method ```FMM1d::evaluate(double x)``` to evalute the total force at position ```x```. 

### main.cpp & RadialBasis.pdf

[main.cpp](main.cpp) contains test to meet the technical requirements in [RadialBasis.pdf](Radial Basis.pdf). It shows that the FMM has a run time complexity of $O(n\log(1/\varepsilon))$ with the maximal error $\varepsilon$, compared to a $O(n^2)$ of the brute force. 

### 18_330_Final.pdf   

[18_330_Final.pdf ](18_330_Final.pdf ) is a detailed report on how the algorithm is implemented and comparative results.

## Acknowledgements

This project is inspired by [Matthias Taus](taus@mit.edu) from MIT Mathematics in the 18.330 class during the fall term in 2017.  