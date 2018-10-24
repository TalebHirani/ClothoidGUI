/*
Copyright (c) 2014, Enrico Bertolazzi

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef CLOTHOID_H
#define CLOTHOID_H
#include <cmath>
#include <sstream>
#include <vector>


template <typename T = double>
std::vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
  *x = val;
  return xs;
};


namespace Clothoid {

  void FresnelCS(double x, double &C, double &S);

  void FresnelCS(int nk, double x, double C[], double S[]);

  void GeneralizedFresnelCS(int nk, double a, double b, double c, double intC[], double intS[]);

  void GeneralizedFresnelCS(double a, double b, double c, double &intC, double &intS);

  int buildClothoid(double x0, double y0, double theta0, double x1, double y1, double theta1, double &k, double &dk, double &L);

  int pointsOnClothoid(double x0, double y0, double theta0, double kappa, double dkappa, double L, uint npts, std::vector<double> &X, std::vector<double> &Y);

}



#endif // CLOTHOID_H
