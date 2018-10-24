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


#include "clothoid.h"

namespace Clothoid {

  using namespace std;

  static const double fn[] = { 0.49999988085884732562,
                                  1.3511177791210715095,
                                  1.3175407836168659241,
                                  1.1861149300293854992,
                                  0.7709627298888346769,
                                  0.4173874338787963957,
                                  0.19044202705272903923,
                                  0.06655998896627697537,
                                  0.022789258616785717418,
                                  0.0040116689358507943804,
                                  0.0012192036851249883877 };

  static const double fd[] = { 1.0,
                                  2.7022305772400260215,
                                  4.2059268151438492767,
                                  4.5221882840107715516,
                                  3.7240352281630359588,
                                  2.4589286254678152943,
                                  1.3125491629443702962,
                                  0.5997685720120932908,
                                  0.20907680750378849485,
                                  0.07159621634657901433,
                                  0.012602969513793714191,
                                  0.0038302423512931250065 };

  static const double gn[] = { 0.50000014392706344801,
                                  0.032346434925349128728,
                                  0.17619325157863254363,
                                  0.038606273170706486252,
                                  0.023693692309257725361,
                                  0.007092018516845033662,
                                  0.0012492123212412087428,
                                  0.00044023040894778468486,
                                 -8.80266827476172521e-6,
                                 -1.4033554916580018648e-8,
                                  2.3509221782155474353e-10 };

  static const double gd[] = { 1.0,
                                  2.0646987497019598937,
                                  2.9109311766948031235,
                                  2.6561936751333032911,
                                  2.0195563983177268073,
                                  1.1167891129189363902,
                                  0.57267874755973172715,
                                  0.19408481169593070798,
                                  0.07634808341431248904,
                                  0.011573247407207865977,
                                  0.0044099273693067311209,
                                 -0.00009070958410429993314 };

  double const pi = 3.14159265358979323846264338328;
  void FresnelCS(double y, double &C, double &S) {

    double const eps = 1E-15;
    double const x   = y > 0 ? y : -y;

    if (x < 1.0) {
      double twofn, fact, denterm, numterm, sum, term;

      double const s = (pi/2)*(x*x);
      double const t = -s*s;

      twofn   =  0.0;
      fact    =  1.0;
      denterm =  1.0;
      numterm =  1.0;
      sum     =  1.0;
      do {
        twofn   += 2.0;
        fact    *= twofn*(twofn-1.0);
        denterm += 4.0;
        numterm *= t;
        term     = numterm/(fact*denterm);
        sum     += term;
      } while (std::abs(term) > eps*std::abs(sum));

      C = x*sum;

      twofn   = 1.0;
      fact    = 1.0;
      denterm = 3.0;
      numterm = 1.0;
      sum     = 1.0/3.0;
      do {
        twofn   += 2.0;
        fact    *= twofn*(twofn-1.0);
        denterm += 4.0;
        numterm *= t;
        term     = numterm/(fact*denterm);
        sum     += term;
      } while (std::abs(term) > eps*std::abs(sum));

      S = (pi/2)*sum*(x*x*x);

    } else if (x < 6.0) {

      double sumn = 0.0;
      double sumd = fd[11];
      for (int k=10; k >= 0; --k) {
        sumn = fn[k] + x*sumn;
        sumd = fd[k] + x*sumd;
      }
      double f = sumn/sumd;

      sumn = 0.0;
      sumd = gd[11];
      for (int k=10; k >= 0; --k) {
        sumn = gn[k] + x*sumn;
        sumd = gd[k] + x*sumd;
      }
      double g = sumn/sumd;

      double U    = (pi/2)*(x*x);
      double SinU = sin(U);
      double CosU = cos(U);
      C = 0.5 + f*SinU - g*CosU;
      S = 0.5 - f*CosU - g*SinU;

    } else {

      double absterm;

      double const s = pi*x*x;
      double const t = -1/(s*s);

      double numterm = -1.0;
      double term    =  1.0;
      double sum     =  1.0;
      double oldterm =  1.0;
      double eps10   =  0.1*eps;

      do {
        numterm += 4.0;
        term    *= numterm*(numterm-2.0)*t;
        sum     += term;
        absterm  = std::abs(term);
        oldterm  = absterm;
      } while (absterm > eps10*std::abs(sum));

      double f = sum/(pi*x);

      numterm = -1.0;
      term    =  1.0;
      sum     =  1.0;
      oldterm =  1.0;

      do {
        numterm += 4.0;
        term    *= numterm*(numterm+2.0)*t;
        sum     += term;
        absterm  = std::abs(term);

        oldterm  = absterm;
      } while (absterm > eps10*std::abs(sum));

      double g = pi*x; g = sum/(g*g*x);

      double U    = (pi/2)*(x*x);
      double SinU = sin(U);
      double CosU = cos(U);
      C = 0.5 + f*SinU - g*CosU;
      S = 0.5 - f*CosU - g*SinU;

    }
    if (y < 0) { C = -C; S = -S; }
  }

  void FresnelCS(int nk, double t, double C[], double S[]) {
    FresnelCS(t,C[0],S[0]);
    if (nk > 1) {
      double tt = (pi/2)*(t*t);
      double ss = sin(tt);
      double cc = cos(tt);
      C[1] = ss*(1/pi);
      S[1] = (1-cc)*(1/pi);
      if (nk > 2) {
        C[2] = (t*ss-S[0])*(1/pi);
        S[2] = (C[0]-t*cc)*(1/pi);
      }
    }
  }

  static void evalXYaLarge(double   a, double   b, double &X, double &Y) {
    double s    = a > 0 ? +1 : -1;
    double absa = std::abs(a);
    double z    = (1/sqrt(pi))*sqrt(absa);
    double ell  = s*b*(1/sqrt(pi))/sqrt(absa);
    double g    = -0.5*s*(b*b)/absa;
    double cg   = cos(g)/z;
    double sg   = sin(g)/z;

    double Cl, Sl, Cz, Sz;
    FresnelCS(ell,   Cl, Sl);
    FresnelCS(ell+z, Cz, Sz);

    double dC0 = Cz - Cl;
    double dS0 = Sz - Sl;

    X = cg*dC0 - s*sg*dS0;
    Y = sg*dC0 + s*cg*dS0;
  }

  static void evalXYaLarge(int nk, double a, double b, double X[], double Y[]) {
    double s    = a > 0 ? +1 : -1;
    double absa = std::abs(a);
    double z    = (1/sqrt(pi))*sqrt(absa);
    double ell  = s*b*(1/sqrt(pi))/sqrt(absa);
    double g    = -0.5*s*(b*b)/absa;
    double cg   = cos(g)/z;
    double sg   = sin(g)/z;


    double Cl[nk], Sl[nk], Cz[nk], Sz[nk];


    FresnelCS(nk, ell,   Cl, Sl);
    FresnelCS(nk, ell+z, Cz, Sz);

    double dC0 = Cz[0] - Cl[0];
    double dS0 = Sz[0] - Sl[0];
    X[0] = cg*dC0 - s*sg*dS0;
    Y[0] = sg*dC0 + s*cg*dS0;
    if (nk > 1) {
      cg /= z;
      sg /= z;
      double dC1 = Cz[1] - Cl[1];
      double dS1 = Sz[1] - Sl[1];
      double DC  = dC1-ell*dC0;
      double DS  = dS1-ell*dS0;
      X[1] = cg*DC - s*sg*DS;
      Y[1] = sg*DC + s*cg*DS;
      if (nk > 2) {
        double dC2 = Cz[2] - Cl[2];
        double dS2 = Sz[2] - Sl[2];
        DC   = dC2+ell*(ell*dC0-2*dC1);
        DS   = dS2+ell*(ell*dS0-2*dS1);
        cg   = cg/z;
        sg   = sg/z;
        X[2] = cg*DC - s*sg*DS;
        Y[2] = sg*DC + s*cg*DS;
      }
    }
  }

  static double LommelReduced(double mu, double nu, double b) {
    double tmp = 1/((mu+nu+1)*(mu-nu+1));
    double res = tmp;
    for (int n = 1; n < 100; ++n) {
      tmp *= (-b/(2*n+mu-nu+1))*(b/(2*n+mu+nu+1));
      res += tmp;
      if (std::abs(tmp) < std::abs(res)*1e-50) break;
    }
    return res;
  }

  static void evalXYazero(int nk, double b, double X[], double Y[]) {
    double sb = sin(b);
    double cb = cos(b);
    double b2 = b*b;
    if (std::abs(b) < 1e-3) {
      X[0] = 1-(b2/6)*(1-(b2/20)*(1-(b2/42)));
      Y[0] = (b/2)*(1-(b2/12)*(1-(b2/30)));
    } else {
      X[0] = sb/b;
      Y[0] = (1-cb)/b;
    }
    double A = b*sb;
    double D = sb-b*cb;
    double B = b*D;
    double C = -b2*sb;
    for (int k = 1; k < nk; ++k) {
      double t1 = LommelReduced(k+0.5,1.5,b);
      double t2 = LommelReduced(k+1.5,0.5,b);
      double t3 = LommelReduced(k+1.5,1.5,b);
      double t4 = LommelReduced(k+0.5,0.5,b);
      X[k] = (k*A*t1 + B*t2 + cb)/(1+k);
      Y[k] = (C*t3 + sb)/(2+k) + D*t4;
    }
  }

  static void evalXYaSmall(double a, double b, int p, double &X,  double &Y) {

    int nkk = 4*p + 3;

    double X0[nkk], Y0[nkk];

    evalXYazero(nkk, b, X0, Y0);

    X = X0[0]-(a/2)*Y0[2];
    Y = Y0[0]+(a/2)*X0[2];

    double t  = 1;
    double aa = -a*a/8 ;
    for (int n=1; n <= p; ++n) {
      t *= aa/(n*(2*n-1));
      double bf = a/(4*n+2);
      int jj = 4*n;
      X += t*(X0[jj]-bf*Y0[jj+2]);
      Y += t*(Y0[jj]+bf*X0[jj+2]);
    }
  }

  static void evalXYaSmall(int nk, double a, double b, int p, double X[], double Y[]) {

    int nkk = nk + 4*p + 2;

    double X0[nkk], Y0[nkk];

    evalXYazero(nkk, b, X0, Y0);

    for (int j=0; j < nk; ++j) {
      X[j] = X0[j]-(a/2)*Y0[j+2];
      Y[j] = Y0[j]+(a/2)*X0[j+2];
    }

    double t  = 1;
    double aa = -a*a/8 ;
    for (int n=1; n <= p; ++n) {
      t *= aa/(n*(2*n-1));
      double bf = a/(4*n+2);
      for (int j = 0; j < nk; ++j) {
        int jj = 4*n+j;
        X[j] += t*(X0[jj]-bf*Y0[jj+2]);
        Y[j] += t*(Y0[jj]+bf*X0[jj+2]);
      }
    }
  }

  void GeneralizedFresnelCS(double a, double b, double c, double &intC, double &intS) {

    double xx, yy;
    if (std::abs(a) < 0.01){
          evalXYaSmall(a, b, 5, xx, yy);
    }
    else {
      evalXYaLarge(a, b, xx, yy);
    }

    double cosc = cos(c);
    double sinc = sin(c);

    intC = xx*cosc - yy*sinc;
    intS = xx*sinc + yy*cosc;
  }

  void GeneralizedFresnelCS(int nk, double a, double b, double c, double intC[], double intS[]) {


    if (std::abs(a) < 0.01){
          evalXYaSmall(nk, a, b, 5, intC, intS);
    }
    else{
          evalXYaLarge(nk, a, b, intC, intS);
    }

    double cosc = cos(c);
    double sinc = sin(c);

    for (int k = 0; k < nk; ++k) {
      double xx = intC[k];
      double yy = intS[k];
      intC[k] = xx*cosc - yy*sinc;
      intS[k] = xx*sinc + yy*cosc;
    }
  }

  static double const CF[] = { 2.989696028701907,  0.716228953608281,
                                 -0.458969738821509, -0.502821153340377,
                                  0.261062141752652, -0.045854475238709 };

  int buildClothoid(double x0, double y0, double theta0, double x1, double y1, double theta1, double &k, double &dk, double &L) {

    double dx  = x1 - x0;
    double dy  = y1 - y0;
    double r   = hypot(dx, dy);
    double phi = atan2(dy, dx);

    double phi0 = theta0 - phi;
    double phi1 = theta1 - phi;

    while (phi0 >  pi) phi0 -= (2*pi);
    while (phi0 < -pi) phi0 += (2*pi);
    while (phi1 >  pi) phi1 -= (2*pi);
    while (phi1 < -pi) phi1 += (2*pi);

    double delta = phi1 - phi0;


    double X  = phi0*(1/pi);
    double Y  = phi1*(1/pi);
    double xy = X*Y;
    Y *= Y; X *= X;
    double A  = (phi0+phi1)*(CF[0]+xy*(CF[1]+xy*CF[2])+(CF[3]+xy*CF[4])*(X+Y)+CF[5]*(X*X+Y*Y));


    double g=0, dg, intC[3], intS[3];
    int niter = 0;
    do {
      GeneralizedFresnelCS(3, 2*A, delta-A, phi0, intC, intS);
      g   = intS[0];
      dg  = intC[2] - intC[1];
      A  -= g/dg;
    } while (++niter <= 10 && std::abs(g) > 1E-12);

    L = r/intC[0];


    k  = (delta-A)/L;
    dk = 2*A/L/L;

    return niter;
  }


  int pointsOnClothoid(double x0, double y0, double theta0, double kappa, double dkappa, double L, uint npts, std::vector<double> &X,std::vector<double> &Y) {

    X.resize(npts);
    Y.resize(npts);

    uint k = 0;
    for (auto t : linspace(0.0, L, npts)) {
      double C[1], S[1];
      GeneralizedFresnelCS(1, dkappa*pow(t,2), kappa*t, theta0, C, S);
      X[k] = x0 + t*C[0];
      Y[k] = y0 + t*S[0];
      k++;
    }

    return 0;
  }

}
