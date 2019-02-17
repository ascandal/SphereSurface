#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include <string.h>
using namespace std;

class Parameterization {
public:
  Parameterization(double radius);
  virtual void PrintPlot3DHeader(ofstream& outfile);
  virtual void PrintPlot3D(ofstream& outfile);
  string GetName();
  virtual double coordinateX(int J, int K, int L, int G);
  virtual double coordinateY(int J, int K, int L, int G);
  virtual double coordinateZ(int J, int K, int L, int G);

  string name;

  double radius;

  // Indices to loop over x, y, z directions and the grid partitions.
  int J, K, L, G;
  int Jmax, Kmax, Lmax, Gmax;

  // Coordinate coordinates.
  double x, y, z;
  double x_lb, x_ub;
  double y_lb, y_ub;
  double z_lb, z_ub;
  double dx, dy, dz;
};


class Rectangular : public Parameterization {
public:
  Rectangular(double radius);
  double coordinateX(int J, int K, int L, int G);
  double coordinateY(int J, int K, int L, int G);
  double coordinateZ(int J, int K, int L, int G);
};


class Spherical : public Parameterization {
public:
  Spherical(double radius);
  double coordinateX(int J, int K, int L, int G);
  double coordinateY(int J, int K, int L, int G);
  double coordinateZ(int J, int K, int L, int G);

  // Angular coordinates.
  double theta1, theta2;
  double theta1_lb, theta1_ub;
  double theta2_lb, theta2_ub;
  double d_theta1, d_theta2;
};


class BoxProjection : public Parameterization {
public:
  BoxProjection(double radius);
  double coordinateX(int J, int K, int L, int G);
  double coordinateY(int J, int K, int L, int G);
  double coordinateZ(int J, int K, int L, int G);
  double projection(double x, double y, double z);
};


class YinYang : public Parameterization {
public:
  YinYang(double radius);
  double coordinateX(int J, int K, int L, int G);
  double coordinateY(int J, int K, int L, int G);
  double coordinateZ(int J, int K, int L, int G);

  // Angular coordinates.
  double theta1, theta2;
  double theta1_lb, theta1_ub;
  double theta2_lb, theta2_ub;
  double d_theta1, d_theta2;
};


class TwoCaps : public Parameterization {
public:
  TwoCaps(double radius);
  double coordinateX(int J, int K, int L, int G);
  double coordinateY(int J, int K, int L, int G);
  double coordinateZ(int J, int K, int L, int G);
  void PrintPlot3DHeader(ofstream& outfile);
  void GridIndexRanges(ofstream& outfile);

  // support variables.
  double rad_B, proj;

  // Angular coordinates.
  double theta1, theta2;
  double theta1_lb, theta1_ub;
  double theta2_lb, theta2_ub;
  double d_theta1, d_theta2;

  int Jmax_cap, Kmax_cap;
  double theta1_cap_lb, theta1_cap_ub;
  double d_theta1_cap, d_theta2_cap;
  double theta2_cap_lb, theta2_cap_ub;
};

#endif
