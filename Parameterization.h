#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include <string.h>
using namespace std;

class Parameterization {
public:
  Parameterization();
  void PrintPlot3DHeader(ofstream& outfile);
  void PrintPlot3D(ofstream& outfile);
  string GetName();

private:
  string name;

  // Indices to loop over x, y, z directions and the grid partitions.
  int J, K, L, G;
  int Jmax, Kmax, Lmax, Gmax;

  // Coordinate coordinates.
  double x, y, z;
  double x_lb, x_ub;
  double y_lb, y_ub;
  double z_lb, z_ub;
  double dx, dy, dz;

  // Angular coordinates.
  double theta1, theta2;
  double theta1_lb, theta1_ub;
  double theta2_lb, theta2_ub;
  double d_theta1, d_theta2;
};


class Rectangular : public Parameterization {
public:
  Rectangular();

private:
};


class Spherical : public Parameterization {
public:
  Spherical();

private:
};


class BoxProjection : public Parameterization {
public:
  BoxProjection();

private:
};


class YinYang : public Parameterization {
public:
  YinYang();

private:
};


class TwoCaps : public Parameterization {
public:
  TwoCaps();
  void PrintPlot3DHeader(ofstream& outfile);

private:
};

#endif
