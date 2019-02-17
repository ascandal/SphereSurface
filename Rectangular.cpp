#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;

Rectangular::Rectangular(double radius) : Parameterization(radius) {
  name = "Rectangular";

  Jmax = 200;
  Kmax = 100;
  Lmax = 1;
  Gmax = 2;

  x_lb = -radius;
  x_ub =  radius;

  dx = (x_ub - x_lb) / (double)(Jmax - 1);
}

double Rectangular::coordinateX(int J, int K, int L, int G) {
  return x_lb + (double)J * dx;
}

double Rectangular::coordinateY(int J, int K, int L, int G) {
  x = x_lb + (double)J * dx;

  y_lb = -sqrt(radius * radius - x * x);
  y_ub = -y_lb;

  dy = (y_ub - y_lb) / (double)(Kmax - 1);

  return y_lb + (double)K * dy;
}

double Rectangular::coordinateZ(int J, int K, int L, int G) {
  x = x_lb + (double)J * dx;

  y_lb = -sqrt(radius * radius - x * x);
  y_ub =  sqrt(radius * radius - x * x);

  dy = (y_ub - y_lb) / (double)(Kmax - 1);

  y = y_lb + (double)K * dy;

  z_lb = -sqrt(max(0.0e0, radius * radius - x * x - y * y));
  z_ub = -z_lb;

  switch (G) {
    case 0:  // Lower z
      z = z_lb;
      break;
    case 1:  // Upper z
      z = z_ub;
      break;
  }

  return z;
}
