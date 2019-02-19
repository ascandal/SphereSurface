#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;


Spherical::Spherical(double radius) : Parameterization(radius) {
  name = "Spherical";

  Jmax = 200;
  Kmax = 100;
  Lmax = 1;
  Gmax = 1;

  theta1_lb =  0.0e0;
  theta1_ub = -2.0e0 * PI;

  d_theta1 = (theta1_ub - theta1_lb) / (double)(Jmax - 1);

  theta2_lb = PI;
  theta2_ub = 0.0e0;

  d_theta2 = (theta2_ub - theta2_lb) / (double)(Kmax - 1);
}

double Spherical::coordinateX(int J, int K, int L, int G) {
  theta1 = theta1_lb + (double)J * d_theta1;
  theta2 = theta2_lb + (double)K * d_theta2;

  return radius * sin(theta2) * sin(theta1);
}

double Spherical::coordinateY(int J, int K, int L, int G) {
  theta1 = theta1_lb + (double)J * d_theta1;
  theta2 = theta2_lb + (double)K * d_theta2;

  return radius * sin(theta2) * cos(theta1);
}

double Spherical::coordinateZ(int J, int K, int L, int G) {
  theta1 = theta1_lb + (double)J * d_theta1;
  theta2 = theta2_lb + (double)K * d_theta2;

  return radius * cos(theta2);
}
