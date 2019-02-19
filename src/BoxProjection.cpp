#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;


BoxProjection::BoxProjection(double radius) : Parameterization(radius) {
  name = "BoxProjection";

  Jmax = 20;
  Kmax = 20;
  Lmax = 1;
  Gmax = 6;

  x_lb = -radius + 0.00001e0;
  x_ub =  radius - 0.00001e0;

  y_lb = -radius + 0.00001e0;
  y_ub =  radius - 0.00001e0;

  z_lb = -radius + 0.00001e0;
  z_ub =  radius - 0.00001e0;

  dx = (x_ub - x_lb) / (double)(Jmax - 1);
  dy = (y_ub - y_lb) / (double)(Jmax - 1);
  dz = (z_ub - z_lb) / (double)(Jmax - 1);
}

double BoxProjection::coordinateX(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Lower z
      x = x_lb + (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_lb;
      break;
    case 1:  // Lower y
      x = x_ub - (double)J * dx;
      y = y_lb;
      z = z_lb + (double)K * dz;
      break;
    case 2:  // Lower x
      x = x_lb;
      y = y_lb + (double)J * dy;
      z = z_lb + (double)K * dz;
      break;
    case 3:  // Upper z
      x = x_ub - (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_ub;
      break;
    case 4:  // Upper y
      x = x_lb + (double)J * dx;
      y = y_ub;
      z = z_lb + (double)K * dz;
      break;
    case 5:  // Upper x
      x = x_ub;
      y = y_ub - (double)J * dy;
      z = z_lb + (double)K * dz;
      break;
  }

  return x * projection(x, y, z);
}

double BoxProjection::coordinateY(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Lower z
      x = x_lb + (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_lb;
      break;
    case 1:  // Lower y
      x = x_ub - (double)J * dx;
      y = y_lb;
      z = z_lb + (double)K * dz;
      break;
    case 2:  // Lower x
      x = x_lb;
      y = y_lb + (double)J * dy;
      z = z_lb + (double)K * dz;
      break;
    case 3:  // Upper z
      x = x_ub - (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_ub;
      break;
    case 4:  // Upper y
      x = x_lb + (double)J * dx;
      y = y_ub;
      z = z_lb + (double)K * dz;
      break;
    case 5:  // Upper x
      x = x_ub;
      y = y_ub - (double)J * dy;
      z = z_lb + (double)K * dz;
      break;
  }

  return y * projection(x, y, z);
}

double BoxProjection::coordinateZ(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Lower z
      x = x_lb + (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_lb;
      break;
    case 1:  // Lower y
      x = x_ub - (double)J * dx;
      y = y_lb;
      z = z_lb + (double)K * dz;
      break;
    case 2:  // Lower x
      x = x_lb;
      y = y_lb + (double)J * dy;
      z = z_lb + (double)K * dz;
      break;
    case 3:  // Upper z
      x = x_ub - (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_ub;
      break;
    case 4:  // Upper y
      x = x_lb + (double)J * dx;
      y = y_ub;
      z = z_lb + (double)K * dz;
      break;
    case 5:  // Upper x
      x = x_ub;
      y = y_ub - (double)J * dy;
      z = z_lb + (double)K * dz;
      break;
  }

  return z * projection(x, y, z);
}

double BoxProjection::projection(double x, double y, double z) {
  return radius / sqrt(x * x + y * y + z * z);
}
