#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;


void TwoCaps::setGridIndexRanges(int G) {
  // Choose index ranges for each grid.
  // Top and bottom caps
  if (G > 0) {
    Jmax = Jmax_cap;
    Kmax = Kmax_cap;
  }
}

/*
 * READ(1) NGRID
 * READ(1) (JD(IG),KD(IG),LD(IG),IG=1,NGRID)
 */
void TwoCaps::PrintPlot3DHeader(ofstream& outfile) {
  outfile << Gmax << endl;
  outfile << Jmax << " " << Kmax << " " << Lmax << endl;
  outfile << Jmax_cap << " " << Kmax_cap << " " << Lmax << endl;
  outfile << Jmax_cap << " " << Kmax_cap << " " << Lmax << endl;
}

TwoCaps::TwoCaps(double radius) : Parameterization(radius) {
  name = "TwoCaps";

  Jmax = 160;
  //Kmax = 36;
  Lmax = 1;
  Gmax = 3;

  // Angle in the xy-plane
  theta1_lb = 0.0e0;
  theta1_ub = 2.0e0 * PI;

  d_theta1 = (theta1_ub - theta1_lb) / (double)(Jmax - 1.0e0);

  // Angle between positive z-axis and radius.
  theta2_lb = (1.0e0 / 5.0e0) * PI;
  theta2_ub = (4.0e0 / 5.0e0) * PI;

  Kmax = (int)((theta2_ub - theta2_lb) / d_theta1) + 1;

  d_theta2 = (theta2_ub - theta2_lb) / (double)(Kmax - 1.0e0);

  //-------------------------------------------
  // FOR THE CAPS
  double cap_overlap = 2.0e0 * d_theta2;      // (1.0e0 / 15.0e0) * PI;

  // Angle in the xy-plane
  theta1_cap_lb = -theta2_lb - cap_overlap;
  theta1_cap_ub =  theta2_lb + cap_overlap;

  Jmax_cap = (int)((theta1_cap_ub - theta1_cap_lb) / d_theta1) + 1;

  d_theta1_cap = (theta1_cap_ub - theta1_cap_lb) / (double)(Jmax_cap - 1);

  // Angle between positive z-axis and radius.
  theta2_cap_lb =  PI / 2.0e0 - theta2_lb - cap_overlap;
  theta2_cap_ub =  PI / 2.0e0 + theta2_lb + cap_overlap;

  Kmax_cap = (int)((theta1_cap_ub - theta1_cap_lb) / d_theta1) + 1;

  d_theta2_cap = (theta2_cap_ub - theta2_cap_lb) / (double)(Kmax_cap - 1);

  //-----------------------------------------------------
  x_lb = -radius - 0.1e0;
  x_ub =  radius + 0.1e0;

  y_lb = -radius - 0.1e0;
  y_ub =  radius + 0.1e0;

  z_lb = -radius;   // + 0.01e0;
  z_ub =  radius;   // - 0.01e0;

  double dx_tmp = sqrt(2.0e0 * radius * radius * (1.0e0 - cos(d_theta1)));

  Kmax_cap = Jmax_cap = (int)((x_ub - x_lb) / dx_tmp) + 1;

  dx = (x_ub - x_lb) / (double)(Jmax_cap - 1);
  dy = (y_ub - y_lb) / (double)(Jmax_cap - 1);
  dz = (z_ub - z_lb) / (double)(Jmax_cap - 1);
}

double TwoCaps::coordinateX(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Low latitude spherical shell
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      x = radius * sin(theta2) * sin(theta1);
      break;
    case 1:  // Top cap
      theta1 = theta1_cap_lb + (double)J * d_theta1_cap;
      theta2 = theta2_cap_lb + (double)K * d_theta2_cap;
      x = radius * sin(theta2) * sin(theta1);

      x = x_lb + (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_ub;

      rad_B = sqrt(x * x + y * y + z * z);
      proj = radius / rad_B;
      x *= proj;
      break;
    case 2:  // Bottom cap
      theta1 = theta1_cap_lb + (double)J * d_theta1_cap;
      theta2 = theta2_cap_lb + (double)K * d_theta2_cap;
      x = radius * sin(theta2) * sin(theta1);

      x = x_ub - (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_lb;

      rad_B = sqrt(x * x + y * y + z * z);
      proj = radius / rad_B;
      x *= proj;
      break;
  }

  return x;
}

double TwoCaps::coordinateY(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Low latitude spherical shell
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      y = radius * sin(theta2) * cos(theta1);
      break;
    case 1:  // Top cap
      theta1 = theta1_cap_lb + (double)J * d_theta1_cap;
      theta2 = theta2_cap_lb + (double)K * d_theta2_cap;
      y = -radius * cos(theta2);

      x = x_lb + (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_ub;

      rad_B = sqrt(x * x + y * y + z * z);
      proj = radius / rad_B;
      y *= proj;
      break;
    case 2:  // Bottom cap
      theta1 = theta1_cap_lb + (double)J * d_theta1_cap;
      theta2 = theta2_cap_lb + (double)K * d_theta2_cap;
      y =  radius * cos(theta2);

      x = x_ub - (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_lb;

      rad_B = sqrt(x * x + y * y + z * z);
      proj = radius / rad_B;
      y *= proj;
      break;
  }

  return y;
}

double TwoCaps::coordinateZ(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Low latitude spherical shell
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      z = radius * cos(theta2);
      break;
    case 1:  // Top cap
      theta1 = theta1_cap_lb + (double)J * d_theta1_cap;
      theta2 = theta2_cap_lb + (double)K * d_theta2_cap;
      z = radius * sin(theta2) * cos(theta1);

      x = x_lb + (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_ub;

      rad_B = sqrt(x * x + y * y + z * z);
      proj = radius / rad_B;
      z *= proj;
      break;
    case 2:  // Bottom cap
      theta1 = theta1_cap_lb + (double)J * d_theta1_cap;
      theta2 = theta2_cap_lb + (double)K * d_theta2_cap;
      z = -radius * sin(theta2) * cos(theta1);

      x = x_ub - (double)J * dx;
      y = y_lb + (double)K * dy;
      z = z_lb;

      rad_B = sqrt(x * x + y * y + z * z);
      proj = radius / rad_B;
      z *= proj;
      break;
  }

  return z;
}
