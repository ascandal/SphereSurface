#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;

Parameterization::Parameterization() {}

string Parameterization::GetName() {
  return name;
}

/*
 * Print the PLOT3D header for the surface parameterization
 * to the output file.
 * READ(1) NGRID
 * READ(1) (JD(IG),KD(IG),LD(IG),IG=1,NGRID)
 */
void Parameterization::PrintPlot3DHeader(ofstream& outfile) {
  outfile << Gmax << endl;
  for (G = 0; G < Gmax; G++) {
    outfile << Jmax << " " << Kmax << " " << Lmax << endl;
  }
}

/*
 * Print the entire PLOT3D formatted surface parameterization to the output file.
 * The PLOT3D format (file extension .xyz) is described by its read method:
 * READ(1) NGRID
 * READ(1) (JD(IG),KD(IG),LD(IG),IG=1,NGRID)
 * DO IG = 1,NGRID
 *   READ(1) (((X(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG)),
 *   &       (((Y(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG)),
 *   &       (((Z(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG)),
 *   &       (((IBLANK(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG))
 * ENDDO
 */
void Parameterization::PrintPlot3D(ofstream& outfile) {

  PrintPlot3DHeader(outfile);

  for (G = 0; G < Gmax; G++) {
    // All x-coordinates
    for (L = 0; L < Lmax; L++) {
    for (K = 0; K < Kmax; K++) {
    for (J = 0; J < Jmax; J++) {
      outfile << coordinateX(J, K, L, G) << " ";
    }
    }
    outfile << endl;
    }

    // All y-coordinates
    for (L = 0; L < Lmax; L++) {
    for (K = 0; K < Kmax; K++) {
    for (J = 0; J < Jmax; J++) {
      outfile << coordinateY(J, K, L, G) << " ";
    }
    }
    outfile << endl;
    }

    // All z-coordinates
    for (L = 0; L < Lmax; L++) {
    for (K = 0; K < Kmax; K++) {
    for (J = 0; J < Jmax; J++) {
      outfile << coordinateZ(J, K, L, G) << " ";
    }
    }
    outfile << endl;
    }

  }  //End: G for-loop
}

double Parameterization::coordinateX(int J, int K, int L, int G) {
  return 0;
}

double Parameterization::coordinateY(int J, int K, int L, int G) {
  return 0;
}

double Parameterization::coordinateZ(int J, int K, int L, int G) {
  return 0;
}




Rectangular::Rectangular() {
  name = "RECTANGULAR";

  Jmax = 200;
  Kmax = 100;
  Lmax = 1;
  Gmax = 1;

  x_lb = -radius + 0.00001e0;
  x_ub =  radius - 0.00001e0;

  dx = (x_ub - x_lb) / (double)(Jmax - 1);
}

double Rectangular::coordinateX(int J, int K, int L, int G) {
  return x_lb + (double)J * dx;
}

double Rectangular::coordinateY(int J, int K, int L, int G) {
  x = x_lb + (double)J * dx;

  y_lb = -sqrt(radius * radius - x * x);
  y_ub =  sqrt(radius * radius - x * x);

  dy = (y_ub - y_lb) / (double)(Kmax - 1);

  return y_lb + (double)K * dy;
}

double Rectangular::coordinateZ(int J, int K, int L, int G) {
  x = x_lb + (double)J * dx;

  y_lb = -sqrt(radius * radius - x * x);
  y_ub =  sqrt(radius * radius - x * x);

  dy = (y_ub - y_lb) / (double)(Kmax - 1);

  y = y_lb + (double)K * dy;

  //z_lb =  0.0e0;
  //z_ub =  sqrt(radius * radius - x * x - y * y);

  //dz = (z_ub - z_lb) / (double)(Lmax - 1);

  return sqrt(max(0.0e0, radius * radius - x * x - y * y));
}


Spherical::Spherical() {
  name = "SPHERICAL";

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



BoxProjection::BoxProjection() {
  name = "BOX_PROJECTION";

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

  rad_B = sqrt(x * x + y * y + z * z);

  proj = radius / rad_B;

  x *= proj;

  return x;
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

  rad_B = sqrt(x * x + y * y + z * z);

  proj = radius / rad_B;

  y *= proj;

  return y;
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

  rad_B = sqrt(x * x + y * y + z * z);

  proj = radius / rad_B;

  z *= proj;

  return z;
}


YinYang::YinYang() {
  name = "YIN_YANG";

  Jmax = 120;
  Kmax = 60;
  Lmax = 1;
  Gmax = 2;

  // Angle in the xy-plane
  theta1_lb = (1.0e0 / 4.0e0) * PI;
  theta1_ub = (7.0e0 / 4.0e0) * PI;
  //theta1_lb = -(3.0e0 / 4.0e0) * PI;
  //theta1_ub = (3.0e0 / 4.0e0) * PI;

  d_theta1 = (theta1_ub - theta1_lb) / (double)(Jmax - 1);

  // Angle between positive z-axis and radius.
  //theta2_lb = (1.0e0 / 4.0e0) * PI;
  //theta2_ub = (3.0e0 / 4.0e0) * PI;
  theta2_lb = (1.0e0 / 5.0e0) * PI;
  theta2_ub = (4.0e0 / 5.0e0) * PI;

  d_theta2 = (theta2_ub - theta2_lb) / (double)(Kmax - 1);
}

double YinYang::coordinateX(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Yin
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      x = radius * sin(theta2) * sin(theta1);
      break;
    case 1:  // Yang
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      //x = radius * sin(theta2) * sin(theta1);
      x = radius * cos(theta2);
      break;
  }

  return x;
}

double YinYang::coordinateY(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Yin
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      y = radius * sin(theta2) * cos(theta1);
      break;
    case 1:  // Yang
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      //y = radius * cos(theta2);
      y = -radius * sin(theta2) * cos(theta1);
      break;
  }

  return y;
}

double YinYang::coordinateZ(int J, int K, int L, int G) {
  switch (G) {
    case 0:  // Yin
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      z = radius * cos(theta2);
      break;
    case 1:  // Yang
      theta1 = theta1_lb + (double)J * d_theta1;
      theta2 = theta2_lb + (double)K * d_theta2;
      //z = radius * sin(theta2) * cos(theta1);
      z = radius * sin(theta2) * sin(theta1);
      break;
  }

  return z;
}



void TwoCaps::GridIndexRanges(ofstream& outfile) {
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

TwoCaps::TwoCaps() {
  name = "TWO_CAPS";

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

  printf("\n( Jmax , Kmax ) = ( %i , %i )\n", Jmax, Kmax);

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