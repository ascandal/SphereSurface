#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;


YinYang::YinYang(double radius) : Parameterization(radius) {
  name = "YinYang";

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
