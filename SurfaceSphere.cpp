/**
 * \file surface-sphere.c
 * \author Angelo L. Scandaliato
 * \brief Outputs the surface mesh of a sphere in PLOT3D format (file extension .xyz).
 * The sphere is centered about the origin (0, 0, 0).
 * There are multiple parameterizations to choose from. Some with overlapping
 * mesh parts to be used with other overset grid methods.
 *
 * - RECTANGULAR
 * - SPHERICAL
 * - BOX_PROJECTION
 * - YIN_YANG
 * - TWO_CAPS
 *
 * To Compile Run:
 *  g++ -std=c++11 SurfaceSphere.cpp Parameterization.cpp -o SurfaceSphere
 *
 * Usage:
 *  ./SurfaceSphere -r 0.002 -p 3
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "Parameterization.h"

using namespace std;

// Parameterization Options
#define RECTANGULAR        (0)
#define SPHERICAL          (1)
#define BOX_PROJECTION     (2)
#define YIN_YANG           (3)
#define TWO_CAPS           (4)

int main(int argc, char **argv)
{
  Parameterization* param;
  double radius = 0.5e0;
  int coordinates = SPHERICAL;
  int c;

  // Parsing command line arguments with getopt
  while ((c = getopt(argc, argv, "r:p:")) != -1) {
    switch (c) {
      case 'r':
        radius = strtod(optarg, NULL);
        break;
      case 'p':
        coordinates = atoi(optarg);
        break;
      case '?':
        if (optopt == 'c')
          cout << "Option -" << optopt << " requires an argument." << endl;
        else if (isprint (optopt))
          cout << "Unknown option -" << optopt << endl;
        else
          cout << "Unknown option character " << optopt << endl;
        return 1;
      default:
        abort();
    }
  }

  switch (coordinates) {
    case RECTANGULAR:
      param = new Rectangular(radius);
      break;
    case SPHERICAL:
      param = new Spherical(radius);
      break;
    case BOX_PROJECTION:
      param = new BoxProjection(radius);
      break;
    case YIN_YANG:
      param = new YinYang(radius);
      break;
    case TWO_CAPS:
      param = new TwoCaps(radius);
      break;
  }

  // Initial log to console.
  cout << ">>> Starting SurfaceSphere <<<\n";
  cout << "Sphere Radius = " << param->radius << endl;
  cout << "Using Parameterization = " << param->GetName() << endl;


  // Open file to write.
  stringstream ss;
  ss << fixed << setprecision(2) << param->radius;
  string fileName = "SphereSurf-" + param->GetName() + "-r" + ss.str() + ".xyz";
  ofstream outfile;
  outfile.open(fileName);

  param->PrintPlot3D(outfile);

  outfile.close();

  cout << "\n>>> SurfaceSphere Mesh File Complete <<<\n";

  return 0;
}
