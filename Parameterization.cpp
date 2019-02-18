#include <math.h>
#include <fstream>
#include <string.h>

#include "Parameterization.h"
#include "GlobalMacros.h"

using namespace std;

Parameterization::Parameterization(double radius) : radius(radius) {}

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

    setGridIndexRanges(G);

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
