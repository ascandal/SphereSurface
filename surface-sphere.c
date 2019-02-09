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
 * Usage:
 *  surface-sphere -d 0.002
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>


#define max(x,y)  (((x) < (y)) ? (y) : (x))
#define min(x,y)  (((x) < (y)) ? (x) : (y))

#define ALLOC_1d_array(type,                                      \
                       array_name,                                \
                       jd)                                        \
{                                                                 \
  array_name = (type *) malloc((jd) * sizeof(type));              \
}

#define ALLOC_2d_array(type,                                      \
                       array_name,                                \
                       jd, kd)                                    \
{                                                                 \
  int ii;                                                         \
                                                                  \
  array_name = (type **) calloc((jd), sizeof(type *));            \
                                                                  \
  for (ii = 0; ii < (jd); ii++                                    \
    (array_name[ii]) = (type *) calloc((kd), sizeof(type));       \
  }                                                               \
}

#define PI   (3.1415926535897932384626433832795028841971693993751)

// Parameterization Options
#define RECTANGULAR        (0)
#define SPHERICAL          (1)
#define BOX_PROJECTION     (2)
#define YIN_YANG           (3)
#define TWO_CAPS           (4)

/*
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
int main(int argc, char **argv)
{
  double radius = 0.5e0;
  int Coordinates = SPHERICAL;
  int c;

  // Parsing command line arguments with getopt
  while ((c = getopt(argc, argv, "r:p:")) != -1) {
    switch (c) {
      case 'r':
        radius = strtod(optarg, NULL);
        break;
      case 'p':
        Coordinates = atoi(optarg);
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  printf("Starting surface-sphere:\n");
  printf("Parameterization : %e\n", param->name);
  printf("Sphere radius : %e\n", radius);

  //------------------------+-------------------------------------------------------------
  // Open File to write     |
  //------------------------+
  FILE *outfile;
  outfile = fopen("SphereSurf.fmt", "w");

  param->PrintPlot3D(outfile);

  fclose(outfile);

  printf("\n\n>>> Surface Mesh File Complete <<<\n\n");

  return 0;
}  // End: main()
