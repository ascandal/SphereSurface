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
  for (ii = 0; ii < (jd); ii++                                   \
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
  if (Coordinates == RECTANGULAR) {
    printf("Parameterization : RECTANGULAR\n");
  } else if (Coordinates == SPHERICAL) {
    printf("Parameterization : SPHERICAL\n");
  } else if (Coordinates == BOX_PROJECTION) {
    printf("Parameterization : BOX_PROJECTION\n");
  } else if (Coordinates == YIN_YANG) {
    printf("Parameterization : YIN_YANG\n");
  } else if (Coordinates == TWO_CAPS) {
    printf("Parameterization : TWO_CAPS\n");
  }
  printf("Sphere radius : %e\n", radius);

  //------------------------+-------------------------------------------------------------
  // Open File to write     |
  //------------------------+
  FILE *outfile;

  outfile = fopen("SphereSurf.fmt", "w");


  // Indices to loop over x, y, and z directions.
  int J, K, L;

  // Maximum values of J, K, L indices.
  int Jmax, Kmax, Lmax;

  int G, Gmax;

  // Coordinate dimensions

  double x, y, z;

  double x_lb, x_ub,
         y_lb, y_ub,
         z_lb, z_ub;

  double dx, dy, dz;


  double theta1, theta2;

  double theta1_lb, theta1_ub;
  double theta2_lb, theta2_ub;
  double d_theta1, d_theta2;

  //---------------------+-----------------------------------------
  // rectangular         |
  //---------------------+
  Jmax = 200;
  Kmax = 100;
  Lmax = 1;
  Gmax = 1;

  x_lb = -radius + 0.00001e0;
  x_ub =  radius - 0.00001e0;

  dx = (x_ub - x_lb) / (double)(Jmax - 1);
  //---------------------+
  // rectangular         |
  //---------------------+-----------------------------------------

  //---------------------+-----------------------------------------
  // spherical           |
  //---------------------+
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
  //---------------------+
  // spherical           |
  //---------------------+-----------------------------------------


  //---------------------+-----------------------------------------
  // BOX_PROJECTION      |
  //---------------------+
  double rad_B, proj;

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
  //---------------------+
  // BOX_PROJECTION      |
  //---------------------+-----------------------------------------


  //---------------------+-----------------------------------------
  // YIN_YANG            |
  //---------------------+
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
  //---------------------+
  // YIN_YANG            |
  //---------------------+-----------------------------------------

  //---------------------+-----------------------------------------
  // TWO_CAPS            |
  //---------------------+
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
  double theta1_cap_lb = -theta2_lb - cap_overlap;
  double theta1_cap_ub =  theta2_lb + cap_overlap;

  int Jmax_cap = (int)((theta1_cap_ub - theta1_cap_lb) / d_theta1) + 1;

  double d_theta1_cap = (theta1_cap_ub - theta1_cap_lb) / (double)(Jmax_cap - 1);

  // Angle between positive z-axis and radius.
  double theta2_cap_lb =  PI / 2.0e0 - theta2_lb - cap_overlap;
  double theta2_cap_ub =  PI / 2.0e0 + theta2_lb + cap_overlap;

  int Kmax_cap = (int)((theta1_cap_ub - theta1_cap_lb) / d_theta1) + 1;

  double d_theta2_cap = (theta2_cap_ub - theta2_cap_lb) / (double)(Kmax_cap - 1);

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

  // FOR THE CAPS
  //-------------------------------------------
  //---------------------+
  // TWO_CAPS            |
  //---------------------+-----------------------------------------


  //  READ(1) NGRID
  fprintf(outfile, "%i\n", Gmax);

  //  READ(1) (JD(IG),KD(IG),LD(IG),IG=1,NGRID)
  if (Coordinates == TWO_CAPS) {
    fprintf(outfile, "%i %i %i\n", Jmax, Kmax, Lmax);
    fprintf(outfile, "%i %i %i\n", Jmax_cap, Kmax_cap, Lmax);
    fprintf(outfile, "%i %i %i\n", Jmax_cap, Kmax_cap, Lmax);
  } else {
    for (G = 0; G < Gmax; G++) {
      fprintf(outfile, "%i %i %i\n", Jmax, Kmax, Lmax);
    }
  }


  for (G = 0; G < Gmax; G++) {
    //---------------------------------------+-------------------------
    // Choose Index Ranges for each grid.    |
    //---------------------------------------+
    if (Coordinates == RECTANGULAR) {
      Jmax = Jmax;
      Kmax = Kmax;
      Lmax = Lmax;

    } else if (Coordinates == SPHERICAL) {
      Jmax = Jmax;
      Kmax = Kmax;
      Lmax = Lmax;

    } else if (Coordinates == BOX_PROJECTION) {
      Jmax = Jmax;
      Kmax = Kmax;
      Lmax = Lmax;

    } else if (Coordinates == YIN_YANG) {
      Jmax = Jmax;
      Kmax = Kmax;
      Lmax = Lmax;

    } else if (Coordinates == TWO_CAPS) {

      switch (G) {
        case 0: // Low latitude spherical shell
          Jmax = Jmax;
          Kmax = Kmax;
          Lmax = Lmax;
          break;
        case 1: // Top cap
          Jmax = Jmax_cap;
          Kmax = Kmax_cap;
          Lmax = Lmax;
          break;
        case 2: // Bottom cap
          Jmax = Jmax_cap;
          Kmax = Kmax_cap;
          Lmax = Lmax;
          break;
      }
    }
    //---------------------------------------+
    // Choose Index Ranges for each grid.    |
    //---------------------------------------+-------------------------





    // All x-coordinates
    for (L = 0; L < Lmax; L++) {
    for (K = 0; K < Kmax; K++) {
    for (J = 0; J < Jmax; J++) {

    //for (J = 0; J < Jmax; J++) {
    //for (K = 0; K < Kmax; K++) {
    //for (L = 0; L < Lmax; L++) {

      if (Coordinates == RECTANGULAR) {
        x = x_lb + (double)J * dx;
      } else if (Coordinates == SPHERICAL) {

        theta1 = theta1_lb + (double)J * d_theta1;
        theta2 = theta2_lb + (double)K * d_theta2;

        x = radius * sin(theta2) * sin(theta1);
      } else if (Coordinates == BOX_PROJECTION) {

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

      } else if (Coordinates == YIN_YANG) {

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
      } else if (Coordinates == TWO_CAPS) {

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
      }

      fprintf(outfile, "%e ", x);

    }
    }
    fprintf(outfile, "\n");
    }



    // All y-coordinates
    for (L = 0; L < Lmax; L++) {
    for (K = 0; K < Kmax; K++) {
    for (J = 0; J < Jmax; J++) {

    //for (J = 0; J < Jmax; J++) {
    //for (K = 0; K < Kmax; K++) {
    //for (L = 0; L < Lmax; L++) {

      if (Coordinates == RECTANGULAR) {
        x = x_lb + (double)J * dx;

        y_lb = -sqrt(radius * radius - x * x);
        y_ub =  sqrt(radius * radius - x * x);

        dy = (y_ub - y_lb) / (double)(Kmax - 1);

        y = y_lb + (double)K * dy;

      } else if (Coordinates == SPHERICAL) {

        theta1 = theta1_lb + (double)J * d_theta1;
        theta2 = theta2_lb + (double)K * d_theta2;

        y = radius * sin(theta2) * cos(theta1);
      } else if (Coordinates == BOX_PROJECTION) {

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
      } else if (Coordinates == YIN_YANG) {

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
      } else if (Coordinates == TWO_CAPS) {

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
      }

      fprintf(outfile, "%e ", y);

    }
    }
    fprintf(outfile, "\n");
    }




    // All z-coordinates
    for (L = 0; L < Lmax; L++) {
    for (K = 0; K < Kmax; K++) {
    for (J = 0; J < Jmax; J++) {

    //for (J = 0; J < Jmax; J++) {
    //for (K = 0; K < Kmax; K++) {
    //for (L = 0; L < Lmax; L++) {

      if (Coordinates == RECTANGULAR) {
        x = x_lb + (double)J * dx;

        y_lb = -sqrt(radius * radius - x * x);
        y_ub =  sqrt(radius * radius - x * x);

        dy = (y_ub - y_lb) / (double)(Kmax - 1);

        y = y_lb + (double)K * dy;

        //z_lb =  0.0e0;
        //z_ub =  sqrt(radius * radius - x * x - y * y);

        //dz = (z_ub - z_lb) / (double)(Lmax - 1);

        z = sqrt(max(0.0e0, radius * radius - x * x - y * y));

      } else if (Coordinates == SPHERICAL) {

        theta1 = theta1_lb + (double)J * d_theta1;
        theta2 = theta2_lb + (double)K * d_theta2;

        z = radius * cos(theta2);
      } else if (Coordinates == BOX_PROJECTION) {

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
      } else if (Coordinates == YIN_YANG) {

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
      } else if (Coordinates == TWO_CAPS) {

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
      }

      fprintf(outfile, "%e ", z);

    }
    }
    fprintf(outfile, "\n");
    }


  }  //End: G for-loop

  fclose(outfile);

  printf("\n\n>>> Surface Mesh File Complete <<<\n\n");

  return 0;
}  // End: main()
