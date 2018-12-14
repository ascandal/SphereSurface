# SurfaceSphere

SurfaceSphere generates the structured surface mesh of a sphere in formatted PLOT3D format (file extension .xyz).

It is designed to be used for high-fidelity computational solvers that need precise control over the surface mesh properties.

* The sphere is centered about the origin (0, 0, 0).
* There are multiple parameterizations to choose from.
* Some of the parameterizations have overlapping mesh parts to be used with overset grid schemes. Mesh joining routines are used as a post process to intersect the mesh parts where they overlap along a seam. In this case for a sphere, final mesh will be water-tight (i.e. have no holes). The joining methods are not included here.

Some surface mesh properies include:

* resolution total number of cells
* ratio of smallest to largest cells
* regularity of cell shape
* presence of singularities

## PLOT3D Format

The PLOT3D formatted file (extension .xyz) is described by its read method.

Structered grids

3D, Whole, formatted, Multi-Block Grid

A link descriping the PLOT3D format can be found here:

https://www.grc.nasa.gov/www/wind/valid/plot3d.html


Written in **Fortran 77** for **3D, Whole, formatted, Multi-Block Grid**

NGRID {int} : the number of grids
X(J,K,L) {double} : the x-coordinate at index (J,K,L)
Y(J,K,L) {double} : the y-coordinate at index (J,K,L)
Z(J,K,L) {double} : the z-coordinate at index (J,K,L)
IBLANK(J,K,L)
JD(IG) {int} : the number of J-parameter indexes for the IG'th grid
KD(IG) {int} : the number of K-parameter indexes for the IG'th grid
LD(IG) {int} : the number of L-parameter indexes for the IG'th grid

```fortran
READ(1) NGRID
READ(1) (JD(IG),KD(IG),LD(IG),IG=1,NGRID)
DO IG = 1,NGRID
  READ(1) (((X(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG)),
  &       (((Y(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG)),
  &       (((Z(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG)),
  &       (((IBLANK(J,K,L),J=1,JD(IG)),K=1,KD(IG)),L=1,LD(IG))
ENDDO
 ```

## Installing SurfaceSphere

gcc -o surface-sphere surface-sphere.c

## Running SurfaceSphere

./surface-sphere -r double -p int

-r {double} the radius of the sphere.
-p {int} the parameterization of the mesh.
      (0) RECTANGULAR
      (1) SPHERICAL
      (2) BOX_PROJECTION
      (3) YIN_YANG
      (4) TWO_CAPS


## SurfaceSphere Parameterizations

* [Rectangular](#rectangular)
* [Spherical](#spherical)
* [Box Projection](#boxprojection)
* [Yin Yang](#yinyang)
* [Two Caps](#twocaps)

### Rectangular

Number of Mesh Parts : 1

### Spherical

Number of Mesh Parts : 1

### Box Projection

Number of Mesh Parts : 6

### Yin Yang

The Yin-Yang mesh avoids singularities at the poles. It has had particular use in atmospherical computations.

Number of Mesh Parts : 2

### Two Caps

The Two Caps mesh avoids singularities at the poles.

Number of Mesh Parts : 3
