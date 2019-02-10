#ifndef GLOBAL_MACROS_H
#define GLOBAL_MACROS_H

#define max(x,y)  (((x) < (y)) ? (y) : (x))
#define min(x,y)  (((x) < (y)) ? (x) : (y))

#define ALLOC_1d_array(type, array_name, jd) {                    \
  array_name = (type *) malloc((jd) * sizeof(type));              \
}

#define ALLOC_2d_array(type, array_name, jd, kd) {                \
  array_name = (type **) calloc((jd), sizeof(type *));            \
  for (int ii = 0; ii < (jd); ii++) {                             \
    (array_name[ii]) = (type *) calloc((kd), sizeof(type));       \
  }                                                               \
}

#define PI   (3.1415926535897932384626433832795028841971693993751)

#endif
