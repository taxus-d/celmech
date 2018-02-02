! useful macros
#ifdef __GFORTRAN__
#define CONCAT(x,y) x/**/y
#else
#define CONCAT(x,y) x ## y
#endif

