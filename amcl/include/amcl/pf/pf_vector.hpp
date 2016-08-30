#pragma once

#include <cstdio>

// The basic vector
struct Pose
{
  double v[3];

  // static functions
  static Pose Sum(const Pose &a, const Pose &b);
  static Pose Sub(const Pose &a, const Pose &b);
  static Pose CoordSum(const Pose &a, const Pose &b);
  static Pose CoordSub(const Pose &a, const Pose &b);

  Pose(const double &x = 0.0, const double &y = 0.0, const double &a = 0.0);
  bool isFinite();
  void dumpToFile(FILE *file, const char *fmt);
};


// The basic matrix
struct Covariance
{
  double m[3][3];

  Covariance();
  bool isFinite();
  void dumpToFile(FILE *file, const char *fmt);
  void unitary(Covariance& r, Covariance& d);
};


// Return a zero vector
//Pose pf_vector_zero();

// Check for NAN or INF in any component
//int pf_vector_finite(Pose a);

// Print a vector
//void pf_vector_fprintf(Pose s, FILE *file, const char *fmt);

// Simple vector addition
//Pose pf_vector_add(Pose a, Pose b);

// Simple vector subtraction
//Pose pf_vector_sub(Pose a, Pose b);

// Transform from local to global coords (a + b)
//Pose pf_vector_coord_add(Pose a, Pose b);

// Transform from global to local coords (a - b)
//Pose pf_vector_coord_sub(Pose a, Pose b);


// Return a zero matrix
//Covariance pf_matrix_zero();

// Check for NAN or INF in any component
//int pf_matrix_finite(Covariance a);

// Print a matrix
//void pf_matrix_fprintf(Covariance s, FILE *file, const char *fmt);

// Compute the matrix inverse.  Will also return the determinant,
// which should be checked for underflow (indicated singular matrix).
//Covariance pf_matrix_inverse(Covariance a, double *det);

// Decompose a covariance matrix [a] into a rotation matrix [r] and a
// diagonal matrix [d] such that a = r * d * r^T.
//void pf_matrix_unitary(Covariance *r, Covariance *d, Covariance a);

