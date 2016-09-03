#include <cmath>
#include <cstdlib>

#include "pf_vector.hpp"
#include "eig3.hpp"


// Return a zero vector
//Pose pf_vector_zero()
Pose::Pose(const double &x, const double &y, const double &a)
{
  v[0] = x;
  v[1] = y;
  v[2] = a;
}


// Check for NAN or INF in any component
//int pf_vector_finite(Pose a)
bool Pose::isFinite()
{
  for (size_t i = 0; i < 3; i++) {
    if (!finite(v[i])) {
      return false;
    }
  }
  return true;
}


// Print a vector
//void pf_vector_fprintf(Pose a, FILE *file, const char *fmt)
void Pose::dumpToFile(FILE *file, const char *fmt)
{
  for (size_t i = 0; i < 3; i++) {
    fprintf(file, fmt, v[i]);
    fprintf(file, " ");
  }
  fprintf(file, "\n");

  return;
}


// Simple vector addition
//Pose pf_vector_add(Pose a, Pose b)
Pose Pose::Sum(const Pose &a, const Pose &b)
{
  Pose c;

  c.v[0] = a.v[0] + b.v[0];
  c.v[1] = a.v[1] + b.v[1];
  c.v[2] = a.v[2] + b.v[2];

  return c;
}


// Simple vector subtraction
//Pose pf_vector_sub(Pose a, Pose b)
Pose Pose::Sub(const Pose &a, const Pose &b)
{
  Pose c;

  c.v[0] = a.v[0] - b.v[0];
  c.v[1] = a.v[1] - b.v[1];
  c.v[2] = a.v[2] - b.v[2];

  return c;
}


// Transform from local to global coords (a + b)
//Pose pf_vector_coord_add(Pose a, Pose b)
Pose Pose::CoordSum(const Pose &a, const Pose &b)
{
  Pose c;

  c.v[0] = b.v[0] + a.v[0] * cos(b.v[2]) - a.v[1] * sin(b.v[2]);
  c.v[1] = b.v[1] + a.v[0] * sin(b.v[2]) + a.v[1] * cos(b.v[2]);
  c.v[2] = b.v[2] + a.v[2];
  c.v[2] = atan2(sin(c.v[2]), cos(c.v[2]));

  return c;
}


// Transform from global to local coords (a - b)
//Pose pf_vector_coord_sub(Pose a, Pose b)
Pose Pose::CoordSub(const Pose &a, const Pose &b)
{
  Pose c;

  c.v[0] = +(a.v[0] - b.v[0]) * cos(b.v[2]) + (a.v[1] - b.v[1]) * sin(b.v[2]);
  c.v[1] = -(a.v[0] - b.v[0]) * sin(b.v[2]) + (a.v[1] - b.v[1]) * cos(b.v[2]);
  c.v[2] = a.v[2] - b.v[2];
  c.v[2] = atan2(sin(c.v[2]), cos(c.v[2]));

  return c;
}


// Return a zero matrix
//Covariance pf_matrix_zero()
Covariance::Covariance()
{
  //int i, j;
  //Covariance c;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      m[i][j] = 0.0; // c.m[i][j] = 0.0;
    }
  }

  //return c;
}


// Check for NAN or INF in any component
//int pf_matrix_finite(Covariance a)
bool Covariance::isFinite()
{
  //int i, j;

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      if (!finite(m[i][j])) {
        return false;
      }
    }
  }

  return true;
}


// Print a matrix
//void pf_matrix_fprintf(Covariance a, FILE *file, const char *fmt)
void Covariance::dumpToFile(FILE *file, const char *fmt)
{
  int i, j;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      fprintf(file, fmt, m[i][j]);
      fprintf(file, " ");
    }
    fprintf(file, "\n");
  }
  return;
}


/*
// Compute the matrix inverse
Covariance pf_matrix_inverse(Covariance a, double *det)
{
  double lndet;
  int signum;
  gsl_permutation *p;
  gsl_matrix_view A, Ai;

  Covariance ai;

  A = gsl_matrix_view_array((double*) a.m, 3, 3);
  Ai = gsl_matrix_view_array((double*) ai.m, 3, 3);

  // Do LU decomposition
  p = gsl_permutation_alloc(3);
  gsl_linalg_LU_decomp(&A.matrix, p, &signum);

  // Check for underflow
  lndet = gsl_linalg_LU_lndet(&A.matrix);
  if (lndet < -1000)
  {
    //printf("underflow in matrix inverse lndet = %f", lndet);
    gsl_matrix_set_zero(&Ai.matrix);
  }
  else
  {
    // Compute inverse
    gsl_linalg_LU_invert(&A.matrix, p, &Ai.matrix);
  }

  gsl_permutation_free(p);

  if (det)
    *det = exp(lndet);

  return ai;
}
*/


// Decompose a covariance matrix [a] into a rotation matrix [r] and a diagonal
// matrix [d] such that a = r d r^T.
//void pf_matrix_unitary(Covariance *r, Covariance *d, Covariance a)
void Covariance::unitary(Covariance &r, Covariance &d)
{
  int i, j;
  /*
  gsl_matrix *aa;
  gsl_vector *eval;
  gsl_matrix *evec;
  gsl_eigen_symmv_workspace *w;

  aa = gsl_matrix_alloc(3, 3);
  eval = gsl_vector_alloc(3);
  evec = gsl_matrix_alloc(3, 3);
  */

  double **aa;
  double eval[3];
  double **evec;

  aa = (double**)malloc(3 * sizeof(double*));
  evec = (double**)malloc(3 * sizeof(double*));
  for (i = 0; i < 3; i++) {
    aa[i] = (double*)malloc(3 * sizeof(double));
    evec[i] = (double*)malloc(3 * sizeof(double));
  }

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      //gsl_matrix_set(aa, i, j, a.m[i][j]);
      aa[i][j] = m[i][j]; //aa[i][j] = a.m[i][j];
    }
  }

  // Compute eigenvectors/values
  /*
  w = gsl_eigen_symmv_alloc(3);
  gsl_eigen_symmv(aa, eval, evec, w);
  gsl_eigen_symmv_free(w);
  */

  eigen_decomposition(aa,evec,eval);

  d = Covariance(); //*d = pf_matrix_zero();
  for (i = 0; i < 3; i++)
  {
    //d->m[i][i] = gsl_vector_get(eval, i);
    d.m[i][i] = eval[i]; //d->m[i][i] = eval[i];
    for (j = 0; j < 3; j++)
    {
      //r->m[i][j] = gsl_matrix_get(evec, i, j);
      r.m[i][j] = evec[i][j]; //r->m[i][j] = evec[i][j];
    }
  }

  //gsl_matrix_free(evec);
  //gsl_vector_free(eval);
  //gsl_matrix_free(aa);

  return;
}

