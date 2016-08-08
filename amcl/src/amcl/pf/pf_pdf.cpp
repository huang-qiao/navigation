#include <cassert>
#include <cmath>
#include <cstdlib>
#include <string>

#include "pf_pdf.hpp"

#include <cstring>
#include <iostream>
#define __FILENAME__ \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define TRACE_FUNC                                                     \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << std::endl;                                            \
  } while (0);

#define TRACE_FUNC_ENTER                                               \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ": [ENTER]" << std::endl;                             \
  } while (0);

#define TRACE_FUNC_EXIT                                                \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ": [EXIT]" << std::endl;                              \
  } while (0);


using namespace amcl::pf;

// Random number generator seed value
static unsigned int pf_pdf_seed;

/**************************************************************************
 * Gaussian
 *************************************************************************/


// Create a gaussian pdf
PdfGaussianPtr PdfGaussian::CreatePdfGaussian(const Pose &x, const Covariance &cx)
{
  PdfGaussianPtr pdf = std::make_shared<PdfGaussian>();
  pdf->x_ = x;
  pdf->cx_ = cx;

  Covariance cd;

  // Decompose the convariance matrix into a rotation matrix and a diagonal matrix.
  MatrixUnitary(pdf->cx_, pdf->cr_, cd);

  pdf->cd_.v[0] = sqrt(cd.m[0][0]);
  pdf->cd_.v[1] = sqrt(cd.m[1][1]);
  pdf->cd_.v[2] = sqrt(cd.m[2][2]);

  // Initialize the random number generator
  srand48(++pf_pdf_seed);

  return pdf;
}

/*
pf_pdf_gaussian_t *pf_pdf_gaussian_alloc(pf_vector_t x, pf_matrix_t cx)
{
  pf_matrix_t cd;
  pf_pdf_gaussian_t *pdf;

  pdf = calloc(1, sizeof(pf_pdf_gaussian_t));

  pdf->x = x;
  pdf->cx = cx;
  //pdf->cxi = pf_matrix_inverse(cx, &pdf->cxdet);

  // Decompose the convariance matrix into a rotation
  // matrix and a diagonal matrix.
  pf_matrix_unitary(&pdf->cr, &cd, pdf->cx);
  pdf->cd.v[0] = sqrt(cd.m[0][0]);
  pdf->cd.v[1] = sqrt(cd.m[1][1]);
  pdf->cd.v[2] = sqrt(cd.m[2][2]);

  // Initialize the random number generator
  //pdf->rng = gsl_rng_alloc(gsl_rng_taus);
  //gsl_rng_set(pdf->rng, ++pf_pdf_seed);
  srand48(++pf_pdf_seed);

  return pdf;
}
*/

/* Destroy the pdf
void pf_pdf_gaussian_free(pf_pdf_gaussian_t *pdf)
{
  //gsl_rng_free(pdf->rng);
  free(pdf);
  return;
}
*/


// Compute the value of the pdf at some point [x].
//double pf_pdf_gaussian_value(pf_pdf_gaussian_t *pdf, pf_vector_t x)
double PdfGaussian::getGaussianValue(const Pose &pose)
{
  int i, j;
  //pf_vector_t z;
  Pose z;
  double zz, p;
  
  //z = pf_vector_sub(x, pdf->x);
  z = pose - x_;

  zz = 0;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      zz += z.v[i] * cxi_.m[i][j] * z.v[j];

  p =  1 / (2 * M_PI * cxdet_) * exp(-zz / 2);
          
  return p;
}

// Generate a sample from the the pdf.
//pf_vector_t pf_pdf_gaussian_sample(pf_pdf_gaussian_t *pdf)
Pose PdfGaussian::genGaussianSample()
{
  //int i, j;
  //pf_vector_t r;
  //pf_vector_t x;
  Pose r, x;

  // Generate a random vector
  for (size_t i = 0; i < 3; i++) {
    r.v[i] = RandomGaussian(cd_.v[i]);
  }

  for (size_t i = 0; i < 3; i++) {
    x.v[i] = x_.v[i];
    for (size_t j = 0; j < 3; j++)
      x.v[i] += cr_.m[i][j] * r.v[j];
  } 
  
  return x;
}

// Draw randomly from a zero-mean Gaussian distribution, with standard deviation sigma.
// We use the polar form of the Box-Muller transformation, explained here:
// http://www.taygeta.com/random/gaussian.html
//double pf_ran_gaussian(double sigma)
double amcl::pf::RandomGaussian(double sigma)
{
  double x1, x2, w, r;

  do {
    do { r = drand48(); } while (r==0.0);
    x1 = 2.0 * r - 1.0;
    do { r = drand48(); } while (r==0.0);
    x2 = 2.0 * r - 1.0;
    w = x1*x1 + x2*x2;
  } while(w > 1.0 || w==0.0);

  //std::cout << "RandomGaussian return: " << (sigma * x2 * sqrt(-2.0*log(w)/w)) << std::endl;
  return(sigma * x2 * sqrt(-2.0*log(w)/w));
}
