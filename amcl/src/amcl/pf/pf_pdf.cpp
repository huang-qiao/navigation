#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#include "pf_pdf.hpp"

// Random number generator seed value
static unsigned int pf_pdf_seed;


/**************************************************************************
 * Gaussian
 *************************************************************************/

// Create a gaussian pdf
PdfGaussianPtr PdfGaussian::CreatePdf(const Pose &x, const Covariance &cx)
{
  Covariance cd;
  //PdfGaussianPtr pdf;

  PdfGaussianPtr pdf = std::make_shared<PdfGaussian>(); //pdf = (PdfGaussian*)malloc(sizeof(PdfGaussian));

  pdf->x = x;
  pdf->cx = cx;
  //pdf->cxi = pf_matrix_inverse(cx, &pdf->cxdet);

  // Decompose the convariance matrix into a rotation
  // matrix and a diagonal matrix.
  pdf->cx.unitary(pdf->cr, cd); //pf_matrix_unitary(&pdf->cr, &cd, pdf->cx);
  pdf->cd.v[0] = sqrt(cd.m[0][0]);
  pdf->cd.v[1] = sqrt(cd.m[1][1]);
  pdf->cd.v[2] = sqrt(cd.m[2][2]);

  // Initialize the random number generator
  //pdf->rng = gsl_rng_alloc(gsl_rng_taus);
  //gsl_rng_set(pdf->rng, ++pf_pdf_seed);
  srand48(++pf_pdf_seed);

  return pdf;
}


/* Destroy the pdf
void pf_pdf_gaussian_free(PdfGaussianPtr pdf)
{
  //gsl_rng_free(pdf->rng);
  free(pdf);
  return;
}
*/

/*
// Compute the value of the pdf at some point [x].
double pf_pdf_gaussian_value(pf_pdf_gaussian_t *pdf, Pose x)
{
  int i, j;
  Pose z;
  double zz, p;

  z = pf_vector_sub(x, pdf->x);

  zz = 0;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      zz += z.v[i] * pdf->cxi.m[i][j] * z.v[j];

  p =  1 / (2 * M_PI * pdf->cxdet) * exp(-zz / 2);

  return p;
}
*/


// Generate a sample from the the pdf.
Pose PdfGaussian::sample()
{
  int i, j;
  Pose r;
  Pose x;

  // Generate a random vector
  for (i = 0; i < 3; i++)
  {
    //r.v[i] = gsl_ran_gaussian(pdf->rng, pdf->cd.v[i]);
    r.v[i] = RandomGaussian(this->cd.v[i]);
  }

  for (i = 0; i < 3; i++)
  {
    x.v[i] = this->x.v[i];
    for (j = 0; j < 3; j++)
      x.v[i] += this->cr.m[i][j] * r.v[j];
  }

  return x;
}

// Draw randomly from a zero-mean Gaussian distribution, with standard
// deviation sigma.
// We use the polar form of the Box-Muller transformation, explained here:
//   http://www.taygeta.com/random/gaussian.html
double RandomGaussian(double sigma)
{
  double x1, x2, w, r;

  do
  {
    do { r = drand48(); } while (r==0.0);
    x1 = 2.0 * r - 1.0;
    do { r = drand48(); } while (r==0.0);
    x2 = 2.0 * r - 1.0;
    w = x1*x1 + x2*x2;
  } while(w > 1.0 || w==0.0);

  return(sigma * x2 * sqrt(-2.0*log(w)/w));
}

#if 0

/**************************************************************************
 * Discrete
 * Note that GSL v1.3 and earlier contains a bug in the discrete
 * random generator.  A patched version of the the generator is included
 * in gsl_discrete.c.
 *************************************************************************/


// Create a discrete pdf
pf_pdf_discrete_t *pf_pdf_discrete_alloc(int count, double *probs)
{
  pf_pdf_discrete_t *pdf;

  pdf = calloc(1, sizeof(pf_pdf_discrete_t));

  pdf->prob_count = count;
  pdf->probs = malloc(count * sizeof(double));
  memcpy(pdf->probs, probs, count * sizeof(double));

  // Initialize the random number generator
  pdf->rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(pdf->rng, ++pf_pdf_seed);

  // Initialize the discrete distribution generator
  pdf->ran = gsl_ran_discrete_preproc(count, probs);

  return pdf;
}


// Destroy the pdf
void pf_pdf_discrete_free(pf_pdf_discrete_t *pdf)
{
  gsl_ran_discrete_free(pdf->ran);
  gsl_rng_free(pdf->rng);
  free(pdf->probs);
  free(pdf);
  return;
}


// Compute the value of the probability of some element [i]
double pf_pdf_discrete_value(pf_pdf_discrete_t *pdf, int i)
{
  return pdf->probs[i];
}


// Generate a sample from the the pdf.
int pf_pdf_discrete_sample(pf_pdf_discrete_t *pdf)
{
  int i;

  i = gsl_ran_discrete(pdf->rng, pdf->ran);
  assert(i >= 0 && i < pdf->prob_count);

  return i;
}

#endif
