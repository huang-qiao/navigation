#pragma once

#include "pf_vector.hpp"

#include <memory>

/**************************************************************************
 * Gaussian
 *************************************************************************/

// Gaussian PDF info
struct PdfGaussian
{
  // Mean, covariance and inverse covariance
  Pose x;
  Covariance cx;
  //Covariance cxi;
  double cxdet;

  // Decomposed covariance matrix (rotation * diagonal)
  Covariance cr;
  Pose cd;

  // A random number generator
  //gsl_rng *rng;

  // Generate a sample from the the pdf.
  Pose sample();

  // Create a gaussian pdf
  static std::shared_ptr<PdfGaussian> CreatePdf(const Pose &x, const Covariance &cx);
};

using PdfGaussianPtr = std::shared_ptr<PdfGaussian>;

//PdfGaussianPtr PdfGaussian::CreatePdf(Pose x, Covariance cx);

// Destroy the pdf
//void pf_pdf_gaussian_free(PdfGaussianPtr pdf);

// Compute the value of the pdf at some point [z].
//double pf_pdf_gaussian_value(pf_pdf_gaussian_t *pdf, pf_vector_t z);

// Draw randomly from a zero-mean Gaussian distribution, with standard
// deviation sigma.
// We use the polar form of the Box-Muller transformation, explained here:
//   http://www.taygeta.com/random/gaussian.html
double RandomGaussian(double sigma);

// Generate a sample from the the pdf.
//Pose PdfGaussian::sample(PdfGaussianPtr pdf);


#if 0

/**************************************************************************
 * Discrete
 *************************************************************************/

// Discrete PDF info
typedef struct
{
  // The list of discrete probs
  int prob_count;
  double *probs;

  // A random number generator
  gsl_rng *rng;

  // The discrete prob generator
  gsl_ran_discrete_t *ran;

} pf_pdf_discrete_t;


// Create a discrete pdf
pf_pdf_discrete_t *pf_pdf_discrete_alloc(int count, double *probs);

// Destroy the pdf
void pf_pdf_discrete_free(pf_pdf_discrete_t *pdf);

// Compute the value of the probability of some element [i]
double pf_pdf_discrete_value(pf_pdf_discrete_t *pdf, int i);

// Generate a sample from the the pdf.
int pf_pdf_discrete_sample(pf_pdf_discrete_t *pdf);
#endif

