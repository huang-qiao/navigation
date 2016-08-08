#pragma once

#include "pf_vector.hpp"

#include <memory>

namespace amcl {
namespace pf {

// Gaussian

// Draw randomly from a zero-mean Gaussian distribution, with standard
// deviation sigma.
// We use the polar form of the Box-Muller transformation, explained here:
//   http://www.taygeta.com/random/gaussian.html
//double pf_ran_gaussian(double sigma);
double RandomGaussian(double sigma);

// Gaussian PDF info
class PdfGaussian
{
public:
  // Generate a sample from the the pdf.
  Pose genGaussianSample();

  // Compute the value of the pdf at some point
  double getGaussianValue(const Pose &pose);

  // factory: create a gaussian pdf
  static std::shared_ptr<PdfGaussian> CreatePdfGaussian(const Pose &x, const Covariance &cx);

private:
  // Mean, covariance and inverse covariance
  Pose x_;
  Covariance cx_;

  Covariance cxi_;
  double cxdet_;

  // Decomposed covariance matrix (rotation * diagonal)
  Covariance cr_;
  Pose cd_;
};

using PdfGaussianPtr = std::shared_ptr<PdfGaussian>;


// Destroy the pdf
//void pf_pdf_gaussian_free(pf_pdf_gaussian_t *pdf);

// Compute the value of the pdf at some point [z].
//double pf_pdf_gaussian_value(pf_pdf_gaussian_t *pdf, pf_vector_t z);


} // namespace pf
} // namespace amcl
