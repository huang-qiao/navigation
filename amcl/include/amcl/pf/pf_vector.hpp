#pragma once

#include <cstdio>
#include "eig3.hpp"

namespace amcl {
namespace pf {

// The basic vector
struct Pose
{
  static Pose CoordAdd(const Pose &a, const Pose &b);
  static Pose CoordSub(const Pose &a, const Pose &b);

  Pose(double x = 0, double y = 0, double th = 0) {
    v[0] = x;
    v[1] = y;
    v[2] = th;
  }
  double x() const { return v[0]; }
  double y() const { return v[1]; }
  double a() const { return v[2]; }

  double v[3];

  bool isFinite();

  inline Pose operator+(const Pose &arg) const {
    Pose p;
    p.v[0] = v[0] + arg.x();
    p.v[1] = v[1] + arg.y();
    p.v[2] = v[2] + arg.a();
    return p;
  }

  inline Pose operator-(const Pose &arg) const {
    Pose p;
    p.v[0] = v[0] - arg.x();
    p.v[1] = v[1] - arg.y();
    p.v[2] = v[2] - arg.a();
    return p;
  }
};

// The basic matrix
struct Covariance
{
  Covariance() {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        m[i][j] = 0.0;
      }
    }
  }
  double m[3][3];

  bool isFinite();
};

// Decompose the convariance matrix into a rotation matrix and a diagonal matrix.
void MatrixUnitary(const Covariance &a, Covariance &rot, Covariance &diag);

} // namespace pf
} // namespace amcl
