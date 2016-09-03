#pragma once

#include <cmath>

#include <cstring>
#include <iostream>
#define __FILENAME__ \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#if 1
#define TRACE_FUNC                                                     \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << std::endl;                                            \
  } while (0);

#define TRACE_FUNC_ENTER                                               \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ":[ENTER]" << std::endl;                              \
  } while (0);

#define TRACE_FUNC_EXIT                                                \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ":[EXIT]" << std::endl;                               \
  } while (0);
#else
#define TRACE_FUNC
#define TRACE_FUNC_ENTER
#define TRACE_FUNC_EXIT
#endif

using namespace std;

inline double normalize(double z) { return atan2(sin(z), cos(z)); }
inline double angle_diff(double a, double b) {
  double d1, d2;
  a = normalize(a);
  b = normalize(b);
  d1 = a - b;
  d2 = 2 * M_PI - fabs(d1);
  if (d1 > 0) d2 *= -1.0;
  if (fabs(d1) < fabs(d2))
    return (d1);
  else
    return (d2);
}
