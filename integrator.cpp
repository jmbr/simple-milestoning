#include "integrator.h"

integrator::integrator(double delta_t,
                       double temperature_)
    : dt(delta_t),
      temperature(temperature_),
      gaussian(0.0, 1.0) {
}

static inline double F(double x) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  return -0.2e1 / 0.625e3 * x3 - x2 / 0.250e3 + 0.4e1 / 0.25e2 * x + 0.2e1 / 0.15e2;
}

double integrator::next(std::mt19937& rng) {
  x += F(x) * dt + sqrt(2.0 * temperature * dt) * gaussian(rng);
  time += dt;

  return x;
}
