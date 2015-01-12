#pragma once

#include <random>

class integrator {
 public:
  integrator(double delta_t, double temperature_);

  double next(std::mt19937& rng);

 public:
  double x;
  double time;

 private:
  const double dt, temperature;
  std::normal_distribution<double> gaussian;
};
