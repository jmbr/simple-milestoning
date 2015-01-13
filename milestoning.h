#pragma once

#include <cstddef>

#include <random>
#include <iostream>

#include <armadillo>

#include "linspace.h"
#include "milestones.h"
#include "integrator.h"

struct milestoning {
  milestoning(unsigned num_milestones, unsigned long max_fragments,
              double delta_t, double temperature, unsigned long random_seed);

  void run();
  void sample_fragments(const milestone* start);
  void compute_results();

  unsigned num_milestones;
  unsigned long max_fragments;
  double delta_t, temperature;
  unsigned long random_seed, max_steps;
  std::mt19937 rng;
  integrator dyna;
  milestone_list milestones;
  arma::mat A, B;
};
