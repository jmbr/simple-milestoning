#include <cmath>

#include <iostream>

#include "milestoning.h"

int main(int argc, char* argv[]) {
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0] << " num-milestones max-fragments "
              << "delta-t temperature random-seed" << std::endl;
    return EXIT_SUCCESS;
  }

  const unsigned num_milestones = atoi(argv[1]);
  const unsigned long max_fragments = atof(argv[2]);
  const double delta_t = atof(argv[3]);
  const double temperature = atof(argv[4]);
  const unsigned long random_seed = atof(argv[5]);

  std::clog << "Using random seed " << random_seed << std::endl;

  milestoning mls(num_milestones, max_fragments,
                  delta_t, temperature, random_seed);

  mls.run();

  return EXIT_SUCCESS;
}

milestoning::milestoning(unsigned num_milestones_,
                         unsigned long max_fragments_,
                         double delta_t_,
                         double temperature_,
                         unsigned long random_seed_)
    : num_milestones(num_milestones_),
      max_fragments(max_fragments_),
      delta_t(delta_t_),
      temperature(temperature_),
      random_seed(random_seed_),
      max_steps(1e14),
      rng(random_seed),
      dyna(delta_t, temperature),
      milestones(linspace<double>(-10.0, 5.0, num_milestones)),
      A(arma::zeros<arma::mat>(num_milestones, num_milestones)),
      B(arma::zeros<arma::mat>(num_milestones, num_milestones)) {
}

void milestoning::run() {
  //  The last milestone is always absorbing, so we don't sample
  //  from it.
  for (unsigned n = 0; n < num_milestones-1; ++n) {
    std::cout << "Sampling from " << milestones[n] << "\n";
    sample_fragments(&milestones[n]);
  }

  compute_results();
}

void milestoning::sample_fragments(const milestone* start) {
  for (unsigned long f = 0; f < max_fragments; ++f) {
    dyna.x = start->x;
    dyna.time = 0.0;

    double xprev;

    unsigned long step;
    for (step = 1; step <= max_steps; ++step) {
      xprev = dyna.x;

      dyna.next(rng);

      // Check to see if we have crossed a milestone and ignore the
      // crossing event unless it happened on a different milestone
      // than the one we started from.
      const milestone* current = milestones.crossed(xprev, dyna.x);
      if (current && current != start) {
        const auto i = start->idx;
        const auto j = current->idx;

        ++A(i, j);
        B(i, j) += dyna.time;

        break;
      }
    }

    if (step == max_steps)
      std::cerr << "Warning: reached maximum number of steps. "
                << "You may want to increase max_steps."
                << std::endl;
  }
}

void milestoning::compute_results() {
  arma::mat K = A;
  arma::mat T = B;

  for (unsigned i = 0; i < num_milestones; ++i) {
    const double nrm = norm(K.row(i), 1);
    if (nrm > 0.0)
      K.row(i) /= nrm;

    for (unsigned j = 0; j < num_milestones; ++j)
      T(i, j) = A(i, j) > 0 ? T(i, j) / A(i, j) : 0.0;
  }

  arma::mat I = arma::eye(num_milestones, num_milestones);
  arma::vec e1 = I.col(0);

  K.row(num_milestones-1) = arma::zeros<arma::rowvec>(num_milestones);

  // By computing the stationary flux in this manner we don't need an
  // eigenvalue solver. In production code, however, it is better to
  // avoid matrix inversions and find the dominant eigenvector
  // directly.
  arma::rowvec q = arma::trans(e1) * arma::inv(I - K);
  q /= norm(q, 1);

  arma::vec t = (K % T) * arma::ones<arma::vec>(num_milestones);

  arma::vec aux = arma::solve(I - K, t);
  const double tau = arma::dot(e1, aux);

  std::cout << "K (transition matrix):\n" << K
            << "T (local MFPT matrix):\n" << T
            << "Estimated stationary flux (q):\n" << q
            << "Estimated MFPT: " << tau << std::endl;
}
