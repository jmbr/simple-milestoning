#pragma once

#include <cstddef>

#include <vector>

struct milestone {
  milestone(double coord, unsigned index) : x(coord), idx(index) {}

  bool crossed(double u, double v) const;

  double x;
  size_t idx;
};

struct milestone_list {
  milestone_list(const std::vector<double>& coordinates);

  const milestone* crossed(double u, double v) const;

  const milestone& operator[](unsigned i) const;

  size_t idx;
  std::vector<milestone> milestones;
};
