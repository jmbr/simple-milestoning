#pragma once

#include <cstddef>

#include <vector>
#include <iostream>

struct milestone {
  milestone(double coord, unsigned index) : x(coord), idx(index) {}

  bool crossed(double u, double v) const;

  double x;
  size_t idx;

  friend std::ostream& operator<<(std::ostream& stream,const milestone& m) {
    return stream << "milestone #" << m.idx << " (x = " << m.x << ")";
  }
};

struct milestone_list {
  milestone_list(const std::vector<double>& coordinates);

  const milestone* crossed(double u, double v) const;

  const milestone& operator[](unsigned i) const;

  size_t idx;
  std::vector<milestone> milestones;
};
