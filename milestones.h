#pragma once

#include <vector>

struct milestone {
  milestone() {}
  milestone(double coord, unsigned index) : x(coord), idx(index) {}

  bool crossed(double u, double v) const;

  double x;
  unsigned idx;
};

struct milestone_list {
  milestone_list(const std::vector<double>& coordinates);

  const milestone* crossed(double u, double v) const;

  const milestone& operator[](unsigned i) const;

  unsigned size() const;

  std::vector<milestone>::iterator begin() { 
    return milestones.begin(); 
  }
  
  std::vector<milestone>::iterator end() { 
    return milestones.end(); 
  }

  unsigned idx;
  std::vector<milestone> milestones;
};
