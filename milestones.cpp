#include <cmath>

#include "milestones.h"

bool milestone::crossed(double u, double v) const {
  return (u < x && x <= v) || (v < x && x <= u);
}

milestone_list::milestone_list(const std::vector<double>& coordinates) {
  idx = 0;
  for (unsigned k = 0; k < coordinates.size(); ++k)
    milestones.push_back(milestone(coordinates[k], idx++));
}

const milestone* milestone_list::crossed(double u, double v) const {
  for (auto& milestone : milestones)
    if (milestone.crossed(u, v))
      return &milestone;
  return nullptr;
}

const milestone& milestone_list::operator[](unsigned i) const {
  return milestones.at(i);
}
