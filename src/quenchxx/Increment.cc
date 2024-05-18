/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "quenchxx/Increment.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

void Increment::ones() {
  this->fields().constantValue(1.0);
}

// -----------------------------------------------------------------------------

oops::LocalIncrement Increment::getLocal(const GeometryIterator & geometryIterator) const {
  int index = 0;
  if (this->geometry()->iteratorDimension() == 2) {
    std::vector<int> variableSizes;
    for (const auto & variableSize : this->geometry()->variableSizes(this->variables())) {
      variableSizes.push_back(static_cast<int>(variableSize));
    }
    size_t valuesSize = std::accumulate(variableSizes.begin(), variableSizes.end(), 0);
    std::vector<double> values(valuesSize);
    for (const auto & var : this->variables().variables()) {
      const auto view = atlas::array::make_view<double, 2>(this->fields().fields().field(var));
      for (size_t jlevel = 0; jlevel < this->geometry()->nlevs(); ++jlevel) {
        values[index] = view(geometryIterator.jnode(), jlevel);
        ++index;
      }
    }
    return oops::LocalIncrement(this->variables(), values, variableSizes);
  } else {
    std::vector<int> variableSizes(this->variables().size(), 1);
    size_t valuesSize = std::accumulate(variableSizes.begin(), variableSizes.end(), 0);
    std::vector<double> values(valuesSize);
    for (const auto & var : this->variables().variables()) {
      const auto view = atlas::array::make_view<double, 2>(this->fields().fields().field(var));
      values[index] = view(geometryIterator.jnode(), geometryIterator.jlevel());
      ++index;
    }
    return oops::LocalIncrement(this->variables(), values, variableSizes);
  }
}

// -----------------------------------------------------------------------------

void Increment::setLocal(const oops::LocalIncrement & localIncrement,
                      const GeometryIterator & geometryIterator) {
  std::vector<double> values = localIncrement.getVals();
  size_t index = 0;
  if (this->geometry()->iteratorDimension() == 2) {
    for (const auto & var : this->variables().variables()) {
      auto view = atlas::array::make_view<double, 2>(this->fields().fields().field(var));
      for (size_t jlevel = 0; jlevel < this->geometry()->nlevs(); ++jlevel) {
        view(geometryIterator.jnode(), jlevel) = values[index];
        ++index;
      }
    }
  } else {
    for (const auto & var : this->variables().variables()) {
      auto view = atlas::array::make_view<double, 2>(this->fields().fields().field(var));
      view(geometryIterator.jnode(), geometryIterator.jlevel()) = values[index];
      ++index;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
