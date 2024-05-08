/*
 * (C) Copyright 2017-2021 UCAR
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "vader/vader.h"

#include "src/VariableChange.h"

namespace oops {
  class Variables;
}

namespace quenchxx {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------

class VariableChange : public quench::VariableChange {
 public:
  static const std::string classname() {return "quenchxx::VariableChange";}

  VariableChange(const eckit::Configuration &,
                 const Geometry &);

  void changeVar(State &,
                 const oops::Variables &) const;
  void changeVarInverse(State &,
                        const oops::Variables &) const;

 private:
  void print(std::ostream & os) const override {os << *vader_;};
  std::vector<eckit::LocalConfiguration> alias_;
  std::unique_ptr<vader::Vader> vader_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
