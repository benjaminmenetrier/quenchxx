/*
 * (C) Copyright 2024-     UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/State.h"
#include "quenchxx/VariablesSwitch.h"

#include "vader/vader.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname()
    {return "quenchxx::VariableChange";}

  // Constructor/destructor
  VariableChange(const eckit::Configuration &,
                 const Geometry &);
  ~VariableChange()
    {}

  // Variable changes: direct and inverse
  void changeVar(State &,
                 const varns::Variables &) const;
  void changeVarInverse(State &,
                        const varns::Variables &) const;

 private:
  // Print
  void print(std::ostream & os) const override
    {os << "VariableChange";};

  // Geometry pointer
  std::shared_ptr<const Geometry> geom_;

  // VADER
  std::unique_ptr<vader::Vader> vader_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
