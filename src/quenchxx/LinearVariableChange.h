/*
 * (C) Copyright 2023  UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "atlas/field.h"

#include "oops/util/Printable.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/Increment.h"
#include "quenchxx/State.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

class LinearVariableChange: public util::Printable {
 public:
  static const std::string classname()
    {return "quenchxx::LinearVariableChange";}

  // Constructor/destructor
  LinearVariableChange(const Geometry &,
                       const eckit::Configuration &);
  ~LinearVariableChange()
    {}

  // Linear variable changes: TL, inverseTL, AD and inverseAD
  void changeVarTL(Increment &,
                   const varns::Variables &) const;
  void changeVarInverseTL(Increment &,
                          const varns::Variables &) const;
  void changeVarAD(Increment &,
                   const varns::Variables &) const;
  void changeVarInverseAD(Increment &,
                          const varns::Variables &) const;

  // Trajectory setup
  void changeVarTraj(const State &,
                     const varns::Variables &)
    {}

 private:
  // Print
  void print(std::ostream & os) const override
    {os << "LinearVariableChange";}

  // Map from output to input variables
  std::map<std::string, std::string> map_;

  // Multiplicative factor
  atlas::FieldSet multiplierFset_;
};
// -----------------------------------------------------------------------------

}  // namespace quenchxx
