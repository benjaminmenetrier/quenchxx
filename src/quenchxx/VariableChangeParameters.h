/*
 * (C) Copyright 2024-     UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/VariableChangeParametersBase.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------
/// VariableChange parameters class

class VariableChangeParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, VariableChangeParametersBase)
 public:
  // VADER parameters
  oops::Parameter<std::map<std::string, std::vector<std::string>>> vaderCustomCookbook{
    "vader custom cookbook", vaderQuenchxxCustomCookbook(), this};
  oops::OptionalParameter<vader::VaderParameters> vaderParam{"vader", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
