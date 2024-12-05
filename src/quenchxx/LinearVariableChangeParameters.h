/*
 * (C) Copyright 2021-2024 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/LinearVariableChangeParametersBase.h"

namespace quenchxx {

// -------------------------------------------------------------------------------------------------
/// LinearVariableChange parameters class

class LinearVariableChangeParameters : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, LinearVariableChangeParametersBase)
 public:
  // ATLAS file (multiplicative factor)
  oops::OptionalParameter<eckit::LocalConfiguration> atlasFile{"atlas file", this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace quenchxx
