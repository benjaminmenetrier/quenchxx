/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <limits>

#include "oops/runs/MakeObsPatched.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"
#include "quenchxx/instantiateQuenchMatrices.h"
#include "quenchxx/Logbook.h"
#include "quenchxx/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::Log::test().setf(std::ios::scientific);
  oops::Log::test().precision(std::numeric_limits<double>::digits10+1);
  quenchxx::instantiateQuenchMatrices();
  oops::MakeObsPatched<quenchxx::Traits> mop;
  quenchxx::Logbook::start();
  run.execute(mop);
  quenchxx::Logbook::stop();
  return 0;
}