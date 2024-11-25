/*
 * (C) Copyright 2024 Meteorologisk Institutt.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/HofX.h"
#include "oops/runs/Run.h"
#include "quenchxx/instantiateQuenchMatrices.h"
#include "quenchxx/Logbook.h"
#include "quenchxx/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  quenchxx::instantiateQuenchMatrices();
  oops::HofX<quenchxx::Traits> hofx;
  quenchxx::Logbook::start();
  run.execute(hofx);
  quenchxx::Logbook::stop();
  return 0;
}