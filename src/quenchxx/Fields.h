/*
 * (C) Copyright 2022 UCAR.
 * (C) Copyright 2023-2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "eckit/serialisation/Stream.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "quenchxx/GeoVaLs.h"
#include "quenchxx/Interpolation.h"
#include "quenchxx/Locations.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {
  class Geometry;

// -----------------------------------------------------------------------------
/// Fields class

class Fields : public util::Printable,
               public util::Serializable,
               private util::ObjectCounter<Fields> {
 public:
  static const std::string classname()
    {return "quenchxx::Fields";}

  // Constructors/destructor
  Fields(const Geometry &,
         const varns::Variables &,
         const util::DateTime &);
  Fields(const Fields &,
         const Geometry &);
  Fields(const Fields &,
         const bool);
  Fields(const Fields &);
  ~Fields()
    {}

  // Basic operators
  void zero();
  void constantValue(const double &);
  void constantValue(const std::vector<double> &);
  void constantValue(const eckit::Configuration &);
  Fields & operator=(const Fields &);
  Fields & operator+=(const Fields &);
  Fields & operator-=(const Fields &);
  Fields & operator*=(const double &);
  void axpy(const double &,
            const Fields &);
  double dot_product_with(const Fields &) const;
  void schur_product_with(const Fields &);
  void dirac(const eckit::Configuration &);
  void random();
  void diff(const Fields &,
            const Fields &);
  double min(const varns::Variables &) const;
  double max(const varns::Variables &) const;
  void interpolate(const Locations &,
                   GeoVaLs &) const;
  void interpolateAD(const Locations &,
                     const GeoVaLs &);
  void forceWith(const Fields &,
                 const varns::Variables &);

  // ATLAS FieldSet
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);
  const atlas::FieldSet & fieldSet() const
    {return fset_;}
  atlas::FieldSet & fieldSet()
    {return fset_;}
  void synchronizeFields();

  // Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  std::shared_ptr<const Geometry> geometry() const
    {return geom_;}
  const varns::Variables & variables() const
    {return vars_;}
  const util::DateTime & time() const
    {return time_;}
  util::DateTime & time()
    {return time_;}
  void updateTime(const util::Duration & dt)
    {time_ += dt;}

  // Serialization
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &,
                   size_t &);
  friend eckit::Stream & operator<<(eckit::Stream &,
                                    const Fields &);
  friend eckit::Stream & operator>>(eckit::Stream &,
                                    Fields &);

  // Grid interpolations
  static std::vector<quenchxx::Interpolation>& interpolations();

  // Duplicate points
  void resetDuplicatePoints();

 private:
  // Print
  void print(std::ostream &) const;

  // Return grid interpolation
  std::vector<quenchxx::Interpolation>::iterator setupGridInterpolation(const Geometry &) const;

  // Return observations interpolation
  std::vector<quenchxx::Interpolation>::iterator setupObsInterpolation(const Locations &) const;

  // Reduce duplicate points
  void reduceDuplicatePoints();

  // Geometry
  std::shared_ptr<const Geometry> geom_;

  // Variables
  varns::Variables vars_;

  // Time
  util::DateTime time_;

  // Fieldset
  mutable atlas::FieldSet fset_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
