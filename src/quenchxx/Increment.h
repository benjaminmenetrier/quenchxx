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

#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "oops/base/LocalIncrement.h"

#include "quenchxx/Fields.h"
#include "quenchxx/GeometryIterator.h"
#include "quenchxx/State.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {
  class Geometry;

// -----------------------------------------------------------------------------
/// Increment class

class Increment : public util::Printable,
                  public util::Serializable,
                  private util::ObjectCounter<Increment> {
 public:
  static const std::string classname()
    {return "quenchxx::Increment";}

  // Constructors/destructor
  Increment(const Geometry &,
            const varns::Variables &,
            const util::DateTime &);
  Increment(const Geometry &,
            const varns::Variables &,
            const util::DateTime &,
            const util::DateTime &);
  Increment(const Geometry &,
            const Increment &);
  Increment(const Increment &,
            const bool);

  // Basic operators
  void diff(const State &,
            const State &);
  void zero()
    {fields_->zero();}
  void zero(const util::DateTime &);
  void ones()
    {this->fields().constantValue(1.0);}
  void dirac(const eckit::Configuration & config)
    {fields_->dirac(config);}
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &,
            const Increment &,
            const bool check = true);
  double dot_product_with(const Increment & dx) const
    {return fields_->dot_product_with(*dx.fields_);}
  void schur_product_with(const Increment & dx)
    {fields_->schur_product_with(*dx.fields_);}
  void random()
    {fields_->random();}
  double max(const varns::Variables & var) const
    {return fields_->max(var);}
  double min(const varns::Variables & var) const
    {return fields_->min(var);}

  // I/O and diagnostics
  void read(const eckit::Configuration & config)
    {fields_->read(config);}
  void write(const eckit::Configuration & config) const
    {fields_->write(config);}
  double norm() const
    {return fields_->norm();}
  const util::DateTime & validTime() const
    {return fields_->time();}
  void updateTime(const util::Duration & dt)
    {fields_->time() += dt;}

  // ATLAS FieldSet accessor
  void toFieldSet(atlas::FieldSet & fset) const
    {fields_->toFieldSet(fset);}
  void fromFieldSet(const atlas::FieldSet & fset)
    {fields_->fromFieldSet(fset);}
  const atlas::FieldSet & fieldSet() const
    {return fields_->fieldSet();}
  atlas::FieldSet & fieldSet()
    {return fields_->fieldSet();}
  void synchronizeFields()
    {fields_->synchronizeFields();}

  // Access to fields
  Fields & fields()
    {return *fields_;}
  const Fields & fields() const
    {return *fields_;}
  std::shared_ptr<const Geometry> geometry() const
    {return fields_->geometry();}

  // Other
  void accumul(const double & zz,
               const State & xx)
    {fields_->axpy(zz, xx.fields());}
  const varns::Variables & variables() const
    {return fields_->variables();}
  void interpolateTL(const Locations & locs,
                     GeoVaLs & gv) const
    {fields_->interpolate(locs, gv);}
  void interpolateAD(const Locations & locs,
                     const GeoVaLs & gv)
    {fields_->interpolateAD(locs, gv);}

  // Serialization
  size_t serialSize() const
    {return fields_->serialSize();}
  void serialize(std::vector<double> & vect) const
    {fields_->serialize(vect);}
  void deserialize(const std::vector<double> & vect,
                   size_t & index)
    {fields_->deserialize(vect, index);}
  friend eckit::Stream & operator<<(eckit::Stream &,
                                    const Increment &);
  friend eckit::Stream & operator>>(eckit::Stream &,
                                    Increment &);

  // Local increment
  oops::LocalIncrement getLocal(const GeometryIterator & geometryIterator) const;
  void setLocal(const oops::LocalIncrement & localIncrement,
                const GeometryIterator & geometryIterator);

 private:
  // Print
  void print(std::ostream &) const;

  // Fields
  std::unique_ptr<Fields> fields_;
};

// -----------------------------------------------------------------------------

}  // namespace quenchxx
