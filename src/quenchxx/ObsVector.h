/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <Eigen/Dense>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quenchxx/ObsSpace.h"
#include "quenchxx/VariablesSwitch.h"

namespace quenchxx {

// -----------------------------------------------------------------------------
/// ObsVector class

class ObsVector : public util::Printable,
                  private util::ObjectCounter<ObsVector> {
 public:
  static const std::string classname()
    {return "quenchxx::ObsVector";}

  explicit ObsVector(const ObsSpace & obsSpace);
  ObsVector(const ObsVector &,
            const bool copy = true);
  ~ObsVector()
    {}

  ObsVector & operator= (const ObsVector &);
  ObsVector & operator*= (const double &);
  ObsVector & operator+= (const ObsVector &);
  ObsVector & operator-= (const ObsVector &);
  ObsVector & operator*= (const ObsVector &);
  ObsVector & operator/= (const ObsVector &);

  void zero();
  void ones();
  void axpy(const double &,
            const ObsVector &);
  void invert();
  void random();
  double dot_product_with(const ObsVector &) const;
  double rms() const;
  void mask(const ObsVector &);
  void sqrt();

  size_t size() const
    {return obsSpace_.sizeGlb();}
  size_t sizeLoc() const
    {return obsSpace_.sizeLoc();}
  size_t nvars() const
    {return vars_.size();}
  void get(const size_t &,
           const size_t &,
           double &) const;
  void set(const size_t &,
           const size_t &,
           const double &);
  const double operator() (const size_t & jvar,
                           const size_t & jo) const
    {double value; this->get(jvar, jo, value); return value;}

  void read(const std::string & name)
    {data_.name() = name; obsSpace_.getdb(data_);}
  void save(const std::string & name) const
    {data_.name() = name; obsSpace_.putdb(data_);}

  Eigen::VectorXd packEigen(const ObsVector &) const;
  size_t packEigenSize(const ObsVector &) const;

 private:
  void fillHalo();
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  const ObsSpace & obsSpace_;
  const varns::Variables & vars_;
  mutable atlas::FieldSet data_;
  const double missing_;
};

//-----------------------------------------------------------------------------

}  // namespace quenchxx
