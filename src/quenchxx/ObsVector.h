/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <Eigen/Dense>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quenchxx/ObsSpace.h"

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

  size_t size() const
    {return obsSpace_.sizeGlb();}
  size_t sizeLoc() const
    {return obsSpace_.sizeLoc();}
  double & operator() (const size_t ii)
    {return data_[ii];}
  const double & operator() (const size_t ii) const
    {return data_[ii];}

  void read(const std::string & name)
    {obsSpace_.getdb(name, data_);}
  void save(const std::string & name) const
    {obsSpace_.putdb(name, data_);}

  Eigen::VectorXd packEigen(const ObsVector &) const;
  size_t packEigenSize(const ObsVector &) const;

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  const ObsSpace & obsSpace_;
  std::vector<double> data_;
  const double missing_;
};

//-----------------------------------------------------------------------------

}  // namespace quenchxx
