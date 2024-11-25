/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#include "quenchxx/ObsVector.h"

#include <math.h>
#include <limits>
#include <random>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/missingValues.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

#include "quenchxx/GeoVaLs.h"

namespace quenchxx {

// -----------------------------------------------------------------------------

ObsVector::ObsVector(const ObsSpace & obsSpace)
  : comm_(obsSpace.getComm()), obsSpace_(obsSpace), vars_(obsSpace_.vars()), data_(),
    missing_(util::missingValue<double>()) {
  oops::Log::trace() << classname() << "::ObsVector starting" << std::endl;

  // Allocate FieldSet
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field(vars_[jvar].name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(obsSpace_.sizeLoc(), 1));
    data_.add(field);
  }

  // Set FieldSet to zero
  zero();

  oops::Log::trace() << classname() << "::ObsVector done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsVector::ObsVector(const ObsVector & other,
                     const bool copy)
  : comm_(other.comm_), obsSpace_(other.obsSpace_), vars_(other.vars_), data_(),
    missing_(util::missingValue<double>()) {
  oops::Log::trace() << classname() << "::ObsVector starting" << std::endl;

  util::copyFieldSet(other.data_, data_);
  if (!copy) {
    zero();
  }

  oops::Log::trace() << classname() << "::ObsVector done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator= starting" << std::endl;

  util::copyFieldSet(rhs.data_, data_);

  oops::Log::trace() << classname() << "::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator*= (const double & zz) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  util::multiplyFieldSet(data_, zz);

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator+= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator+= starting" << std::endl;

  util::addFieldSets(data_, rhs.data_);

  oops::Log::trace() << classname() << "::operator+= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator-= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator-= starting" << std::endl;

  util::subtractFieldSets(data_, rhs.data_);

  oops::Log::trace() << classname() << "::operator-= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator*= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator*= starting" << std::endl;

  util::multiplyFieldSets(data_, rhs.data_);

  oops::Log::trace() << classname() << "::operator*= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

ObsVector & ObsVector::operator/= (const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::operator/= starting" << std::endl;

  util::divideFieldSets(data_, rhs.data_);

  oops::Log::trace() << classname() << "::operator/= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsVector::zero() {
  oops::Log::trace() << classname() << "::zero starting" << std::endl;

  util::zeroFieldSet(data_);

  oops::Log::trace() << classname() << "::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::ones() {
  oops::Log::trace() << classname() << "::ones starting" << std::endl;

  for (auto & field : data_) {
    auto view = atlas::array::make_view<double, 2>(field);
    view.assign(1.0);
  }
  data_.set_dirty(false);

  oops::Log::trace() << classname() << "::ones done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::invert() {
  oops::Log::trace() << classname() << "::invert starting" << std::endl;

  atlas::FieldSet tmp;
  util::copyFieldSet(data_, tmp);
  ones();
  util::divideFieldSets(data_, tmp);

  oops::Log::trace() << classname() << "::invert done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::axpy(const double & zz,
                     const ObsVector & rhs) {
  oops::Log::trace() << classname() << "::axpy starting" << std::endl;

  atlas::FieldSet tmp;
  util::copyFieldSet(rhs.data_, tmp);
  util::multiplyFieldSet(tmp, zz);
  util::addFieldSets(data_, tmp);

  oops::Log::trace() << classname() << "::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::random() {
  oops::Log::trace() << classname() << "::random starting" << std::endl;

  // Global vector
  std::vector<double> randomPertGlb;

  if (comm_.rank() == 0) {
    // Random numbers generator
    static std::mt19937 generator(2);
    static std::normal_distribution<double> randomDistrib(0.0, 1.0);

    // Generate perturbations
    std::vector<double> randomPertTmp(vars_.size()*obsSpace_.sizeGlb());
    for (size_t jj = 0; jj < vars_.size()*obsSpace_.sizeGlb(); ++jj) {
      randomPertTmp[jj] = randomDistrib(generator);
    }

    // Reorder perturbations
    randomPertGlb.resize(vars_.size()*obsSpace_.sizeGlb());
    for (size_t jo = 0; jo < obsSpace_.sizeGlb(); ++jo) {
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        randomPertGlb[jo*vars_.size()+jvar] =
          randomPertTmp[obsSpace_.order()[jo]*vars_.size()+jvar];
      }
    }
  }

  // Define counts and displacements
  std::vector<int> dataCounts;
  std::vector<int> dataDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    dataCounts.push_back(vars_.size()*obsSpace_.sizeVec()[jt]);
  }
  dataDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    dataDispls.push_back(dataDispls[jt]+dataCounts[jt]);
  }

  // Global vector( 
  std::vector<double> randomPertLoc;

  // Scatter perturbation
  comm_.scatterv(randomPertGlb.begin(), randomPertGlb.end(), dataCounts, dataDispls,
    randomPertLoc.begin(), randomPertLoc.end(), 0);

  // Format data
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field = data_[vars_[jvar].name()];
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < obsSpace_.sizeLoc(); ++jo) {
      view(jo, 0) = randomPertLoc[jo*vars_.size()+jvar];
    }
  }

  oops::Log::trace() << classname() << "::random done" << std::endl;
}

// -----------------------------------------------------------------------------

double ObsVector::dot_product_with(const ObsVector & other) const {
  oops::Log::trace() << classname() << "::dot_product_with starting" << std::endl;

  double zz = util::dotProductFieldSets(data_, other.data_, vars_.variables(), comm_);

  oops::Log::trace() << classname() << "::dot_product_with done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

double ObsVector::rms() const {
  oops::Log::trace() << classname() << "::rms starting" << std::endl;

  double zz = util::dotProductFieldSets(data_, data_, vars_.variables(), comm_);
  if (obsSpace_.sizeGlb() > 0) {
    comm_.allReduceInPlace(zz, eckit::mpi::sum());
    zz = std::sqrt(zz/static_cast<double>(obsSpace_.sizeGlb()));
  }

  oops::Log::trace() << classname() << "::rms done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsVector::mask(const ObsVector & mask) {
  oops::Log::trace() << classname() << "::mask starting" << std::endl;

  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    const atlas::Field maskField = mask.data_[vars_[jvar].name()];
    const auto maskView = atlas::array::make_view<double, 2>(maskField);
    atlas::Field field = data_[vars_[jvar].name()];
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < obsSpace_.sizeLoc(); ++jo) {
      if (maskView(jo, 0) == missing_) view(jo, 0) = missing_;
    }
  }

  oops::Log::trace() << classname() << "::mask done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::sqrt() {
  oops::Log::trace() << classname() << "::sqrt starting" << std::endl;

  util::sqrtFieldSet(data_);

  oops::Log::trace() << classname() << "::sqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsVector::get(const size_t & jvar,
                    const size_t & ii,
                    double & value) const {
  const atlas::Field field = data_[vars_[jvar].name()];
  ASSERT(ii < field.shape(0));
  const auto view = atlas::array::make_view<double, 2>(field);
  value = view(ii, 0);
}

// -----------------------------------------------------------------------------

void ObsVector::set(const size_t & jvar,
                    const size_t & ii,
                    const double & value) {
  atlas::Field field = data_[vars_[jvar].name()];
  ASSERT(ii < field.shape(0));
  auto view = atlas::array::make_view<double, 2>(field);
  view(ii, 0) = value;
}

// -----------------------------------------------------------------------------

Eigen::VectorXd ObsVector::packEigen(const ObsVector & mask) const {
  oops::Log::trace() << classname() << "::packEigen starting" << std::endl;

  Eigen::VectorXd vec(packEigenSize(mask));
  size_t ii = 0;
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    const atlas::Field maskField = mask.data_[vars_[jvar].name()];
    const auto maskView = atlas::array::make_view<double, 2>(maskField);
    atlas::Field field = data_[vars_[jvar].name()];
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < obsSpace_.sizeLoc(); ++jo) {
      if ((view(jo, 0) != missing_) && (maskView(jo, 0) != missing_)) {
        vec(ii++) = view(jo, 0);
      }
    }
  }

  oops::Log::trace() << classname() << "::packEigen done" << std::endl;
  return vec;
}

// -----------------------------------------------------------------------------

size_t ObsVector::packEigenSize(const ObsVector & mask) const {
  oops::Log::trace() << classname() << "::packEigenSize starting" << std::endl;

  size_t ii = 0;
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    const atlas::Field maskField = mask.data_[vars_[jvar].name()];
    const auto maskView = atlas::array::make_view<double, 2>(maskField);
    atlas::Field field = data_[vars_[jvar].name()];
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < obsSpace_.sizeLoc(); ++jo) {
      if ((view(jo, 0) != missing_) && (maskView(jo, 0) != missing_)) {
        ii++;
      }
    }
  }

  oops::Log::trace() << classname() << "::packEigenSize done" << std::endl;
  return ii;
}

// -----------------------------------------------------------------------------

void ObsVector::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  std::vector<double> zmin(vars_.size(), std::numeric_limits<double>::max());
  std::vector<double> zmax(vars_.size(), -std::numeric_limits<double>::max());
  std::vector<double> zavg(vars_.size(), 0.0);

  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    const atlas::Field field = data_[vars_[jvar].name()];
    const auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < obsSpace_.sizeLoc(); ++jo) {
      if (view(jo, 0) < zmin[jvar]) zmin[jvar] = view(jo, 0);
      if (view(jo, 0) > zmax[jvar]) zmax[jvar] = view(jo, 0);
      zavg[jvar] += view(jo, 0);
    }
  }

  comm_.allReduceInPlace(zmin.begin(), zmin.end(), eckit::mpi::min());
  comm_.allReduceInPlace(zmax.begin(), zmax.end(), eckit::mpi::max());
  if (obsSpace_.sizeGlb() > 0) {
    comm_.allReduceInPlace(zavg.begin(), zavg.end(), eckit::mpi::sum());
    const double norm = 1.0/static_cast<double>(obsSpace_.sizeGlb());
    for (auto & item : zavg) {
      item *= norm;
    }
  }

  os << "quenchxx[" << obsSpace_.sizeGlb() << "]: Min=" << zmin[0] << ", Max=" << zmax[0]
    << ", Average=" << zavg[0];
  // TODO(Benjamin): print all variables

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx