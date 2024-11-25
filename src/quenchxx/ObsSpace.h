/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/util/Point.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#ifdef ECSABER
#include "quenchxx/Variables.h"
namespace varns = quenchxx;
#else
#include "oops/base/Variables.h"
namespace varns = oops;
#endif

namespace eckit {
  class Configuration;
}

namespace quenchxx {
  class GeoVaLs;
  class ObsVector;
  class Geometry;

// -----------------------------------------------------------------------------
/// ObsSpace class

class ObsSpace : public util::Printable,
                 private util::ObjectCounter<ObsSpace> {
 public:
  static const std::string classname()
    {return "quenchxx::ObsSpace";}

  enum Write : bool {
    Screened = true,
    Original = false
  };

  ObsSpace(const eckit::Configuration &,
           const Geometry &,
           const util::DateTime &,
           const util::DateTime &,
           const bool lscreened = false);
  ~ObsSpace();

  const eckit::mpi::Comm & getComm() const
    {return comm_;}

  void putdb(const atlas::FieldSet &) const;
  void getdb(atlas::FieldSet &) const;

  std::vector<atlas::Point3> locations(const util::DateTime &,
                                       const util::DateTime &) const;
  std::vector<size_t> timeSelect(const util::DateTime &,
                                 const util::DateTime &) const;
  void generateDistribution(const eckit::Configuration &);
  void printJo(const ObsVector &,
               const ObsVector &);
  void screenObservations(const ObsVector &,
                          const GeoVaLs &) const;
  void saveObservations() const {}


  const size_t & sizeGlb() const
    {return nobsGlb_;}
  const size_t & sizeLoc() const
    {return nobsLoc_;}
  const std::vector<size_t> & sizeVec() const
    {return nobsLocVec_;}
  const std::vector<int> & order() const
    {return order_;}
  const varns::Variables & vars() const
    {return vars_;}
  std::vector<atlas::Point3> & locations() const // TODO(Benjamin): to remove
    {return locs_;}

 private:
  void print(std::ostream &) const;
  void read(const std::string &);
  void write(const std::string &,
             const bool &) const;

  const util::DateTime winbgn_;
  const util::DateTime winend_;
  const bool lscreened_;
  const eckit::mpi::Comm & comm_;
  const std::shared_ptr<const Geometry> geom_;
  mutable std::vector<util::DateTime> times_;
  mutable std::vector<atlas::Point3> locs_;
  mutable std::vector<atlas::FieldSet> data_;
  mutable std::vector<util::DateTime> screenedTimes_;
  mutable std::vector<atlas::Point3> screenedLocations_;
  mutable std::vector<atlas::FieldSet> screenedData_;
  std::string nameIn_;
  std::string nameOut_;
  size_t nobsLoc_;
  size_t nobsGlb_;
  const varns::Variables vars_;
  std::vector<size_t> nobsLocVec_;
  std::vector<int> order_;
};
// -----------------------------------------------------------------------------
}  // namespace quenchxx