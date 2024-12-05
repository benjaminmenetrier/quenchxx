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
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/VariablesSwitch.h"

namespace eckit {
  class Configuration;
}

namespace quenchxx {
  class GeoVaLs;
  class ObsVector;

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
  void saveObservations() const
    {}

  const size_t & sizeOwn() const
    {return nobsOwn_;}
  const size_t & sizeOwn(const size_t & it) const
    {return nobsOwnVec_[it];}
  const size_t & sizeLoc() const
    {return nobsLoc_;}
  const size_t & sizeGlb() const
    {return nobsGlb_;}
  const std::vector<int> & order() const
    {return order_;}
  const varns::Variables & vars() const
    {return vars_;}
  std::vector<atlas::Point3> & locations() const
    {return locs_;}
  void fillHalo(atlas::FieldSet & fset) const;

 private:
  void print(std::ostream &) const;
  void read(const std::string &);
  void write(const std::string &,
             const bool &) const;
  void fillHalo();
  template <typename T>
  void splitObservations(const atlas::grid::Distribution &,
                         const std::vector<T> &,
                         const std::vector<T> &,
                         std::vector<int> &);
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
  size_t nobsOwn_;
  size_t nobsLoc_;
  size_t nobsGlb_;
  const varns::Variables vars_;
  std::vector<size_t> nobsOwnVec_;
  std::vector<int> order_;
  eckit::LocalConfiguration distribution_;
  std::vector<int> sendBufIndex_;
  size_t nSend_;
  size_t nRecv_;
  std::vector<int> dataSendCounts_;
  std::vector<int> dataRecvCounts_;
  std::vector<int> dataSendDispls_;
  std::vector<int> dataRecvDispls_;
};

// -----------------------------------------------------------------------------

template <typename T>
void ObsSpace::splitObservations(const atlas::grid::Distribution & distribution,
                                 const std::vector<T> & longitude,
                                 const std::vector<T> & latitude,
                                 std::vector<int> & partition) {
  oops::Log::trace() << classname() << "::splitObservations starting" << std::endl;

  // TODO(Benjamin): other distributions
  if (true) {
    // Using nearest neighbor
    atlas::util::IndexKDTree search;
    search.reserve(geom_->grid().size());
    size_t jnode = 0;
    for (const auto & lonLat : geom_->grid().lonlat()) {
      atlas::PointLonLat pointLonLat(lonLat);
      pointLonLat.normalise();
      atlas::PointXY point(pointLonLat);
      search.insert(point, jnode);
      ++jnode;
    }
    search.build();

    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      // Find MPI task
      atlas::PointLonLat pointLonLat(longitude[jo], latitude[jo]);
      pointLonLat.normalise();

      // Search nearest neighborass
      atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);

      // Define partition
      partition[jo] = distribution.partition(neighbor[0].payload());
      ++nobsOwnVec_[partition[jo]];
    }
  }

  oops::Log::trace() << classname() << "::splitObservations done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
