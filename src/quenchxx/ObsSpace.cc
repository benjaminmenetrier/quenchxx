/*
 * (C) Copyright 2023 Meteorologisk Institutt
 *
 */

#include "quenchxx/ObsSpace.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/KDTree.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "quenchxx/Geometry.h"
#include "quenchxx/GeoVaLs.h"
#include "quenchxx/ObsVector.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace quenchxx {

// -----------------------------------------------------------------------------

ObsSpace::ObsSpace(const eckit::Configuration & config,
                   const Geometry & geom,
                   const util::DateTime & bgn,
                   const util::DateTime & end,
                   const bool lscreened)
  : winbgn_(bgn), winend_(end), lscreened_(lscreened), comm_(geom.getComm()),
    geom_(new Geometry(geom)), nobsLoc_(0), nobsGlb_(0),
    vars_(config.getStringVector("variables")) {
  oops::Log::trace() << classname() << "::ObsSpace starting" << std::endl;

  nameIn_.clear();
  nameOut_.clear();
  if (config.has("ObsData")) {
    const eckit::LocalConfiguration dataConfig(config, "ObsData");
    if (lscreened) {
      if (dataConfig.has("ObsDataInScreened")) {
        nameIn_ = dataConfig.getString("ObsDataInScreened.filepath");
        oops::Log::trace() << classname() << "::ObsSpace reading screened observations from "
          << nameIn_ << std::endl;
        read(nameIn_);
      }
      if (dataConfig.has("ObsDataOutScreened")) {
        nameOut_ = dataConfig.getString("ObsDataOutScreened.filepath");
      }
    } else {
      if (dataConfig.has("ObsDataIn")) {
        nameIn_ = dataConfig.getString("ObsDataIn.filepath");
        oops::Log::trace() << classname() << "::ObsSpace reading observations from " << nameIn_
          << std::endl;
        read(nameIn_);
      }
      if (dataConfig.has("ObsDataOut")) {
        nameOut_ = dataConfig.getString("ObsDataOut.filepath");
      }
    }
  }

  oops::Log::trace() << classname() << "::ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSpace::~ObsSpace() {
  oops::Log::trace() << classname() << "::~ObsSpace starting" << std::endl;

  if (!nameOut_.empty()) {
    oops::Log::trace() << classname() << "::~ObsSpace saving nameOut = " << nameOut_ << std::endl;
    write(nameOut_, Write::Original);
  }

  oops::Log::trace() << classname() << "::~ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::putdb(const atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::putdb starting" << std::endl;

  // Check Fields size
  for (const auto & field : fset) {
    ASSERT(field.shape(0) == nobsLoc_);
    ASSERT(field.shape(1) == 1);
  }

  bool existingFieldSet = false;
  for (auto & dataFset : data_) {
    if (dataFset.name() == fset.name()) {
      // Existing FieldSet, clear and reset fields
      existingFieldSet = true;
      dataFset.clear();
      for (const auto & field : fset) {
        dataFset.add(field);
      }
    }
  }

  if (!existingFieldSet) {
    // FieldSet name not found, inserting new FieldSet
    data_.push_back(fset);
  }

  oops::Log::trace() << classname() << "::putdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::getdb(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::getdb starting" << std::endl;

  fset.clear();
  bool existingFieldSet = false;
  for (const auto & dataFset : data_) {
    if (dataFset.name() == fset.name()) {
      // FieldSet found, add fields
      existingFieldSet = true;
      for (const auto & field : dataFset) {
        fset.add(field);
      }
    }
  }

  if (!existingFieldSet) {
    // FieldSet name not found
    std::string message = "Cannot find group " + fset.name() + " in observation database,"
      + "existing groups are: ";
    for (const auto & dataFset : data_) {
      message += dataFset.name() + " ";
    }
    throw eckit::Exception(message, Here());
  }

  oops::Log::trace() << classname() << "::getdb done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<atlas::Point3> ObsSpace::locations(const util::DateTime & t1,
                                               const util::DateTime & t2) const {
  oops::Log::trace() << classname() << "::locations starting" << std::endl;

  std::vector<size_t> olist = timeSelect(t1, t2);
  const size_t nobs = olist.size();
  std::vector<atlas::Point3> locs(nobs);
  for (size_t jo = 0; jo < nobs; ++jo) {
    locs[jo] = locs_[olist[jo]];
  }

  oops::Log::trace() << classname() << "::locations done" << std::endl;
  return locs;
}

// -----------------------------------------------------------------------------

std::vector<size_t> ObsSpace::timeSelect(const util::DateTime & t1,
                                         const util::DateTime & t2) const {
  oops::Log::trace() << classname() << "::timeSelect starting" << std::endl;

  std::vector<size_t> mask;
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    if (t1 == t2) {
      if (times_[jo] == t1) {
        mask.push_back(jo);
      }
    } else {
      if (times_[jo] > t1 && times_[jo] <= t2) {
        mask.push_back(jo);
      }
    }
  }

  oops::Log::trace() << classname() << "::timeSelect done" << std::endl;
  return mask;
}

// -----------------------------------------------------------------------------

void ObsSpace::generateDistribution(const eckit::Configuration & config) {
  oops::Log::trace() << classname() << "::generateDistribution starting" << std::endl;

  // No forecast: assimilation window start time must be equal to final time
  ASSERT(winend_ == winbgn_);

  // Parameters
  nobsGlb_ = config.getInt("density");
  const std::vector<double> vert_coord_avg = geom_->vert_coord_avg(config.getString("variable"));
  const double bottom = *std::min_element(vert_coord_avg.begin(), vert_coord_avg.end());
  const double top = *std::max_element(vert_coord_avg.begin(), vert_coord_avg.end());

  // Random number generator
  static std::mt19937 generator(winbgn_.timestamp());
  static std::uniform_real_distribution<double> randomDistrib(0.0, 1.0);

  // Global vectors
  std::vector<std::string> colNames;
  std::vector<double> locsGlb;

  // Define source grid horizontal distribution
  const atlas::grid::Distribution distribution = geom_->partitioner().partition(geom_->grid());

  // Local number of observations
  nobsLocVec_.resize(comm_.size());

  if (comm_.rank() == 0) {
    // Resize vectors
    std::vector<double> locsGlbTmp(3*nobsGlb_);
    std::vector<size_t> orderGlbTmp(nobsGlb_);

    // Generate locations
    size_t iobs = 0;
    while (iobs < nobsGlb_) {
      // Generate location
      const double lon = 360.0*randomDistrib(generator);
      const double lat = 90.0-std::acos(1.0-2.0*randomDistrib(generator))*180.0/M_PI;

      // Check point validity
      bool validPoint = true;
      if (!geom_->grid().domain().global()) {
        const atlas::functionspace::StructuredColumns fs(geom_->functionSpace());
        const atlas::RegularGrid grid(fs.grid());
        atlas::PointLonLat p({lon, lat});
        grid.projection().lonlat2xy(p);
        validPoint = grid.domain().contains(p);
      }

      // Process valid point
      if (validPoint) {
        locsGlbTmp[3*iobs] = lon;
        locsGlbTmp[3*iobs+1] = lat;
        locsGlbTmp[3*iobs+2] = bottom+(top-bottom)*randomDistrib(generator);
        orderGlbTmp[iobs] = iobs;
        ++iobs;
      }
    }

    // Split observations between tasks
    std::vector<int> partition(nobsGlb_);
    if (true) {  // TODO(Benjamin): other distributions
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
        atlas::PointLonLat pointLonLat(locsGlbTmp[3*jo], locsGlbTmp[3*jo+1]);
        pointLonLat.normalise();

        // Search nearest neighbor
        atlas::util::IndexKDTree::ValueList neighbor = search.closestPoints(pointLonLat, 1);

        // Define partition
        partition[jo] = distribution.partition(neighbor[0].payload());
        ++nobsLocVec_[partition[jo]];
      }
    }

    // Reorder vectors
    locsGlb.resize(3*nobsGlb_);
    order_.resize(nobsGlb_);
    std::vector<int> iobsLocVec(comm_.size(), 0);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      size_t offset = 0;
      for (int jt = 0; jt < partition[jo]; ++jt) {
        offset += nobsLocVec_[jt];
      }
      offset += iobsLocVec[partition[jo]];
      locsGlb[3*offset+0] = locsGlbTmp[3*jo+0];
      locsGlb[3*offset+1] = locsGlbTmp[3*jo+1];
      locsGlb[3*offset+2] = locsGlbTmp[3*jo+2];
      order_[offset] = orderGlbTmp[jo];
      ++iobsLocVec[partition[jo]];
    }
  }

  // Broadcast number of observations
  comm_.broadcast(nobsLocVec_, 0);
  nobsLoc_ = nobsLocVec_[comm_.rank()];

  // Allocation of local vector
  std::vector<double> locsLoc(3*nobsLoc_);

  // Define counts and displacements
  std::vector<int> locsCounts;
  std::vector<int> locsDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    locsCounts.push_back(3*nobsLocVec_[jt]);
  }
  locsDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    locsDispls.push_back(locsDispls[jt]+locsCounts[jt]);
  }

  // Scatter locations
  comm_.scatterv(locsGlb.begin(), locsGlb.end(), locsCounts, locsDispls,
    locsLoc.begin(), locsLoc.end(), 0);

  // Format data
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    locs_.push_back(atlas::Point3(locsLoc[3*jo], locsLoc[3*jo+1], locsLoc[3*jo+2]));
  }

  // Generate times
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    times_.push_back(winbgn_);
  }

  // Generate observations error
  const std::vector<double> err = config.getDoubleVector("error");
  atlas::FieldSet fset;
  fset.name() = config.getString("obserror");
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    atlas::Field field(vars_[jvar].name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nobsLoc_, 1));
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t jo = 0; jo < nobsLoc_; ++jo) {
      view(jo, 0) = err[jvar];
    }
    fset.add(field);
  }
  this->putdb(fset);

  oops::Log::trace() << classname() << "::generateDistribution done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::printJo(const ObsVector &,
                       const ObsVector &) {
  oops::Log::trace() << classname() << "::printJo starting" << std::endl;

  oops::Log::info() << "ObsSpace::printJo not implemented" << std::endl;

  oops::Log::trace() << classname() << "::printJo starting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::screenObservations(const ObsVector & dep,
                                  const GeoVaLs & gv) const {
  oops::Log::trace() << classname() << "::screenObservations starting" << std::endl;

  // Define whether observations are screened or not
  // TODO(Benjamin): more clever screening
  std::vector<bool> validObs(nobsLoc_);
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    validObs[jo] = true;
  }
  const int nobsLocScr = std::count(validObs.cbegin(), validObs.cend(), true);

  // Apply screening on time, location and columns
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    if (validObs[jo]) {
      screenedTimes_.push_back(times_[jo]);
      screenedLocations_.push_back(locs_[jo]);
    }
  }
  for (const auto & dataFset : data_) {
    atlas::FieldSet fset;
    fset.name() = dataFset.name();
    for (const auto & dataField : dataFset) {
      atlas::Field field(dataField.name(), atlas::array::make_datatype<double>(),
        atlas::array::make_shape(nobsLocScr, 1));
      const auto dataView = atlas::array::make_view<double, 2>(dataField);
      auto view = atlas::array::make_view<double, 2>(field);
      size_t joScr = 0;
      for (size_t jo = 0; jo < nobsLoc_; ++jo) {
        if (validObs[jo]) {
          view(joScr, 0) = dataView(jo, 0);
          ++joScr;
        }
      }
      fset.add(field);
    }
    screenedData_.push_back(fset);
  }

  // Write screened observations
  write(nameOut_, Write::Screened);

  oops::Log::trace() << classname() << "::screenObservations done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::read(const std::string & filePath) {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Global size and vectors
  int ngrp;
  std::vector<int> partition;
  std::vector<int> timesGlb;
  std::vector<double> locsGlb;
  std::vector<double> dataGlb;

  // Counts and displacements
  nobsLocVec_.resize(comm_.size());

  // Define source grid horizontal distribution
  const atlas::grid::Distribution distribution = geom_->partitioner().partition(geom_->grid());

  // NetCDF IDs
  int retval, ncid, nobs_id, dateTime_id, order_id, longitude_id, latitude_id, height_id, data_id;
  std::vector<int> group_ids;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    std::string ncFilePath = filePath + ".nc";
    oops::Log::info() << "Reading file: " << ncFilePath << std::endl;
    if (retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid)) ERR(retval);

    // Get dimension
    if (retval = nc_inq_dimid(ncid, "Location", &nobs_id)) ERR(retval);
    if (retval = nc_inq_dimlen(ncid, nobs_id, &nobsGlb_)) ERR(retval);

    // Get groups list
    if (retval = nc_inq_grps(ncid, &ngrp, NULL)) ERR(retval);
    group_ids.resize(ngrp);
    if (retval = nc_inq_grps(ncid, NULL, group_ids.data())) ERR(retval);

    // Get order if available
    order_.resize(nobsGlb_);
    retval = nc_inq_varid(ncid, "order", &order_id);
    if (retval == NC_NOERR) {
      if (retval = nc_get_var_int(ncid, order_id, order_.data())) ERR(retval);
    } else {
      for (size_t jo = 0; jo < nobsGlb_; ++jo) {
        order_[jo] = jo;
      }
    }

    // Get MetaData
    int meta_group_id;
    if (retval = nc_inq_grp_ncid(ncid, "MetaData", &meta_group_id)) ERR(retval);

    // Get dateTime from MetaData
    std::vector<int64_t> dateTime(nobsGlb_);
    if (retval = nc_inq_varid(meta_group_id, "dateTime", &dateTime_id)) ERR(retval);
    if (retval = nc_get_var_long(meta_group_id, dateTime_id, dateTime.data())) ERR(retval);
    const std::string dateTime_units_key = "units";
    size_t attlen;
    if (retval = nc_inq_attlen(meta_group_id, dateTime_id, dateTime_units_key.c_str(),
      &attlen)) ERR(retval);
    char **dateTime_units_char = reinterpret_cast<char**>(malloc(attlen*sizeof(char*)));
    memset(dateTime_units_char, 0, attlen*sizeof(char*));
    if (retval = nc_get_att_string(meta_group_id, dateTime_id, dateTime_units_key.c_str(),
      dateTime_units_char)) ERR(retval);
    std::string dateTime_units_value(*dateTime_units_char);
    util::DateTime start(dateTime_units_value.substr(14, 20));

    // Get longitude, latitude and height from MetaData
    std::vector<float> longitude(nobsGlb_);
    std::vector<float> latitude(nobsGlb_);
    std::vector<float> height(nobsGlb_);
    if (retval = nc_inq_varid(meta_group_id, "longitude", &longitude_id)) ERR(retval);
    if (retval = nc_get_var_float(meta_group_id, longitude_id, longitude.data()))
      ERR(retval);
    if (retval = nc_inq_varid(meta_group_id, "latitude", &latitude_id)) ERR(retval);
    if (retval = nc_get_var_float(meta_group_id, latitude_id, latitude.data())) ERR(retval);
    if (retval = nc_inq_varid(meta_group_id, "height", &height_id)) ERR(retval);
    if (retval = nc_get_var_float(meta_group_id, height_id, height.data())) ERR(retval);

    // Split observations between tasks
    // TODO(Benjamin): other distributions
    partition.resize(nobsGlb_);
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
        ++nobsLocVec_[partition[jo]];
      }
    }

    // Get ordered time and location
    timesGlb.resize(6*nobsGlb_);
    locsGlb.resize(3*nobsGlb_);
    std::vector<int> iobsLocVec(comm_.size(), 0);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      size_t offset = 0;
      for (int jt = 0; jt < partition[jo]; ++jt) {
        offset += nobsLocVec_[jt];
      }
      offset += iobsLocVec[partition[jo]];
      util::DateTime startTest = start + util::Duration(dateTime[jo]);
      startTest.toYYYYMMDDhhmmss(timesGlb[6*offset+0], timesGlb[6*offset+1], timesGlb[6*offset+2],
        timesGlb[6*offset+3], timesGlb[6*offset+4], timesGlb[6*offset+5]);
      locsGlb[3*offset+0] = static_cast<double>(longitude[jo]);
      locsGlb[3*offset+1] = static_cast<double>(latitude[jo]);
      locsGlb[3*offset+2] = static_cast<double>(height[jo]);
      ++iobsLocVec[partition[jo]];
    }
  }

  // Broadcast number of observations for each task
  comm_.broadcast(nobsLocVec_, 0);
  nobsLoc_ = nobsLocVec_[comm_.rank()];
  nobsGlb_ = 0;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    nobsGlb_ += nobsLocVec_[jt];
  }

  // Allocation of local vectors
  std::vector<int> timesLoc(6*nobsLoc_);
  std::vector<double> locsLoc(3*nobsLoc_);

  // Define counts and displacements
  std::vector<int> timesCounts;
  std::vector<int> locsCounts;
  std::vector<int> dataCounts;
  std::vector<int> timesDispls;
  std::vector<int> locsDispls;
  std::vector<int> dataDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    timesCounts.push_back(6*nobsLocVec_[jt]);
    locsCounts.push_back(3*nobsLocVec_[jt]);
    dataCounts.push_back(vars_.size()*nobsLocVec_[jt]);
  }
  timesDispls.push_back(0);
  locsDispls.push_back(0);
  dataDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    timesDispls.push_back(timesDispls[jt]+timesCounts[jt]);
    locsDispls.push_back(locsDispls[jt]+locsCounts[jt]);
    dataDispls.push_back(dataDispls[jt]+dataCounts[jt]);
  }

  // Scatter times and locations
  comm_.scatterv(timesGlb.begin(), timesGlb.end(), timesCounts, timesDispls,
    timesLoc.begin(), timesLoc.end(), 0);
  comm_.scatterv(locsGlb.begin(), locsGlb.end(), locsCounts, locsDispls,
    locsLoc.begin(), locsLoc.end(), 0);

  // Format local times and locations
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    times_.push_back(util::DateTime(timesLoc[6*jo], timesLoc[6*jo+1], timesLoc[6*jo+2],
      timesLoc[6*jo+3], timesLoc[6*jo+4], timesLoc[6*jo+5]));
    locs_.push_back(atlas::Point3(locsLoc[3*jo], locsLoc[3*jo+1], locsLoc[3*jo+2]));
  }

  // Broadcast number of groups
  comm_.broadcast(ngrp, 0);

  for (int jgrp = 0; jgrp < ngrp; ++jgrp) {
    // Get group name
    std::string grpName;
    size_t grpNameLen = 0;
    if (comm_.rank() == 0) {
      if (retval = nc_inq_grpname_len(group_ids[jgrp], &grpNameLen)) ERR(retval);
      grpName.resize(grpNameLen);
      if (retval = nc_inq_grpname(group_ids[jgrp], &grpName[0])) ERR(retval);
    }
    oops::mpi::broadcastString(comm_, grpName, 0);

    if (grpName.substr(0, grpNameLen-1) != "MetaData") {
      // Non-MetaData group
      if (comm_.rank() == 0) {
        std::vector<float> dataVar(nobsGlb_);
        dataGlb.resize(vars_.size()*nobsGlb_);
        for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
          // Get data
          if (retval = nc_inq_varid(group_ids[jgrp], vars_[jvar].name().c_str(), &data_id))
            ERR(retval);
          if (retval = nc_get_var_float(group_ids[jgrp], data_id, dataVar.data())) ERR(retval);

          // Get ordered data
          std::vector<int> iobsLocVec(comm_.size(), 0);
          for (size_t jo = 0; jo < nobsGlb_; ++jo) {
            size_t offset = 0;
            for (int jt = 0; jt < partition[jo]; ++jt) {
              offset += nobsLocVec_[jt];
            }
            offset += iobsLocVec[partition[jo]];
            dataGlb[vars_.size()*offset+jvar] = static_cast<double>(dataVar[jo]);
            ++iobsLocVec[partition[jo]];
          }
        }
      }
    
      // Allocation of local vector
      std::vector<double> dataLoc(vars_.size()*nobsLoc_);

      // Scatter data
      comm_.scatterv(dataGlb.begin(), dataGlb.end(), dataCounts, dataDispls,
        dataLoc.begin(), dataLoc.end(), 0);

      // Format data
      atlas::FieldSet fset;
      fset.name() = grpName.substr(0, grpNameLen-1);
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        atlas::Field field(vars_[jvar].name(), atlas::array::make_datatype<double>(),
          atlas::array::make_shape(nobsLoc_, 1));
        auto view = atlas::array::make_view<double, 2>(field);
        for (size_t jo = 0; jo < nobsLoc_; ++jo) {
          view(jo, 0) = dataLoc[vars_.size()*jo+jvar];
        }
        fset.add(field);
      }
      data_.push_back(fset);
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if (retval = nc_close(ncid)) ERR(retval);
  }

  // Add one second if time_ = winbgn_ and  winend_ > winbgn_
  if (winend_ > winbgn_) {
    for (auto & time : times_) {
      if (time == winbgn_) {
        time += util::Duration("PT1S");
      }
    }
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::write(const std::string & filePath,
                     const bool & writeScreened) const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // TODO NOW

/*
  std::vector<util::DateTime> * times = 0;
  std::vector<atlas::Point3> * locs = 0;
  std::vector<atlas::FieldSet> * data = 0;

  // Set pointers
  if (writeScreened) {
    times = &screenedTimes_;
    locs = &screenedLocations_;
    data = &screenedData_;
  } else {
    times = &times_;
    locs = &locs_;
    data = &data_;
  }
  ASSERT(times);
  ASSERT(locs);
  ASSERT(data);

  // Number of columns
  const size_t ncol = data->size();

  // Global vectors
  std::vector<std::string> colNames;
  std::vector<int> timesGlb;
  std::vector<double> locsGlb;
  std::vector<double> colsGlb;

  // Allocation of local vectors
  std::vector<int> timesLoc(6*nobsLoc_);
  std::vector<double> locsLoc(3*nobsLoc_);
  std::vector<double> colsLoc(ncol*nobsLoc_);

  // Format data
  for (size_t jo = 0; jo < nobsLoc_; ++jo) {
    (*times)[jo].toYYYYMMDDhhmmss(timesLoc[6*jo], timesLoc[6*jo+1], timesLoc[6*jo+2],
      timesLoc[6*jo+3], timesLoc[6*jo+4], timesLoc[6*jo+5]);
    locsLoc[3*jo+0] = (*locs)[jo][0];
    locsLoc[3*jo+1] = (*locs)[jo][1];
    locsLoc[3*jo+2] = (*locs)[jo][2];
    size_t jc = 0;
    for (auto const & vec : *data) {
      for (size_t jo = 0; jo < nobsLoc_; ++jo) {
        colsLoc[ncol*jo+jc] = vec.second[jo];
      }
      ++jc;
    }
  }

  // Global size check
  ASSERT(nobsLocVec_.size() == comm_.size());

  if (comm_.rank() == 0) {
    // Resize vectors
    timesGlb.resize(6*nobsGlb_);
    locsGlb.resize(3*nobsGlb_);
    colsGlb.resize(ncol*nobsGlb_);
  }

  // Define counts and displacements
  std::vector<int> timesCounts;
  std::vector<int> locsCounts;
  std::vector<int> colsCounts;
  std::vector<int> timesDispls;
  std::vector<int> locsDispls;
  std::vector<int> colsDispls;
  for (size_t jt = 0; jt < comm_.size(); ++jt) {
    timesCounts.push_back(6*nobsLocVec_[jt]);
    locsCounts.push_back(3*nobsLocVec_[jt]);
    colsCounts.push_back(ncol*nobsLocVec_[jt]);
  }
  timesDispls.push_back(0);
  locsDispls.push_back(0);
  colsDispls.push_back(0);
  for (size_t jt = 0; jt < comm_.size()-1; ++jt) {
    timesDispls.push_back(timesDispls[jt]+timesCounts[jt]);
    locsDispls.push_back(locsDispls[jt]+locsCounts[jt]);
    colsDispls.push_back(colsDispls[jt]+colsCounts[jt]);
  }

  // Gather times, locations and data
  comm_.gatherv(timesLoc, timesGlb, timesCounts, timesDispls, 0);
  comm_.gatherv(locsLoc, locsGlb, locsCounts, locsDispls, 0);
  comm_.gatherv(colsLoc, colsGlb, colsCounts, colsDispls, 0);

  if (comm_.rank() == 0) {
    // Define locations vector
    std::vector<int> Locations(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      Locations[jo] = jo;
    }

    // Find start dateTime
    util::DateTime start(timesGlb[0], timesGlb[1], timesGlb[2], timesGlb[3], timesGlb[4],
      timesGlb[5]);
    for (size_t jo = 1; jo < nobsGlb_; ++jo) {
      util::DateTime startTest(timesGlb[6*jo], timesGlb[6*jo+1], timesGlb[6*jo+2], timesGlb[6*jo+3],
        timesGlb[6*jo+4], timesGlb[6*jo+5]);
      if (startTest < start) {
        start = startTest;
      }
    }

    // Define dateTime
    std::vector<int64_t> dateTime(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      util::DateTime startTest(timesGlb[6*jo], timesGlb[6*jo+1], timesGlb[6*jo+2], timesGlb[6*jo+3],
        timesGlb[6*jo+4], timesGlb[6*jo+5]);
      dateTime[jo] = (startTest-start).toSeconds();
    }

    // Define longitude
    std::vector<float> longitude(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      longitude[jo] = locsGlb[3*jo];
    }

    // Define latitude
    std::vector<float> latitude(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      latitude[jo] = locsGlb[3*jo+1];
    }

    // Define height
    std::vector<float> height(nobsGlb_);
    for (size_t jo = 0; jo < nobsGlb_; ++jo) {
      height[jo] = locsGlb[3*jo+2];
    }

    // NetCDF IDs
    int retval, ncid, nobs_id, d_id[1], Locations_id, order_id, metaData_id, dateTime_id,
      longitude_id, latitude_id, height_id;
    std::vector<int> group_ids;
    std::vector<int> groupVar_ids;

    // Create NetCDF file
    const std::string ncFilePath = writeScreened ? filePath + "_screened.nc" : filePath + ".nc";
    if (retval = nc_create(ncFilePath.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid)) ERR(retval);

    // Global attributes
    const std::string ioda_layout_key = "_ioda_layout";
    const std::string ioda_layout_value = "ObsGroup";
    const char *ioda_layout_char[1] = {ioda_layout_value.c_str()};
    if (retval = nc_put_att_string(ncid, NC_GLOBAL, ioda_layout_key.c_str(), 1, ioda_layout_char))
      ERR(retval);
    const std::string ioda_layout_version_key = "_ioda_layout_version";
    const int64_t ioda_layout_version_value = 0;
    if (retval = nc_put_att_long(ncid, NC_GLOBAL, ioda_layout_version_key.c_str(), NC_INT64, 1,
      &ioda_layout_version_value)) ERR(retval);

    // Create dimension
    if (retval = nc_def_dim(ncid, "Location", NC_UNLIMITED, &nobs_id)) ERR(retval);
    d_id[0] = nobs_id;

    // Missing value
    const std::string fillValue_key = "_FillValue";

    // Define global variables
    if (retval = nc_def_var(ncid, "Location", NC_INT, 1, d_id, &Locations_id)) ERR(retval);
    if (retval = nc_put_att_int(ncid, Locations_id, fillValue_key.c_str(), NC_INT, 1,
      &util::missingValue<int>())) ERR(retval);
    if (retval = nc_def_var(ncid, "order", NC_INT, 1, d_id, &order_id)) ERR(retval);
    if (retval = nc_put_att_int(ncid, order_id, fillValue_key.c_str(), NC_INT, 1,
      &util::missingValue<int>())) ERR(retval);

    // Define metadata group
    if (retval = nc_def_grp(ncid, "MetaData", &metaData_id)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "dateTime", NC_INT64, 1, d_id, &dateTime_id))
      ERR(retval);
    if (retval = nc_put_att_long(metaData_id, dateTime_id, fillValue_key.c_str(), NC_INT64, 1,
      &util::missingValue<int64_t>())) ERR(retval);
    const std::string dateTime_units_key = "units";
    const std::string dateTime_units_value = "seconds since " + start.toString();
    const char *dateTime_units_char[1] = {dateTime_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, dateTime_id, dateTime_units_key.c_str(),
      1, dateTime_units_char)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "longitude", NC_FLOAT, 1, d_id, &longitude_id))
      ERR(retval);
    if (retval = nc_put_att_float(metaData_id, longitude_id, fillValue_key.c_str(), NC_FLOAT, 1,
      &util::missingValue<float>())) ERR(retval);
    const std::string longitude_units_key = "units";
    const std::string longitude_units_value = "degrees_east";
    const char *longitude_units_char[1] = {longitude_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, longitude_id, longitude_units_key.c_str(), 1,
      longitude_units_char)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "latitude", NC_FLOAT, 1, d_id, &latitude_id))
      ERR(retval);
    if (retval = nc_put_att_float(metaData_id, latitude_id, fillValue_key.c_str(), NC_FLOAT, 1,
      &util::missingValue<float>())) ERR(retval);
    const std::string latitude_units_key = "units";
    const std::string latitude_units_value = "degrees_north";
    const char *latitude_units_char[1] = {latitude_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, latitude_id, latitude_units_key.c_str(), 1,
      latitude_units_char)) ERR(retval);
    if (retval = nc_def_var(metaData_id, "height", NC_FLOAT, 1, d_id, &height_id))
      ERR(retval);
    if (retval = nc_put_att_float(metaData_id, height_id, fillValue_key.c_str(), NC_FLOAT, 1,
      &util::missingValue<float>())) ERR(retval);
    const std::string height_units_key = "units";
    const std::string height_units_value = "m";
    const char *height_units_char[1] = {height_units_value.c_str()};
    if (retval = nc_put_att_string(metaData_id, height_id, height_units_key.c_str(), 1,
      height_units_char)) ERR(retval);

    // Define data groups
    for (auto const & vec : *data) {
      int group_id;
      if (retval = nc_def_grp(ncid, vec.first.c_str(), &group_id)) ERR(retval);
      group_ids.push_back(group_id);
      int groupVar_id;
      if (retval = nc_def_var(group_id, vars_[0].name().c_str(), NC_FLOAT, 1, d_id, &groupVar_id))
        ERR(retval);
      if (retval = nc_put_att_float(group_id, groupVar_id, fillValue_key.c_str(), NC_FLOAT, 1,
        &util::missingValue<float>())) ERR(retval);
      groupVar_ids.push_back(groupVar_id);
    }

    // End definition mode
    if (retval = nc_enddef(ncid)) ERR(retval);

    // Write data
    const size_t init = 0;
    const size_t nobs = Locations.size();
    if (retval = nc_put_vara_int(ncid, Locations_id, &init, &nobs, Locations.data())) ERR(retval);
    if (retval = nc_put_vara_int(ncid, order_id, &init, &nobs, order_.data())) ERR(retval);
    if (retval = nc_put_vara_long(metaData_id, dateTime_id, &init, &nobs, dateTime.data()))
      ERR(retval);
    if (retval = nc_put_vara_float(metaData_id, longitude_id, &init, &nobs, longitude.data()))
      ERR(retval);
    if (retval = nc_put_vara_float(metaData_id, latitude_id, &init, &nobs, latitude.data()))
      ERR(retval);
    if (retval = nc_put_vara_float(metaData_id, height_id, &init, &nobs, height.data()))
      ERR(retval);
    size_t ncol = group_ids.size();
    for (size_t jc = 0; jc < ncol; ++jc) {
      std::vector<float> vec(nobsGlb_);
      for (size_t jo = 0; jo < nobsGlb_; ++jo) {
        vec[jo] = static_cast<float>(colsGlb[jo*ncol+jc]);
      }
      if (retval = nc_put_var_float(group_ids[jc], groupVar_ids[jc], vec.data())) ERR(retval);
    }

    // Close file
    if (retval = nc_close(ncid)) ERR(retval);
  }

  // Reset pointers
  times = 0;
  locs = 0;
  data = 0;
*/
  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpace::print(std::ostream & os) const {
  oops::Log::trace() << classname() << "::print starting" << std::endl;

  os << "ObsSpace: assimilation window = " << winbgn_ << " to " << winend_ << std::endl;
  os << "ObsSpace: file in = " << nameIn_ << ", file out = " << nameOut_ << std::endl;

  oops::Log::trace() << classname() << "::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
