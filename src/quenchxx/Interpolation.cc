/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "quenchxx/Interpolation.h"

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

#include "quenchxx/VariablesSwitch.h"

// -----------------------------------------------------------------------------

namespace quenchxx {

// -----------------------------------------------------------------------------

Interpolation::Interpolation(const Geometry & geom,
                             const std::string & srcUid,
                             const atlas::Grid & dstGrid,
                             const atlas::FunctionSpace & dstFspace,
                             const std::string & dstUid)
  : srcUid_(srcUid), dstUid_(dstUid), dstFspace_(dstFspace) {
  oops::Log::trace() << classname() << "::Interpolation starting" << std::endl;

  // Get interpolation type
  const std::string type = geom.interpolation().getString("interpolation type");

  // Setup interpolation
  if (type == "atlas interpolation wrapper") {
    atlasInterpWrapper_ = std::make_shared<saber::interpolation::AtlasInterpWrapper>(
      geom.partitioner(), geom.functionSpace(), dstGrid, dstFspace);
  } else if (type == "regional") {
    regionalInterp_ = std::make_shared<atlas::Interpolation>(
      atlas::util::Config("type", "regional-linear-2d"),
      geom.functionSpace(), dstFspace);
  } else if (type == "unstructured") {
    // Get ghost points
    const auto ghostView = atlas::array::make_view<int, 1>(dstFspace.ghost());
    for (atlas::idx_t jnode = 0; jnode < dstFspace.ghost().shape(0); ++jnode) {
      ghostVector_.push_back(ghostView(jnode));
    }

    // Get longitudes/latitudes
    std::vector<double> lons;
    std::vector<double> lats;
    const auto lonLatField = dstFspace.lonlat();
    const auto lonLatView = atlas::array::make_view<double, 2>(lonLatField);
    for (atlas::idx_t jnode = 0; jnode < lonLatField.shape(0); ++jnode) {
      ghostVector_.push_back(ghostView(jnode));
      if (ghostVector_[jnode] == 0) {
        lons.push_back(lonLatView(jnode, 0));
        lats.push_back(lonLatView(jnode, 1));
      }
    }

    // Setup unstructured interpolator
    unstructuredInterp_ = std::make_shared<oops::UnstructuredInterpolator>(geom.interpolation(),
      geom.generic(), lats, lons);
  } else {
    throw eckit::Exception("wrong interpolation type", Here());
  }

  oops::Log::trace() << classname() << "::Interpolation done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::execute(const atlas::FieldSet & srcFieldSet,
                            atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;

  if (atlasInterpWrapper_) {
    atlasInterpWrapper_->execute(srcFieldSet, tgtFieldSet);
  }
  if (regionalInterp_) {
    regionalInterp_->execute(srcFieldSet, tgtFieldSet);
  }
  if (unstructuredInterp_) {
    // Exchange FieldSet halo
    atlas::FieldSet fset = util::copyFieldSet(srcFieldSet);
    fset.haloExchange();

    // Apply unstructured interpolator
    const varns::Variables vars(fset.field_names());
    std::vector<double> vals;
    unstructuredInterp_->apply(vars, fset, vals);

    // Format data
    size_t index = 0;
    for (auto & tgtField : tgtFieldSet) {
      auto tgtView = atlas::array::make_view<double, 2>(tgtField);
      for (atlas::idx_t jlevel = 0; jlevel < tgtView.shape(1); ++jlevel) {
        for (atlas::idx_t jnode = 0; jnode < tgtView.shape(0); ++jnode) {
          if (ghostVector_[jnode] == 0) {
            tgtView(jnode, jlevel) = vals[index];
            ++index;
          }
        }
      }
    }
    tgtFieldSet.haloExchange();
  }


  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeAdjoint(atlas::FieldSet & srcFieldSet,
                                   const atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::executeAdjoint starting" << std::endl;

  if (atlasInterpWrapper_) {
    atlasInterpWrapper_->executeAdjoint(srcFieldSet, tgtFieldSet);
  }
  if (regionalInterp_) {
    regionalInterp_->execute_adjoint(srcFieldSet, tgtFieldSet);
  }
  if (unstructuredInterp_) {
    // Format data
    std::vector<double> vals;
    for (const auto & tgtField : tgtFieldSet) {
      const auto tgtView = atlas::array::make_view<double, 2>(tgtField);
      for (atlas::idx_t jlevel = 0; jlevel < tgtField.shape(1); ++jlevel) {
        for (atlas::idx_t jnode = 0; jnode < tgtField.shape(0); ++jnode) {
          if (ghostVector_[jnode] == 0) {
            vals.push_back(tgtView(jnode, jlevel));
          }
        }
      }
    }

    // Apply unstructured interpolator, adjoint
    const varns::Variables vars(tgtFieldSet.field_names());
    unstructuredInterp_->applyAD(vars, srcFieldSet, vals);

    // Exchange FieldSet halo, adjoint
    srcFieldSet.adjointHaloExchange();
    srcFieldSet.set_dirty();
  }

  oops::Log::trace() << classname() << "::executeAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------


void Interpolation::insertVerticalInterpolation(const std::string & var,
                                                const std::vector<std::array<size_t, 2>> & stencil,
                                                const std::vector<std::array<double, 2>> & weights,
                                                const std::vector<size_t> & stencilSize) {
  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;

  if (verStencil_.find(var) != verStencil_.end()) {
    throw eckit::Exception("vertical interpolation already computed for this variables");
  }
  ASSERT(stencil.size() == stencilSize.size());
  ASSERT(weights.size() == stencilSize.size());
  verStencil_.insert({var, stencil});
  verWeights_.insert({var, weights});
  verStencilSize_.insert({var, stencilSize});

  oops::Log::trace() << classname() << "::insertVerticalInterpolation starting" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeVertical(const atlas::FieldSet & srcFieldSet,
                                    atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::executeVertical starting" << std::endl;

  for (auto & tgtField : tgtFieldSet) {
    const std::string var = tgtField.name();
    const auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[var]);
    auto tgtView = atlas::array::make_view<double, 2>(tgtField);
    tgtView.assign(0.0);
    for (size_t jo = 0; jo < verStencilSize_.at(var).size(); ++jo) {
      for (size_t jj = 0; jj < verStencilSize_.at(var)[jo]; ++jj) {
        tgtView(jo, 0) += verWeights_.at(var)[jo][jj]*srcView(jo, verStencil_.at(var)[jo][jj]);
      }
    }
  }

  oops::Log::trace() << classname() << "::executeVertical done" << std::endl;
}

// -----------------------------------------------------------------------------

void Interpolation::executeVerticalAdjoint(atlas::FieldSet & srcFieldSet,
                                           const atlas::FieldSet & tgtFieldSet) const {
  oops::Log::trace() << classname() << "::executeVerticalAdjoint starting" << std::endl;

  for (const auto & tgtField : tgtFieldSet) {
    const std::string var = tgtField.name();
    auto srcView = atlas::array::make_view<double, 2>(srcFieldSet[var]);
    const auto tgtView = atlas::array::make_view<double, 2>(tgtField);
    srcView.assign(0.0);
    for (size_t jo = 0; jo < verStencilSize_.at(var).size(); ++jo) {
      for (size_t jj = 0; jj < verStencilSize_.at(var)[jo]; ++jj) {
        srcView(jo, verStencil_.at(var)[jo][jj]) += verWeights_.at(var)[jo][jj]*tgtView(jo, 0);
      }
    }
  }

  oops::Log::trace() << classname() << "::executeVerticalAdjoint done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace quenchxx
