// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2016 LSST/AURA
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef LSST_MEAS_ASTROM_TanSipFitter_INCLUDED
#define LSST_MEAS_ASTROM_TanSipFitter_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/geom/LinearTransform.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/astrom/PolynomialTransform.h"
#include "lsst/meas/astrom/SipTransform.h"

namespace lsst { namespace meas { namespace astrom {


/**
 *  Control object for outlier rejection in ScaledPolynomialTransformFitter.
 *
 *  See ScaledPolynomialTransformFitter::rejectOutliers for more information.
 */
class OutlierRejectionControl {
public:

    OutlierRejectionControl() :
        nSigma(3), nClipMin(1), nClipMax(5)
    {}

    LSST_CONTROL_FIELD(nSigma, double, "Number of sigma to clip at.");

    LSST_CONTROL_FIELD(nClipMin, int, "Always clip at least this many matches.");

    LSST_CONTROL_FIELD(nClipMax, int, "Never clip more than this many matches.");

};


/**
 *  A fitter class for scaled polynomial transforms
 *
 *  This class allows for iteration between actual fitting, outlier rejection,
 *  and estimation of intrinsic scatter.  It also provides access to the
 *  current model for debugging via an afw::table::BaseCatalog that contains
 *  the input values, output values, and the best-fit-transformed input values.
 *
 *  ScaledPolynomialTransformFitter has two public construction methods:
 *
 *   - fromMatches fits a transform that maps intermediate world coordinates
 *     to pixel coordinates, using as inputs a match of reference coordinates
 *     to measured source positions.  This sets up the fitter to do outlier
 *     rejection and intrinsic scatter estimation.
 *
 *   - fromGrid fits an inverse to an existing transform, using the
 *     to-be-inverted transform to populate a grid of positions.  Because the
 *     to-be-inverted transform is assumed to be known exactly, fitters
 *     initialized with this method have no need for outlier rejection or
 *     scatter estimation.
 *
 *  In either case, the fitter creates affine transforms that map the input
 *  and output data points onto [-1, 1].  It then fits a polynomial
 *  transform that, when composed with the input scaling transform and the
 *  inverse of the output scaling transform, maps the input data points to the
 *  output data points.
 *
 *  The fitter can be used in an outlier-rejection loop with the following
 *  pattern (with user-defined convergence criteria):
 *  @code
 *  while (true) {
 *      fitter.fit();
 *      if (converged) break;
 *      fitter.updateModel();
 *      fitter.computeIntrinsicScatter();
 *      fitter.rejectOutliers();
 *  }
 *  @endcode
 *  This pattern fits a model, uses that model to transform the input
 *  points, estimates intrinsic scatter from the differences between
 *  model-transformed points and output data points, and finally rejects
 *  outliers by sigma-clipping those differences before repeating.
 *
 *  ScaledPolynomialTransformFitter instances should be confined to a single thread.
 */
class ScaledPolynomialTransformFitter {
public:

    /**
     *  Initialize a fit from intermediate world coordinates to pixels using
     *  source/reference matches.
     *
     *  @param[in] maxOrder         Maximum polynomial order for the fit.
     *
     *  @param[in] matches          Vector of source-reference matches.
     *                              The centroids and centroid errors from the
     *                              sources are used, along with the coords
     *                              from the reference objects.
     *
     *  @param[in] initialWcs       Initial WCS, used to transform reference
     *                              positions to intermediate world coordinates and
     *                              populate the "ref" field in the data catalog.
     *
     *  @param[in] intrinsicScatter Initial intrinsic scatter to be added in
     *                              quadrature with source measurement errors in
     *                              setting the uncertainty on each match.
     *
     *  This initializes the data catalog with the following fields:
     *   - ref_id:             reference object ID
     *   - src_id:             source ID from the match
     *   - src_{x,y}:          source position in pixels
     *   - src_{xSigma,ySigma,x_y_Cov}:  uncertainty on source positions,
     *                         including measurement errors and the inferred
     *                         intrinsic scatter.
     *   - intermediate_{x,y}: reference positions in intermediate world
     *                         coordinates.
     *   - initial_{x,y}:      reference positions transformed by initialWcs.
     *   - model_{x,y}:        reference positions transformed by current
     *                         best-fit distortion.
     *   - valid:              Flag that is set for objects that have not been
     *                         rejected as outliers.
     */
    static ScaledPolynomialTransformFitter fromMatches(
        int maxOrder,
        afw::table::ReferenceMatchVector const & matches,
        afw::image::Wcs const & initialWcs,
        double intrinsicScatter
    );

    /**
     *  Initialize a fit that inverts an existing transform by evaluating and
     *  fitting to points on a grid.
     *
     *  @param[in] maxOrder   Maximum polynomial order for the fit.
     *
     *  @param[in] bbox       Bounding box for the grid in the input
     *                        coordinates of the toInvert transform (or
     *                        equivalently the output coordinates of the
     *                        transform to be fit).
     *
     *  @param[in] nGridX     Number of grid points in the X direction.
     *
     *  @param[in] nGridY     Number of grid points in the Y direction.
     *
     *  @param[in] toInvert   Transform to invert
     *
     *  This initializes the data catalog with the following fields:
     *   - output:  grid positions passed as input to the toInvert transform,
     *              treated as output data points when fitting its inverse.
     *   - input:   grid positions generated as the output of the toInvert
     *              transform, treated as input data points when fitting its
     *              inverse.
     *   - model:   best-fit-transformed input points.
     *
     *  The updateIntrinsicScatter and rejectOutliers methods cannot be used
     *  on fitters initialized using this method.
     */
    static ScaledPolynomialTransformFitter fromGrid(
        int maxOrder,
        afw::geom::Box2D const & bbox,
        int nGridX, int nGridY,
        ScaledPolynomialTransform const & toInvert
    );

    /**
     *  Perform a linear least-squares fit of the polynomial coefficients.
     *
     *  @param[in]  order    The maximum order of the polynomial transform.
     *                       If negative (the default) the maxOrder from
     *                       construction is used.
     *
     *  After fitting, updateModel() should be called to apply the model
     *  transform to the input points if the user wants to make use of the
     *  data catalog for diagnostics.  Calling fit() alone is sufficient if
     *  the user is only interested in getting the model transform itself.
     */
    void fit(int order=-1);

    /**
     *  Update the 'model' field in the data catalog using the current best-
     *  fit transform.
     *
     *  Should be called after fit() and before updateIntrinsicScatter() if
     *  the fitter is being used in an outlier-rejection loop.
     */
    void updateModel();

    /**
     *  Infer the intrinsic scatter in the offset between the data points and
     *  the best-fit model, and update the uncertainties accordingly.
     *
     *  The scatter is calculated as the RMS scatter that maximizes the
     *  likelihood of the current positional offsets (assumed fixed) assuming
     *  the per-data-point uncertainties are the quadrature sum of the
     *  intrinsic scatter and that source's centroid measurement uncertainty.
     *
     *  Should be called after updateModel() and before rejectOutliers() if
     *  the fitter is being used in an outlier-rejection loop.
     */
    double updateIntrinsicScatter();

    /**
     *  Return the current estimate of the intrinsic scatter.
     *
     *  Because the intrinsic scatter is included in the uncertainty used in
     *  both the fit and outlier rejection, this may be called either before
     *  updateIntrinsicScatter() to obtain the scatter used in the last such
     *  operation or after updateIntrinsicScatter() to obtain the scatter to
     *  be used in the next operation.
     */
    double getIntrinsicScatter() const { return _intrinsicScatter; }

    /**
     *  Mark outliers in the data catalog using sigma clipping.
     *
     *  A data point is declare an outlier if:
     *   - the offset between the model prediction and the data position,
     *     weighted by the uncertainty (including inferred intrinsic scatter)
     *     of that point exceeds ctrl.nSigma
     *   - the number of points already with larger weighted offsets is less
     *     than ctrl.nClipMax.
     *  In addition, the worst (in weighted offset) ctrl.nClipMin data points
     *  are always rejected.
     *
     *  @return A pair of the smallest rejected weighted offset (units of sigma)
     *          and the number of point rejected.
     *
     *  Should be called after updateIntrinsicScatter() and before fit() if
     *  the fitter is being used in an outlier rejection loop.
     */
    std::pair<double,size_t> rejectOutliers(OutlierRejectionControl const & ctrl);

    /**
     *  Return a catalog of data points and model values for diagnostic purposes.
     *
     *  The values in the returned catalog should not be modified by the user.
     *
     *  For information about the schema, either introspect it programmatically
     *  or see fromMatches and fromGrid.
     */
    afw::table::BaseCatalog const & getData() const { return _data; }

    /**
     *  Return the best-fit transform
     *
     *  If f is an instance of this class, then for an
     *  arbitrary point p:
     *  @code
     *  f.getTransform()(p) == f.getOutputScaling().invert()(f.getPoly()(f.getInputScaling()(p)))
     *  @endcode
     */
    ScaledPolynomialTransform const & getTransform() const { return _transform; }

    /**
     *  Return the polynomial part of the best-fit transform.
     */
    PolynomialTransform const & getPoly() const { return _transform.getPoly(); }

    /**
     *  Return the input scaling transform that maps input data points to [-1, 1].
     */
    afw::geom::AffineTransform const & getInputScaling() const { return _transform.getInputScaling(); }

    /**
     *  Return the output scaling transform that maps output data points to [-1, 1].
     */
    afw::geom::AffineTransform const & getOutputScaling() const { return _outputScaling; }

private:

    class Keys;

    ScaledPolynomialTransformFitter(
        afw::table::BaseCatalog const & data,
        Keys const & keys,
        int maxOrder,
        double intrinsicScatter,
        afw::geom::AffineTransform const & inputScaling,
        afw::geom::AffineTransform const & outputScaling
    );

    double computeIntrinsicScatter() const;

    // Normally it's not safe to use a reference as a data member because the
    // class holding it can't control when the referenced object gets
    // destroyed, but this points to one of two singletons (which never get
    // destroyed).
    Keys const & _keys;
    double _intrinsicScatter;
    afw::table::BaseCatalog _data;
    afw::geom::AffineTransform _outputScaling;
    ScaledPolynomialTransform _transform;
    // 2-d generalization of the Vandermonde matrix: evaluates polynomial at
    // all data points when multiplied by a vector of packed polynomial
    // coefficients.
    Eigen::MatrixXd _vandermonde;
};

}}} // namespace lsst::meas::astrom

#endif // !LSST_MEAS_ASTROM_TanSipFitter_INCLUDEDtr