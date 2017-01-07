// -*- lsst-C++ -*-

#include <sstream>
#include <utility>
#include <vector>

#include "lsst/pex/exceptions.h"
#include "lsst/meas/astrom/astrometry_net.h"

namespace lsst {
namespace meas {
namespace astrom {

namespace {

struct timer_baton {
    solver_t* s;
    double timelimit;
};

static time_t timer_callback(void* baton) {
    struct timer_baton* tt = static_cast<struct timer_baton*>(baton);
    solver_t* solver = tt->s;
    if (solver->timeused > tt->timelimit)
        solver->quit_now = 1;
    return 1;
}

}  // namespace lsst::meas::astrom::<anonymous>

MultiIndex::MultiIndex(std::string const & filepath) : _multiindex(multiindex_new(filepath.c_str())) {
    if (!_multiindex) {
        std::ostringstream os;
        os << "Could not read multi-index star file " << filepath;
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, os.str());
    }
}

void MultiIndex::addIndex(std::string const & filepath, bool metadataOnly) {
    int const flags = metadataOnly ? INDEX_ONLY_LOAD_METADATA : 0;
    if (multiindex_add_index(_multiindex.get(), filepath.c_str(), flags)) {
        std::ostringstream os;
        os << "Failed to read multiindex from \"" << filepath << "\""
            << " with meatadaOnly=" << metadataOnly;
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, os.str());
    }
}

int MultiIndex::isWithinRange(double ra, double dec, double radius_deg) {
    return index_is_within_range(multiindex_get(_multiindex.get(), 0), ra, dec, radius_deg);
}

void MultiIndex::reload() {
    if (multiindex_reload_starkd(_multiindex.get())) {
        std::ostringstream os;
        os << "Failed to reload multi-index star file " << getName();
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, os.str());
    }
}


Solver::Solver() : _solver(solver_new()) {}

Solver::~Solver() {
    // Working around a bug in Astrometry.net: doesn't take ownership of the field.
    // unseemly familiarity with the innards... but valgrind-clean.
    starxy_free(_solver->fieldxy);
    _solver->fieldxy = NULL;
}

lsst::afw::table::SimpleCatalog Solver::getCatalog(
    std::vector<index_t*> inds,
    lsst::afw::coord::Coord const &ctrCoord,
    lsst::afw::geom::Angle const &radius,
    const char* idCol,
    std::vector<std::string> const& filterNameList,
    std::vector<std::string> const& magColList,
    std::vector<std::string> const& magErrColList,
    const char* starGalCol,
    const char* varCol,
    bool uniqueIds)
{
    if ((filterNameList.size() != magColList.size()) || (filterNameList.size() != magErrColList.size())) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterError,
            "Filter name, mag column, and mag error column vectors must be the same length.");
    }
    std::vector<lsst::meas::astrom::detail::MagColInfo> magColInfoList;
    for (size_t i=0; i<filterNameList.size(); ++i) {
        lsst::meas::astrom::detail::MagColInfo mc;
        mc.filterName = filterNameList[i];
        mc.magCol = magColList[i];
        mc.magErrCol = magErrColList[i];
        magColInfoList.push_back(mc);
    }
    return lsst::meas::astrom::detail::getCatalogImpl(inds, ctrCoord, radius,
        idCol, magColInfoList, starGalCol, varCol, uniqueIds);
}

std::shared_ptr<lsst::daf::base::PropertyList> Solver::getSolveStats() const {
    // Gather solve stats...
    auto qa = std::make_shared<daf::base::PropertyList>();
    // FIXME -- Ticket #1875 prevents dotted-names from working with toString().
    qa->set("meas_astrom*an*n_tried", _solver->numtries);
    qa->set("meas_astrom*an*n_matched", _solver->nummatches);
    qa->set("meas_astrom*an*n_scaleok", _solver->numscaleok);
    qa->set("meas_astrom*an*n_cxdxcut", _solver->num_cxdx_skipped);
    qa->set("meas_astrom*an*n_meanxcut", _solver->num_meanx_skipped);
    qa->set("meas_astrom*an*n_radeccut", _solver->num_radec_skipped);
    qa->set("meas_astrom*an*n_scalecut", _solver->num_abscale_skipped);
    qa->set("meas_astrom*an*n_verified", _solver->num_verified);
    qa->set("meas_astrom*an*time_used", _solver->timeused);
    qa->set("meas_astrom*an*best_logodds", _solver->best_logodds);
    if (_solver->best_index) {
        index_t* ind = _solver->best_index;
        qa->set("meas_astrom*an*best_index*id", ind->indexid);
        qa->set("meas_astrom*an*best_index*hp", ind->healpix);
        qa->set("meas_astrom*an*best_index*nside", ind->hpnside);
        qa->set("meas_astrom*an*best_index*name", std::string(ind->indexname));
    }
    if (_solver->have_best_match) {
        MatchObj* mo = &(_solver->best_match);
        std::string s = boost::str(boost::format("%i") % mo->star[0]);
        for (int i=1; i<mo->dimquads; i++)
            s = s + boost::str(boost::format(", %i") % mo->star[i]);
        qa->set("meas_astrom*an*best_match*starinds", s);
        qa->set("meas_astrom*an*best_match*coderr", std::sqrt(mo->code_err));
        qa->set("meas_astrom*an*best_match*nmatch", mo->nmatch);
        qa->set("meas_astrom*an*best_match*ndistract", mo->ndistractor);
        qa->set("meas_astrom*an*best_match*nconflict", mo->nconflict);
        qa->set("meas_astrom*an*best_match*nfield", mo->nfield);
        qa->set("meas_astrom*an*best_match*nindex", mo->nindex);
        qa->set("meas_astrom*an*best_match*nbest", mo->nbest);
        qa->set("meas_astrom*an*best_match*logodds", mo->logodds);
        qa->set("meas_astrom*an*best_match*parity", mo->parity ? 0 : 1);
        qa->set("meas_astrom*an*best_match*nobjs", mo->objs_tried);
    }
    return qa;
}

std::shared_ptr<lsst::afw::image::Wcs> Solver::getWcs() {
    MatchObj* match = solver_get_best_match(_solver.get());
    if (!match)
        return std::shared_ptr<afw::image::Wcs>();
    tan_t* wcs = &(match->wcstan);

    afw::geom::Point2D crpix(wcs->crpix[0], wcs->crpix[1]);
    std::shared_ptr<afw::coord::Coord const> crval
        (new afw::coord::Coord(wcs->crval[0] * afw::geom::degrees,
                             wcs->crval[1] * afw::geom::degrees));
    return afw::image::makeWcs(*crval, crpix,
                             wcs->cd[0][0], wcs->cd[0][1],
                             wcs->cd[1][0], wcs->cd[1][1]);
}

void Solver::run(double cpulimit) {
    solver_log_params(_solver.get());
    struct timer_baton tt;
    if (cpulimit > 0.) {
        tt.s = _solver.get();
        tt.timelimit = cpulimit;
        _solver->userdata = &tt;
        _solver->timer_callback = timer_callback;
    }
    solver_run(_solver.get());
    if (cpulimit > 0.) {
        _solver->timer_callback = NULL;
        _solver->userdata = NULL;
    }
}

/**
 * Add indices to the solver
 *
 * The indices are bare pointers whose memory is managed by the caller.
 * Typically the indices are owned by a MultiIndex object owned by the caller.
 */
void Solver::addIndices(std::vector<index_t*> inds) {
    for (std::vector<index_t*>::iterator pind = inds.begin();
         pind != inds.end(); ++pind) {
        lsst::meas::astrom::detail::IndexManager man(*pind);
//            printf("Checking index \"%s\"\n", man.index->indexname);
        if (_solver->use_radec) {
            double ra,dec,radius;
            xyzarr2radecdeg(_solver->centerxyz, &ra, &dec);
            radius = distsq2deg(_solver->r2);
            if (!index_is_within_range(man.index, ra, dec, radius)) {
                                 //printf("Not within RA,Dec range\n");
                continue;
            }
        }
        // qlo,qhi in arcsec
        double qlo, qhi;
        solver_get_quad_size_range_arcsec(_solver.get(), &qlo, &qhi);
        if (!index_overlaps_scale_range(man.index, qlo, qhi)) {
//                printf("Not within quad scale range\n");
            continue;
        }
//            printf("Adding index.\n");
        if (index_reload(man.index)) {
            throw LSST_EXCEPT(lsst::pex::exceptions::IoError,
                              "Failed to index_reload() an astrometry_net_data index file -- out of file descriptors?");
        }

        solver_add_index(_solver.get(), man.index);
    }
}

void Solver::setImageSize(int width, int height) {
    solver_set_field_bounds(_solver.get(), 0, width, 0, height);
    double hi = hypot(width, height);
    double lo = 0.1 * std::min(width, height);
    solver_set_quad_size_range(_solver.get(), lo, hi);
}

void Solver::setStars(lsst::afw::table::SourceCatalog const & srcs, int x0, int y0) {
    // convert to Astrometry.net "starxy_t"
    starxy_free(_solver->fieldxy);
    const size_t N = srcs.size();
    starxy_t *starxy = starxy_new(N, true, false);
    for (size_t i=0; i<N; ++i) {
        double const x    = srcs[i].getX();
        double const y    = srcs[i].getY();
        double const flux = srcs[i].getPsfFlux();
        starxy_set(starxy, i, x - x0, y - y0);
        starxy_set_flux(starxy, i, flux);
    }
    // Sort the array
    starxy_sort_by_flux(starxy);

    starxy_free(solver_get_field(_solver.get()));
    solver_free_field(_solver.get());
    solver_set_field(_solver.get(), starxy);
    solver_reset_field_size(_solver.get());
    // Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver.get());
}


lsst::afw::geom::Angle healpixDistance(int hp, int nside, lsst::afw::coord::Coord const& coord) {
    lsst::afw::coord::IcrsCoord icrs = coord.toIcrs();
    return lsst::afw::geom::Angle(healpix_distance_to_radec(hp, nside, icrs.getLongitude().asDegrees(),
                                                            icrs.getLatitude().asDegrees(), NULL),
                                  lsst::afw::geom::degrees);
}

}}}  // namespace lsst::meas::astrom