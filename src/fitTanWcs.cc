#include <cmath>

#include "boost/make_shared.hpp"
#include "boost/shared_array.hpp"
#include "boost/scoped_array.hpp"
#include "Eigen/Dense"
#include "fitsio.h"
#include "gsl/gsl_linalg.h"

#include "lsst/afw/coord.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image.h"
#include "lsst/meas/astrom/matchOptimisticB.h"
#include "lsst/meas/astrom/fitTanWcs.h"

double const DegPerRad = 180./M_PI;

using std::sin;
using std::cos;
using std::pow;

namespace lsst { namespace meas { namespace astrom {

namespace {

double calXi(double a, double d, double A, double D) {
    return cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*cos(d)*sin(a-A)*sin(a-A)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
           -cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_D(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.);
}

double calEta(double a, double d, double A, double D) {
    return (cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
           -sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_D(double a, double d, double A, double D) {
    return -pow(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)-1.;
}

} // sub-namespace anonymous

PTR(lsst::afw::image::Wcs) fitTanWcs(
    lsst::afw::table::ReferenceMatchVector const &matches,
    bool verbose
) {
    int npair = matches.size();
    ProxyVector img;
    ProxyVector cat;
    boost::scoped_array<double> x(new double[npair]);
    boost::scoped_array<double> y(new double[npair]);
    boost::scoped_array<double> u(new double[npair]);
    boost::scoped_array<double> v(new double[npair]);

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    for (int i = 0; i < npair; i++) {
        lsst::afw::geom::Point3D v = matches[i].first->getCoord().getVector();
        cx += v[0];
        cy += v[1];
        cz += v[2];
        Sx += matches[i].second->getX();
        Sy += matches[i].second->getY();
    }
    cx /= npair;
    cy /= npair;
    cz /= npair;
    lsst::afw::coord::Coord cmean(lsst::afw::geom::Point3D(cx, cy, cz));
    lsst::afw::geom::PointD crval = cmean.getPosition(lsst::afw::geom::radians);
    lsst::afw::geom::PointD crpix = lsst::afw::geom::Point2D(Sx/npair, Sy/npair);

    int order = 1;
    int ncoeff = (order+1)*(order+2)/2 - 1;
    int ndim = ncoeff * 2 + 2;

    boost::scoped_array<int> xorder(new int[ncoeff]);
    boost::scoped_array<int> yorder(new int[ncoeff]);

    int n = 0;
    for (int i = 1; i <= order; i++) {
        for (int k = 0; k <= i; k++) {
            int j = i - k;
            xorder[n] = j;
            yorder[n] = k;
            n++;
        }
    }

    boost::scoped_array<double> a_data(new double[ndim*ndim]);
    boost::scoped_array<double> b_data(new double[ndim]);

    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            a_data[i*ndim+j] = 0.0;
        }
        b_data[i] = 0.0;
    }

    int iexp = 0; int nexp = 1;
    double w1 = 1.0;
    for (int i = 0; i < npair; i++) {
        double ra = matches[i].first->getRa().asRadians();
        double dec = matches[i].first->getDec().asRadians();
        double xi    = calXi  (ra, dec, crval[0], crval[1]);
        double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
        double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
        double eta   = calEta  (ra, dec, crval[0], crval[1]);
        double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
        double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
        double u = matches[i].second->getX() - crpix[0];
        double v = matches[i].second->getY() - crpix[1];
        if (verbose) {
            std::cout << u << "," << v << " --> " << xi << "," << eta << std::endl;
        }
        int i0 = ncoeff * 2 * iexp;
        for (int k = 0; k < ncoeff; k++) {
            for (int l = 0; l < ncoeff; l++) {
                a_data[(i0+k)*ndim+i0+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
                    pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
                a_data[(i0+ncoeff+k)*ndim+i0+ncoeff+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
                    pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
            }
            b_data[i0+k] += xi * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            b_data[i0+ncoeff+k] += eta * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
        }
        int j0 = ncoeff * 2 * nexp;
        for (int k = 0; k < ncoeff; k++) {
            a_data[(i0+k)*ndim+j0+iexp*2]          += -xi_A  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(i0+k)*ndim+j0+iexp*2+1]        += -xi_D  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(i0+ncoeff+k)*ndim+j0+iexp*2]   += -eta_A * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1] += -eta_D * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(j0+iexp*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+j0+iexp*2];
            a_data[(j0+iexp*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+j0+iexp*2+1];
            a_data[(j0+iexp*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2];
            a_data[(j0+iexp*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1];
        }
        a_data[(j0+iexp*2)*ndim+j0+iexp*2]     += (xi_A * xi_A + eta_A * eta_A) * w1;
        a_data[(j0+iexp*2)*ndim+j0+iexp*2+1]   += (xi_D * xi_A + eta_D * eta_A) * w1;
        a_data[(j0+iexp*2+1)*ndim+j0+iexp*2]   += (xi_A * xi_D + eta_A * eta_D) * w1;
        a_data[(j0+iexp*2+1)*ndim+j0+iexp*2+1] += (xi_D * xi_D + eta_D * eta_D) * w1;
        
        b_data[j0+iexp*2]   += -(xi * xi_A + eta * eta_A) * w1;
        b_data[j0+iexp*2+1] += -(xi * xi_D + eta * eta_D) * w1;
    }
    
    gsl_matrix_view a = gsl_matrix_view_array(a_data.get(), ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data.get(), ndim);
    
    boost::shared_ptr<gsl_vector> c(gsl_vector_alloc(ndim), gsl_vector_free);

    if (verbose) {
        for (int i = 0; i < ndim; ++i) {
            for (int j = 0; j < ndim; ++j) {
                std::cout << a_data[i*ndim+j] << " ";
            }
            std::cout << std::endl;
        }
        for (int i = 0; i < ndim; ++i) {
            std::cout << b_data[i] << " ";
        }
        std::cout << std::endl;
    }

    //int s;
    /*
    boost::shared_ptr<gsl_permutation> p(gsl_permutation_alloc(ndim), gsl_permutation_free);

    gsl_linalg_LU_decomp(&a.matrix, p.get(), &s);
    gsl_linalg_LU_solve(&a.matrix, p.get(), &b.vector, c.get());
    */
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, c.get());

    if (verbose) {
        for (int i = 0; i < ncoeff; i++) {
            printf("%2d %12.5e %12.5e\n", i, c->data[i], c->data[ncoeff+i]);
        }
        printf("\n");
        printf("   %12.5e %12.5e\n", c->data[ncoeff*2], c->data[ncoeff*2+1]);
        printf("\n");
    }

    crval[0] += c->data[ncoeff*2];
    crval[1] += c->data[ncoeff*2+1];

    Eigen::Matrix2d cd; cd << c->data[0], c->data[1], c->data[ncoeff], c->data[ncoeff+1];
    
    crval[0] *= DegPerRad;
    crval[1] *= DegPerRad;
    cd *= DegPerRad;
    return boost::make_shared<lsst::afw::image::TanWcs>(crval, crpix, cd);
}

}}} // namespace lsst::meas::astrom
