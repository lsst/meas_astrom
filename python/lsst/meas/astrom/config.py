import math
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom


class AstromNetDataConfig(pexConfig.Config):
    from lsst.pex.config import Field

    idColumn = Field(
        str,
        '''Column name of the ID number''',
        default='id')

    defaultFilterColumn = Field(
        str,
        '''Column name of the default filter''',
        default='mag')

    starGalaxyColumn = Field(
        str,
        '''Column name of the star/galaxy flag''',
        default=None)

    variableColumn = Field(
        str,
        '''Column name of the star variability flag''',
        default=None)

    magErrorMap = Field(
        dict,
        '''Mapping from mag to mag error column names''',
        default=None)

    # ?
    availableFilters = Field(
        list,
        '''Filter column names''',
        default=[])

class AstromConfig(pexConfig.Config):
    # matchThreshold
    # numBrightStars (?)
    # pixelScaleUncertainty
    # [force]blindSolve
    # distanceForCatalogueMatchinArcsec
    # cleaningParameter [ugh]
    # calculateSip
    # sipOrder
    from lsst.pex.config import Field

    matchThreshold = Field(
        float,
        '''Matching threshold for Astrometry.net solver (log-odds)''',
        default=math.log(1e12), check=lambda x: x > math.log(1e6))
    pixelScaleUncertainty = Field(
        float,
        '''Range of scales around the estimated WCS to allow when solving''',
        default = 1.1, check=lambda x: x>1.)
    numBrightStars = Field(
        int,
        '''Number of stars to use in Astrometry.net solving''',
        default=50, check = lambda x: x>=10)
    forceBlindSolve = Field(
        bool,
        '''Ignore any input WCS and do a blind solve?''',
        default=False)
    catalogMatchDist = Field(
        #afwGeom.Angle,
        float,
        '''When matching image to reference catalog stars, how big should
        the matching radius be?''',
        default=1.,#* afwGeom.arcseconds,
        check=lambda x: x>0)
    cleaningParameter = Field(
        float,
        '''Sigma-clipping parameter in sip/cleanBadPoints.py''',
        default=3., check=lambda x: x>0)
    calculateSip = Field(
        bool,
        '''Compute polynomial SIP distortion terms?''',
        default=True)
    sipOrder = Field(
        int,
        '''Polynomial order of SIP distortion terms''',
        default=4, check=lambda x: x>=2)
    
    # defaultFilterName = Field(
    #     str,
    #     '''Default column name for mag in astrometry_net_data files'''
    #     default='mag')
    # defaultIdName = Field(
    #     str,
    #     '''Default column name for ID num in astrometry_net_data files'''
    #     default='id')
