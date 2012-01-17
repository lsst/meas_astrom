import math
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom


class AstromNetDataConfig(pexConfig.Config):
    from lsst.pex.config import Field, ListField

    idColumn = Field(
        str,
        '''Column name of the ID number''',
        default='id')

    defaultMagColumn = Field(
        str,
        '''Default mag column name''',
        default='mag')

    starGalaxyColumn = Field(
        str,
        '''Column name of the star/galaxy flag''',
        default=None)

    variableColumn = Field(
        str,
        '''Column name of the star variability flag''',
        default=None)

    magErrorColumnMap = Field(
        dict,
        '''Mapping from LSST filter name to mag error column name''',
        default=None)

    magColumnMap = Field(
        dict,
        '''Mapping from LSST filter name to mag column name''',
        default=None)

    # ?
    availableFilters = Field(
        list,
        '''Filter column names''',
        default=[])

    indexFiles = ListField(str, default=[],
                           doc='''Astrometry.net index filenames''')

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

    # maxCpuTime

    matchThreshold = Field(
        float,
        '''Matching threshold for Astrometry.net solver (log-odds)''',
        default=math.log(1e12), check=lambda x: x > math.log(1e6))
    maxStars = Field(
        int,
        '''Maximum number of stars to use in Astrometry.net solving''',
        default=50, check = lambda x: x>=10)

    #forceBlindSolve = Field(
    #    bool,
    #    '''Ignore any input WCS and do a blind solve?''',
    #    default=False)

    # useImageWcs = ChoiceField(
    #     str,
    #     '''Use the WCS header information in the input exposure?  We
    #     can take advantage of the pixel scale and RA,Dec center to
    #     restrict the search and make astrometric matching faster.
    #     '''
    #     default='all',
    #     allowed={'all': 'Use all available WCS information',
    #              'pixelscale': 'Use the pixel scale only',
    #              'radec': 'Use the RA,Dec center only',
    #              })

    useWcsPixelScale = Field(
        bool,
        '''Use the pixel scale from the input exposure\'s WCS headers?''',
        default=True)

    useWcsRaDecCenter = Field(
        bool, default=True,
        doc='''Use the RA,Dec center information from the input exposure\'s WCS headers?''')

    useWcsParity = Field(
        bool, default=True,
        doc='''Use the parity (flip / handedness) of the image from the input exposure\'s WCS headers?''')

    raDecSearchRadius = Field(
        float,
        '''When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure\'s WCS to search for a solution.''',
        default=1.)

    pixelScaleUncertainty = Field(
        float,
        '''Range of pixel scales, around the value in the WCS header, to search''',
        default = 1.1, check=lambda x: x>1.)

    # forceParity?
    # pixelScale?
    # pixelScaleRange?
    # forceRaDecCenter?
    # doTrim ?

    # forceImageSize = Field(
    #     tuple,
    #     '''Ignore the size of the input exposure and assume this
    #     image size instead.''',
    #     optional=True,
    #     check=lambda x: (len(x) == 2 and
    #                      type(x[0]) is int and type(x[1]) is int))


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
