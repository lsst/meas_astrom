import math
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom


class AstromNetDataConfig(pexConfig.Config):
    from lsst.pex.config import Field, ListField

    idColumn = Field(
        '''Column name of the ID number''',
        str,
        default='id')

    defaultMagColumn = Field(
        '''Default mag column name''',
        str,
        default='mag')

    starGalaxyColumn = Field(
        '''Column name of the star/galaxy flag''',
        str,
        default=None)

    variableColumn = Field(
        '''Column name of the star variability flag''',
        str,
        default=None)

    magErrorColumnMap = Field(
        '''Mapping from LSST filter name to mag error column name''',
        dict,
        default=None)

    magColumnMap = Field(
        '''Mapping from LSST filter name to mag column name''',
        dict,
        default=None)

    # ?
    # availableFilters = Field(
    #     list,
    #     '''Filter column names''',
    #     default=[])

    indexFiles = ListField(dtype=str, default=[],
                           doc='''Astrometry.net index filenames''')

class AstromConfig(pexConfig.Config):
    from lsst.pex.config import Field

    maxCpuTime = Field(
        '''Maximum CPU time to spend solving, in seconds''',
        float,
        default=0., check=lambda x: x>=0.)

    matchThreshold = Field(
        '''Matching threshold for Astrometry.net solver (log-odds)''',
        float,
        default=math.log(1e12), check=lambda x: x > math.log(1e6))

    maxStars = Field(
        '''Maximum number of stars to use in Astrometry.net solving''',
        int,
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
        '''Use the pixel scale from the input exposure\'s WCS headers?''',
        bool,
        default=True)

    useWcsRaDecCenter = Field(
        dtype=bool, default=True,
        doc='''Use the RA,Dec center information from the input exposure\'s WCS headers?''')

    useWcsParity = Field(
        dtype=bool, default=True,
        doc='''Use the parity (flip / handedness) of the image from the input exposure\'s WCS headers?''')

    raDecSearchRadius = Field(
        '''When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure\'s WCS to search for a solution.''',
        float,
        default=1.)

    pixelScaleUncertainty = Field(
        '''Range of pixel scales, around the value in the WCS header, to search''',
        float,
        default = 1.1, check=lambda x: x>1.)

    # forceParity?
    # forcePixelScale?
    # forceRaDecCenter?
    # forcePixelScaleRange?
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
        '''When matching image to reference catalog stars, how big should
        the matching radius be?''',
        float,
        default=1.,#* afwGeom.arcseconds,
        check=lambda x: x>0)

    cleaningParameter = Field(
        '''Sigma-clipping parameter in sip/cleanBadPoints.py''',
        float,
        default=3., check=lambda x: x>0)

    calculateSip = Field(
        '''Compute polynomial SIP distortion terms?''',
        bool,
        default=True)

    sipOrder = Field(
        '''Polynomial order of SIP distortion terms''',
        int,
        default=4, check=lambda x: x>=2)

