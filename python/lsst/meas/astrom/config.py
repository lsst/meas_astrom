import math
import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom


class AstrometryNetDataConfig(pexConfig.Config):
    from lsst.pex.config import Field, ListField, DictField

    idColumn = Field(
        '''Column name of the ID number''',
        str,
        default='id')

    defaultMagColumn = Field(
        '''Default mag column name''',
        str,
        default='mag')

    defaultMagErrorColumn = Field(
        '''Default mag error column name''',
        str,
        default='')

    starGalaxyColumn = Field(
        '''Column name of the star/galaxy flag''',
        str,
        default=None)

    variableColumn = Field(
        '''Column name of the star variability flag''',
        str,
        default=None)

    magErrorColumnMap = DictField(
        doc='''Mapping from LSST filter name to mag error column name''',
        keytype=str,
        itemtype=str,
        default={})

    magColumnMap = DictField(
        doc='''Mapping from LSST filter name to mag column name''',
        keytype=str,
        itemtype=str,
        default={})

    # ?
    # availableFilters = Field(
    #     list,
    #     '''Filter column names''',
    #     default=[])

    indexFiles = ListField(dtype=str, default=[],
                           doc='''Astrometry.net index filenames''')

class MeasAstromConfig(pexConfig.Config):
    from lsst.pex.config import Field, RangeField, DictField

    maxCpuTime = RangeField(
        '''Maximum CPU time to spend solving, in seconds''',
        float,
        default=0., min=0.)

    matchThreshold = RangeField(
        '''Matching threshold for Astrometry.net solver (log-odds)''',
        float,
        default=math.log(1e12), min=math.log(1e6))

    maxStars = RangeField(
        '''Maximum number of stars to use in Astrometry.net solving''',
        int,
        default=50, min=10)

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

    raDecSearchRadius = RangeField(
        '''When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure\'s WCS to search for a solution.''',
        float,
        default=1., min=0.)

    pixelScaleUncertainty = RangeField(
        '''Range of pixel scales, around the value in the WCS header, to search.  If the value of this field is X and the nominal scale is S, the range searched will be  S/X to S*X''',
        float,
        default = 1.1, min=1.001)

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

    catalogMatchDist = RangeField(
        #afwGeom.Angle,
        '''When matching image to reference catalog stars, how big should
        the matching radius be?''',
        float,
        default=1.,#* afwGeom.arcseconds,
        min=0.)

    cleaningParameter = RangeField(
        '''Sigma-clipping parameter in sip/cleanBadPoints.py''',
        float,
        default=3., min=0.)

    calculateSip = Field(
        '''Compute polynomial SIP distortion terms?''',
        bool,
        default=True)

    sipOrder = RangeField(
        '''Polynomial order of SIP distortion terms''',
        int,
        default=4, min=2)

    filterMap = DictField(
        doc="Mapping from input filter to catalogue filter",
        keytype=str, itemtype=str,
        default={},
        optional=True)
