import matplotlib
matplotlib.use('Agg')
from pylab import *


'''
matches
----------------------
= det::matchRaDec(_catSet, _imgSet, _distInArcsec);
in afw/src/detection/SourceMatch.cc
is a vector<det::SourceMatch>
SourceMatch: Source first, Source second, double distance
Source.{X_astrom, Y_astrom, ..., RA_object, DEC_object}

refsources
----------------------
solver.getCatalogue
is an afw::detection::SourceSet
Source.getRa(), getDec()

--- both img sources and ref sources should have x,y,ra,dec set
 by meas_astrom:sip:MatchSrcToCatalogue.cc : findMatches()

'''
def wcsPlots(wcs, imgsources, refsources, matches, W, H, prefix):
	print 'WCS plots'
	clf()

	ix = [s.getXAstrom() for s in imgsources]
	iy = [s.getYAstrom() for s in imgsources]
	plot(ix, iy, 'r.')

	rx = [s.getXAstrom() for s in refsources]
	ry = [s.getYAstrom() for s in refsources]
	plot(rx, ry, 'bo')#, mec='b', mfc='none')

	print 'ref source 0:', dir(refsources[0])

	print 'match 0:', dir(matches[0])

	'''
	ref source 0: ['SWIGSharedPtrUpcast', '__class__', '__del__', '__delattr__', '__dict__', '__doc__', '__eq__', '__getattr__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__str__', '__swig_destroy__', '__swig_getmethods__', '__swig_setmethods__', '__weakref__', '_s', 'getAmpExposureId', 'getApDia', 'getApFlux', 'getApFluxErr', 'getAtmCorrFlux', 'getAtmCorrFluxErr', 'getChi2', 'getDec', 'getDecAstrom', 'getDecAstromErr', 'getDecErrForDetection', 'getDecErrForWcs', 'getDecFlux', 'getDecFluxErr', 'getDecObject', 'getDecPeak', 'getFilterId', 'getFlagForAssociation', 'getFlagForDetection', 'getFlagForWcs', 'getId', 'getInstFlux', 'getInstFluxErr', 'getIxx', 'getIxxErr', 'getIxy', 'getIxyErr', 'getIyy', 'getIyyErr', 'getModelFlux', 'getModelFluxErr', 'getMovingObjectId', 'getNonGrayCorrFlux', 'getNonGrayCorrFluxErr', 'getObjectId', 'getPetroFlux', 'getPetroFluxErr', 'getProcHistoryId', 'getPsfFlux', 'getPsfFluxErr', 'getRa', 'getRaAstrom', 'getRaAstromErr', 'getRaErrForDetection', 'getRaErrForWcs', 'getRaFlux', 'getRaFluxErr', 'getRaObject', 'getRaPeak', 'getSky', 'getSkyErr', 'getSnr', 'getSourceId', 'getTaiMidPoint', 'getTaiRange', 'getXAstrom', 'getXAstromErr', 'getXFlux', 'getXFluxErr', 'getXPeak', 'getYAstrom', 'getYAstromErr', 'getYFlux', 'getYFluxErr', 'getYPeak', 'isNull', 'setAmpExposureId', 'setApDia', 'setApFlux', 'setApFluxErr', 'setAtmCorrFlux', 'setAtmCorrFluxErr', 'setChi2', 'setDec', 'setDecAstrom', 'setDecAstromErr', 'setDecErrForDetection', 'setDecErrForWcs', 'setDecFlux', 'setDecFluxErr', 'setDecObject', 'setDecPeak', 'setFilterId', 'setFlagForAssociation', 'setFlagForDetection', 'setFlagForWcs', 'setId', 'setInstFlux', 'setInstFluxErr', 'setIxx', 'setIxxErr', 'setIxy', 'setIxyErr', 'setIyy', 'setIyyErr', 'setModelFlux', 'setModelFluxErr', 'setMovingObjectId', 'setNonGrayCorrFlux', 'setNonGrayCorrFluxErr', 'setNotNull', 'setNull', 'setObjectId', 'setPetroFlux', 'setPetroFluxErr', 'setProcHistoryId', 'setPsfFlux', 'setPsfFluxErr', 'setRa', 'setRaAstrom', 'setRaAstromErr', 'setRaErrForDetection', 'setRaErrForWcs', 'setRaFlux', 'setRaFluxErr', 'setRaObject', 'setRaPeak', 'setSky', 'setSkyErr', 'setSnr', 'setSourceId', 'setTaiMidPoint', 'setTaiRange', 'setXAstrom', 'setXAstromErr', 'setXFlux', 'setXFluxErr', 'setXPeak', 'setYAstrom', 'setYAstromErr', 'setYFlux', 'setYFluxErr', 'setYPeak', 'this', 'toString']
	match 0: ['__class__', '__del__', '__delattr__', '__dict__', '__doc__', '__getattr__', '__getattribute__', '__getitem__', '__hash__', '__init__', '__len__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setitem__', '__str__', '__swig_destroy__', '__swig_getmethods__', '__swig_setmethods__', '__weakref__', 'clone', 'distance', 'first', 'second', 'this']
	'''

	# Only getRa() is non-zero.
	#r0 = refsources[0]
	#print r0.getRa(), r0.getRaAstrom(), r0.getRaObject()

	#xy = array([wcs.skyToPixel(r.getRa(), r.getDec()) for r in refsources])
	rx,ry = [],[]
	for r in refsources:
		xy = wcs.skyToPixel(r.getRa(), r.getDec())
		rx.append(xy[0])
		ry.append(xy[1])
	plot(rx, ry, 'bo', mec='b', mfc='none')

	for m in matches:
		plot([m.first.getXAstrom(), m.second.getXAstrom()],
			 [m.first.getYAstrom(), m.second.getYAstrom()], 'g.-')

	#print rx, ry
	#print refsources
	#print imgsources
	#print matches

	axis('scaled')
	axis([0, W, 0, H])

	fn = prefix + '-matches.png'
	print 'Saving', fn
	savefig(fn)
	
