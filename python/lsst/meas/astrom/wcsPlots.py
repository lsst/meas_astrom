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
	clf()

	x = [s.getXAstrom() for s in imgsources]
	y = [s.getYAstrom() for s in imgsources]
	plot(x, y, 'r.')

	axis([0, W, 0, H])
	savefig(prefix + '-matches.png')
	
