'''
We used to have AstrometryNetDataConfig() use the pex_config
mechanism, but we need nested lists, so we do this home-brew version
instead.
'''


def _checkMagMap(magmap):
    '''
    Checks the validity of a magnitude column map in AstrometryNetDataConfig.
    '''
    if not isinstance(magmap, dict):
        raise RuntimeError('Mag maps must be dicts')
    for k,v in magmap.items():
        if not isinstance(k, str):
            raise RuntimeError('Mag maps must be dicts mapping str->str: got bad key \"%s\"' % str(k))
        if not isinstance(v, str):
            raise RuntimeError('Mag maps must be dicts mapping str->str: got bad value \"%s\"' % str(v))
        if not (len(k) > 0 and len(v) > 0):
            raise RuntimeError('Mag maps items must be non-empty: got bad values \"%s\" -> \"%s\"' % (str(k), str(v)))
        
def _checkIndexList(indexList):
    '''
    Checks the validity of an index list in AstrometryNetDataConfig.
    '''
    if not isinstance(indexList, list):
        raise RuntimeError('indexList config item must be a list')
    for k in indexList:
        if not isinstance(k, str):
            raise RuntimeError('indexList config items must be strings: got bad one \"%s\"' % str(k))
        if len(k) == 0:
            raise RuntimeError('indexList config items must be non-empty strings')

def _checkMultiIndexList(multiIndexList):
    '''
    Checks the validity of a multi_index list in AstrometryNetDataConfig.
    '''
    if not isinstance(multiIndexList, list):
        raise RuntimeError('multiIndexList config item must be a list')
    for k in multiIndexList:
        if not isinstance(k, list):
            raise RuntimeError('multiIndexList config items must be lists: got bad one \"%s\"' % str(k))
        if len(k) == 0:
            raise RuntimeError('multiIndexList config items must be non-empty lists')
        for kk in k:
            if not isinstance(kk, str):
                raise RuntimeError('multiIndexList config items must be strings: got bad one \"%s\"' % str(kk))
            if len(kk) == 0:
                raise RuntimeError('multiIndexList config items must be non-empty strings')

class AstrometryNetDataConfig(object):
    '''
    Astrometry.net data config object.  This is a plain-python config
    structure similar to pexConfig.

    For examples of use, see tests/astrometry_net_data/photocal/andConfig*.py

    '''
    fields = [
        ('idColumn', str, 'id', None,
         'Column name (in the index files) of the ID number of reference sources'),
        ('defaultMagColumn', str, 'mag', None,
         'Default column name (in the index files) of the reference source mag'),
        ('defaultMagErrorColumn', str, '', None,
         'Default column name (in the index files) of the reference source mag error'),
        ('starGalaxyColumn', str, None, None,
         'Column name of the star/galaxy flag'),
        ('variableColumn', str, None, None,
         'Column name of the star variability flag'),
        ('magErrorColumnMap', dict, {}, _checkMagMap,
         'Mapping from LSST filter name to mag error column name'),
        ('magColumnMap', dict, {}, _checkMagMap,
         'Mapping from LSST filter name to mag column name'),
        ('indexFiles', list, [], _checkIndexList,
          'List of Astrometry.net index filenames'),
        ('multiIndexFiles', list, [], _checkMultiIndexList,
         'Astrometry.net multi-index filename lists.  Each item in this list must itself be a list of filenames.  The first filename is the file that contains the star kd-tree and tag-along tables.  Subsequent filenames must be files containing just the non-star index parts (quads and code kd-tree).  Note that this means you may need to repeat the first filename if it contains a star kd-tree and the first index.'),
        ('allowCache', bool, True, None,
         'Allow use of cache for reading index file regions?'),
         ]

    def load(self, fn):
        # Hold on to your socks!
        loc = dict(root=self)
        execfile(fn, globals(), loc)
    
    def __init__(self, **kwargs):
        self.setDefaults()
        for k,v in kwargs.items():
            self.set(k, v)

    def setDefaults(self):
        for nm,typ,deef,check,doc in AstrometryNetDataConfig.fields:
            self.set(nm, deef)

    def set(self, k, v):
        setattr(self, k, v)

    def __setattr__(self, k, v):
        for nm,typ,deef,check,doc in AstrometryNetDataConfig.fields:
            if k != nm:
                continue
            if typ is not None:
                if v is None:
                    pass
                elif not isinstance(v, typ):
                    raise RuntimeError(('Attempted to set AstrometryNetDataConfig'
                                        ' field \"%s\" to type %s, but need type %s') %
                                        (nm, str(typ), str(type(v))))
            if check is not None:
                check(v)
            # Looks ok; set it!
            object.__setattr__(self, nm, v)
            return

        raise RuntimeError('Attempted to set invalid AstrometryNetDataConfig'
                           ' field \"%s\"' % k)
