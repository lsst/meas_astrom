'''
HDU #14  Binary Table:  5 columns x 14193 rows
COL NAME             FORMAT
1 mag              1E
2 id               1K
3 mag_err          1E
4 starnotgal       1L
5 variable         1L
'''

root.starGalaxyColumn = 'starnotgal'
root.variableColumn = 'variable'
filters = ('mag')
root.magColumnMap = dict([(f,f) for f in filters])
root.magErrorColumnMap = dict([(f, f + '_err') for f in filters])
root.indexFiles = ['index-2033.fits',]
