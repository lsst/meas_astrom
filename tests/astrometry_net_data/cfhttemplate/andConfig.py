root.defaultMagColumn = 'r'
root.defaultMagErrorColumn = 'r_err'
root.starGalaxyColumn = 'star'
#root.variableColumn = 'variable'
filters = ('u','g','r','i','z')
root.magColumnMap = dict([(f,f) for f in filters])
root.magErrorColumnMap = dict([(f, f + '_err') for f in filters])
root.indexFiles = ['D1.index']#, 'D2.index', 'D3.index', 'D4.index']
