#!/usr/bin/env python

from optparse import OptionParser
import sys
import os

if __name__ == '__main__':
    parser = OptionParser(usage='%(prog)s <image pattern>')
    
    opt,args = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        print 'Need image filename pattern (including %(raft)s and %(sensor)s)'
        print 'eg, "imsim-v85471052-r%(raft)s-s%(sensor)s-photom.png"'
        sys.exit(-1)

    imgpat = args[0]
    W,H = 50,50

    print '<html><body><table>'

    for rr in [['',   '01','02','03',     '',],
               ['10', '11', '12', '13', '14',],
               ['20', '21', '22', '23', '24',],
               ['30', '31', '32', '33', '34',],
               ['',   '41', '42', '43',   '',]]:
        print '  <tr>'
        for r in rr:
            if len(r):
                print '    <th>raft ',r,'</th>'
            else:
                print '    <th />'
                
        print '  </tr>'
        print '  <tr>'
        for r in rr:
            if not len(r):
                print '    <td />'
                continue
                
            print '    <td><table>'
            for sr in ['%i' % i for i in range(3)]:
                print '      <tr>'
                for sc in ['%i' % i for i in range(3)]:
                    s = sr + sc
                    imgfn = imgpat % (dict(raft=r, sensor=s))
                    thumbfn = 'thumb-' + imgfn
                    if not os.path.exists(thumbfn) and os.path.exists(imgfn):
                        cmd = 'pngtopnm %s | pnmscale -width %i -height %i | pnmtopng > %s' % (imgfn, W, H, thumbfn)
                        print >>sys.stderr, 'Running:', cmd
                        os.system(cmd)
                    print '        <td><a href="%s"><img width="%i" height="%i" border="0" src="%s"></a></td>' % (imgfn, W, H, thumbfn)
                print '      </tr>'
            print '    </table></td>'
        print '  </tr>'
    print '</table>'
        
