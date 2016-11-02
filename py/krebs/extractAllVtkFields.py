#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2016  Michael Welter and Thierry Fredrich

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
#!/usr/bin/env python2
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))


from krebs.extractVtkFields import Extractor
import h5py
import os, sys
from os.path import basename, splitext
import myutils


if __name__ == '__main__':
  fn = sys.argv[1]
  fn, _ = myutils.splitH5PathsFromFilename(fn)
  grpnames = [ sys.argv[2] ]
  prefix = sys.argv[3] if len(sys.argv)>=4 else ""
  f = h5py.File(fn, 'r')
  dirs = myutils.walkh5(f['.'], grpnames[0])
  print 'found items %s to search pattern %s' % (str(dirs), grpnames[0])  
  for d in dirs:  
    e = Extractor(f, [d], recursive = True) 
    print 'found datasets:', e.getDatasetPaths(), 'in', d
    e.write("%s%s%s.vtk" % (prefix, splitext(basename(fn))[0], d.replace('/','-') if d!='/' else ''))
