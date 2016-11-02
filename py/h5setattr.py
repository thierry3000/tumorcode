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
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import os,sys
import optparse  #Note: Deprecated since version 2.7: The optparse module is deprecated and will not be developed further; development will continue with the argparse module.
import numpy as np
import myutils


if __name__ == '__main__':
  parser = optparse.OptionParser(usage = 'Usage: %prog FILENAME ... H5PATH ATTRIBUTE-NAME VALUE\n\nNote:H5PATH supports wildcards. VALUE is evaluated as python expression.')
  options, args = parser.parse_args()
  filenames, pattern, attrname, attrdatastr = args[:-3], args[-3], args[-2], args[-1]

  data = eval(attrdatastr, {'os': os, 'np' : np, 'sys' : sys}, {})
    
  for fn in filenames:
    with h5py.File(fn, 'r+') as f:
      items = myutils.walkh5(f, pattern)
      items = [f[i] for i in items]
      for item in items:
        print '%s:%s.attrs["%s"] = %s' % (fn, item.name, attrname, str(data))
        a = item.attrs
        try:
          del a[attrname]
        except KeyError: 
          pass
        a[attrname] = data