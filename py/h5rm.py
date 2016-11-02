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
import optparse #Note: Deprecated since version 2.7. Use argparse instead
import myutils

def check_allow_delete(name, opts):
  if opts.force:
    return True
  ok = raw_input(' [y|.]')
  return ok=='y'


def delete_items(fn, pattern, opts):
  with h5py.File(fn, 'r+') as f:
    items = myutils.walkh5(f, pattern)
    for item in items:
      print '%s:%s' % (fn, item),
      if not opts.dry and check_allow_delete(item, opts):
        del f[item]
      print '\n',

if __name__ == '__main__':
  parser = optparse.OptionParser()
  parser.add_option("-f",dest="force", help="no user interaction", default=False, action="store_true")
  parser.add_option("-n",dest="dry", help="simulate, don't delete anything", default=False, action="store_true")
  options, args = parser.parse_args()
  filenames, pattern = args[:-1], args[-1]
  for fn in filenames:
    delete_items(fn, pattern, options)