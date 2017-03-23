#!/usr/bin/python2.7
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

from __future__ import print_function
import platform
import socket

# symbols which should be visible to importers of this module
__all__ = [ 'name' ] 
name = 'unknown'


def getname():
  # now that does return the full name
  fullname = socket.getfqdn()
  releasename = platform.release()
  
  #fullname, _, _ = socket.gethostbyaddr(hostname) # gethostname already returns the full name
  #n = os.uname()[1] ?? what was wrong with that???
  if 'sleipnir' in fullname:
    return 'sleipnir'
  if 'durga' in fullname:
    return 'durga'
  durga = [ ('arm%i' % i) for i in range(0, 33) ]
  sleipnir = [ ('leg%02i' % i) for i in range(0, 25) ]
  snowden = [ ('leak%02i' % i) for i in range(10,65)]
  snowden2 = [ ('leak%01i' % i) for i in range(1,10)]
  if any((s in fullname) for s in snowden):
    return 'snowden'
  if any((s in fullname) for s in snowden2):
   return 'snowden'
  if any((s in fullname) for s in durga):
    return 'durga'
  elif any((s in fullname) for s in sleipnir):
    return 'sleipnir'
  if 'snowden.lusi' in fullname:
      return 'snowden'
  if 'dema' in fullname:
    if 'gentoo' in releasename:
      #assume we are on olafs old system
      return 'lusi-gentoo'
    elif not 'gentoo' in releasename and 'ARCH' in releasename:
      return 'lusi-arch'
    else:
      return 'is there is new system?'

# executed on import (i.e. always)
name = getname()

# when executed as script
if __name__ == '__main__':
  print(name)
  
