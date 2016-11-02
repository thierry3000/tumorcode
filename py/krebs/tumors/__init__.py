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
# -*- coding: utf-8 -*-
from os.path import basename, dirname, join, splitext
import os,sys
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..'))

import krebsutils # import of this must come in front of import of detailedo2 libs because some required initialization stuff on the c++ side (see mainboost.cpp)
sys.path.append(os.path.join(os.path.dirname(__file__),'../../../lib'))

import qsub


if sys.flags.debug:
  tumor = __import__('libtumor_d', globals(), locals())
else:
  tumor = __import__('libtumor_', globals(), locals())
  
def run_faketum(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor.run_faketum_(configstr)

def run_bulktissue_no_vessels(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor.run_bulktissue_no_vessels_(configstr)  

def run_bulktissue_w_vessels(configstr):
  if qsub.is_client:
    qsub.printClientInfo()
  return tumor.run_bulktissue_with_vessels_(configstr)
