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
import os, sys
from os.path import join, basename, dirname, splitext
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))

import subprocess
import Queue

import krebsjobs.parameters.parameterSetsAdaption as parameterSetsAdaption

if __name__ == '__main__':
  factory = getattr(parameterSetsAdaption, 'value_list')
  num_of_sets = len(factory)
  q = Queue.LifoQueue()
  for i in range(num_of_sets):
    q.put(i)
    
  while not q.empty():
    i=q.get()
    subprocess.call(["echo", "starting set: %i"%i])
    #command = 'python2'
    #arguments = '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py value_list --listindex=%i'
    #arguments = '-m scoop -vvv --hostfile=hosts.txt adaption_deap.py value_list --listindex=%i'
    #single_command = 'python2 -m scoop -vvv --hostfile=hosts.txt /localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py value_list --listindex=%i'
    #print('running command: %s' % arguments % i)
    #subprocess.check_call(['python2',  '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', 'value_list', '--listindex=%i'% i])
    #subprocess.call(['python2',  '-m scoop -vvv --hostfile=hosts.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', 'value_list', '--listindex=%i'% i], shell=True)
    subprocess.check_call(["/localdisk/thierry/tc_install/utils/bash_wrapper.sh %i" %i ],shell=True)