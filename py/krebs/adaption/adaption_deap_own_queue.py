#!/usr/bin python2
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
  paramsetNameToConsider = 'value_list_pressure_2'
  vesselFileNameToConsider = '/home/usersHR/thierry/mychosen/vessels-default-typeE-11x15L130-sample00.h5'
  factory = getattr(parameterSetsAdaption, paramsetNameToConsider)
  num_of_sets = len(factory)
  q = Queue.LifoQueue()
  for i in range(num_of_sets):
    q.put(i)
  
  while not q.empty():
    i=q.get()
    #subprocess.call(["echo", "starting set: %i"%i])
    #command = 'python2'
    #arguments = '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py value_list --listindex=%i'
    #arguments = '-m scoop -vvv --hostfile=hosts.txt adaption_deap.py value_list --listindex=%i'
    #single_command = 'python2 -m scoop -vvv --hostfile=hosts.txt /localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py value_list --listindex=%i'
    #print('running command: %s' % arguments % i)
    #subprocess.check_call(['python2',  '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', 'value_list', '--listindex=%i'% i])
    #subprocess.call(['python2',  '-m scoop -vvv --hostfile=hosts.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', 'value_list', '--listindex=%i'% i], shell=True)
    #subprocess.check_call(['/localdisk/thierry/tc_install/utils/bash_wrapper.sh %i' %i ], shell=True)
    #subprocess.check_call(["/localdisk/thierry/tc_install/utils/bash_wrapper.sh", "%i"%i ])
    #subprocess.check_call(['python2', '-c', 'import time; time.sleep(3)'])
    ''' amazingly this line does the job!'''
    #subprocess.check_call(['python2', '-m', 'scoop', '--nice', '19', '--tunnel', '--hostfile', '/home/usersHR/thierry/oldHostsDeap.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', 'value_list', '--listindex=%i'% i])
    if 0: #this works on ichthys
      process = subprocess.Popen(['python2', '-m', 'scoop', '--nice', '19','--hostfile', '/home/usersHR/thierry/oldHostsDeap.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
      returncode = process.wait()
      print('returncode : %i' % returncode)
      print('std.out from subprocess %s' % process.stdout.read())
    if 0: #this works on snowden, note nice is on different possition!
      # processes are left, when strg + c this
      process = subprocess.Popen(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/oldHostsDeapFew.txt','--nice 19', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
      returncode = process.wait()
      print('returncode : %i' % returncode)
      print('std.out from subprocess %s' % process.stdout.read())
    ''' snowden '''
    if 0: #INTERRUPT WORKS WELL WITH THIS COMMAND ON SNOWDEN
      #subprocess.check_call(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/oldHostsDeap.txt','--nice 19', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i])
      #subprocess.check_call(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/hostsIchthys.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i])
      subprocess.check_call(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/armsDeapDurga.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i])
    ''' ichthys '''
    if 1:
      process = subprocess.check_call(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/hostsIchthys.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i,])
      #process = subprocess.Popen(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/hostsIchthys.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
      #returncode = process.wait()
      #print('returncode : %i' % returncode)
      #print('std.out from subprocess %s' % process.stdout.read())
      
      #subprocess.check_call(['python2', '-m', 'scoop', '--hostfile', '/home/usersHR/thierry/hostsIchthys.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i, '--outputFileFolder', '%s' % '/home/userHR/thierry/output_ichthys'])
    
    #subprocess.call(['python2', '-m', 'scoop', '--nice', '19', '--hostfile', '/home/usersHR/thierry/oldHostsDeap.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
    #subprocess.check_call(['/localdisk/thierry/tc_install/utils/bash_wrapper.sh %i'%i ],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    #process = subprocess.check_call(['python2', '-m', 'scoop', '--nice', '19', '--hostfile', '/home/usersHR/thierry/oldHostsDeap.txt', '/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py', '--fileName', '%s'%vesselFileNameToConsider ,'%s' % paramsetNameToConsider, '--listindex=%i'% i])
