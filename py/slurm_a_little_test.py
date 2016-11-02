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
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 09:30:02 2016

@author: thierry
"""
from __future__ import print_function
import pyslurm
import time
def display(node_dict):
    
    if node_dict:

        date_fields = [ 'boot_time', 'slurmd_start_time', 'last_update', 'reason_time' ]

        print('{0:*^80}'.format(''))
        for key, value in node_dict.iteritems():

            print("{0} :".format(key))
            for part_key in sorted(value.iterkeys()):

                if part_key in date_fields:
                    ddate = value[part_key]
                    if ddate == 0:
                        print("\t{0:<17} : N/A".format(part_key))
                    else:
                        ddate = pyslurm.epoch2date(ddate)
                        print("\t{0:<17} : {1}".format(part_key, ddate))
                elif ('reason_uid' in part_key and value['reason'] is None):
                    print("\t{0:<17} : ".format(part_key))
                else: 
                    print("\t{0:<17} : {1}".format(part_key, value[part_key]))

            print('{0:*^80}'.format(''))

def display_nodes_from_github():
    try:

        Nodes = pyslurm.node()
        node_dict = Nodes.get()

        if len(node_dict) > 0:

            display(node_dict)

            print()
            print("Node IDs - {0}".format(Nodes.ids()))

        else:
    
            print("No Nodes found !")

    except ValueError as e:
        print("Error - {0}".format(e.args[0]))

def job_test():
  try:
    a = pyslurm.job()

    jobs =  a.get()
    print(jobs)
  except ValueError as e:
      print("Job list error - {0}".format(e.args[0]))
def stats():
  try:
        stats = pyslurm.statistics()
        s = stats.get()
        #display(s)
        print(s)
  except ValueError as e:
      print("Error - {0}".format(e.args[0]))
if __name__=='__main__':
  #display_nodes_from_github()
  #job_test()
  stats()