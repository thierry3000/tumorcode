#!/bin/bash
# -*- coding: utf-8 -*-
#"""
#Created on Thu Jun 22 15:19:23 2017
#
#@author: thierry
#"""

# This file is part of tumorcode project.
# (http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)
# 
# Copyright (C) 2016  Michael Welter and Thierry Fredrich
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#python2 -m scoop -vvv --hostfile=hosts.txt /localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py value_list --listindex=$1
#python2 -m scoop --nice 19 --hostfile="/home/usersHR/thierry/oldHostsDeap.txt" "/localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py" "value_list3" --listindex=$1 --fileName="/home/usersHR/thierry/mychosen/vessels-default-typeE-11x15L130-sample00.h5"

# Unfortunatelly I was not able to get the loop working in python so I had to do this around it!
echo "hello from bash"
echo "Argument 1: $1"
#python2 -m scoop /tmp/try.py

python2 -m scoop -vvv --nice 19 --hostfile="/home/usersHR/thierry/oldHostsDeap.txt" /localdisk/thierry/tc_install/py/krebs/adaption/adaption_deap.py "value_list3" --listindex=1 --pythonpath=/localdisk/thierry/local/lib64/python2.7/site-packages

