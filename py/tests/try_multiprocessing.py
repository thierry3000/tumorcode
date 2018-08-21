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
from __future__ import print_function
from scoop import futures
import scoop
import numpy as np

threshold = 0.5
def helloWorld(value_around):
    array = scoop.shared.getConst('numbers_to_proceed')
#    if(array[value]<threshold):
#    #print(array[value])
#      return value
    return_list=[]
    for (i,aNumber) in enumerate(array):
      if abs(aNumber-value_around)/value_around <0.1:
        return_list.append(i)
    return return_list

if __name__ == "__main__":
    the_numbers = np.random.randn(2000000)
    the_points = [0.2, 0.25,0.3,0.4,0.5, 0.6, 0.7,0.8,0.9]
    scoop.shared.setConst(numbers_to_proceed = the_numbers)
    list_of_lists = list(futures.map(helloWorld, the_points))
    b_print = False;
    if b_print:
      print(list_of_lists)
    for (aList, aNumber) in zip(list_of_lists, the_points):
      if b_print:
        print("value +- 10percent: %f" % aNumber)
      for aEntry in aList:
        if b_print:
          print("%i : %f" % (aEntry, the_numbers[aEntry]))