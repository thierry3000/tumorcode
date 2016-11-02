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
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar

def calc_rc(r_a,alpha):
  a = np.power(r_a,alpha)
  a = 2*a
  print(a)
  return np.power(2*np.power(r_a,alpha),1/float(alpha))

def func(alpha, *data ):
  return np.power(data[0],alpha)-2*np.power(data[1], alpha)
  
if __name__=="__main__":
  print(calc_rc(13.1,1))
  def func2(alpha):
    return(func(6.3,5,alpha))
  data = (17.45,8.079)
  alpha = fsolve(func,4.5,args = data,full_output=True)
  alpha2 = minimize_scalar(func,bounds=(1,6), method='bounded',args = data)
  print(alpha)  
  print(alpha2.x)