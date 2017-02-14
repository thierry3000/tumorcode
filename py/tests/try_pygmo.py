# -*- coding: utf-8 -*-
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
import sys
#from os.path import  basename, dirname
  
import identifycluster
if identifycluster.name == 'snowden':
  import matplotlib
  matplotlib.use('Agg') 
#import qsub
#import dicttoinfo
#import myutils
#import h5py
import h5files
import krebs.adaption
#import itertools
#import time
#import string
#import numpy as np

#sys.path.append(join(dirname(__file__),'/localdisk/thierry/tc_install/py/krebs/adaptionAnalysis'))
#from krebs.adaption import getVesselTypes



#import krebs.adaption

from krebsjobs.parameters import parameterSetsAdaption
#import krebsjobs.submitAdaption
#import krebsjobs.parameters.tumorParameters as parameterSets_tum

#import krebsutils
#if sys.flags.debug: # library search path made known by initialization of krebsutils
#  adaption_cpp = __import__('libadaption_d', globals(), locals())
#else:
#  adaption_cpp = __import__('libadaption_', globals(), locals())

#globals
#calculated_mean_cap_flow = 0.0

# velocity in capillary about 1mm/s = 1000mu/s
# typical diameter 5mu e.g r=2.5mu
# 3.1415*2.5*2.5*1000 about 2e4
suggested_cap_flow = 20000.;

if 0:
  from PyGMO import *
  class my_problem(problem.base):
      """
      De Jong (sphere) function implemented purely in Python.
  
      USAGE: my_problem(dim=10)
  
      * dim problem dimension
      """
  
      def __init__(self, dim=4):
          # First we call the constructor of the base class telling PyGMO
          # what kind of problem to expect ('dim' dimensions, 1 objective, 0 contraints etc.)
          super(my_problem,self).__init__(dim)
  
          # We set the problem bounds (in this case equal for all components)
          self.set_bounds(5.12, 6.12)
  
      # Reimplement the virtual method that defines the objective function.
      def _objfun_impl(self, x):
  
          # Compute the sphere function
          f = sum([x[i] ** 2 for i in range(self.dimension)])
  
          # Note that we return a tuple with one element only. In PyGMO the objective functions
          # return tuples so that multi-objective optimization is also possible.
          return (f, )
  
      # Finally we also reimplement a virtual method that adds some output to the __repr__ method
      def human_readable_extra(self):
          return "\n\t Problem dimension: " + str(self.__dim)
def own():
  prob = my_problem()  # Create a 10-dimensional problem
  #algo = algorithm.bee_colony(gen=500)  # 500 generations of bee_colony algorithm
  algo = algorithm.pso(gen=2000)
  isl = island(algo, prob, 20)  # Instantiate population with 20 individuals
  archi = archipelago(algo,prob,8,20)
  #isl.evolve(1)  # Evolve the island once
  #isl.join()
  #print(isl.population.champion.f)
  
  print min([isl.population.champion.f for isl in archi])
  archi.evolve(10)
  print min([isl.population.champion.f for isl in archi])
if 0:
  class my_tricky_problem(problem.base):
      """
      De Jong (sphere) function implemented purely in Python.
  
      USAGE: my_problem(dim=10)
  
      * dim problem dimension
      """
      def __init__(self, dim=4, astring='bar', other=3.3):
          # First we call the constructor of the base class telling PyGMO
          # what kind of problem to expect ('dim' dimensions, 1 objective, 0 contraints etc.)
          super(my_tricky_problem,self).__init__(dim)
          #super(my_tricky_problem,self).__init__(astring)
          # We set the problem bounds (in this case equal for all components)
          self.set_bounds(5.12, 6.12)
          self.astring=astring
          self.other=other
          print('Init: ' + self.astring)
      def initialize(self,other,astring):
        self.astring = astring
        self.other = other
      # Reimplement the virtual method that defines the objective function.
      def _objfun_impl(self, x):
          print('obj: ' + self.astring + ' other: %s' % self.other)
          # Compute the sphere function
          f = sum([x[i] ** 2 +self.other for i in range(self.dimension)])
  
          # Note that we return a tuple with one element only. In PyGMO the objective functions
          # return tuples so that multi-objective optimization is also possible.
          return (f, )
  
      # Finally we also reimplement a virtual method that adds some output to the __repr__ method
      def human_readable_extra(self):
          return "\n\t Problem dimension: " + str(self.__dim)
def own2():
  prob = my_tricky_problem(5,'foo',100.0)  # Create a 10-dimensional problem
  #prob.initialize(5.5,'bar')  
  #algo = algorithm.bee_colony(gen=500)  # 500 generations of bee_colony algorithm
  algo = algorithm.pso(gen=2)
  isl = island(algo, prob, 1)  # Instantiate population with 20 individuals
  archi = archipelago(algo,prob,2,1)
  #isl.evolve(1)  # Evolve the island once
  #isl.join()
  #print(isl.population.champion.f)
  
  print min([isl.population.champion.f for isl in archi])
  archi.evolve(10)
  print min([isl.population.champion.f for isl in archi])
  
def pase_to_cpp():
  vesselgroup = h5files.open('/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5', 'r', search = False)['adaption/vessels_after_adaption']
  factory = getattr(parameterSetsAdaption, 'pyGmo')
  krebs.adaption.adaption_cpp.doAdaptionOptimization(vesselgroup, factory['adaption'], factory['calcflow'])

def serial():
  prob = problem.schwefel(dim = 200)
  algo = algorithm.de(gen = 500)
  isl = island(algo,prob,20)
  print isl.population.champion.f
  isl.evolve(10)
  print isl.population.champion.f
  
def parallel():
  prob = problem.schwefel(dim = 200)
  algo = algorithm.de(gen = 500)
  archi = archipelago(algo,prob,8,20)
  print min([isl.population.champion.f for isl in archi])
  archi.evolve(10)
  print min([isl.population.champion.f for isl in archi])
  
def other():
  prob = problem.lennard_jones(5)
  algo = algorithm.bee_colony(gen = 10) #instantiates artificial bee colony setting 10 generations for each algorithmic call
  topo = topology.ring() #defines an empty ring topology
  archi = archipelago(algo,prob,8,20,topology=topo) #connects 8 islands in a ring topology
  archi.topology.draw() #requires networkx to be installed

if __name__=='__main__':
  #abc = problem.py_example()
  if 0:
    serial();
  if 0:
    parallel();
  if 0:
    own2();
  if 0:
    other();
  if 1:
    pase_to_cpp();