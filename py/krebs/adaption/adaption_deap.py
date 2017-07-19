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
import datetime
import identifycluster
if identifycluster.name == 'snowden':
  import matplotlib
  matplotlib.use('Agg') 

import adaption
import krebsjobs.parameters.parameterSetsAdaption as parameterSetsAdaption

import operator
import random
import copy
from scoop import futures, shared

import numpy

from deap import base
from deap import creator
from deap import tools

import h5py

creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Particle", list, fitness=creator.FitnessMax, speed=list, 
    smin=None, smax=None, pmin=None, pmax= None, best=None, adaptionParameters=None)
n = 560
GEN = 10

def generate(size, pmin, pmax, smin, smax):
    part = creator.Particle(random.uniform(pmin, pmax) for _ in range(size)) 
    part.speed = [random.uniform(smin, smax) for _ in range(size)]
    part.smin = smin
    part.smax = smax
    part.pmax = pmax
    part.pmin = pmin
    factory = getattr(parameterSetsAdaption, shared.getConst('adaptionParams'))
    myInt = int(shared.getConst('parameterListIndex'))
    if( shared.getConst('parameterListIndex') is not None ):
      factory = factory[myInt]
    vfile_name = shared.getConst('vesselInputFile')
    grp_name = shared.getConst('vessel_grp')
    factory['adaption'].update(
      vesselFileName = vfile_name,
      vesselGroupName = grp_name,
      )
    part.adaptionParameters=factory
    return part

def updateParticle(part, best, phi1, phi2):
    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = futures.map(operator.mul, u1, map(operator.sub, part.best, part))
    v_u2 = futures.map(operator.mul, u2, map(operator.sub, best, part))
    part.speed = list(futures.map(operator.add, part.speed, map(operator.add, v_u1, v_u2)))
    for i, speed in enumerate(part.speed):
        if speed < part.smin:
            part.speed[i] = part.smin
        elif speed > part.smax:
            part.speed[i] = part.smax
    part[:] = list(map(operator.add, part, part.speed))
    for (ind,v) in enumerate(part[:]):
      if v<part.pmin:
        part[ind]=part.pmin
      if v>part.pmax:
        part[ind]=part.pmax
        
toolbox = base.Toolbox()
toolbox.register("particle", generate, size=3, pmin=0.5, pmax=4, smin=-1, smax=1)
toolbox.register("population", tools.initRepeat, list, toolbox.particle)
toolbox.register("update", updateParticle, phi1=2.0, phi2=2.0)
toolbox.register("evaluate", adaption.doit_optimize_deap)

def main():
  
  pop = toolbox.population(n=n)
  stats = tools.Statistics(lambda ind: ind.fitness.values)
  stats.register("avg", numpy.mean)
  stats.register("std", numpy.std)
  stats.register("min", numpy.min)
  stats.register("max", numpy.max)

  logbook = tools.Logbook()
  logbook.header = ["gen", "evals"] + stats.fields

  
  best = None
  print("starting loop")
  for g in range(GEN):
      fitnesses = list(futures.map(toolbox.evaluate,pop))
      print("generation: %i " % g)
      for ind, fit in zip(pop,fitnesses):
        ind.fitness.values = fit
      for part in pop:  
        if not part.best or part.best.fitness < part.fitness:
          part.best = creator.Particle(part)
          part.best.fitness.values = part.fitness.values
        if not best or best.fitness < part.fitness:
          best = creator.Particle(part)
          best.fitness.values = part.fitness.values
      for part in pop:
        toolbox.update(part, best)
       #Gather all the fitnesses in one list and print the stats
      #logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
      #print(logbook.stream)
      
  print("best fitness.values :")
  print(best.fitness.values)
  print("best")
  print(best)
  return best

if __name__ == "__main__":
  ''' this uses a sytem queueing system'''
  import argparse
  parser = argparse.ArgumentParser(description='Compute adaption see Secomb model', formatter_class=argparse.ArgumentDefaultsHelpFormatter)  
  parser.add_argument('AdaptionParamSet')
  parser.add_argument('--listindex', type=int, help="index of value list" )
  parser.add_argument('--outputFileFolder', type=str, help="where to store the output, default is working directory")
  parser.add_argument('--fileName', type=str,help="name of vesselfile to optimize")
  goodArguments, otherArguments = parser.parse_known_args()
  print("running with %s" % goodArguments.AdaptionParamSet)
  
  if( str(goodArguments.AdaptionParamSet).startswith('value_li' ) ):
    bla = int(goodArguments.listindex)
    print("list index: %i chosen" % bla)
    shared.setConst(parameterListIndex= bla)
  #adaption parameters
  shared.setConst(adaptionParams=goodArguments.AdaptionParamSet)
  #file
  if( goodArguments.fileName ):
    vfile_name=goodArguments.fileName;
    vessel_grp ="vessels" #default place for vesselgenerator
  else:
    if 0: # 2d set
        vfile_name = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5"
    if 1: # 3d set
      vfile_name = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-9x11L600-sample13_adption_p_human_guess.h5"
    vessel_grp ="adaption/vessels_after_adaption"
  shared.setConst(vesselInputFile=vfile_name)
  #group
  shared.setConst(vessel_grp=vessel_grp)
  best = main()
  print(os.getcwd())
  with h5py.File('deap_results_%s.h5' % str(goodArguments.AdaptionParamSet)) as f:
  #f=h5py.File('deap_results_%s.h5' % str(goodArguments.AdaptionParamSet))
    currentTime = str(datetime.datetime.now().time())
    myGroup = f.create_group(str(goodArguments.listindex))
    myGroup.attrs.create("time", data=currentTime)
    myGroup.create_dataset('best' , data=best)
    myGroup.create_dataset('fitness values', data=best.fitness.values)
    myGroup.attrs.create("GEN", data=GEN)
    myGroup.attrs.create("n", data=n)
    myGroup.attrs.create("vfile", data=vfile_name)
    myGroup.attrs.create("vessel_grp", data=vessel_grp)
    myGroup.attrs.create("params", data=goodArguments.AdaptionParamSet)
    myGroup.attrs.create("paramListIndex", data = goodArguments.listindex)
    myGroup.attrs.create("type", data="variance opt")
  #f.close()