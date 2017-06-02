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
import identifycluster
if identifycluster.name == 'snowden':
  import matplotlib
  matplotlib.use('Agg') 

import adaption
import krebsjobs.parameters.parameterSetsAdaption as parameterSetsAdaption

import operator
import random
from scoop import futures

import numpy

from deap import base
from deap import creator
from deap import tools

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Particle", list, fitness=creator.FitnessMax, speed=list, 
    smin=None, smax=None, best=None, adaptionParameters=None)


def generate(size, pmin, pmax, smin, smax):
    part = creator.Particle(random.uniform(pmin, pmax) for _ in range(size)) 
    part.speed = [random.uniform(smin, smax) for _ in range(size)]
    part.smin = smin
    part.smax = smax
    factory = getattr(parameterSetsAdaption, "deap_test")
    vfile_name = "/localdisk/thierry/vessel_trees_better/my_chosen/PSO_data_vessels-large_2d-typeE-17x1L600-sample05_adption_p_human_guess.h5"
    grp_name = "adaption/vessels_after_adaption"
    factory['adaption'].update(
      vesselFileName = vfile_name,
      vesselGroupName = grp_name,
      )
    part.adaptionParameters=factory
    return part

def updateParticle(part, best, phi1, phi2):
    u1 = (random.uniform(0, phi1) for _ in range(len(part)))
    u2 = (random.uniform(0, phi2) for _ in range(len(part)))
    v_u1 = futures.map(operator.mul, u1, futures.map(operator.sub, part.best, part))
    v_u2 = futures.map(operator.mul, u2, futures.map(operator.sub, best, part))
    part.speed = list(futures.map(operator.add, part.speed, futures.map(operator.add, v_u1, v_u2)))
    for i, speed in enumerate(part.speed):
        if speed < part.smin:
            part.speed[i] = part.smin
        elif speed > part.smax:
            part.speed[i] = part.smax
    part[:] = list(futures.map(operator.add, part, part.speed))

def main():
  toolbox = base.Toolbox()
  toolbox.register("particle", generate, size=3, pmin=0.5, pmax=4, smin=-1, smax=1)
  toolbox.register("population", tools.initRepeat, list, toolbox.particle)
  toolbox.register("update", updateParticle, phi1=2.0, phi2=2.0)
  toolbox.register("evaluate", adaption.doit_optimize_deap)
  
  
  pop = toolbox.population(n=100)
  stats = tools.Statistics(lambda ind: ind.fitness.values)
  stats.register("avg", numpy.mean)
  stats.register("std", numpy.std)
  stats.register("min", numpy.min)
  stats.register("max", numpy.max)

  logbook = tools.Logbook()
  logbook.header = ["gen", "evals"] + stats.fields

  GEN = 2
  best = None
  print("starting loop")
  for g in range(GEN):
      fitnesses = list(futures.map(toolbox.evaluate,pop))
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
      
  print("best pop:")
  print(pop)
  print("best")
  print(best)
  return pop, logbook, best

if __name__ == "__main__":
  main()
  
  