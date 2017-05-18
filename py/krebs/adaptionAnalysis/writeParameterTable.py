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
if __name__ == '__main__':
  import os.path, sys
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../..'))
  sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))

import krebsjobs.parameters.parameterSetsBulkTissueTumor as parameterSets_tum
from krebsjobs.submitAdaption import create_auto_dicts
from krebsjobs.submitAdaption import typelist
#from submitAdaption import 
#create_auto_dicts(param_group)
#typelist = 'typeA typeB typeC typeD typeE typeF typeG typeH typeI'
type_to_label_dict={"typeA": "RC1","typeB": "RC2","typeC": "RC3","typeD": "RC4","typeE": "RC5","typeF": "RC6","typeG": "RC7","typeH": "RC8","typeI": "RC9"}

def create_main_adaption_table():
  #### input
  parameter_set_name = 'p3d_bigger_H2'
  
  adaptionParamsList = []
  ### change here for different types    
  #parameter_set_name_from_pipe = parameter_set_name #begins with auto_
  #parameter_set_name =  parameter_set_name_from_pipe[5:]
  print('Found param identifiyer: %s' % parameter_set_name)
  type_to_paramset = create_auto_dicts(parameter_set_name+'_')
#  \newcommand{\wallstress}{\tau_w}
#\newcommand{\refflow}{Q_{ref}}
#\newcommand{\conducRef}{S_0}
#\newcommand{\metabConst}{k_m}
#\newcommand{\conducConst}{k_c}
#\newcommand{\shrinkConst}{k_s}
  with open('/localdisk/thierry/high_end_research_docs/adaption/tables/adaption_parameters.tex','wb') as texfile:
    #texfile.write('\\begin{tabular}{lllllll}\n ');
    texfile.write('\\begin{tabular}{ccccccc}\n ');
    texfile.write(' & $k_m$ & $\\refflow$ & $L$ & $\\conducConst$ & $\\conducRef$ & $\\shrinkConst$ \\\\ \n \\hline \n');
    for t in typelist:
      print(t)
      print(type_to_paramset[t])
      adaptionParams_this_type = getattr(parameterSets_tum, type_to_paramset[t])
      #print(adaptionParams_this_type['adaption'])
      #texfile.write('%s \\\\ %s\\midrule \n' % ('blub', 'bla'))
      texfile.write('%s & %0.2f & %0.0f & %0.f & %0.2f & %0.0f &%0.2f \\\\ \n' % (type_to_label_dict[t], adaptionParams_this_type['adaption']['k_m'], adaptionParams_this_type['adaption']['Q_refdot'],adaptionParams_this_type['adaption']['cond_length'],adaptionParams_this_type['adaption']['k_c'],adaptionParams_this_type['adaption']['S_0'],adaptionParams_this_type['adaption']['k_s']))

    texfile.write('\\hline\n \\end{tabular}\n')  


        
if __name__ == '__main__':
  create_main_adaption_table()
