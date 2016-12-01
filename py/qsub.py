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
#! /usr/bin/env python2
# -*- coding: utf-8 -*-

import os,sys
import cStringIO
import base64
import cPickle
import subprocess
import math
import time
import re

def printClientInfo():
  import socket
  print('run this on client: %s' % socket.gethostname())
  print('invoked by: %s' % sys._getframe(1).f_code.co_name)
  print('cwd: %s' % os.getcwd())

__all__ = [ 'parse_args', 'submit', 'exe', 'func', 'is_client']

''' globals '''
defaultMemory = '1MB'
global goodArgumentsQueue
goodArgumentsQueue = {} #empty namespace
installed_queue_system = 'foo'
is_client = False


def parse_args(argv):
  import argparse
  parserQueue = argparse.ArgumentParser(prog='qsub',description='Queing system parser.')
  memory_option = parserQueue.add_argument('-m', '--memory', help= 'Memory assigned by the queing system', type=str, default = '2GB')
  global defaultMemory
  defaultMemory = memory_option.default
  parserQueue.add_argument('--q-local', help= ' Do not submit to queue, even if queuing system is pressent', default=False, action='store_true')
  parserQueue.add_argument('--q-dry', help= 'Do not run but print configuration to be submitted', default=False, action='store_true')
  parserQueue.add_argument('--q-verbose', help= 'more output', default=False, action='store_true')
  localgoodArgumentsQueue, otherArgumentsQueue = parserQueue.parse_known_args()  
  global goodArgumentsQueue
  goodArgumentsQueue = localgoodArgumentsQueue



def identify_installed_submission_program_():
  qsys = os.environ.get('QUEUE_SUBMISSION_PROGRAM', None)
  if qsys:
    if qsys in ('qsub','sbatch'):
      return qsys
    else:
      print('qsub.py WARNING: env. var. QUEUE_SUBMISSION_PROGRAM is set to unsupported "%s"' % qsys)
  # check if qsub is there
  try:
    subprocess.Popen(['qsub', '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
  except OSError, e:
    if goodArgumentsQueue.q_verbose:
      print ('qsub.py: we are not using qsub.')
  else:
    return 'qsub'
  #check if slurm works
  try:
    subprocess.Popen(['srun','--version'],stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
  except OSError, e:
    if goodArgumentsQueue.q_verbose:
      print('qsub.py: we are not using slurm')
  else:
    return 'sbatch'
  if goodArgumentsQueue.q_verbose: #opts_['run_verbose']:
    print('qsub.py: no supported queue system found.')
  return None

  
def determine_submission_program_():
  returnstring = 'return'
  if goodArgumentsQueue.q_local:
    returnstring = 'run_locally'
  else:
    returnstring = identify_installed_submission_program_()
    if returnstring is None:
      print('Warning: no supported queueing system found -> run locally')
      returnstring = 'run_locally'
  return returnstring



def fmtDate_(days, hours):
    '''turn any number of days and hours into hours h which are
       0 <= d < 24. Also handle fractional days and hours. Fractional
       hours are simply rounded up to full hours.'''
    if days is None:
        days = 0
    if hours is None:
        hours = 0
    i_days = int(days)
    h_days = (days-i_days)*24.0
    d_hours = int(hours/24)
    h_hours = hours - d_hours*24.0
    days  = i_days + d_hours
    hours = int(math.ceil(h_days + h_hours))
    return days, hours
    

def write_directives_qsub_(f,name=None, days=None, hours=None, num_cpus=1, outdir=None, export_env=False, jobfiledir=None, change_cwd=False, dependsOnJob = None):
  mem =goodArgumentsQueue.memory  
  print >>f, '#PBS -j oe'
  if jobfiledir and not outdir: #DEPRECATED
    outdir = jobfiledir
  if outdir is not None:
    print >>f, '#PBS -o %s' % (outdir)
  if name:
    print >>f, '#PBS -N %s' % name
  if num_cpus > 1:
    print >>f, '#PBS -l nodes=1:ppn=%i' % num_cpus
  if days or hours:
    days, hours = fmtDate_(days, hours)
    print >>f, '#PBS -l walltime=%i:00:00' % max(1, hours + days*24)
  mem =goodArgumentsQueue.memory
  if mem:
    if re.match(r'^\d+(kB|MB|GB)$', mem) is None:
      raise RuntimeError('mem argument needs integer number plus one of kB, MB, GB')
    print >>f, '#PBS -l mem=%s' % mem
  if export_env:
    print >>f, '#PBS -V'
  if dependsOnJob:
    print >>f, '#PBS -W depend=afterok:%s' % dependsOnJob
    
    
def write_directives_slurm_(f, name=None, days=None, hours=None, num_cpus=1, outdir=None, export_env=False, jobfiledir=None, change_cwd=False):
  #print >>f, '#PBS -j oe'
  #if jobfiledir and not outdir: #DEPRECATED
  #  outdir = jobfiledir
  #if outdir is not None:
  #  print >>f, '#PBS -o %s' % (outdir)
  #print >>f, 'cd $SLURM_SUBMIT_DIR'
  
  mem =goodArgumentsQueue.memory  
  if name:
    print >>f, '#SBATCH --job-name=%s' % name
  if num_cpus == 1:
    print >>f, '#SBATCH --cpus-per-task=1'
    print >>f, '#SBATCH --ntasks=1'
    print >>f, '#SBATCH --partition=onenode'
  if num_cpus > 1:
    print >>f, '#SBATCH --cpus-per-task=%i' % num_cpus
    print >>f, '#SBATCH --ntasks=1'
    print >>f, '#SBATCH --nodes=1'
    print >>f, '#SBATCH --partition=onenode'
    print >>f, '#SBATCH --ntasks-per-node=1'
    #print >>f, 'export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE'
    #print >>f, '#SBATCH --resv-ports'
  if days or hours:
    days, hours = fmtDate_(days, hours)
    print >>f, '#SBATCH --time=%i-%i:00:00' % (days, hours)
  if mem:
    if re.match(r'^\d+(kB|MB|GB)$', mem) is None:
      raise RuntimeError('mem argument needs integer number plus one of kB, MB, GB')
    print >>f, '#SBATCH --mem=%s' % mem
  #if export_env:
  #  print >>f, '#PBS -V'


def submit_(interpreter, submission_program, script):
  global opts_
  # determine how to run
  if submission_program == 'run_locally':
    submission_program = interpreter
  # verbose output
  if goodArgumentsQueue.q_verbose or goodArgumentsQueue.q_dry:
    print (' submission to %s ' % submission_program).center(30,'-')
    print script
    print ''.center(30,'-')
  # run stuff
  if not goodArgumentsQueue.q_dry:
    time.sleep(1.1)
    if submission_program == 'python': # running python script with python locally?!! We can do it like so
      exec script in dict(), dict()
    else: # all the other cases go like so!
      subprocess.call("%s <<EOFQSUB\n%s\nEOFQSUB" % (submission_program, script), shell=True)
      #subprocess.check_output("%s <<EOFQSUB\n%s\nEOFQSUB" % (submission_program, script), shell=True)

pyfuncscript_ = """\
import imp as imp__
import cPickle as cPickle__
import base64 as base64__
import sys
sys.path.append('%s')
qsub = __import__('qsub', globals(), locals())
setattr(qsub, 'is_client', True)
%s
from your_module__ import *
args__, kwargs__ = cPickle__.loads(base64__.urlsafe_b64decode('%s'))
%s(*args__, **kwargs__)
"""
pychangecwd_ = """\
import os as os__
try: os__.chdir('%s')
except KeyError: pass
"""


class Func(object):
  '''
    A Wrapper around a function call. This allows you to submit python function calls to qsub.

    It does this by submitting a small 'outer' script which imports the module in which the
    function resides and calls the function by its name with the supplied parameters.
    Therefore your script must be file on the filesystem. It does not work in interactive mode.

    WARNING: Your script is NOT COPIED! It is loaded as is at the time of execution on the
    cluster node. This is inconsistent with the normal behavior of qsub!

    You give an object of this kind to the submit function in order to run it.
  '''
  interpreter = 'python'
  def __init__(self, func, *args, **kwargs):
    self.func = func
    self.args = (args, kwargs)
    self.name_hint = func.__name__

  def generate_script(self, qsubopts):
    func, (args, kwargs) = self.func, self.args
    # generate the python script
    a = base64.urlsafe_b64encode(cPickle.dumps((args,kwargs)))
    # doesn't work so well because if the module is not the main module,
    # the __module__ variable can point to something else than the file 
    # where func is defined. Example: /py/krebsjobs/__init__.pyc instead
    # of /krebsjobs/submitDetailedO2.py!
    #   m = __import__(func.__module__)
    #   filename = os.path.abspath(m.__file__)
    # next try:
    filename = func.__globals__['__file__']  # seems simple enough?! '__globals__' just contains the global variables that the function sees.
    filename = os.path.abspath(filename) # don't want to get in trouble because of changing work dir
    functionname = func.__name__
    if filename.endswith('.pyc'):
      load_string = "your_module__ = imp__.load_compiled('your_module__', '%s')" % filename
    else:
      load_string = "your_module__ = imp__.load_source('your_module__', '%s')" % filename
    s = pyfuncscript_ % (os.path.abspath(os.path.dirname(__file__)), load_string, a, functionname)
    if qsubopts.get('change_cwd', False):
      s = (pychangecwd_ % os.getcwd()) + s
    return s

func = Func




class Exe(object):
  '''
    A wrapper around a system call to execute a program.

    You give an object of this kind to the submit function in order to run it.

    cmd - string or sequence of objects which can be converted to strings
  '''
  interpreter = 'sh'
  def __init__(self, *cmds):
    for i, cmd in enumerate(cmds):
      if isinstance(cmd, str) or not hasattr(cmd, '__iter__'):
        cmds[i] = [cmd]
    self.cmds = cmds
    self.name_hint = None
  def generate_script(self, qsubopts):
    lines = [
      ' '.join(list(str(q) for q in cmd)) for cmd in self.cmds
      ]
    if qsubopts.get('change_cwd', False):
      lines = [ 'cd %s' % os.getcwd() ] + lines
    return '\n'.join(lines)

exe = Exe

class SrunExe(object):
  '''
    A wrapper around a system call to execute a program.

    You give an object of this kind to the submit function in order to run it.

    cmd - string or sequence of objects which can be converted to strings
  '''
  interpreter = 'sh'
  def __init__(self, *cmds):
    for i, cmd in enumerate(cmds):
      if isinstance(cmd, str) or not hasattr(cmd, '__iter__'):
        cmds[i] = [cmd]
    self.cmds = cmds
    self.name_hint = None
  def generate_script(self, qsubopts):
    lines = [
      ' '.join(list(str(q) for q in cmd)) for cmd in self.cmds
      ]
    if qsubopts.get('change_cwd', False):
      lines = [ 'cd %s' % os.getcwd() ] + lines
    return '\n'.join(lines)

sexe = SrunExe    


def submit_qsub(obj, submission_program, **qsubopts):
  '''
    Run stuff with qsub. See Exe and Func.
    Qsub parameters are special. WARNING: parameters might change!

    qsubopts can contain:
      name -> job name (string)
      days -> runtime in days (any number)
      num_cpus -> processors per job (int)
      mem -> memory usage (string)
      outdir -> directory where stdout/err are written to (string)
      export_env -> export current environment to launched script (bool)
  '''
  # default name
  if not 'name' in qsubopts:
    if obj.name_hint:
      qsubopts['name'] = obj.name_hint
    else:
      qsubopts['name'] = 'unnamed'
  # interpreter string
  cases = {
    'sh' : '#!/bin/sh',
    'python' : ('#!/usr/bin/env python%i'  % sys.version_info.major)
  }
  first_line = cases[obj.interpreter]
  # add qsub stuff + python script
  f = cStringIO.StringIO()
  print >>f, first_line
  write_directives_qsub_(f, **qsubopts)
  print >>f, obj.generate_script(qsubopts)
  return submit_(obj.interpreter, submission_program, f.getvalue())
  

def submit_slurm(obj, submission_program, **slurmopts):
  '''
    Run stuff with slurm. See Exe and Func.
    Slurm parameters are special. WARNING: parameters might change!
  '''
  # default name
  if not 'name' in slurmopts:
    if obj.name_hint:
      slurmopts['name'] = obj.name_hint
    else:
      slurmopts['name'] = 'unnamed'
  # interpreter string
  cases = {
    'sh' : '#!/bin/sh',
    'python' : ('#!/usr/bin/env python%i'  % sys.version_info.major)
    #'python' : '#!/bin/sh'
  }
  first_line = cases[obj.interpreter]
  #print first_line
  # add qsub stuff + python script
  f = cStringIO.StringIO()
  print >>f, first_line
  write_directives_slurm_(f, **slurmopts)
  old=True
  if old:
    print >>f,obj.generate_script(slurmopts)
    return submit_(obj.interpreter,submission_program, f.getvalue())
  else:
    print >>f, 'echo "Maybe a stupid idea, but I will create a .py file now"'
    string_abuffer_name = 'buffer_%s.py' % time.time()
    if 'thierry_slurm_buffer' not in os.listdir('/localdisk/tmp'):
      os.mkdir('/localdisk/tmp/thierry_slurm_buffer')
      mypyscriptfile = open('/localdisk/tmp/thierry_slurm_buffer/%s' % string_abuffer_name,"w+")
      mypyscriptfile.write(obj.generate_script(slurmopts))
      mypyscriptfile.close()
      if obj.interpreter == 'python':
        #print >>f, ('srun --resv-port mpirun /usr/bin/python2 %s' % mypyscriptfile.name)
        print >>f, ('/usr/bin/python2 %s' % mypyscriptfile.name)
      if obj.interpreter == 'sh':
        #print >>f, ('srun --resv-port mpirun /bin/sh %s' % mypyscriptfile.name)
        print >>f, ('/bin/sh %s' % mypyscriptfile.name)
        print(f.getvalue())
        #print(f.getvalue())
        #print >>f, obj.generate_script(slurmopts)
        #print(f.getvalue())
        submit_('sh', submission_program, f.getvalue())


def submit(obj, **qsubopts):
    if 'mem' in qsubopts and goodArgumentsQueue.memory == defaultMemory:
      print('Memory setting provided by program')
      global goodArgumentsQueue
      goodArgumentsQueue.memory = qsubopts.pop('mem')
    if not goodArgumentsQueue.memory == defaultMemory:
      if 'mem' in qsubopts:
        print('OVERRIDE Memory setting provided by program')
        qsubopts.pop('mem') #pop and not storing!!!
    print(goodArgumentsQueue.memory)
    print(defaultMemory)
    
    prog = determine_submission_program_()
    if prog == 'sbatch':
      submit_slurm(obj, prog, **qsubopts)
    else:
      submit_qsub(obj, prog, **qsubopts)
