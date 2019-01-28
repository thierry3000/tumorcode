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
#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
from os.path import join, dirname
sys.path.append(join(dirname(__file__),'.'))
import stat

import qsub


def ClientWorker(*args, **kwargs):
  print ("hello, i'm running on the cluster. my args are %s and %s" % (str(args), str(kwargs)))
  print ("my cwd is: %s" % os.getcwd())

if not qsub.is_client and __name__=='__main__':
  argv = qsub.parse_args(sys.argv)

  jobIDinPy = qsub.submit(qsub.func(ClientWorker, 'Test Argument', a_kw_arg = 9001),
                name = 'qsub-script-python-test',
                num_cpus = 1,
                days = 0,
                hours = 0.1,
                mem = '100MB',
                change_cwd = True)
  
  test_program_name = 'client-test-script.sh'
  with open(test_program_name, 'w') as f:
    f.write('''#! /bin/sh\n''')
    f.write('''echo "hello, i'm running on the cluster"\n''')
    f.write('''echo "my cwd is:"\n''')
    f.write('''pwd\n''')
    f.write('''echo "my job id was %i"''' %jobIDinPy)
  os.chmod(test_program_name, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)
  
  qsub.submit(qsub.exe(['ulimit -c unlimited'],
                       [join('.',test_program_name)]),
              name = 'qsub-script-executable-submission-test',
              mem = '100MB',
              days = 0,
              hours = 0.1,
              num_cpus = 1,
              change_cwd = True)
