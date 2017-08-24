##############################################################################
# Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of Spack.
# Created by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647188
#
# For details, see https://github.com/llnl/spack
# Please also see the NOTICE and LICENSE files for our notice and the LGPL.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (as
# published by the Free Software Foundation) version 2.1, February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
# conditions of the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##############################################################################

#     SPACK support by Thierry Fredrich 8.8.17
#
# T.F.
# on my maschine this file resides at: 
# /localdisk/thierry/local/spack/var/spack/repos/builtin/packages/vbl-master/package.py
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install vbl-master
#
# You can edit this file again by typing:
#
#     spack edit vbl-master
#
# See the Spack documentation for more information on packaging.
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
from spack import *


class Tumorcode(CMakePackage):
    """FIXME: Put a proper description of your package here."""
    homepage = "http://www.uni-saarland.de/fak7/rieger/homepage/research/biological_physics/tumor/tumor.html"
    url      = "https://github.com/thierry3000/tumorcode/archive/master.zip"

    # version('2014-10-08', git='https://github.com/example-project/example.git',commit='9d38cd4e2c94c3cea97d0e2924814acc')
    version('develop', git='https://github.com/thierry3000/tumorcode.git',branch='spack')
    # python
    depends_on('py-h5py')
    depends_on('py-scipy')
    depends_on('py-numpy')
    depends_on('py-matplotlib')
    # others
    depends_on('boost+python')
    depends_on('trilinos')
    depends_on('hdf5')
    depends_on('vbl')

    def cmake_args(self):
        args = ['-DCMAKE_BUILD_TYPE=Release',
                '-DUSE_ADAPTION=OFF',
                '-DUSE_MILOTTI_MTS=ON',
                ]
        return args
    #def install(self, spec, prefix):
        # FIXME: Unknown build system
        # make()
    #    make('install')
