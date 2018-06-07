/**
This file is part of tumorcode project.
(http://www.uni-saarland.de/fak7/rieger/homepage/research/tumor/tumor.html)

Copyright (C) 2017 Thierry Fredrich

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
*/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_SPAWNS 2

int main( int argc, char *argv[] )
{
  int np = NUM_SPAWNS;
  int errcodes[NUM_SPAWNS];
  MPI_Comm parentcomm, intercomm;

  MPI_Init( &argc, &argv );
  MPI_Comm_get_parent( &parentcomm );
  if (parentcomm == MPI_COMM_NULL) {
    MPI_Comm_spawn( "hello_mpi", MPI_ARGV_NULL, np, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &intercomm, errcodes );
    printf("I'm the parent.\n");
  } else {
    printf("I'm the spawned.\n");
  }
  fflush(stdout);
  MPI_Finalize();
  return 0;
}
