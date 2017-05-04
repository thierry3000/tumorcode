#ifndef TRILINOS_HELPER_H
#define TRILINOS_HELPER_H

#include <bits/types/time_t.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack_Utils.h>
#include <fstream>

void vector_print(Epetra_Vector &v);
void sparse_matrix_print(Epetra_CrsMatrix &mat);

#endif //TRILINOS_HELPER_H
