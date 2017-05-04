#include "trilinos_helper.h"

void vector_print(Epetra_Vector &v)
{
  // thierrys sparse output
  time_t rawtime;
  struct tm*timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  std::ofstream myMatrixFile;
  std::string myString ("Epetra_Matrix");
  myString = myString + asctime(timeinfo) +".txt";
  myMatrixFile.open(myString.c_str());
  //flowsys.sys->Print(myMatrixFile);
  v.Print(myMatrixFile);
//     std::string myString2 ("Epetra_Sparsity_Pattern");
//     myString2 = myString2 + asctime(timeinfo) + ".ps";
//     Ifpack_PrintSparsity(mat,myString2.c_str());
  myMatrixFile.close();
}

void sparse_matrix_print(Epetra_CrsMatrix &mat)
{
  // thierrys sparse output
  time_t rawtime;
  struct tm*timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  std::ofstream myMatrixFile;
  std::string myString ("Epetra_Matrix");
  myString = myString + asctime(timeinfo) +".txt";
  myMatrixFile.open(myString.c_str());
  //flowsys.sys->Print(myMatrixFile);
  mat.Print(myMatrixFile);
  std::string myString2 ("Epetra_Sparsity_Pattern");
  myString2 = myString2 + asctime(timeinfo) + ".ps";
  Ifpack_PrintSparsity(mat,myString2.c_str());
  myMatrixFile.close();
}
