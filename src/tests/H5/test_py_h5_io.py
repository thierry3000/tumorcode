import h5py
import py_h5_tumorcode_io

if __name__ == '__main__':
  print(py_h5_tumorcode_io.greet())
  print("is mpi on?")
  print(py_h5_tumorcode_io.is_mpi_on())
  py_h5_tumorcode_io.run_file_test()
