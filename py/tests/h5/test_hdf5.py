import h5py
import numpy as np

abc = np.random.randn(100);

def safe():
	with h5py.File('test.h5', 'w') as f:
		f.create_dataset('aTestSet', data=abc);
def read():
	with h5py.File('test.h5', 'r') as f:
		bla=f['aTestSet']
		print(bla)
		bla_as_numpy = np.asarray(bla)
		print(bla.shape)

if __name__ == '__main__':
	safe()
	read()
