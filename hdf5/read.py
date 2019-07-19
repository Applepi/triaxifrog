import h5py
import numpy as np

with h5py.File('random.hdf5', 'r') as f:
   data = f['default']
   print(min(data))
   print(max(data))
   print(data[:15])
