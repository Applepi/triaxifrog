import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt


filename1 = sys.argv[1]
filename2 = sys.argv[2]
command = sys.argv[3]


testdata = np.loadtxt(filename1, delimiter="\t")
testdata = testdata.astype(int)
posdata = np.loadtxt(filename2, delimiter="\t")

with h5py.File('ifrog.hdf5', 'w') as f:
    dset = f.create_dataset(filename1, data=testdata)
    dset = f.create_dataset(filename2, data=posdata)

hf = h5py.File('ifrog.hdf5', 'r')
    
with h5py.File('ifrog.hdf5', 'r') as f:
   data = f[filename1]
   stepdata = f[filename2]

print(command)

if command in ['plot']:
   plt.matshow(testdata, aspect='auto', vmin=1000, vmax=16000)
   plt.show()
    

#stepdata = hf.get('20190717_880_bbo2_pos')

#hf = h5py.File('test.hdf5', 'r')
#hf.close()

### Needs work to get working
#for i in stepdata:
#    print(stepdata[i])
#    value = stepdata[i]
#    newvalue = stepdata[i] * 10e-6 * 2 / 3 * 10e8
#    stepdata[i] = newvalue



##stepvalue(um)*10^-6*2/3*10^8
