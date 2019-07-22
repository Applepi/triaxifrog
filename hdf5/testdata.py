import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

if len(sys.argv) == 4:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    command = sys.argv[3]
elif len(sys.argv) == 3:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
else:
    print("Please provide 2D spectrum data and 1d position data as aruments")



if len(sys.argv) == 4:
    print(command + ' starting')
    if command in ['plot']:
        
        plt.matshow(testdata, aspect='auto', vmin=1000, vmax=16000)
        plt.show()
    if command in ['poscoro']:
        for i in stepdata:
            print(i)
        hf.close()
    if command in ['import']:
        testdata = np.loadtxt(filename1, delimiter="\t")
        testdata = testdata.astype(int)
        posdata = np.loadtxt(filename2, delimiter="\t")
        with h5py.File('ifrog.hdf5', 'w') as f:
            dset = f.create_dataset(filename1, data=testdata)
            dset = f.create_dataset(filename2, data=posdata)
        print(command + ' complete')




#hfr = h5py.File('ifrog.hdf5', 'r')
#hfw = h5py.File('ifrog.hdf5', 'w')
    
#with h5py.File('ifrog.hdf5', 'r') as f:
#   data = f[filename1]
#   stepdata = f[filename2]

        
#stepdata = hf.get('20190717_880_bbo2_pos')

#hf = h5py.File('test.hdf5', 'r')
#hf.close()

### Needs work to get working
#for i in stepdata:
#    print(i)
#    value = stepdata[i]
#    newvalue = stepdata[i] * 10e-6 * 2 / 3 * 10e8
#    stepdata[i] = newvalue
    
##stepvalue(um)*10^-6*2/3*10^8
