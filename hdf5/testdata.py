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
        f = h5py.File('ifrog.hdf5', 'r')
        plt.matshow(f[filename1], aspect='auto', vmin=1000, vmax=16000)
        plt.show()
        f.close()
    if command in ['poscoro']:
        #print(len(f[filename2]))
        index = -1 
        f = h5py.File('ifrog.hdf5', 'r+')
        corodata = f[filename2]
        for i in f[filename2]:
            #print(i)
            index = index + 1 
            #print(index)
            #Magic number to set data around t=0
            magicnum = 9858.6325
            i = (i + magicnum) * (10e-6) * (2/(3*10e8))
            corodata[index] = i
            ##stepvalue(um)*10^-6*2/3*10^8
            #-9858.6325 = t=0
            #print(i)
        index = -1
        for i in f[filename1]:
            index = index + 1
            #print(index)
            wl = -0.81967 * index + 440
            print(wl)
            #pix num to WL
            #y=-0.081967x+440
           # #x = pixnum
            #y = wl (nm)
        f.close()
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
