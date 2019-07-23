import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import pypret

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
        corodata = f[filename2+'coro']
        corowl = f['corowl']
        plt.matshow(f[filename1], aspect='auto', vmin=1000, vmax=16000)
        plt.show()
        f.close()
    if command in ['mesh']:
        f = h5py.File('ifrog.hdf5', 'r')
        corodata = f[filename2+'coro'][::]
        corowl = f['corowl']
        coroinds = corodata[::].argsort()
        sorted_coro = corodata[coroinds[:-1]]
        sorted_grid = f[filename1][::,::][:,coroinds[:-1]]
        

        ft = pypret.FourierTransform(256, dt=2.5e-15)
        pulse = pypret.Pulse(ft, 400e-9)
        pulse.spectrum = pypret.lib.gaussian(pulse.wl, x0=800e-9, sigma=20e-9)
        # print the accurate FWHM of the temporal intensity envelope
        # propagate it through 1cm of BK7 (remove first ord)
        phase = np.exp(1.0j * pypret.material.BK7.k(pulse.wl) * 0.01)
        pulse.spectrum = pulse.spectrum * phase
        insertion = np.linspace(-0.025, 0.025, 128)  # insertion in m
        pnps = pypret.PNPS(pulse, "ifrog", "shg")
        # calculate the measurement trace
        pnps.calculate(pulse.spectrum, sorted_coro)
        original_spectrum = pulse.spectrum
        # and plot it
        #mesh = pypret.mesh_data.MeshData(f[filename1][::],  corowl[::], corodata[::])
        
        mesh = pypret.mesh_data.MeshData(sorted_grid[:,10000:],  corowl[::], sorted_coro[10000:])
        pypret.graphics.MeshDataPlot(mesh, plot=True)
        #ret = pypret.retrieval.retriever.BaseRetriever._retrieve_begin(mesh, pulse.specturm ,weights=None)
        
        ret = pypret.Retriever(pnps, "copra", verbose=True, maxiter=300)
        ret.retrieve(mesh, pulse.spectrum, weights=None)
        pypret.random_gaussian(pulse, 50e-15, phase_max=0.0)
        pypret.graphics.MeshDataPlot(mesh, plot=True)

        #pypret.mesh_data.MeshData(f[filename1][::], corodata[::]) 
        f.close()
    if command in ['poscoro']:
        index = -1 
        f = h5py.File('ifrog.hdf5', 'r+')
        corodata = f[filename2+'coro']
        corowl = f['corowl']
        for i in f[filename2]:
            index = index + 1 
            #Magic number to set data around t=0
            magicnum = 9858.6325
            i = (i + magicnum) * (10e-6) * (2/(3*10e8))
            corodata[index] = i
            ##stepvalue(um)*10^-6*2/3*10^8
            #-9858.6325 = t=0
        index = -1
        for i in f[filename1]:
            index = index + 1
            wl = -0.81967 * index + 440
            corowl[index] = wl +100
            #pix num to WL
            #y=-0.081967x+440
            #x = pixnum
            #y = wl (nm)
        f.close()
    if command in ['import']:
        testdata = np.loadtxt(filename1, delimiter="\t")
        testdata = testdata.astype(int)
        posdata = np.loadtxt(filename2, delimiter="\t")
        with h5py.File('ifrog.hdf5', 'w') as f:
            f.create_dataset(filename1, data=testdata[:,:-1])
            f.create_dataset(filename2, data=posdata)
            f.create_dataset(filename2+'coro', data=posdata)
            f.create_dataset('corowl', data=np.zeros((1021)))
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
