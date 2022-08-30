###################################
## Utilities for utilityautomata ##
###################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Generate random density field

def randfield(nx,ny,cutoff):
    
    # return a smooth, randomly varying field between zero and 1. As cutoff
    # increases, the structure of the field becomes smaller wavelength. 
    f = np.fft.fftshift(np.fft.fft2(np.random.rand(ny,nx)))
    [Nx,Ny] = np.meshgrid(range(0,nx),range(0,ny))
    f[np.sqrt((Nx-nx/2)**2+(Ny-ny/2)**2)>cutoff]=0
    f = np.real(np.fft.ifft2(np.fft.ifftshift(f)))
    f = f-np.amin(f)
    f = f/np.amax(f)
    
    return f

