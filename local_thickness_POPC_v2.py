import MDAnalysis as mda
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import cm
import glob

import warnings
warnings.filterwarnings("ignore")

#pdbs = glob.glob("./*.pdb")
pdbs = ['system.pdb']

for pdb in pdbs[:]:
    u = mda.Universe(pdb)
    p = u.select_atoms('name PO4').positions
    p_z_mean = p[:, 2].mean()
    p_up = u.select_atoms('name PO4 and prop z > ' + str(p_z_mean)).positions
    p_down = u.select_atoms('name PO4 and prop z < ' + str(p_z_mean)).positions
    print('Lipids in the upper leaflet: ', p_up.shape[0])
    print('Lipids in the lower leaflet: ', p_down.shape[0])

    popc_xy = u.select_atoms('name PO4 and resname POPC').positions[:, :2]
    
    popc_up_xy = u.select_atoms('name PO4 and resname POPC and prop z > ' + str(p_z_mean)).positions[:, :2]
    popc_down_xy = u.select_atoms('name PO4 and resname POPC and prop z < ' + str(p_z_mean)).positions[:, :2]
    
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()

    grid_x, grid_y = np.mgrid[x_min:x_max:0.5, y_min:y_max:0.5]

    grid_up_z = griddata(p_up[:, :2], p_up[:, 2], (grid_x, grid_y), method='nearest')
    grid_down_z = griddata(p_down[:, :2], p_down[:, 2], (grid_x, grid_y), method='nearest')
    grid_tn = grid_up_z - grid_down_z
    
    print('Min./max. thickness: ', grid_tn.min(), grid_tn.max())
    print('Av. thickness: ', grid_tn.mean(), ' ± ', grid_tn.std())
        
    fig = plt.figure(figsize = (8, 4))
    
    ax1 = fig.add_subplot(121)
    im1 = ax1.imshow(grid_tn, extent=(x_min, x_max, y_min, y_max), origin='lower', interpolation='none', vmin = 30, vmax = 50, cmap = 'Greys')
    ax1.plot(popc_up_xy[:,0], popc_up_xy[:,1], 'r.', ms=7, alpha = 0.7)
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])
    fig.colorbar(im1)
    ax1.set_title(pdb.split('/')[-1][:-4] + ' UP')
    
    ax2 = fig.add_subplot(122)
    im2 = ax2.imshow(grid_tn, extent=(x_min, x_max, y_min, y_max), origin='lower', interpolation='none', vmin = 30, vmax = 50, cmap = 'Greys')
    ax2.plot(popc_down_xy[:,0], popc_down_xy[:,1], 'r.', ms=7, alpha = 0.7)
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])
    fig.colorbar(im2)
    ax2.set_title(pdb.split('/')[-1][:-4] + ' DOWN')

    fig.tight_layout()

    plt.show()