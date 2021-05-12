import MDAnalysis as mda
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import cm
import glob

import warnings
warnings.filterwarnings("ignore")

pdbs = glob.glob("./*.pdb")

violin_data = []
violin_title = ['']

for pdb in pdbs[:]:
    u = mda.Universe(pdb)
    p = u.select_atoms('name PO4').positions
    p_z_mean = p[:, 2].mean()
    p_up = u.select_atoms('name PO4 and prop z > ' + str(p_z_mean)).positions
    p_down = u.select_atoms('name PO4 and prop z < ' + str(p_z_mean)).positions
    print('Lipids in the upper leaflet: ', p_up.shape[0])
    print('Lipids in the lower leaflet: ', p_down.shape[0])

    pope_xy = u.select_atoms('name PO4 and resname POPE').positions[:, :2]
    popg_xy = u.select_atoms('name PO4 and resname POPG').positions[:, :2]
    
    popg_up_xy = u.select_atoms('name PO4 and resname POPG and prop z > ' + str(p_z_mean)).positions[:, :2]
    popg_down_xy = u.select_atoms('name PO4 and resname POPG and prop z < ' + str(p_z_mean)).positions[:, :2]
    if len(u.select_atoms('not resname POPE POPG PW ION CL NA K')) > 0:
        anti_up_xy = np.vstack([i.atoms.center_of_geometry() for i in u.select_atoms('not resname POPE POPG PW ION CL NA K and prop z > ' + str(p_z_mean)).residues])[:, :2]
        anti_down_xy = np.vstack([i.atoms.center_of_geometry() for i in u.select_atoms('not resname POPE POPG PW ION CL NA K and prop z < ' + str(p_z_mean)).residues])[:, :2]
        anti_xy = np.vstack([i.atoms.center_of_geometry() for i in u.select_atoms('not resname POPE POPG PW ION CL NA K').residues])[:, :2]
    
    x_min, x_max = p[:, 0].min(), p[:, 0].max()
    y_min, y_max = p[:, 1].min(), p[:, 1].max()

    grid_x, grid_y = np.mgrid[x_min:x_max:0.5, y_min:y_max:0.5]

    grid_up_z = griddata(p_up[:, :2], p_up[:, 2], (grid_x, grid_y), method='nearest')
    grid_down_z = griddata(p_down[:, :2], p_down[:, 2], (grid_x, grid_y), method='nearest')
    grid_tn = grid_up_z - grid_down_z
    
    print('Min./max. thickness: ', grid_tn.min(), grid_tn.max())
    print('Av. thickness: ', grid_tn.mean(), ' ± ', grid_tn.std())
    
    pope_tn = griddata(p_up[:, :2], p_up[:, 2], pope_xy, method='nearest') - griddata(p_down[:, :2], p_down[:, 2], pope_xy, method='nearest')
    popg_tn = griddata(p_up[:, :2], p_up[:, 2], popg_xy, method='nearest') - griddata(p_down[:, :2], p_down[:, 2], popg_xy, method='nearest')
    if len(u.select_atoms('not resname POPE POPG PW ION CL NA K')) > 0:
        anti_tn = griddata(p_up[:, :2], p_up[:, 2], anti_xy, method='nearest') - griddata(p_down[:, :2], p_down[:, 2], anti_xy, method='nearest')
    
    print('POPE thickness ', pope_tn.mean(), ' ± ', pope_tn.std())
    violin_data.append(pope_tn)
    violin_title.append(pdb.split('/')[-1][:-4] + ' POPE')
    print('POPG thickness ', popg_tn.mean(), ' ± ', popg_tn.std())
    violin_data.append(popg_tn)
    violin_title.append(pdb.split('/')[-1][:-4] + ' POPG')
    if len(u.select_atoms('not resname POPE POPG PW ION CL NA K')) > 0:
        print('Antiseptic thickness ', anti_tn.mean(), ' ± ', anti_tn.std())
        violin_data.append(anti_tn)
        violin_title.append(pdb.split('/')[-1][:-4] + ' Antisept.')
    
    fig = plt.figure(figsize = (8, 4))
    
    ax1 = fig.add_subplot(121)
    im1 = ax1.imshow(grid_tn, extent=(x_min, x_max, y_min, y_max), origin='lower', interpolation='none', vmin = 30, vmax = 50, cmap = 'Greys')
    ax1.plot(popg_up_xy[:,0], popg_up_xy[:,1], 'r.', ms=7, alpha = 0.7)
    if len(u.select_atoms('not resname POPE POPG PW ION CL NA K')) > 0:
        ax1.plot(anti_up_xy[:,0], anti_up_xy[:,1], 'b.', ms=7, alpha = 0.7)
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])
    fig.colorbar(im1)
    ax1.set_title(pdb.split('/')[-1][:-4] + ' UP')
    
    ax2 = fig.add_subplot(122)
    im2 = ax2.imshow(grid_tn, extent=(x_min, x_max, y_min, y_max), origin='lower', interpolation='none', vmin = 30, vmax = 50, cmap = 'Greys')
    ax2.plot(popg_down_xy[:,0], popg_down_xy[:,1], 'r.', ms=7, alpha = 0.7)
    if len(u.select_atoms('not resname POPE POPG PW ION CL NA K')) > 0:
        ax2.plot(anti_down_xy[:,0], anti_down_xy[:,1], 'b.', ms=7, alpha = 0.7)
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])
    fig.colorbar(im2)
    ax2.set_title(pdb.split('/')[-1][:-4] + ' DOWN')

    fig.tight_layout()

    plt.show()

fig_viol = plt.figure(figsize = (10,4))
ax = fig_viol.add_axes([0,0,1,1])

parts = ax.violinplot(violin_data, showmeans = False, showextrema=False)
for pc in parts['bodies'][::1]:
    pc.set_color('deepskyblue')
    pc.set_alpha(0.5)

for pc in parts['bodies'][::2]:
    pc.set_color('orangered')
    pc.set_alpha(0.5)

for pc in parts['bodies'][::3]:
    pc.set_color('black')
    pc.set_alpha(0.5)

# set style for the axes
ax.set_xticks(np.arange(len(violin_title)) )
ax.set_xticklabels(violin_title)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(15)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(15)
plt.xticks(rotation=90)
plt.ylabel(r'Thickness, $\AA$',  fontsize=15)
plt.grid()
fig_viol.tight_layout()
plt.show()