"""
This plots figure used in publication for ball model
Plots are saved in figure folder

Author: Arnab Majumdar
Date: 24.06.2025
"""
import sys
import tarfile
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import ListedColormap
import h5py
import paths

# ignore warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

###################################
# (DE) Funktionsbereich
# (ENG) Functions
###################################

# analytical SAS function
def J1(x):
    """
    function desc:
    spherical bessel function
    """
    if x==0:
        y = 0
    else:
        y = (np.sin(x)-x*np.cos(x))/x**2
    return y

def ball (qmax,qmin,Npts,scale,bg,sld,sld_sol,rad):
    """
    function desc:
    analytical SAS pattern of sphere
    """
    vol=(4/3)*np.pi*rad**3
    # SLD unit 10^-5 \AA^-2
    del_rho=sld-sld_sol
    q_arr=np.linspace(qmin,qmax,Npts)
    FormFactor=np.zeros(len(q_arr))
    for q_idx, q_cur in enumerate(q_arr):
        if q_cur==0:
            # Form factor unit 10^-5 \AA
            FormFactor[q_idx]=vol*del_rho
        else:
            # Form factor unit 10^-5 \AA
            FormFactor[q_idx]=3*vol*del_rho*J1(q_cur*rad)/(q_cur*rad)
    # Intensity unit 10^-10 \AA^2
    Iq_arr = (scale)*np.abs(FormFactor)**2+bg
    return Iq_arr, q_arr

##############################################
# (DE) Henunterläde input Daten von tar ordner
# (ENG) Download input data from tar file
##############################################

# open tarfile
tar_path = paths.data / "dataset.tar.xz"
# extract data folder
with tarfile.open(tar_path, "r:xz") as tar:
    for member in tar.getmembers():
        if member.name.startswith("dataset/data/"):
            tar.extract(member, path=paths.data)
# working dir (or data folder) path
working_dir = paths.data / "dataset"

##################
# (DE) Input Daten
# (ENG) Input data
##################

### struct gen ###
# box side lengths (float values)
length_a= 40.
length_b=40.
length_c=40.
# number of cells in each direction (int values)
nx= 40
ny= 40
nz= 40
# element details
el_type='lagrangian'
el_order=1
# calculate mid point of structure (simulation box)
mid_point=np.array([length_a/2, length_b/2, length_c/2])

### sim gen ###
# model name
sim_model= 'ball'
# simulation parameters
dt= 1.
t_end= 0.
n_ensem=1
# calculate time array
t_arr=np.arange(0,t_end+dt, dt)

### model param ###
# model params
ball_rad=10
ball_sld=2
qclean_sld=0
# dir name for model param
model_param_dir_name = ('rad' + '_' + str(ball_rad) +
                         '_' + 'sld' + '_' + str(ball_sld)
                           + '_' + 'qclean_sld' + '_' + str(qclean_sld)
                           ).replace('.', 'p')

### scatt cal ###
# decreitization params
# number of categories and method of categorization
num_cats= [3, 101, 3]
method_cats=['direct', 'direct', 'extend']
# scatt_cal params
sig_file='signal.h5'
scan_vec=np.array([1, 0, 0])
Q_range=np.array([0., 1.])
num_points=100
num_orientation=10

"""
read folder structure
"""

# folder structure
## mother folder name (working dir + )
### save length values as strings
### decimal points are replaced with p
length_a_str=str(length_a).replace('.','p')
length_b_str=str(length_a).replace('.','p')
length_c_str=str(length_a).replace('.','p')

### save num_cell values as strings
nx_str=str(nx)
ny_str=str(ny)
nz_str=str(nz)

### save element_order values as strings
el_order_str='lagrangian_' + str(el_order)


mother_dir_name = length_a_str+'_' + length_b_str+'_' + length_c_str\
      + '_' + nx_str + '_' + ny_str+'_' + nz_str + '_' + el_order_str
data_dir=os.path.join(working_dir, 'data')
mother_dir = os.path.join(data_dir, mother_dir_name)

# read structure info
data_file=os.path.join(mother_dir, 'structure/struct.h5')

# folder name for simulation
sim_dir=os.path.join(mother_dir, 'simulation')

# folder name for model
model_dir_name= (sim_model + '_tend_' + str(t_end) + '_dt_' + str(dt) \
    + '_ensem_' + str(n_ensem)).replace('.','p')
model_dir=os.path.join(sim_dir,model_dir_name)

# folder name for model with particular run param
model_param_dir=os.path.join(model_dir,model_param_dir_name)

# create folder for figure (one level up from data folder)
# figure_dir=os.path.join(mother_dir, '../../figure/')
figure_dir=paths.figures
os.makedirs(figure_dir, exist_ok=True)
## folder for this suit of figures
plot_dir= figure_dir #os.path.join(figure_dir, sim_model +'_paper')
# os.makedirs(plot_dir, exist_ok=True)
# plott setting
# plot according to categorization methods
fig, ax = plt.subplots(figsize=(7, 5))
markers=['o', '^', 's']
colors=['m', 'r', 'k']
ms_arr=[3, 5, 3]

for i,t in enumerate(t_arr):
    # time_dir name
    t_dir_name=f't{i:0>3}'
    t_dir=os.path.join(model_param_dir, t_dir_name)
    for idx, (method_cat, num_cat) in enumerate(zip(method_cats, num_cats)):
        # method_cat=method_cats[idx]
        # num_cat=num_cats[idx]

        # retrieve scattering setting for this particular case
        scatt_settings='cat_' + method_cat + '_' + str(num_cat) + 'Q_' \
            + str(Q_range[0]) + '_' + str(Q_range[1]) + '_' + 'orien_' + '_' + str(num_orientation)
        scatt_settings=scatt_settings.replace('.', 'p')

        # ball geometry
        vol_ball=(4/3)*np.pi*ball_rad**3
        vol_box=length_a*length_b*length_c

        # volume for normalization
        vol_norm=vol_box

        # numerical intensity
        ## read I vs Q signal file
        Iq_data_file_name='Iq_{0}.h5'.format(scatt_settings)
        Iq_data_file=os.path.join(t_dir,Iq_data_file_name)
        Iq_data=h5py.File(Iq_data_file,'r')
        Iq=Iq_data['Iq'][:] # unit fm^2
        Iq=Iq/vol_norm*10**2 # unit (fm^2 / \AA^3) * 10^2 = cm^-1
        q=Iq_data['Q'][:]
        Iq_data.close()

        # loglog plot
        ax.loglog(q, Iq,color=colors[idx], linestyle='', marker=markers[idx],
                   markersize=ms_arr[idx],
                     label='Cat method: ' + method_cat + ', ' + str(num_cat) + ' categories' )

    # ananlytical intensity
    ## Intensity unit 10^-10 \AA^2
    Iq_ana,q_ana=ball(qmax=np.max(q),qmin=np.min(q),Npts=100,
                scale=1,bg=0,sld=ball_sld,sld_sol=0,rad=ball_rad)
    ## Normalize by volume
    ## (Before * 10**2) Intensity unit 10^-10 \AA^-1 = 10 ^-2 cm^-1
    ## (after * 10**2) Intensity unit cm^-1
    Iq_ana = (Iq_ana / vol_norm) * 10**2

    # loglog plot
    ax.loglog(q_ana, Iq_ana, 'b', label='Analytical calculation', zorder=-10)
    # ax.loglog(q, Iq,'r', linestyle='', marker='o', markersize=3, label='Numerical calculation')

    # plot formatting
    ## legend
    ax.legend(fontsize=14)
    ## labels
    ax.set_xlabel(r'Q [$\mathrm{\AA}^{-1}$]', fontsize=14)
    ax.set_ylabel(r'I(Q) [$\mathrm{cm}^{-1}$]', fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ## SANS upper boundary Q=1 \AA^-1
    ax.set_xlim(right=1)

    # read node sld
    # save pseudo atom info
    ensem_dir_name='ensem{0:0>3}'.format(0)
    ensem_dir=os.path.join(t_dir, ensem_dir_name)
    scatt_cal_dir_name='scatt_cal_' + scatt_settings
    scatt_cal_dir=os.path.join(ensem_dir, scatt_cal_dir_name)
    scatt_cal_data_file_name='scatt_cal.h5'
    scatt_cal_data_file=os.path.join(scatt_cal_dir, scatt_cal_data_file_name)
    scatt_cal_data=h5py.File(scatt_cal_data_file,'r')
    node_pos=scatt_cal_data['node_pos'][:]
    node_sld=scatt_cal_data['node_sld'][:]
    pseudo_pos=scatt_cal_data['pseudo_pos'][:]
    pseudo_b=scatt_cal_data['pseudo_b'][:]
    pseudo_b_cat_val=scatt_cal_data['pseudo_b_cat_val'][:]
    pseudo_b_cat_idx=scatt_cal_data['pseudo_b_cat_idx'][:]
    scatt_cal_data.close()

    # determine sld min and max for plotting
    sld_min=np.min(node_sld,0)
    sld_max=np.max(node_sld,0)

    if el_type=='lagrangian':
        num_node_x=el_order*nx+1
        num_node_y=el_order*ny+1
        num_node_z=el_order*nz+1
    # plotting node SLD
    ## cutting at z = cut_frac * length_z
    cut_frac=0.5
    node_pos_3d=node_pos.reshape(num_node_x, num_node_y, num_node_z, 3)
    z_idx= np.floor(cut_frac*(num_node_z)).astype(int)
    z_val=node_pos_3d[0, 0, z_idx , 2]
    ## figure specification
    plot_file_name='Iq_ball.pdf'
    plot_file=os.path.join(plot_dir,plot_file_name)
    #fig_ins, ax_ins = plt.subplots(figsize=(5, 5))
    ax_ins = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0, 0.365, 0.45, 0.45),
                         bbox_transform=ax.transAxes)
    ## image plot
    ### .T is required to exchange x and y axis
    ### origin is 'lower' to put it in lower left corner
    node_sld_3d=node_sld.reshape(num_node_x, num_node_y, num_node_z)
    cmap = ListedColormap(['white', '#440154'])
    img = ax_ins.imshow(node_sld_3d[:,:,z_idx].T,
                         extent=[0, length_a, 0, length_b],
                           origin='lower', vmin=sld_min, vmax=sld_max,
                             interpolation='bilinear', cmap=cmap)

    ax_ins.axis('off')

    ## add mesh
    cell_x=length_a/nx
    cell_y=length_a/ny
    for idx1 in range(nx):
        for idx2 in range(ny):
            rect_center_x=idx1*cell_x
            rect_center_y=idx2*cell_y
            ax_ins.add_patch(Rectangle((rect_center_x, rect_center_y), cell_x, cell_y,
                                        edgecolor='k', facecolor='none', linewidth=0.5))

    plt.savefig(plot_file, format='pdf')
    plt.close(fig)



"""
This plots figure used in publication for box model
Plots are saved in figure folder

Author: Arnab Majumdar
Date: 24.06.2025
"""
import sys
import os
import tarfile
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py
import paths

# ignore warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

# analytical SAS function
def arrange_order(a,b,c):
    zusammen=[a,b,c]
    zusammen.sort()
    return zusammen[0], zusammen[1], zusammen[2]

def gauss_legendre_double_integrate(func, domain1, domain2, deg):
    x, w = np.polynomial.legendre.leggauss(deg)
    xgrid, ygrid=np.meshgrid(x,x)
    x=xgrid.reshape((1,np.size(xgrid)))
    y=ygrid.reshape((1,np.size(ygrid)))
    wx, wy=np.meshgrid(w,w)
    w1=wx.reshape((1,np.size(wx)))
    w2=wy.reshape((1,np.size(wy)))
    s1 = (domain1[1] - domain1[0])/2
    a1 = (domain1[1] + domain1[0])/2
    s2 = (domain2[1] - domain2[0])/2
    a2 = (domain2[1] + domain2[0])/2
    return np.sum(s1*s2*w1*w2*func(s1*x + a1,s2*y + a2))

def box (qmax,qmin,Npts,scale,bg,sld,sld_sol,length_a,length_b,length_c):
    length_a, length_b, length_c= arrange_order(length_a,length_b,length_c)
    vol_box=length_a*length_b*length_c
    # SLD unit 10^-5 \AA^-2
    del_rho_box=sld-sld_sol
    q_arr=np.linspace(qmin,qmax,Npts) 
    Aq_arr=np.zeros(len(q_arr))
    for i in range(len(q_arr)):
        q=q_arr[i]
        func=lambda alpha, psi:\
            (del_rho_box*vol_box*(np.sinc((1/np.pi)*q*length_a/2*np.sin(alpha)*np.sin(psi)))*\
            (np.sinc((1/np.pi)*q*length_b/2*np.sin(alpha)*np.cos(psi)))*\
            (np.sinc((1/np.pi)*q*length_c/2*np.cos(alpha))))**2*\
            np.sin(alpha)
            
        psi_lim=np.pi
        alpha_lim=np.pi/2 
        # Amplitude unit (\AA^3 * 10^-5 \AA^-2)^2 = 10^-10 \AA^2
        Aq_arr[i]=(1/psi_lim)*gauss_legendre_double_integrate(func,[0, alpha_lim],[0, psi_lim],76)
    Iq_arr = scale*Aq_arr # Intensity unit 10^-10 \AA^2
    return Iq_arr, q_arr

##############################################
# (DE) Henunterläde input Daten von tar ordner
# (ENG) Download input data from tar file
##############################################

# open tarfile
tar_path = paths.data / "dataset.tar.xz"
# extract data folder
with tarfile.open(tar_path, "r:xz") as tar:
    for member in tar.getmembers():
        if member.name.startswith("dataset/data/"):
            tar.extract(member, path=paths.data)
# working dir (or data folder) path
working_dir = paths.data / "dataset"

"""
Input data
"""
### struct gen ###
# box side lengths (float values) 
length_a=40. 
length_b=40. 
length_c=40.
# number of cells in each direction (int values)
nx=40 
ny=40 
nz=40 
# element details
el_type='lagrangian'
el_order=1
# calculate mid point of structure (simulation box)
mid_point=np.array([length_a/2, length_b/2, length_c/2])


### sim gen ###
sim_model='box'
# simulation parameters
dt=1.
t_end=0.
n_ensem=1
# calculate time array
t_arr=np.arange(0,t_end+dt, dt)

### model param ###
# model params
box_sld=2
qclean_sld=0
# dir name for model param
model_param_dir_name = ('sld' + '_' + str(box_sld) + 
                        '_' + 'qclean_sld' + '_' + str(qclean_sld)
                        ).replace('.', 'p')

### scatt_cal ###
# decreitization params
# number of categories and method of categorization
num_cat=3
method_cat='extend'
sig_file='signal.h5'
scan_vec=np.array([1, 0, 0])
Q_range=np.array([0., 1.]) # mention integers as float
num_points=100
num_orientation_arr=[10, 50, 200]
# scatt settings will be defined later

"""
read folder structure
"""

# folder structure
## mother folder name
### save length values as strings
### decimal points are replaced with p
length_a_str=str(length_a).replace('.','p')
length_b_str=str(length_a).replace('.','p')
length_c_str=str(length_a).replace('.','p')

### save num_cell values as strings
nx_str=str(nx)
ny_str=str(ny)
nz_str=str(nz)

### save element_order values as strings
el_order_str='lagrangian_' + str(el_order)

mother_dir_name = length_a_str+'_' + length_b_str+'_' + length_c_str\
      + '_' + nx_str + '_' + ny_str + '_' + nz_str + '_' + el_order_str
data_dir=os.path.join(working_dir, 'data')
mother_dir = os.path.join(data_dir, mother_dir_name)

# read structure info
data_file=os.path.join(mother_dir, 'structure/struct.h5')

# folder name for simulation
sim_dir=os.path.join(mother_dir, 'simulation')

# folder name for model
model_dir_name= (sim_model + '_tend_' + str(t_end) + '_dt_' + str(dt) \
    + '_ensem_' + str(n_ensem)).replace('.','p')
model_dir=os.path.join(sim_dir,model_dir_name)

# folder name for model with particular run param
model_param_dir=os.path.join(model_dir,model_param_dir_name)

# create folder for figure (one level up from data folder)
# figure_dir=os.path.join(mother_dir, '../../figure/')
# os.makedirs(figure_dir, exist_ok=True)
## folder for this suit of figures
plot_dir= paths.figures # os.path.join(figure_dir, sim_model + '_paper')
os.makedirs(plot_dir, exist_ok=True)
# plott setting 
# plotting according to number of orientations
fig, ax = plt.subplots(figsize=(7, 5))
markers=['^', 's', 'o']
colors=['g', 'b', 'r']
ms_arr=[3, 3, 3]

for i in range(len(t_arr)):
    t=t_arr[i]
    # time_dir name
    t_dir_name='t{0:0>3}'.format(i)
    t_dir=os.path.join(model_param_dir, t_dir_name)

    for num_orien_idx, num_orientation in enumerate(num_orientation_arr):
        # retrieve scattering setting for this particular case
        scatt_settings='cat_' + method_cat + '_' + str(num_cat) + 'Q_' \
            + str(Q_range[0]) + '_' + str(Q_range[1]) + '_' + 'orien_' + '_' + str(num_orientation)
        scatt_settings=scatt_settings.replace('.', 'p')

        # box geometry
        vol_box=length_a*length_b*length_c

        # volume for normalization
        vol_norm=vol_box

        # numerical intensity
        ## read I vs Q signal file
        Iq_data_file_name='Iq_{0}.h5'.format(scatt_settings) 
        Iq_data_file=os.path.join(t_dir,Iq_data_file_name)
        Iq_data=h5py.File(Iq_data_file,'r')
        Iq=Iq_data['Iq'][:] # unit fm^2
        Iq=Iq/vol_norm*10**2 # unit (fm^2 / \AA^3) * 10^2 = cm^-1
        q=Iq_data['Q'][:]
        Iq_data.close()
    
        # loglog plot
        ax.loglog(q, Iq,color=colors[num_orien_idx], linestyle='', 
                  marker=markers[num_orien_idx], markersize=ms_arr[num_orien_idx],
                   label= 'Number of orientations = {0}'.format(num_orientation))

    # ananlytical intensity 
    ## Intensity unit 10^-10 \AA^2
    Iq_ana,q_ana= box(qmax=np.max(q),qmin=np.min(q),Npts=100,
                scale=1,bg=0,sld=box_sld,sld_sol=0,
                length_a=length_a, length_b=length_b, length_c=length_c)
    ## Normalize by volume
    ## (Before * 10**2) Intensity unit 10^-10 \AA^-1 = 10 ^-2 cm^-1
    ## (after * 10**2) Intensity unit cm^-1
    Iq_ana = (Iq_ana / vol_norm) * 10**2
    
    ax.loglog(q_ana, Iq_ana, 'gray', label= 'Analytical calculation')

    # plot formatting
    ## legend
    ax.legend()
    ## labels
    ax.set_xlabel(r'Q [$\mathrm{\AA}^{-1}$]')
    ax.set_ylabel(r'I(Q) [$\mathrm{cm}^{-1}$]')
    ## SANS upper boundary Q=1 \AA^-1
    ax.set_xlim(right=1)

    # read node sld
    # save pseudo atom info
    ensem_dir_name='ensem{0:0>3}'.format(0)
    ensem_dir=os.path.join(t_dir, ensem_dir_name)
    scatt_cal_dir_name='scatt_cal_' + scatt_settings
    scatt_cal_dir=os.path.join(ensem_dir, scatt_cal_dir_name)
    scatt_cal_data_file_name='scatt_cal.h5'
    scatt_cal_data_file=os.path.join(scatt_cal_dir, scatt_cal_data_file_name)
    scatt_cal_data=h5py.File(scatt_cal_data_file,'r')
    node_pos=scatt_cal_data['node_pos'][:]
    node_sld=scatt_cal_data['node_sld'][:]
    pseudo_pos=scatt_cal_data['pseudo_pos'][:]
    pseudo_b=scatt_cal_data['pseudo_b'][:]
    pseudo_b_cat_val=scatt_cal_data['pseudo_b_cat_val'][:]
    pseudo_b_cat_idx=scatt_cal_data['pseudo_b_cat_idx'][:]
    scatt_cal_data.close()

    # determine sld min and max for plotting
    sld_min=np.min(node_sld,0)
    sld_max=np.max(node_sld,0)

    if el_type=='lagrangian':
        num_node_x=el_order*nx+1
        num_node_y=el_order*ny+1
        num_node_z=el_order*nz+1
    # plotting node SLD
    ## cutting at z = cut_frac * length_z
    cut_frac=0.5
    node_pos_3d=node_pos.reshape(num_node_x, num_node_y, num_node_z, 3)
    z_idx= np.floor(cut_frac*(num_node_z)).astype(int)
    z_val=node_pos_3d[0, 0, z_idx , 2]
    ## figure specification
    plot_file_name='Iq_box_orien.pdf'
    plot_file=os.path.join(plot_dir,plot_file_name)
    #fig_ins, ax_ins = plt.subplots(figsize=(5, 5))
    ax_ins = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0, 0.35, 0.5, 0.5),
                        bbox_transform=ax.transAxes)
    ## image plot
    ### .T is required to exchange x and y axis 
    ### origin is 'lower' to put it in lower left corner 
    node_sld_3d=node_sld.reshape(num_node_x, num_node_y, num_node_z)
    cmap = ListedColormap(['white', '#440154'])
    img = ax_ins.imshow(node_sld_3d[:,:,z_idx].T, 
                    extent=[0, length_a, 0, length_b], 
                    origin='lower', vmin=0, vmax=sld_max, 
                    interpolation='bilinear', cmap=cmap)
    
    ax_ins.axis('off')

    ## add mesh
    cell_x=length_a/nx
    cell_y=length_a/ny
    for idx1 in range(nx):
        for idx2 in range(ny):
            rect_center_x=idx1*cell_x
            rect_center_y=idx2*cell_y
            ax_ins.add_patch(Rectangle((rect_center_x, rect_center_y), cell_x, cell_y, 
                                    edgecolor='k', facecolor='none', linewidth=0.5))
    plt.savefig(plot_file, format='pdf')
    plt.close(fig)

"""
This plots figure used in publication for ball in box model
Plots are saved in figure folder

Author: Arnab Majumdar
Date: 24.06.2025
"""
import sys
import tarfile
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py
import fitz
import paths

# ignore warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)

# analytical SAS function
def J1(x):
    if x==0:
        return 0
    else:
        return (np.sin(x)-x*np.cos(x))/x**2
    
def arrange_order(a,b,c):
    zusammen=[a,b,c]
    zusammen.sort()
    return zusammen[0], zusammen[1], zusammen[2]

def gauss_legendre_double_integrate(func, domain1, domain2, deg):
    x, w = np.polynomial.legendre.leggauss(deg)
    xgrid, ygrid=np.meshgrid(x,x)
    x=xgrid.reshape((1,np.size(xgrid)))
    y=ygrid.reshape((1,np.size(ygrid)))
    wx, wy=np.meshgrid(w,w)
    w1=wx.reshape((1,np.size(wx)))
    w2=wy.reshape((1,np.size(wy)))
    s1 = (domain1[1] - domain1[0])/2
    a1 = (domain1[1] + domain1[0])/2
    s2 = (domain2[1] - domain2[0])/2
    a2 = (domain2[1] + domain2[0])/2
    return np.sum(s1*s2*w1*w2*func(s1*x + a1,s2*y + a2))

def ball_in_box_ecc(qmax,qmin,Npts,scale,scale2,bg,sld_box, sld_ball,sld_sol,
                    length_a,length_b,length_c,radius, origin_shift):
    length_a, length_b, length_c=arrange_order(length_a,length_b,length_c)
    vol_box=length_a*length_b*length_c
    vol_ball=(4/3)*np.pi*radius**3
    # SLD unit 10^-5 \AA^-2
    del_rho_box=sld_box-sld_sol
    del_rho_ball=sld_ball-sld_sol
    q_arr=np.linspace(qmin,qmax,Npts) 
    Aq_arr=np.zeros(len(q_arr))
    for i in range(len(q_arr)):
        q=q_arr[i]
        if q==0:
            func=lambda alpha, psi:\
                np.abs(del_rho_box*vol_box*(np.sinc((1/np.pi)*q*length_a/2*np.sin(alpha)*np.sin(psi)))*\
                (np.sinc((1/np.pi)*q*length_b/2*np.sin(alpha)*np.cos(psi)))*\
                (np.sinc((1/np.pi)*q*length_c/2*np.cos(alpha))))**2*\
                    np.sin(alpha)
        else:
            func=lambda alpha, psi:\
                np.abs(del_rho_box*vol_box*(np.sinc((1/np.pi)*q*length_a/2*np.sin(alpha)*np.sin(psi)))*\
                (np.sinc((1/np.pi)*q*length_b/2*np.sin(alpha)*np.cos(psi)))*\
                (np.sinc((1/np.pi)*q*length_c/2*np.cos(alpha)))-\
                scale2*3*vol_ball*(del_rho_box-del_rho_ball)*\
                    J1(q*radius)/(q*radius)*\
                    np.vectorize(complex)(np.cos(q * np.sin(alpha) * np.sin(psi) * origin_shift[0]\
                                    + q * np.sin(alpha) * np.cos(psi) * origin_shift[1]\
                                        + q *  np.cos(psi) * origin_shift[2]),
                            np.sin(q * np.sin(alpha) * np.sin(psi) * origin_shift[0]\
                                    + q * np.sin(alpha) * np.cos(psi) * origin_shift[1]\
                                        + q *  np.cos(psi) * origin_shift[2])))**2*\
                    np.sin(alpha)
        
        psi_lim=np.pi
        alpha_lim=np.pi/2
        # Amplitude unit (\AA^3 * 10^-5 \AA^-2)^2 = 10^-10 \AA^2 
        Aq_arr[i]=(1/psi_lim)*gauss_legendre_double_integrate(func,[0, alpha_lim],[0, psi_lim],76)
    Iq_arr = scale*Aq_arr + bg # Intensity unit 10^-10 \AA^2    
    return Iq_arr, q_arr

##############################################
# (DE) Henunterläde input Daten von tar ordner
# (ENG) Download input data from tar file
##############################################

# open tarfile
tar_path = paths.data / "dataset.tar.xz"
# extract data folder
with tarfile.open(tar_path, "r:xz") as tar:
    for member in tar.getmembers():
        if member.name.startswith("dataset/data/"):
            tar.extract(member, path=paths.data)
# extract ppt folder
with tarfile.open(tar_path, "r:xz") as tar:
    for member in tar.getmembers():
        if member.name.startswith("dataset/ppt/"):
            tar.extract(member, path=paths.data)
# working dir (or data folder) path
working_dir = paths.data / "dataset"

##################
# (DE) Input Daten
# (ENG) Input data
##################

### struct gen ###
# box side lengths (float values)
length_a=40.
length_b=40.
length_c=40.
# number of cells in each direction (int values)
nx=40
ny=40
nz=40
# element details
el_type='lagrangian'
el_order=1
# calculate mid point of structure (simulation box)
mid_point=np.array([length_a/2, length_b/2, length_c/2])

### sim_gen ###
sim_model='bib_ecc'
# simulation parameters
dt=1.
t_end=0.
n_ensem=1
# calculate time array
t_arr=np.arange(0,t_end+dt, dt)

### model_param ###
# model params
ball_rad=10
sld_in=2
sld_out=1
ecc_vec_arr=[np.array([5, 5, 5]), np.array([0, 0, 0])]
qclean_sld=0
# dir name for model param
# will be defined later

### scatt_cal ###
num_cat=3
method_cat='extend'
sig_file='signal.h5'
scan_vec=np.array([1, 0, 0])
Q_range=np.array([0., 1.])
num_points=100
num_orientation=300
# scatt settengs
scatt_settings='cat_' + method_cat + '_' + str(num_cat) + 'Q_' \
    + str(Q_range[0]) + '_' + str(Q_range[1]) + '_' + 'orien__' + str(num_orientation)
scatt_settings=scatt_settings.replace('.', 'p')

"""
read folder structure
"""
# folder structure
## mother folder name
### save length values as strings
### decimal points are replaced with p
length_a_str=str(length_a).replace('.','p')
length_b_str=str(length_a).replace('.','p')
length_c_str=str(length_a).replace('.','p')

### save num_cell values as strings
nx_str=str(nx)
ny_str=str(ny)
nz_str=str(nz)

### save element_order values as strings
el_order_str='lagrangian_' + str(el_order)

mother_dir_name = length_a_str+'_' + length_b_str+'_' + length_c_str\
      + '_' + nx_str + '_' + ny_str + '_' + nz_str + '_' + el_order_str
data_dir=os.path.join(working_dir, 'data')
mother_dir = os.path.join(data_dir, mother_dir_name)

# read structure info
data_file=os.path.join(mother_dir, 'structure/struct.h5')

# folder name for simulation
sim_dir=os.path.join(mother_dir, 'simulation')

# folder name for model
model_dir_name= (sim_model + '_tend_' + str(t_end) + '_dt_' + str(dt) \
    + '_ensem_' + str(n_ensem)).replace('.','p')
model_dir=os.path.join(sim_dir,model_dir_name)

# folder name for model with particular run param
# will be defined later

# create folder for figure (one level up from data folder)
figure_dir=paths.figures #os.path.join(mother_dir, '../../figure/')
# os.makedirs(figure_dir, exist_ok=True)
## folder for this suit of figures
plot_dir=figure_dir # os.path.join(figure_dir, sim_model + '_paper')
# os.makedirs(plot_dir, exist_ok=True)
# plott setting 
# plotting according to ecc vector
fig, ax = plt.subplots(figsize=(7, 5))
markers=['o', '^']
num_colors=['r', 'b']
ms_arr=[3, 3]
ana_linestyles=['-', '--']

for idx in range(len(ecc_vec_arr)):   
    # dir name for model param
    ecc_vec=ecc_vec_arr[idx]
    model_param_dir_name = ('rad' + '_' + str(ball_rad) + '_' +
                            'sld_in' + '_' + str(sld_in) + '_' +
                            'sld_out' + '_' + str(sld_out) + '_' +
                            'x' + '_' + str(ecc_vec[0]) + '_' + 
                            'y' + '_' + str(ecc_vec[1]) + '_' + 
                            'z' + '_' + str(ecc_vec[2]) + '_' + 
                            'qclean_sld' + '_' + str(qclean_sld)).replace('.', 'p')

    # folder name for model with particular run param
    model_param_dir=os.path.join(model_dir,model_param_dir_name)
    for i in range(len(t_arr)):
        t=t_arr[i]
        # time_dir name
        t_dir_name='t{0:0>3}'.format(i)
        t_dir=os.path.join(model_param_dir, t_dir_name)

        # box geometry
        vol_box=length_a*length_b*length_c

        # volume for normalization
        vol_norm=vol_box

        # numerical intensity
        ## read I vs Q signal file
        Iq_data_file_name='Iq_{0}.h5'.format(scatt_settings) 
        Iq_data_file=os.path.join(t_dir,Iq_data_file_name)
        Iq_data=h5py.File(Iq_data_file,'r')
        Iq=Iq_data['Iq'][:] # unit fm^2
        Iq=Iq/vol_norm*10**2 # unit (fm^2 / \AA^3) * 10^2 = cm^-1
        q=Iq_data['Q'][:]
        Iq_data.close()

        # ananlytical intensity 
        ## Intensity unit 10^-10 \AA^2
        Iq_ana,q_ana= ball_in_box_ecc(qmax=np.max(q),qmin=np.min(q),Npts=100,
                    scale=1,scale2=1,bg=0,sld_box=sld_out, sld_ball=sld_in, sld_sol=0,
                    length_a=length_a, length_b=length_b, length_c=length_c, radius=ball_rad, 
                    origin_shift=ecc_vec)
        
        ## Normalize by volume
        ## (Before * 10**2) Intensity unit 10^-10 \AA^-1 = 10 ^-2 cm^-1
        ## (after * 10**2) Intensity unit cm^-1
        Iq_ana = (Iq_ana / vol_norm) * 10**2
        
        # loglog plot
        ax.loglog(q, Iq,color=num_colors[idx], linestyle='', marker=markers[idx], markersize=ms_arr[idx], 
                  label= r'Numerical calculation: $ \overrightarrow{R}_{\mathrm{ecc}} = $' + '({0}, {1}, {2})'.format(ecc_vec[0], ecc_vec[1], ecc_vec[2]))
        ax.loglog(q_ana, Iq_ana, 'gray', linestyle=ana_linestyles[idx], 
                  label= r'Analytical calculation: $ \overrightarrow{R}_{\mathrm{ecc}} = $' + '({0}, {1}, {2})'.format(ecc_vec[0], ecc_vec[1], ecc_vec[2]))

        # plot formatting
        ## legend
        ax.legend(loc='lower left')
        ## labels
        ax.set_xlabel(r'Q [$\mathrm{\AA}^{-1}$]')
        ax.set_ylabel(r'I(Q) [$\mathrm{cm}^{-1}$]')
        ## SANS upper boundary Q=1 \AA^-1
        ax.set_xlim(right=1)

ax_ins1 = inset_axes(ax, width="100%", height="100%", bbox_to_anchor=(0, 0.4, 0.5, 0.5),
        bbox_transform=ax.transAxes)

ax_ins1.axis('off')

# Load the PDF and select the first page
doc = fitz.open(working_dir / 'ppt/ballinboxEccentric.pdf')
page = doc.load_page(0)

# Define crop area (x0, y0, x1, y1) in points
crop_rect = fitz.Rect(250, 50, 550, 350)  # adjust as needed

# Render the cropped area to a pixmap (image)
pix = page.get_pixmap(clip=crop_rect, dpi=200)
img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)

# Plot with matplotlib
ax_ins1.imshow(img)

# save plot
plot_file_name='Iq_bib_ecc.pdf'
plot_file=os.path.join(plot_dir,plot_file_name)
plt.savefig(plot_file, format='pdf')
plt.close(fig)