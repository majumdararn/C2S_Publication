import numpy as np
from scipy.integrate import quad, simpson
from scipy.special import spherical_jn
import matplotlib.pyplot as plt
import h5py
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import paths
import tarfile

# funktionsbereich

# Gaussian f(k)
def f(k):
    """
    Input frequency distrobution f(k)
    1. normal distribution is assumed
    """
    return (
        np.exp(-(k - k0)**2 / (2 * sigma**2))
        / (sigma * np.sqrt(2 * np.pi))
    )

# g(r)
def g(r):
    """
    Calculation of g(r)
    """
    def integrand(k):
        return (
            4 * np.pi
            * k**2
            * f(k)
            * spherical_jn(0, k * r)
        )

    return quad(integrand, 0, np.inf)[0]

########################
# analytical calculation
########################

# Parameters (analytical calculation)
# input
k0 = 0.32
sigma = 0.02
phi_1=0.5
phi_2=1-phi_1
del_rho = 1
# CRW calculation
r_arr = np.linspace(0, 250, 200)
## step 1
g_arr = np.array([g(r) for r in r_arr])
g_arr /= g_arr[0]
## step 2
Gamma_arr = (2 / np.pi) * np.arcsin(g_arr)
## step 3
Q_arr = np.linspace(0.001, 1, 200)
Iq = np.zeros(len(Q_arr))
for i, Q in enumerate(Q_arr):
    integrand = (
        phi_1 * phi_2 * del_rho**2
        * 4 * np.pi
        * r_arr**2
        * Gamma_arr
        * spherical_jn(0, Q * r_arr)
    )
Iq[i] = simpson(integrand, x=r_arr)

###############################
# numerical calcularion
###############################

# Step 1. read Cont2Sas output
tar_path = paths.data / "crw.tar.gz"
# extract data folder
with tarfile.open(tar_path, "r:gz") as tar:
    tar.extractall(path=paths.data)
# with tarfile.open(tar_path, "r:gz") as tar:
#     for member in tar.getmembers():
#         if member.name.startswith("dataset/data/"):
#             tar.extract(member, path=paths.data)

# smaller box
h5_file_1=paths.data / 'Iq_box_250.h5'
file_1=h5py.File(h5_file_1,'r')
# vol_norm=628.3185**3
# Iq_num=file['Iq'][:]/vol_norm*10**2
Iq_num_1=file_1['Iq'][:]/1.55e14
q_num_1=file_1['Q'][:]
file_1.close()

# smaller box
h5_file_2=paths.data / 'Iq_box_600.h5'
file_2=h5py.File(h5_file_2,'r')
# vol_norm=628.3185**3
# Iq_num=file['Iq'][:]/vol_norm*10**2
Iq_num_2=file_2['Iq'][:]/1.55e14
q_num_2=file_2['Q'][:]
file_2.close()

# sld for smaller box
h5_file_3=paths.data / 'sld_box_250.h5'
file_3=h5py.File(h5_file_3,'r')
# Iq_num=file['Iq'][:]/vol_norm*10**2
sld=file_3['sld'][:]
file_3.close()
# reshape and undo r clean
sld_3d=sld.reshape((125+1, 125+1, 125+1))+1.5


# Step 2. plotting
# plot Iq (right)
fig, ax1 = plt.subplots(1,2,figsize=(12, 5))
ax1[1].plot(Q_arr, Iq, label='analytical')
ax1[1].plot(q_num_1, Iq_num_1, '*', ms=4, label='Smaller box')
ax1[1].plot(q_num_2, Iq_num_2, '^', ms=4, label='Bigger box')
ax1[1].axvline(x=k0, color='black', linestyle='--', label=f'Input Peak k0 ({k0})')
ax1[1].set_xlabel("Q")
ax1[1].set_ylabel("I(Q)")
ax1[1].set_xscale("log")
ax1[1].set_xlim([0.09,0.7])
ax1[1].set_ylim([0,1e-5])
ax1[1].legend()
ax1[1].grid()

# plot sld (left)
img = ax1[0].imshow(sld_3d[:,:,0].T,
                    extent=[0, 250, 0, 250],
                    origin='lower',
                    vmin=1, vmax=2,
                    interpolation='none')
# Get the colors corresponding to 1 and 2
cmap = img.cmap
norm = img.norm
legend_elements = [Patch(facecolor=cmap(norm(1)), edgecolor='black', label=r'$0^{\mathrm{th}}$ phase (SLD = $\beta_0$)'),
                   Patch(facecolor=cmap(norm(2)), edgecolor='black', label=r'$1^{\mathrm{st}}$ phase (SLD = $\beta_1$)'),
                   ]
# create legend
ax1[0].legend(handles=legend_elements, 
              loc='upper right',
              frameon=True, facecolor='white',
              framealpha=1.0)
ax1[0].set_xlabel(r'X $[\mathrm{\AA}]$')
ax1[0].set_ylabel(r'Y $[\mathrm{\AA}]$')
plt.savefig(paths.figures / 'crw.pdf', format='pdf')
