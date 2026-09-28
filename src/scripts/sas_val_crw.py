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
del_rho = 1
# CRW calculation
r_arr = np.linspace(0, 500, 500)
## step 1
g_arr = np.array([g(r) for r in r_arr])
g_arr /= g_arr[0]
## step 2
Gamma_arr = (2 / np.pi) * np.arcsin(g_arr)
## step 3
Q_arr = np.linspace(0.001, 1, 200)
Iq = np.zeros(len(Q_arr))
for i, Q in enumerate(Q_arr):
    # 0.5 * 0.5 signifies phi_1 * phi_2 (phi_1 = phi_2)
    integrand = (
        0.5 * 0.5 * del_rho**2
        * 4 * np.pi
        * r_arr**2
        * Gamma_arr
        * spherical_jn(0, Q * r_arr)
    )
    
    Iq[i] = simpson(integrand, x=r_arr)

# # Parameters
# k0 = 0.32
# sigma = 0.02


# # Gaussian f(k)
# def f(k):
#     return (
#         np.exp(-(k - k0)**2 / (2 * sigma**2))
#         / (sigma * np.sqrt(2 * np.pi))
#     )


# # g(r)
# def g(r):
#     def integrand(k):
#         return (
#             4 * np.pi
#             * k**2
#             * f(k)
#             * spherical_jn(0, k * r)
#         )

#     return quad(integrand, 0, np.inf)[0]


# # r grid
# r_arr = np.linspace(0, 250, 200)

# # Calculate g(r)
# g_arr = np.array([g(r) for r in r_arr])

# g_arr /= g_arr[0]

# # Gamma(r)
# Gamma_arr = (2 / np.pi) * np.arcsin(g_arr)


# # Q grid
# Q_arr = np.linspace(0.001, 1, 200)


# # I(Q)
# Iq = np.zeros(len(Q_arr))

# for i, Q in enumerate(Q_arr):

#     integrand = (
#         0.5 * 0.5 * 1**2
#         * 4 * np.pi
#         * r_arr**2
#         * Gamma_arr
#         * spherical_jn(0, Q * r_arr)
#     )
    
#     Iq[i] = simpson(integrand, x=r_arr)

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
vol_norm=250**3
# Iq_num=file['Iq'][:]/vol_norm*10**2 (1.3e7)
Iq_num_1=file_1['Iq'][:]/vol_norm
q_num_1=file_1['Q'][:]
file_1.close()

# bigger box
h5_file_2=paths.data / 'Iq_box_600.h5'
file_2=h5py.File(h5_file_2,'r')
vol_norm=600**3
# Iq_num=file['Iq'][:]/vol_norm*10**2 (2e8)
Iq_num_2=file_2['Iq'][:]/vol_norm
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
fs=10
unit_factor=1e-2 # converting to cm-1
fig, ax1 = plt.subplots(1,2,figsize=(15, 5))
ax1[1].plot(Q_arr, Iq * unit_factor, label='Analytical formula')
ax1[1].plot(q_num_1, Iq_num_1 * unit_factor, '*', ms=4, label=r'Numerical (L = 250 $\mathrm{\AA}$)')
ax1[1].plot(q_num_2, Iq_num_2 * unit_factor, '^', ms=4, label=r'Numerical (L = 600 $\mathrm{\AA}$)')
ax1[1].axvline(x=k0, color='black', linestyle='--', label=fr'$Q_{{\mathrm{{peak}}}} = {k0}\,\mathrm{{\AA}}^{{-1}}$')
ax1[1].set_xlabel(r"Q $[\mathrm{\AA}^{-1}]$", fontsize=fs)
ax1[1].set_ylabel(r"I(Q) $[\mathrm{cm}^{-1}]$", fontsize=fs)
ax1[1].set_xscale("log")
ax1[1].set_xlim([0.09,0.6])
ax1[1].set_ylim([0,7])
ax1[1].legend(fontsize=fs)
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
              framealpha=1.0, fontsize=fs)
ax1[0].set_xlabel(r'X $[\mathrm{\AA}]$', fontsize=fs)
ax1[0].set_ylabel(r'Y $[\mathrm{\AA}]$', fontsize=fs)
plt.savefig(paths.figures / 'crw.pdf', format='pdf')

#
q_mask = (q_num_1 >= 0.09) & (q_num_1 <= 0.6)

q_num_1_cut = q_num_1[q_mask]
Iq_num_1_cut = Iq_num_1[q_mask]
# f_cut = f[mask]
q_p_1=q_num_1_cut[np.argmax(Iq_num_1_cut)]
print(Iq_num_1_cut)
# print(np.max(Iq))
# q_p_1_l=q_num_1_cut[np.argmax(Iq_num_1_cut)-1]
# q_p_1_r=q_num_1_cut[np.argmax(Iq_num_1_cut)+1]
# q_p_err=(q_p_1_r-q_p_1_l)*0.5

q_num_2_cut = q_num_2[q_mask]
Iq_num_2_cut = Iq_num_2[q_mask]

q_p_2=q_num_2_cut[np.argmax(Iq_num_2_cut)]

file_11=paths.output / 'qp1.txt'
file_21=paths.output / 'qp2.txt'
# file_12=paths.output / 'qperr1.txt'

# with open(file_11, "w") as f:
#     f.write(str(q_p_1))

# with open(file_12, "w") as f:
#     f.write(str(q_p_err))

np.savetxt(file_11, [q_p_1], fmt="%.2f")
np.savetxt(file_21, [q_p_2], fmt="%.2f")
# np.savetxt(file_12, [q_p_err], fmt="%.2f")