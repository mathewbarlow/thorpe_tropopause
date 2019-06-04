# coded in Python 3.6
# Version 1.0 finalized 1 June 2019
# by Matt Barlow, please send comments to Mathew_Barlow@uml.edu

# Python code for the axysmmetric PV-inversion model from
# Thorpe, A., 1985: Diagnosis of balanced vortex structure using potential
# vorticity.  J. Atmos. Sci., 42, 397-406.

# The version is set for Fig. 4 in T85.  The results only match the figure
# if the sine in eq. 24 is squared and that change is used here.

# The results appear to match very closely with those shown in T85.

# Poorly coded in basic Fortran style for ease of debugging by the author.
# May cause nausea, cramps, and uncontrollable weeping.

# phi and q are on half-levels (-dZ/2, dZ/2, 3*dZ/2, ...) and theta is on
# whole levels (0, dZ, ...).  All variables are converted to whole levels
# for plotting

# T85 seems to imply that phi and q are on whole levels but that doesn't seem
# to make much difference other than making the boundary conditions more
# complicated, so I have used the easier approach here

# variables that end in '_w' are on whole levels
# variables are, by default, non-dimensional.  The dimensional version of a
# variable ends in '_dim'

# the program has five main parts:
# 1. set up
# 2. calculating phi at lateral boundary (R=L)
# 3. calculating phi everywhere else
# 4. calculating other variables from phi
# 5. plotting

# Phi is calculated via relaxation


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ***** 1. Set up parameters and domain ***********************************

# length and height of domain (non-dimensional)
L = 3.0
H = 1.0

# resolution
dR = 0.06
dZ = 0.025

# over-relaxation parameter
beta = 1.7

# PV in troposphere
qt = 1.0

# potential temperature boundary conditions
theta_c = 0.4
R0 = 0.5

# constants for dimensionalizing
H_const = 10*1000  # in m
theta0_const = 300  # in K
N2_const = 1E-4  # in s^-2
f_const = 5E-5  # in s^-1
g_const = 10  # in m/s^2
L_const = np.sqrt(N2_const)*H_const/f_const  # in m

# 1D length and height variables
# height at half-levels, extending above and below vertical boundaries
# -dZ/2, +dZ/2, ..., H - dZ/2, H + dZ/2
Z = np.arange(-dZ/2, H + dZ, dZ)
R = np.arange(0, L + dR, dR)

# height at whole (regular levels)
# 0, dZ, ..., H - dZ, H
Z_w = np.arange(0, H + dZ, dZ)
nz_w = Z_w.size

# note: nz is number of half-levels
nz = Z.size
nr = R.size

# definte pot temp at bottom and top
theta_bot_final = np.zeros(nr)
theta_top_final = np.zeros(nr)

theta_bot_final += theta_c/(1 + (R/R0)**2)
theta_top_final += theta_bot_final + 1

# define PV
q = np.zeros((nz, nr))
i = 0
while i < nr:
    j = 0
    while j < nz:
        q[j, i] = 1 + 10*np.sin(np.pi*Z[j])**2/((1 + (R[i]/R0)**2)**2)
        j = j + 1
    i = i + 1

q_final = np.copy(q)

# ***** 2. Calculate phi at lateral boundary (R = L) **********************
# lateral boundary variables end in '0'
# assuming phi0(Z=0)=0

phi0 = np.zeros(nz)
phi0_new = np.zeros(nz)
q0 = np.zeros(nz)
theta_w0 = np.zeros(nz-1)

theta_bot = np.zeros_like(theta_bot_final)
theta_bot += theta_bot_final[-1]
theta_top = np.zeros_like(theta_top_final)
theta_top += theta_top_final[-1]

q0[:] = q[:, -1]

# analytic solution at lateral boundary
phi_lytic = np.zeros(nz)
b_t = theta_bot[-1]
c_t = qt/2
j = 0
while j < nz:
    phi_lytic[j] = b_t*Z[j] + c_t*Z[j]**2
    j = j + 1

phi_lytic = phi_lytic - (phi_lytic[0] + phi_lytic[1])/2

phi0 = np.copy(phi_lytic)
phi0_new = np.copy(phi_lytic)

# numerical solution at lateral boundary
change = 1
count = 1
while (change > 1E-4):

    j = 1
    while j < nz-1:

        dZZ_phi0 = phi0[j+1] + phi0_new[j-1] - 2*phi0[j]

        phi0_new[j] = (phi0[j] +
                       beta*((dZZ_phi0*(dR/dZ)**2 - q0[j]*dR**2) /
                        (2*q0[j] + 2*(dR/dZ)**2)))

        j = j + 1

    phi0_new[0] = phi0_new[1] - dZ*theta_bot[-1]
    phi0_new[-1] = phi0_new[-2] + dZ*theta_top[-1]

    phi0_change = phi0 - phi0_new

    change = np.amax(np.abs(phi0_new-phi0))

    phi0 = np.copy(phi0_new)

    # calculate theta from phi
    j = 1
    while j < nz-2:
        theta_w0[j] = (phi0[j+1] - phi0[j])/dZ
        j = j + 1

    theta_w0[0] = theta_bot_final[-1]
    theta_w0[-1] = theta_top_final[-1]

    phi0 = phi0 - (phi0[0] + phi0[1])/2

    print(change)
    count = count + 1



# ***** 3. Calculate phi for the rest of the domain **********************

phi = np.zeros((nz, nr))
phi_new = np.zeros((nz, nr))
p = np.zeros((nz, nr))

phi = phi + phi0[:, None]
phi_start = np.copy(phi)
phi_new = np.copy(phi_start)

theta_w = np.zeros((nz_w, nr))
phi_w = np.zeros((nz_w, nr))
q_w = np.zeros((nz_w, nr))

out = np.zeros_like(phi)
qsmooth = np.zeros_like(phi)

count = 1
change = 1
perchange = 1

#while not (perchange < 1E-4 and count > 500):
while count < 1800:

    # gradually lower the theta BCs toward final values
    cnum = 500  # how many steps to final values
    if count <= cnum:
        theta_bot = (count*(theta_bot_final - theta_bot_final[-1])/cnum +
                     theta_bot_final[-1])
        theta_top = (count*(theta_top_final - theta_top_final[-1])/cnum +
                     theta_top_final[-1])
        q = (count*(q_final - q_final[:, -1][:, None])/cnum +
                     q_final[:, -1][:, None])
    else:
        theta_bot = np.copy(theta_bot_final)
        theta_top = np.copy(theta_top_final)
        q = np.copy(q_final)

    # R=0: use L'Hopital's rule (T85 eq 20)
    j = 1
    while j < nz-1:
        dRR_phi = 2*phi[j, 1] - 2*phi[j, 0]
        dZZ_phi = phi[j+1, 0] + phi_new[j-1, 0] - 2*phi[j, 0]

        phi_new[j, 0] = phi[j, 0] + (beta * (2*dRR_phi + dR*dR) *
               ((2*dRR_phi/(dR*dR) + 1)*dZZ_phi/(dZ*dZ*q[j, 0]) - 1)/
               (2*(1 + (((2*dRR_phi/(dR*dR) + 1)**2)/q[j, 0])*
               (dR*dR / (dZ*dZ)))))
        j = j + 1

    phi_new[0, 0] = phi_new[1, 0] - dZ*theta_bot[0]
    phi_new[-1, 0] = phi_new[-2, 0] + dZ*theta_top[0]

    # R > 0 & R < L
    i = 1
    while i < nr-1:
        j = 1
        while j < nz-1:
            dRR_phi = phi[j, i+1] + phi_new[j, i-1] - 2*phi[j, i]
            dZZ_phi = phi[j+1, i] + phi_new[j-1, i] - 2*phi[j, i]
            dR_phi = phi[j, i+1] - phi_new[j, i-1]

            p[j, i] = dRR_phi + (((1 + dR_phi / (i*dR*dR))**2) *
                      (dZZ_phi / q[j, i]) * (dR*dR/(dZ*dZ))
                      - 3*dR_phi / (2*i) - dR*dR)

            phi_new[j, i] = phi[j, i] + (beta*p[j, i] /
                            (2*(1 + ((1 + dR_phi / (i*dR*dR))**2) *
                            (dR*dR/(dZ*dZ)) / q[j, i])))         
            j = j + 1

        phi_new[0, i] = phi_new[1, i] - dZ*theta_bot[i]
        phi_new[-1, i] = phi_new[-2, i] + dZ*theta_top[i]
        i = i + 1

    # R = L, lateral boundardy
    phi_new[:, -1] = phi0       # numerically calculated
#    phi_new[:, -1] = phi_lytic  # analytically calculated
#    phi_new[:, -1] = phi_new[:, -2]  # no slope
#    phi_new[:, -1] = 2*phi_new[:, -2] - phi_new[:, -3]  # fixed slope

    phi_change = phi - phi_new

    change = np.amax(np.abs(phi_new - phi))  # max change abs val
    indchange = np.argmax(np.abs(phi_new - phi))
    newval = np.ndarray.flatten(np.abs(phi_new))[indchange]
    perchange = change/newval  # max change as percent
    phi = np.copy(phi_new)
    count = count + 1

    # calculate theta
    i = 0
    while i < nr:
        j = 0
        while j < nz-1:
            theta_w[j, i] = (phi[j+1, i] - phi[j, i])/dZ
            j = j + 1
        i = i + 1

    print(count, change, perchange)


# ***** 4. Calculate other variables from phi **********************

# phi and q on whole levels
i = 0
while i < nr:
    j = 0
    while j < nz-1:
        phi_w[j, i] = (phi[j+1, i] + phi[j, i])/2
        q_w[j, i] = (q[j+1, i] + q[j, i])/2
        j = j + 1
    i = i + 1

phi_w = phi_w - phi_w[0, -1]

# theta already calculated
theta_dim_w = theta_w*theta0_const*N2_const*H_const/g_const

# calculate v
v_w = np.zeros((nz_w, nr))
i = 1
while i < nr-1:
    j = 0
    while j < nz_w:
        dR_phi = phi_w[j, i+1] - phi_w[j, i-1]
        v_w[j, i] = (dR_phi/(2*dR)) / (np.sqrt(1.0 + dR_phi/(R[i]*dR)))
        j = j + 1
    i = i + 1

v_w[:, -1] = 2*v_w[:, -2] - v_w[:, -3]
v_dim_w = v_w*np.sqrt(N2_const)*H_const

# calculate geopotential

geo_w = phi_w - v_w*v_w/2.0

# calculate r
r_w = np.zeros_like(phi_w)

i = 1
while i < nr-1:
    j = 0
    while j < nz_w:
        dR_phi = phi_w[j, i+1] - phi_w[j, i-1]
        r_w[j, i] = R[i] / (np.sqrt(1.0 + dR_phi/(R[i]*dR)))
        j = j + 1
    i = i + 1

r_w[:, 0] = 0.0
r_w[:, -1] = R[-1]

# calculate zeta/f (zof)

zof_w = np.zeros_like(phi_w)

i = 1
while i < nr-1:
    j = 0
    while j < nz_w:
        dR_phi = phi_w[j, i+1] - phi_w[j, i-1]
        dRR_phi = phi_w[j, i+1] + phi_w[j, i-1] - 2*phi_w[j, i]
        zof_w[j, i] = ((1 + dR_phi/(R[i]*dR))**2 / (3*dR_phi/(2*R[i]*dR)
                      + 1 - dRR_phi/(dR*dR)))
        j = j + 1
    i = i + 1

i = 0
j = 0
while j < nz_w:
    dRR_phi = 2*phi_w[j, i+1] - 2*phi_w[j, i]
    zof_w[j, i] = 2*dRR_phi/(dR*dR) + 1
    j = j + 1

zof_w[:, -1] = 1.0

# calculate inertial stability

s_w = np.zeros_like(phi_w)

i = 1
while i < nr-1:
    j = 0
    while j < nz_w:
        dR_phi = phi_w[j, i+1] - phi_w[j, i-1]
        s_w[j, i] = (1 + dR_phi/(R[i]*dR))**2
        j = j + 1
    i = i + 1

i = 0
j = 0
while j < nz_w:
    s_w[j, i] = zof_w[j, i]**2
    j = j + 1


# calcuate PV (for comparison to specified PV)

q_calc_w = np.zeros_like(q_w)

i = 0
while i < nr:
    j = 1
    while j < nz_w-1:
        dZZ_phi = phi_w[j+1, i] + phi_w[j-1, i] - 2*phi_w[j, i]
        q_calc_w[j, i] = zof_w[j, i]*dZZ_phi/(dZ*dZ)
        j = j + 1
    i = i + 1

q_calc_w[0, :] = q_w[1, :]
q_calc_w[-1, :] = q_w[-2, :]


# ***** 5. Plot **********************

plt.close(fig='all')

fig = plt.figure(figsize=(4, 8))

fig.subplots_adjust(hspace=0.564)

plt.subplot(3, 1, 1, adjustable='box', aspect=1.5)


clevs = np.arange(0, 90+3, 3)
vlevs = np.arange(-12, 12+4, 4)

# figure out positions for contour labels for theta
clevs_lab = np.array([6, 12, 18, 24])
clevs_lab_pos = np.zeros((clevs_lab.size, 2))
i = 0
while i < clevs_lab.size:
    ix = 45
    iy = np.argmin(np.abs(theta_dim_w[:, ix] - clevs_lab[i]))
    x, y = R[ix], Z_w[iy]
    clevs_lab_pos[i] = x, y
    i = i + 1

# figure out positions for contour labels for v
vlevs = np.array([0, 4, 8, 12, 16, 20, 24, 28, 32])
vlevs_lab = np.array([12, 16, 20, 24, 28, 32])

vlevs_lab_pos = np.zeros((vlevs_lab.size, 2))
i = 0
while i < vlevs_lab.size:
    ix = 10
    iy = np.argmin(np.abs(v_dim_w[:, ix] - vlevs_lab[i]))
    x, y = R[ix], Z_w[iy]
    vlevs_lab_pos[i] = x, y
    i = i + 1

cs = plt.contour(R, Z_w, theta_dim_w, clevs, colors='0.4',
                 linestyles='dashed', linewidths=1)

rcParams['contour.negative_linestyle'] = 'dotted'

cw = plt.contour(R, Z_w, v_dim_w, vlevs, colors='black', linewidths=1)

for c in cs.collections:
    c.set_dashes([(0, (6.0, 3.0))])

plt.clabel(cs, clevs_lab, inline=1, fontsize=6, fmt='%2.0f',
           manual=clevs_lab_pos, inline_spacing=15)
plt.clabel(cw, vlevs_lab, inline=1, fontsize=6, fmt='%2.0f',
           inline_spacing=10)
plt.xlabel(r'$\overline{R}$',  fontsize=8, labelpad=2)
plt.ylabel(r'$\overline{Z}$', fontsize=8, rotation=0, labelpad=5)
plt.xticks(np.array([0.0, 0.5, 1, 1.5, 2.0, 2.5, 3.0]), fontsize=6)
plt.yticks(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize=6)
plt.title('Wind Speed & Pot Temp (dashed)', fontsize=8, fontweight='normal')



plt.subplot(3, 1, 2, adjustable='box', aspect=1.5)

# figure out positions for contour labels for r
clevs = np.array([0, 1.2, 2.4, 3.6, 4.8, 6, 7.2, 8.4, 9.6])
clevs_lab = np.array([1.2, 2.4])
i = 0
while i < clevs_lab.size:
    iy = 25
    ix = np.argmin(np.abs(r_w[iy, :] - clevs_lab[i]))
    x, y = R[ix], Z_w[iy]
    clevs_lab_pos[i] = x, y
    i = i + 1

# figure out positions for contour labels for zof
vlevs = np.array([1.2, 2.4, 3.6, 4.8, 6, 7.2])
vlevs_lab = np.array([1.2, 2.4])

cs = plt.contour(R, Z_w, q_w, clevs, colors='black', linestyles='solid',
                 linewidths=1)

cw = plt.contour(R, Z_w, zof_w, vlevs, colors='0.4', linestyles='dashed',
                 linewidths=1)


rlabels = plt.clabel(cs, clevs_lab, inline=1, fontsize=6, fmt='%1.1f',
                     inline_spacing=10)
for l in rlabels:
    l.set_rotation(0)
plt.clabel(cw, vlevs_lab, inline=1, fontsize=6, fmt='%1.1f',
           inline_spacing=10)

plt.xlabel(r'$\overline{R}$',  fontsize=8, labelpad=2)
plt.ylabel(r'$\overline{Z}$', fontsize=8, rotation=0, labelpad=5)
plt.xticks(np.array([0.0, 0.5, 1, 1.5, 2.0, 2.5, 3.0]), fontsize=6)
plt.yticks(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize=6)
plt.title('(Rel Vort)/f & Radius (dashed)', fontsize=8, fontweight='normal')



plt.subplot(3, 1, 3, adjustable='box', aspect=1.5)

# figure out positions for contour labels for geo
clevs = np.array([-1.2, -1.1,
                  -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                  0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
clevs_lab = np.array([-0.2, 0.0, 0.2, 0.4])
clevs_lab_pos = np.zeros((clevs_lab.size, 2))
i = 0
while i < clevs_lab.size:
    ix = 45
    iy = np.argmin(np.abs(geo_w[:, ix] - clevs_lab[i]))
    x, y = R[ix], Z_w[iy]
    clevs_lab_pos[i] = x, y
    i = i + 1

# figure out positions for contour labels for s
vlevs = np.array([12, 24, 36, 48, 60, 72])
vlevs_lab = np.array([12, 24, 36])

cs = plt.contour(R, Z_w, geo_w, clevs, colors='0.4',
                 linestyles='dashed', linewidths=1)

cw = plt.contour(R, Z_w, s_w, vlevs, colors='black', linewidths=1)

for c in cs.collections:
    c.set_dashes([(0, (6.0, 3.0))])

rlabels = plt.clabel(cs, clevs_lab, inline=1, fontsize=6, fmt='%1.1f',
                     inline_spacing=15)
for l in rlabels:
    l.set_rotation(0)
plt.clabel(cw, vlevs_lab, inline=1, fontsize=6, fmt='%1.0f',
           inline_spacing=10)

plt.xlabel(r'$\overline{R}$',  fontsize=8, labelpad=2)
plt.ylabel(r'$\overline{Z}$', fontsize=8, rotation=0, labelpad=5)
plt.xticks(np.array([0.0, 0.5, 1, 1.5, 2.0, 2.5, 3.0]), fontsize=6)
plt.yticks(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize=6)
plt.title('Stability & Geo Hgt (dashed)', fontsize=8,
          fontweight='normal')

plt.suptitle('Thorpe (1985): Fig. 5', fontsize=10)

plt.savefig('thorpe_1985_fig_5_calc.png', bbox_inches='tight', dpi=300)


plt.figure(2)

clevs = np.array([0, 1.2, 2.4, 3.6, 4.8, 6, 7.2, 8.4, 9.6])
clevs_lab = np.array([1.2, 2.4])


cs = plt.contour(R, Z_w, q_w, clevs, colors='black', linestyles='solid',
                 linewidths=1)

cs2 = plt.contour(R, Z_w, q_calc_w, clevs, colors='red', linestyles='solid',
                  linewidths=1)

rlabels = plt.clabel(cs, clevs_lab, inline=1, fontsize=6, fmt='%1.1f',
                     inline_spacing=10)


plt.xlabel(r'$\overline{R}$',  fontsize=10, labelpad=2)
plt.ylabel(r'$\overline{Z}$', fontsize=10, rotation=0, labelpad=5)
plt.xticks(np.array([0.0, 0.5, 1, 1.5, 2.0, 2.5, 3.0]), fontsize=10)
plt.yticks(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize=10)
plt.title('Specified q (black) & Calculated q (red)', fontsize=14,
          fontweight='normal')

plt.show()
