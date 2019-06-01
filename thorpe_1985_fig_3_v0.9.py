# coded in Python 3.6
# beta version finalized 1 June 2019
# by Matt Barlow, please send comments to Mathew_Barlow@uml.edu

# Python code for the axysmmetric PV-inversion model from
# Thorpe, A., 1985: Diagnosis of balanced vortex structure using potential
# vorticity.  J. Atmos. Sci., 42, 397-406.

# The version is set for the TS1 case, which is Fig. 3 in T85

# NOTE:  This version is not able to exactly reproduce the results from T85 -
# the solution comes close but then moves on while converging.

# Results are very sensitive to all the different parameters and choices for
# calculating tropopause height - I have tried several different variants
# from what is used here but haven't been able to exactly reproduce Fig. 3 in
# T85.

# The core of the code closely reproduces Fig. 5 from T85, so the problem
# is probably associated with how tropopause height is calculated or some
# related calculation.


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

# Phi is calculated via relaxation, and tropopause height is re-calculated
# every step using the specified potential temperature of the tropopause.


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
beta = 1.0

# height of tropopause at lateral boundary
ht = 0.6

# PV in stratosphere and troposphere
qs = 3.0
qt = 1.0

# potential temperature boundary conditions
theta_c = 0.24
del_theta_t = -0.24
R0 = 0.5

# constants for dimensionalizing
H_const = 16.7*1000  # in m
theta0_const = 300  # in K
N2_const = 1E-4  # in s^-2
f_const = 5E-5  # in s^-1
g_const = 10  # in m/s^2

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

# definte pot temp at bottom, tropopause, and top
theta_bot_final = np.zeros(nr)
theta_trop_final = np.zeros(nr)
theta_top_final = np.zeros(nr)

theta_bot_final += theta_c/(1.0 + (R/R0)**2)
theta_trop_final += del_theta_t/(1.0 + (R/R0)**2) + ht
theta_top_final += ht + qs*(1 - ht)



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
theta_trop = np.zeros_like(theta_trop_final)
theta_trop += theta_trop_final[-1]

theta_trop_calc = np.zeros_like(theta_trop)

j = 0
while j < nz:
    if Z[j] <= ht:
        q0[j] = qt
    else:
        q0[j] = qs
    j = j + 1

ht_calc0 = ht

# analytic solution at lateral boundary
phi_lytic = np.zeros(nz)
b_t = theta_bot[-1]
c_t = qt/2
b_s = theta_top[-1] - qs
c_s = qs/2
z_tp = (theta_trop[-1] - theta_bot[-1])/qt
a_s = (b_t - b_s)*z_tp + (c_t - c_s)*z_tp**2
j = 0
while j < nz:
    if Z[j] < ht_calc0:
        phi_lytic[j] = b_t*Z[j] + c_t*Z[j]**2
    else:
        phi_lytic[j] = a_s + b_s*Z[j] + c_s*Z[j]**2
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

    # from theta, determine tropopause height
    j = 0
    while j < nz-3:
        if theta_w0[j] < theta_trop[-1]:
            if theta_w0[j+1] > theta_trop[-1]:
                ht_calc0 = Z_w[j] + ((Z_w[j+1] - Z_w[j])*
                              (theta_trop[-1] - theta_w0[j])
                              /(theta_w0[j+1] - theta_w0[j]))     
        j = j + 1

    # from tropopause height, recalculate q
    j = 1
    while j < nz-2:
        dist = Z[j] - ht_calc0
        if np.abs(dist) <= dZ/2:
            a = (dist + dZ/2)/dZ
            b = 1 - a
            q0[j] = a*qs + b*qt
        else:
            if Z[j] < ht_calc0:
                q0[j] = qt
            if Z[j] > ht_calc0:
                q0[j] = qs
        j = j + 1

        q0[0] = qt
        q0[-1] = qs

    phi0 = phi0 - (phi0[0] + phi0[1])/2

    print(change, ht_calc0)
    count = count + 1



# ***** 3. Calculate phi for the rest of the domain **********************

phi = np.zeros((nz, nr))
phi_new = np.zeros((nz, nr))
q = np.zeros((nz, nr))
p = np.zeros((nz, nr))

phi = phi + phi0[:, None]
q = q + q0[:, None]
phi_start = np.copy(phi)
phi_new = np.copy(phi_start)

theta_w = np.zeros((nz_w, nr))
phi_w = np.zeros((nz_w, nr))
q_w = np.zeros((nz_w, nr))

out = np.zeros_like(phi)
qsmooth = np.zeros_like(phi)

ht_calc = np.zeros(nr)
ht_calc += ht_calc0

count = 1
change = 1
perchange = 1

#while not (perchange < 1E-4 and count > 100):
while count < 10000:

    # gradually lower the theta BCs toward final values
    cnum = 100  # how many steps to final values
    if count <= cnum:
        theta_bot = (count*(theta_bot_final - theta_bot_final[-1])/cnum +
                     theta_bot_final[-1])
        theta_top = (count*(theta_top_final - theta_top_final[-1])/cnum +
                     theta_top_final[-1])
        theta_trop = (count*(theta_trop_final - theta_trop_final[-1])/cnum +
                      theta_trop_final[-1])
    else:
        theta_bot = np.copy(theta_bot_final)
        theta_top = np.copy(theta_top_final)
        theta_trop = np.copy(theta_trop_final)

    # R=0: use L'Hopital's rule (T85 eq 20)
    j = 1
    while j < nz-1:
        dRR_phi = 2*phi[j, 1] - 2*phi[j, 0]
        dZZ_phi = phi[j+1, 0] + phi_new[j-1, 0] - 2*phi[j, 0]

        phi_new[j, 0] = phi[j, 0] + (beta*(2*dRR_phi + dR*dR)*
               ((2*dRR_phi/(dR*dR) + 1)*dZZ_phi/(dZ*dZ*q[j, 0]) - 1)/
               (2*(1 +  (((2*dRR_phi/(dR*dR) + 1)**2)/q[j, 0])*
               (dR*dR/(dZ*dZ)))))
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
            
            p[j, i] = dRR_phi + (((1 + dR_phi/(i*dR*dR))**2)*
                      (dZZ_phi/q[j, i])*(dR*dR/(dZ*dZ))
                      -3*dR_phi/(2*i) - dR*dR)

            phi_new[j, i] = phi[j, i] + (beta*p[j, i]/
                            (2*(1 + ((1 + dR_phi/(i*dR*dR))**2)*
                            (dR*dR/(dZ*dZ))/q[j, i])))         
            j = j + 1

        phi_new[0, i] = phi_new[1, i] - dZ*theta_bot[i]
        phi_new[-1, i] = phi_new[-2, i] + dZ*theta_top[i]
        i = i + 1

    # R = L, lateral boundardy
#    phi_new[:, -1] = phi0       # numerically calculated
#    phi_new[:, -1] = phi_lytic  # analytically calculated
#    phi_new[:, -1] = phi_new[:, -2]  # no slope
    phi_new[:, -1] = 2*phi_new[:, -2] - phi_new[:, -3]  # fixed slope

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

    # from theta, determine tropopause height
    i = 0
    while i < nr:
        j = 0
        while j < nz-1:
            if theta_w[j, i] < theta_trop[i]:
                if theta_w[j+1, i] > theta_trop[i]:

                    ht_calc[i] = Z_w[j] + ((Z_w[j+1] - Z_w[j])*
                              (theta_trop[i] - theta_w[j, i])
                              /(theta_w[j+1, i] - theta_w[j, i]))

            j = j + 1
        i = i + 1

    # for a smoother tropopause, if desired
#    smooth = np.copy(ht_calc)
#    i = 1
#    while i < nr-1:
#        smooth[i] = (ht_calc[i-1] + 2*ht_calc[i] + ht_calc[i+1])/4
#        i = i + 1
#    ht_calc = np.copy(smooth)

    # from tropopause height, recalculate q
    i = 0
    while i < nr:
        j = 1
        while j < nz-1:
            if Z[j+1] < ht_calc[i]:
                q[j, i] = qt
            if Z[j-1] > ht_calc[i]:
                q[j, i] = qs
            dist = Z[j] - ht_calc[i]
            if np.abs(dist) < dZ/2:
                a = (dist + dZ/2)/dZ
                b = 1 - a
                q[j, i] = a*qs + b*qt
            j = j + 1
        i = i + 1
    q[0, :] = qt

    print(count,
          np.format_float_scientific(change, precision=2, exp_digits=1),
          np.format_float_scientific(perchange, precision=2, exp_digits=1),
          np.around(ht_calc[0], decimals=2))


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

q_calc_w[0, :] = 2*q_calc_w[1, :] - q_calc_w[2, :]
q_calc_w[-1, :] = 2*q_calc_w[-2, :] - q_calc_w[-3, :]


# calculate tropopause theta (for comparison to given tropopause theta)

theta_trop_calc = np.copy(theta_trop)
i = 0
while i < nr-1:
    theta_trop_calc[i] = np.interp(ht_calc[i], Z_w, theta_w[:, i])
    i = i + 1


# ***** 5. Plot **********************

plt.close(fig='all')

fig = plt.figure(figsize=(4, 8))

fig.subplots_adjust(hspace=0.6)

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
vlevs = np.array([4, 8, 12, 16, 18])
vlevs_lab = np.array([4, 8, 12, 16, 18])
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
           manual=vlevs_lab_pos, inline_spacing=10)
plt.xlabel(r'$\overline{R}$',  fontsize=8, labelpad=2)
plt.ylabel(r'$\overline{Z}$', fontsize=8, rotation=0, labelpad=5)
plt.xticks(np.array([0.0, 0.5, 1, 1.5, 2.0, 2.5, 3.0]), fontsize=6)
plt.yticks(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize=6)
plt.title('Wind Speed & Pot Temp (dashed)', fontsize=8, fontweight='normal')



plt.subplot(3, 1, 2, adjustable='box', aspect=1.5)

# figure out positions for contour labels for r
clevs_lab = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
clevs_lab_pos = np.zeros((clevs_lab.size, 2))
i = 0
while i < clevs_lab.size:
    iy = 25
    ix = np.argmin(np.abs(r_w[iy, :] - clevs_lab[i]))
    x, y = R[ix], Z_w[iy]
    clevs_lab_pos[i] = x, y
    i = i + 1

# figure out positions for contour labels for zof
vlevs_lab = np.array([1.2, 1.8, 2.4])
vlevs = np.array([0.6, 1.2, 1.8, 2.4, 3.0])

cs = plt.contour(R, Z_w, r_w, clevs_lab, colors='0.4',
                 linestyles='dashed', linewidths=1)

cw = plt.contour(R, Z_w, zof_w, vlevs, colors='black', linewidths=1)

for c in cs.collections:
    c.set_dashes([(0, (6.0, 3.0))])

rlabels = plt.clabel(cs, clevs_lab, inline=1, fontsize=6, fmt='%1.1f',
                     manual=clevs_lab_pos, inline_spacing=10)
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
clevs = np.array([-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25,
                  0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6])
clevs_lab = np.array([0.0, 0.2, 0.4])
clevs_lab_pos = np.zeros((clevs_lab.size, 2))
i = 0
while i < clevs_lab.size:
    ix = 45
    iy = np.argmin(np.abs(geo_w[:, ix] - clevs_lab[i]))
    x, y = R[ix], Z_w[iy]
    clevs_lab_pos[i] = x, y
    i = i + 1

# figure out positions for contour labels for s
vlevs = np.array([1.2, 2.4, 3.6, 4.8, 6.0, 7.2, 8.4])
vlevs_lab = np.array([1.2, 2.4, 3.6])

cs = plt.contour(R, Z_w, geo_w, clevs, colors='0.4',
                 linestyles='dashed', linewidths=1)

cw = plt.contour(R, Z_w, s_w, vlevs, colors='black', linewidths=1)

for c in cs.collections:
    c.set_dashes([(0, (6.0, 3.0))])

rlabels = plt.clabel(cs, clevs_lab, inline=1, fontsize=6, fmt='%1.1f',
                     manual=clevs_lab_pos, inline_spacing=15)
for l in rlabels:
    l.set_rotation(0)
plt.clabel(cw, vlevs_lab, inline=1, fontsize=6, fmt='%1.1f',
           inline_spacing=10)

plt.xlabel(r'$\overline{R}$',  fontsize=8, labelpad=2)
plt.ylabel(r'$\overline{Z}$', fontsize=8, rotation=0, labelpad=5)
plt.xticks(np.array([0.0, 0.5, 1, 1.5, 2.0, 2.5, 3.0]), fontsize=6)
plt.yticks(np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize=6)
plt.title('Stability & Geo Hgt (dashed)', fontsize=8,
          fontweight='normal')

plt.suptitle('Thorpe (1985): TS1 Case - Fig. 3', fontsize=10)

plt.savefig('thorpe_1985_fig_3_calc.png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots(1, 1)
plt.imshow(q, origin='lower', cmap='Reds')
ax.set_aspect(1/3)
plt.title('Gridded PV on half-levels')

fig, ax = plt.subplots(1, 1)
ax.set_aspect(1)
clevs = np.arange(-2, 2+0.1, 0.1)
cd = plt.contour(R, Z_w, q_w - q_calc_w, clevs)
plt.clabel(cd, clevs, fmt='%2.1f')
plt.title('PV: Specified - Calculated')

fig, ax = plt.subplots(1, 1)
plt.plot(q_w[:, 0], Z_w, color='k', marker='o')
plt.plot(q_calc_w[:, 0], Z_w, color='r', marker='o')
plt.title('PV: Specified (black) & Calculated (red) at R=0')

fig, ax = plt.subplots(1, 1)
plt.plot(R, theta_trop, color='k', marker='o')
plt.plot(R, theta_trop_calc, color='r', marker='o')
plt.title('Tropopause Theta: Specified (black) & Calculated (red)')

plt.show()