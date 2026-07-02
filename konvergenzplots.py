import matplotlib.pyplot as plt

h = [5, 10, 50, 100, 200]
#k=0
L2_u_k0 = [2.3646e-5,9.4993e-5,2.8990e-3,1.1711e-2,1.1711e-2]
maxDiff_u_k0 = [4.7580e-7,2.4780e-6,5.5333e-5,1.6113e-4,1.6113e-4]
L2_sigma_k0 = [816.5050,1616.6914,8494.0280,11638.3170,11638.3170]
maxDiff_sigma_k0 = [229.1400,695.3303,13328.7057,17502.6486,17502.6486]

#k=2
L2_u_k2 = [5.0867e-9,5.7667e-8,5.8839e-5,6.3816e-4,5.9148e-4]
maxDiff_u_k2 = [1.0000e-10,5.6000e-9,1.3112e-6,1.3302e-5,1.2619e-5]
L2_sigma_k2 = [8.4522e-1,4.3032,357.0260,1719.3667,1626.5437]
maxDiff_sigma_k2 = [1.0054e-2,1.7091e-1,200.0349,6043.4148,5970.5905]

#k=4
L2_u_k4 = [ 0.0, 5.0000e-10, 6.0459e-7, 1.6516e-5, 1.4634e-5]
maxDiff_u_k4 = [0.0, 1.0000e-10, 1.7000e-8, 4.3630e-7, 3.9490e-7]
L2_sigma_k4 = [ 6.8370e-2, 3.0309e-1, 17.7566, 93.6906,87.6881]
maxDiff_sigma_k4 = [ 8.1108e-5, 1.3670e-3, 8.7338e-1, 14.6612, 34.7003]

#L2-fehler Durchbiegung
plt.figure(figsize=(7,5))

plt.loglog(h, L2_u_k0, 'o-', label='k = 0')
plt.loglog(h, L2_u_k2, 's-', label='k = 2')
plt.loglog(h, L2_u_k4, '^-', label='k = 4')

plt.gca().invert_xaxis()

plt.xlabel("Gitterweite h [mm]")
plt.ylabel(r"$\|u_h-u_{ref}\|_{L^2(\Omega)}$")
plt.title("L2-Fehler der Durchbiegung")
plt.legend()


#L2-Fehler Momententensor

plt.figure(figsize=(7,5))

plt.loglog(h, L2_sigma_k0, 'o-', label='k = 0')
plt.loglog(h, L2_sigma_k2, 's-', label='k = 2')
plt.loglog(h, L2_sigma_k4, '^-', label='k = 4')

plt.gca().invert_xaxis()

plt.xlabel("Gitterweite h [mm]")
plt.ylabel(r"$\|\sigma_h-\sigma_{ref}\|_{L^2(\Omega)}$")
plt.title("L2-Fehler des Momententensors")
plt.legend()

plt.show()

#maxFehler Durchbiegung
plt.figure(figsize=(7,5))

plt.loglog(h, maxDiff_u_k0, 'o-', label='k = 0')
plt.loglog(h, maxDiff_u_k2, 's-', label='k = 2')
plt.loglog(h, maxDiff_u_k4, '^-', label='k = 4')

plt.gca().invert_xaxis()

plt.xlabel("Gitterweite h [mm]")
plt.ylabel(r"max. Fehler u")
plt.title("Maximaler Fehler der Durchbiegung")
plt.legend()

#maxFehler Momentensor
plt.figure(figsize=(7,5))

plt.loglog(h, maxDiff_sigma_k0, 'o-', label='k = 0')
plt.loglog(h, maxDiff_sigma_k2, 's-', label='k = 2')
plt.loglog(h, maxDiff_sigma_k4, '^-', label='k = 4')

plt.gca().invert_xaxis()

plt.xlabel("Gitterweite h [mm]")
plt.ylabel(r"max. Fehler $\sigma$")
plt.title("Maximaler des Momententensors")
plt.legend()

plt.show()