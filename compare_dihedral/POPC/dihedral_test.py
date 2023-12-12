import numpy as np
import matplotlib.pyplot as plt

# Double Bond

# *******************************
# ***** OpenFF calculations *****
# *******************************
# (Type 2)

phi_rad = np.linspace(-np.pi, np.pi, num=1000)
U = 13.197694603020*(1+np.cos(2*phi_rad-np.pi)) + -1.745204449790*(1+np.cos(1*phi_rad-np.pi))

# U = k*(1+np.cos(n*phi-phi0))

# peak at cos = 1

# *******************************
# ***** MacRog calculations *****
# *******************************
# (Type 3/9: Ryckaert-Bellemans function)

constants = [44.4498, -2.39164, -81.248, 8.87284, 34.797, -4.50675]
V = 0
n = 0

for C in constants:
        V += C*(np.cos(phi_rad-np.pi))**n
        n += 1

# HCCH
# C0 = 44.4498 | C1 = -2.39164 | C2 = -81.248 | C3 = 8.87284 | C4 = 34.797 | C5 = -4.50675

phi = np.linspace(-180, 180, num=1000)
plt.figure(1)
plt.plot(phi, U)
plt.plot(phi, V)
plt.xlabel("Degrees")
plt.ylabel("$V_d$ ($kJ mol^{-1}$)")
plt.title("Double Bond C-C Dihedral Potential | CCCC")
plt.plot(phi,U, color='red', label='OpenFF')
plt.plot(phi,V, color='blue', label='MacRog')
plt.legend()
plt.xlim(-180, 180)
plt.xticks([-180,-120,-60,0,60,120,180]);
plt.show()