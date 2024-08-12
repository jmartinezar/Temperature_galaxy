from matplotlib import pyplot as plt
from include import *

x = np.arange(0, 3e04, 10)

plt.plot(x, fun_OII(x), label='OII, T={}'.format(T))
plt.plot(x, fun_SII(x), label='SII,T={}'.format(T))
plt.axhline(y=1.5, color='red', linestyle='--', label='1.5')
plt.axhline(y=0.34, color='red', linestyle='--', label='0.34')
plt.axhline(y=0.44, color='red', linestyle='--', label='0.44')
plt.xscale('log')
plt.xlabel(r'$N_{e}$', fontsize='18')
plt.grid(True, which='both')
plt.legend()
plt.savefig("1.pdf")

plt.plot(x, fun_NII(x), label='NII, N_e={}'.format(Ne))
plt.plot(x, fun_OIII(x), label='OIII, N_e={}'.format(Ne))
plt.xlabel(r'$T[K]$', fontsize='18')
plt.ylabel('$j_{DP}/j_{SD}$', fontsize='18')
plt.legend()
plt.yscale('log')
plt.grid(True, which='both')
plt.savefig("2.pdf")
