import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

########################### Variables #########################################

T = 1e4 # Temperature

########################### SII #########################################
A_6731 = 8.8e-04
A_6716 = 2.6e-04

g_S = 4 # Electronic density of S band
g_D_32 = 4 # Electronic density of D_32 band
g_D_52 = 6 # Electronic density of D_52 band

Omega_SD_32 = 2.76
Omega_SD_52 = 4.14
Omega_DD = 7.47

########################### OII #########################################
A_3729 = 3.6e-05
A_3726 = 1.6e-04

Omega_SD_32_ = 0.536
Omega_SD_52_ = 0.804
Omega_DD_ = 1.17

C = 8.6e-06/np.sqrt(T)

########################### OIII #########################################
Ne = 1e2

A_SD = 1.6
A_SP = 2.3e-01
A_DP_2 = 6.8e-03
A_DP_1 = 2.0e-02
A_DP = A_DP_2 + A_DP_1

g_S = 1
g_D = 5

Omega_SD = 0.58
Omega_SP = 0.29
Omega_DP_1 = 2.29
Omega_DP_2 = 0.29
Omega_DP = Omega_DP_2 + Omega_DP_1

l_SD_OIII = 4.363e-07
l_SP_OIII = 2.321e-7
l_DP2_OIII = 5.007e-7
l_DP1_OIII = 4.959e-7

c = 3e08 #[m/s] light velocity
h = 6.626e-34 #[J*s] Planck constante
kB = 1.38e-23 #[JK^{-1}] Boltzmann constant

########################### NII #########################################
Ne =1e2

A_SD_ = 1
A_SP_ = 3.3e-2
A_DP_1_ = 9.8e-04
A_DP_2_ = 3.0e-03
A_DP_ = A_DP_2_ +A_DP_1_

g_S = 1
g_D = 5

Omega_SD_ = 0.83
Omega_SP_ = 0.29
Omega_DP_1_ = 2.64
Omega_DP_2_ = 0.29
Omega_DP_ = Omega_DP_2_ + Omega_DP_1_

l_SD_NII = 5.755e-07
l_SP_NII = 3.063e-7
l_DP2_NII = 6.583e-7
l_DP1_NII = 6.548e-7

########################### Functions #########################################

def fun_SII(N):
    A = 1 + Omega_DD/Omega_SD_32 + Omega_DD/Omega_SD_52
    B1 = (A_6731*g_D_32)/(Omega_SD_32)
    B2 = (A_6716*g_D_52)/(Omega_SD_52)
    D = (A_6716*g_D_52)/(A_6731*g_D_32)
    return D*((C*A*N+B1)/(C*A*N+B2))

def fun_OII(x):
    A = 1 + Omega_DD_/Omega_SD_32_ + Omega_DD_/Omega_SD_52_
    B1 = (A_3726*g_D_32)/(Omega_SD_32_)
    B2 = (A_3729*g_D_52)/(Omega_SD_52_)
    D = (A_3729*g_D_52)/(A_3726*g_D_32)
    return D*((C*A*x+B1)/(C*A*x+B2))

def E(l):
    return h*c/l

def exp_E_SD(T):
    return np.exp(-E(l_SD_OIII)/(kB*T))

def fun_OIII(x):
    gamma = (Omega_SD + Omega_SP)/Omega_SP
    A1 = (1 + Omega_SD*exp_E_SD(x)/(gamma*Omega_DP))*Ne*8.6e-6
    A2 = (A_SD + A_SP)*g_S/(Omega_SP*gamma) + A_SD*exp_E_SD(x)*g_S/(gamma*Omega_DP)
    B1 = (1 + Omega_SD*exp_E_SD(x)/(Omega_DP*gamma))*Ne*8.6e-06
    B2 = A_DP*g_D/(Omega_DP*gamma)

    D = g_S*exp_E_SD(x)/g_D
    D1 = (A_DP_1/l_DP1_OIII+A_DP_2/l_DP2_OIII)/(A_SD/l_SD_OIII)

    return D1*(A1/np.sqrt(x) + A2)/(D*(B1/np.sqrt(x) + B2))

def exp_E_SD2(T):
    return np.exp(-E(l_SD_NII)/(kB*T))

def fun_NII(x):
    gamma = (Omega_SD_ + Omega_SP_)/Omega_SP_
    A1 = (1 + Omega_SD_*exp_E_SD2(x)/(gamma*Omega_DP_))*Ne*8.6e-6
    A2 = (A_SD_ + A_SP_)*g_S/(Omega_SP_*gamma) + A_SD_*exp_E_SD2(x)*g_S/(gamma*Omega_DP_)
    B1 = (1 + Omega_SD_*exp_E_SD2(x)/(Omega_DP_*gamma))*Ne*8.6e-06
    B2 = A_DP_*g_D/(Omega_DP_*gamma)

    D = g_S*exp_E_SD2(x)/g_D
    D1 = (A_DP_1_/l_DP1_NII+A_DP_2_/l_DP2_NII)/(A_SD_/l_SD_NII)

    return D1*(A1/np.sqrt(x) + A2)/(D*(B1/np.sqrt(x) + B2))

def div_err(d1, d1_err, d2, d2_err):
  # Función que implementa la propagación de errores para una división
  # Argumentos recibidos:
  # d1: Numerador de la división
  # d1_err: Error del numerador
  # d2: Denominador
  # d2_err: Error del denominador
  # Retorna: Error de la fracción
  return d1*(d1_err/d1 + d2_err/d2)/d2

def sol(N_e, Temp, f1, f2):
  # Se implementa una función que calcule la densidad electrónica y la temperatura
  # Inputs:
  # N_e: Densidad electrónica tentativa
  # Temp: Temperatura tentativa
  # f1: Razón de flujo para la ecuación de OII
  # f2: Razón de flujo para la ecuación de OIII
  # Output: densidad electrónica y temperatura

  #Se guardan los valores para poder compararlos al hacer el cambio
  auxN = N_e
  auxT = Temp

  #Se definen las funciones para hallar las raíces
  func1 = lambda x: fun_OII(x) - f1
  func2 = lambda x: fun_OIII(x) - f2

  #Se cambia el N_e
  N_e = fsolve(func1, N_e)
  Ne = N_e

  #Se cambia el T
  Temp = fsolve(func2, Temp)
  C = 8.6e-06/np.sqrt(Temp)

  # El while inicia si la diferencia entre los valores tentativos y los obtenidos no es mayor que 10^(-6)
  # Cada iteración repite el proceso
  while((abs(auxN - N_e)>1e-6) or (abs(auxT - Temp)>1e-6)):
    auxN = N_e
    auxT = Temp

    N_e = fsolve(func1, N_e)
    Ne = N_e
    Temp = fsolve(func2, Temp)
    C = 8.6e-06/np.sqrt(Temp)

  return N_e, Temp
