import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

from astroquery.sdss import SDSS
import pandas as pd
from include import *

query = open('query.sql', 'r').read()

data = SDSS.query_sql(query)

flux = []

for ii in range(50):
  a = (data['oii_3729_flux'][ii + 100])/(data['oii_3726_flux'][ii + 100])
  a_err = div_err(data['oii_3729_flux'][ii + 100], data['oii_3729_flux_err'][ii + 100], data['oii_3726_flux'][ii + 100], data['oii_3726_flux_err'][ii + 100])
  b = ((data['oiii_4959_flux'][ii + 100]) + (data['oiii_5007_flux'][ii + 100]))/(data['oiii_4363_flux'][ii + 100])
  b_err = div_err(((data['oiii_4959_flux'][ii + 100]) + (data['oiii_5007_flux'][ii + 100])), ((data['oiii_4959_flux_err'][ii + 100]) + (data['oiii_5007_flux_err'][ii + 100])), data['oiii_4363_flux'][ii + 100], data['oiii_4363_flux_err'][ii + 100])

  if(a >= 0.34 and a <= 1.5 and b > 100):
    aux = [ii+100, data['specObjid'][ii+100], data['class'][ii+100], a, a_err, b, b_err]
    flux.append(aux)

flux = np.asarray(flux)
flux

sols = []

# Valores tentativos
N_e = 100
T = 10000

# Se hallan los valores para cada objeto en nuestros datos
for ii in range(flux.shape[0] - 1):
  a = float(flux[ii][3])
  a_err = float(flux[ii][4])
  b = float(flux[ii][5])
  b_err = float(flux[ii][6])

  Ne = 100
  C = 8.6e-8
  solutions = sol(N_e, T, a, b)

  Ne = 100
  C = 8.6e-8
  solaux = sol(N_e, T, a+a_err, b+b_err)

  aux = [flux[ii][0], flux[ii][1], flux[ii][2], solutions[0][0], abs(solaux[0][0] - solutions[0][0]), solutions[1][0], abs(solaux[1][0] - solutions[1][0])]
  sols.append(aux)

sols = np.asarray(sols)
sols

index_values = ['first', 'second', 'third', 'fourth', 'fifth']

column_values = ['Number', 'ObjID', 'Class', 'N_e', 'N_e_err', 'T[K]', 'T_err[K]']

df = pd.DataFrame(data = sols, index = index_values, columns = column_values)

print(df)
