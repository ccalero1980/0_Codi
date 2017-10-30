import numpy as np
import sys
# PARAMETRES: Introdueixo la frequencia i la magnitud del camp magnetic exteriorment

namein = './inconf_circular30  : file name\n'
#namein = './inconf_10  : file name\n'

freq_field = float(sys.argv[1]) #Hz
B = float(sys.argv[2]) # Tesla
chi = float(sys.argv[3]) # chi
# ----------------------------------------

long_filament = 1.e-5 # m
radi_filament = 3.e-7  #m
radi_bola = 5.e-6
densitat = 8000. # kg/m^3
massa_filament = densitat*3.1416*radi_filament**2.*long_filament # kg
viscositat = 0.001 #Pa*s
gamma_filament = 6*3.1416*viscositat*(long_filament*radi_filament)**0.5
bend_stiffness = 5.e-25 # J*m (ordre de magnitud, en el cas estudiat pel Gauger es de 4.5e-22)
mu0 = 4*3.1416*1.e-7 # T^2*m^3/J = T^2*m*s^2/kg
Bz = B*0.01 # Tesla

Ncycles = 10 # Nombre de cicles del camp magnetic a simular
moment_magnetic_filament = 3.72e-17 # J/T
susceptibilitat_bola = chi*(4*3.1416/3.)*radi_bola**3./mu0 #J/T^2

t0 = 1.e-6 # segons
tinercial = massa_filament/gamma_filament

oseen = 0
blake = 1

epsD = mu0*susceptibilitat_bola*B*moment_magnetic_filament/(4.*3.1416*bend_stiffness*long_filament**2.)
epsBferro = B*moment_magnetic_filament*long_filament/bend_stiffness
spermnum = (viscositat*freq_field*long_filament**4./bend_stiffness)**0.25

fout1= open('eps','w')
tot = 'epsD'+str(epsD)+'_epsBferro'+str(epsBferro)+'_spermnum'+str(spermnum)
fout1.write('epsD%.2f_epsBferro%.2f_spermnum%.2f\n'%(epsD, epsBferro, spermnum))
fout1.close()


print 'NON-DIMENSIONAL PARAMETERS:\n'
print 'epsD = %.5f, epsBferro = %.5f, spermnum = %.5f\n'%(epsD, epsBferro, spermnum)

fpars = open('parameters.dat','w')
fpars.write( '\n\t PARAMETRES DEL PROBLEMA:\n')
fpars.write( '\nlong_filament (m) = %le\nmassa_filament (kg) = %le\nt0 = %le\n'%(long_filament, massa_filament, t0))
fpars.write( 'radi_filament (m)= %le \nradi_bola (m)= %le \n'%(radi_filament, radi_bola))
fpars.write( 'viscositat (en kg/(m*s)) = %le\ngamma_filament (en kg/s) = %le\n'%(viscositat, gamma_filament))
fpars.write( 'bend_stiffness (en m^3*kg/s^2) = %le\n'%bend_stiffness)
fpars.write( 'Bfield (T) = %le\nfreq_field (1/s) = %le\n'%(B, freq_field))
fpars.write( 'moment_magnetic_filament (en kg*m^2/(Tesla s^2) = %le\nsusceptibilitat_bola (en J/T^2) = %le\n'%(moment_magnetic_filament, susceptibilitat_bola))
fpars.write( 'tinercial = %le (segons)\n'%tinercial)


print '\nlong_filament (m) = %le, massa_filament (kg) = %le, t0 = %le'%(long_filament, massa_filament, t0)
print 'radi_filament (m)= %le, radi_bola (m)= %le'%(radi_filament, radi_bola)
print 'viscositat (en kg/(m*s)) = %le, gamma_filament (en kg/s) = %le'%(viscositat, gamma_filament)
print 'bend_stiffness (en m^3*kg/s^2) = %le'%bend_stiffness
print 'Bfield (T) = %le, freq_field (1/s) = %le'%(B, freq_field)
print 'moment_magnetic_filament (en kg*m^2/(Tesla s^2) = %le, susceptibilitat_bola (en J/T^2) = %le'%(moment_magnetic_filament, susceptibilitat_bola)
print 'tinercial = %le (segons)'%tinercial

# CANVI UNITATS:
# Unitats del problema: long_filament, massa_filament, t0 = 10^-6s (de manera que viscositat = 1 en unitats del problema).

radi_filament = radi_filament/long_filament  # en long_filament
viscositat = viscositat*(long_filament/massa_filament)*t0 # en massa_filament/(long_filament*t0)
#gamma_filament = 6*3.1416*viscositat*radi_filament**0.5 #en massa_filament/t0
gamma_filament = gamma_filament*t0/massa_filament
bend_stiffness = (bend_stiffness/(long_filament**3.*massa_filament))*t0**2 #  en long_filament**3 * massa_filament/t0^2

B = B # Tesla
freq_field = freq_field*t0 #1/t0
moment_magnetic_filament = (moment_magnetic_filament/(massa_filament*long_filament**2.))*t0**2. # en massa_filament*long_filament**2/(t0^2*Tesla)
#susceptibilitat_bola = susceptibilitat_bola* 4*3.1416*1.e-7*(1./long_filament**3.) # mu0*chi, en long_filament^3
susceptibilitat_bola = susceptibilitat_bola*(1./long_filament**2.)*(1./massa_filament)*t0**2. # chi, en long_filament^2*massa_filament/(t0^2 * Tesla^2)
mu0 = mu0*(massa_filament/(long_filament*t0**2)) # mu0 en T^2*long_filament*t0^2/massa_filament

massa_filament = massa_filament/massa_filament
long_filament = long_filament/long_filament

fpars.write( '\n\nCANVI UNITATS:\n')
fpars.write( 'long_filament (lf) = %f\nradi_filament (en lf)= %f\nmassa_filament (mf) = %f\n'%(long_filament, radi_filament, massa_filament))
fpars.write( 'radi_bola (en lf)= %f\n\n'%(radi_bola))
fpars.write( 'viscositat (en mf/(lf*t0)) = %f\ngamma_filament (en mf/t0) = %f\n\n'%(viscositat, gamma_filament))
fpars.write( 'bend_stiffness (en lf^3*mf/t0^2) = %f\n'%bend_stiffness)
fpars.write( 'Bfield (T) = %f\nfreq_field (1/t0) = %f\n'%(B, freq_field))
fpars.write( 'moment_magnetic_filament (en mf*lf^2/(Tesla t0^2) = %f\nsusceptibilitat_bola (en long_filament^3) = %f\n'%(moment_magnetic_filament, susceptibilitat_bola))

print '\nCANVI UNITATS:'
print 'long_filament (lf) = %f, radi_filament (en lf)= %f, massa_filament (mf) = %f'%(long_filament, radi_filament, massa_filament)
print 'radi_bola (en lf)= %f'%(radi_bola)
print 'viscositat (en mf/(lf*t0)) = %f, gamma_filament (en mf/t0) = %f'%(viscositat, gamma_filament)
print 'bend_stiffness (en lf^3*mf/t0^2) = %f'%bend_stiffness
print 'Bfield (T) = %f, freq_field (1/t0) = %f'%(B, freq_field)
print 'moment_magnetic_filament (en mf*lf^2/(Tesla t0^2) = %f, susceptibilitat_bola (en long_filament^3) = %f'%(moment_magnetic_filament, susceptibilitat_bola)
# QUANTITATS BASIQUES

tinercial = (massa_filament/gamma_filament) # en t0

print 'tinercial = %.8f (t0)'%tinercial
if tinercial > 0.1*(1./freq_field): print 'ERROR tinercial'


dt = tinercial/1000.
Tsim = Ncycles*(2.*np.pi/(freq_field))
Ntimesteps = int(Tsim/dt)

fpars.write( '\n\t SIMULACIO:\n')
fpars.write( 'timestep = %le (t0)\nSimulation_Time = %le (t0)\nNtimesteps = %i\n'%(dt, Tsim, Ntimesteps))
if oseen==1:fpars.write('OSEEN\n')
elif blake==1:fpars.write('BLAKE\n')
print 'timestep = %le (t0), Simulation_Time = %le (t0), Ntimesteps = %i'%(dt, Tsim, Ntimesteps)
#
if freq_field <100*t0:
  steps_in_cycle = 500000
elif freq_field >= 100*t0 and freq_field<=1000*t0:
  steps_in_cycle = 50000
elif freq_field >1000*t0:
  steps_in_cycle = 5000
number_cycles = int(float(Ntimesteps)/float(steps_in_cycle))


fout = open('init', 'w')
fout.write(namein)
fout.write('./       : output dir\n')
fout.write('%i	     : number of cycles\n'%number_cycles)
#fout.write('%i	     : number of cycles\n'%1)
fout.write('%f       : timestep\n'%dt)
fout.write('%le       : bending modulus\n'%bend_stiffness)
fout.write('%f       : mass filament\n'%massa_filament)
fout.write('%f       : viscosity\n'%viscositat)
fout.write('%f       : gamma filament\n'%gamma_filament)
fout.write('0.3           : hydr radius (in bond lengths)\n')
fout.write('1.0e-8        : constraint tolerance\n')
fout.write('%i	     : steps per cycle\n'%steps_in_cycle)
#fout.write('%i	     : steps per cycle\n'%1)
fout.write('%f       : Fx\n'%0.0)
fout.write('%f       : Fz\n'%0.0)
fout.write('%f       : Bx\n'%B)
fout.write('%f       : Bz\n'%Bz)
fout.write('%f       : freq_field \n'%freq_field)
fout.write('%f       : susceptibilitat_bola\n'%susceptibilitat_bola)
fout.write('%f       : moment_magnetic_filament\n'%moment_magnetic_filament)
fout.write('%f       : mu0 permeability vacuum \n'%mu0)
fout.write('%i       : Oseen\n'%oseen)
fout.write('%i       : Blake\n'%blake)




epsD = mu0*susceptibilitat_bola*B*moment_magnetic_filament/(4.*3.1416*bend_stiffness*long_filament**2.)
epsBferro = B*moment_magnetic_filament*long_filament/bend_stiffness
spermnum = (viscositat*freq_field*long_filament**4./bend_stiffness)**0.25


print 'NON-DIMENSIONAL PARAMETERS:\n'
print 'epsD = %.5f, epsBferro = %.5f, spermnum = %.5f\n'%(epsD, epsBferro, spermnum)




fout.close()
