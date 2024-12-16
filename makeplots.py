import numpy as np
import matplotlib.pyplot as plt

plot_auto_corr, auto_correlation = np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/autocorrelation.csv', delimiter=',')

stand_err, energy_values = np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/error_auto.csv', delimiter=',')

beta_energy, beta_energy_std, beta_values = np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/beta_data.csv', delimiter=',')

S_energys_big, pro_ens_big, pot_ens_big, S_array_big=np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/proton_sep_data_big_10.csv', delimiter=',')

S_energys_small, pro_ens_small, pot_ens_small, S_array_small=np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/proton_sep_data_small.csv', delimiter=',')

S_energys_small_2, pro_ens_small_2, pot_ens_small_2, S_array_small_2, std_small_2=np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/proton_sep_data_small_2.csv', delimiter=',')

binding_energy_array=np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/im_time_data_mult.csv', delimiter=',')

Energys_im_time=np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/im_time_data_1_run.csv', delimiter=',')

helium_beta_energy, helium_beta_energy_std, helium_beta_values =np.genfromtxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/beta_data_helium.csv', delimiter=',')

'plot autocorrelation'
n_combin=len(stand_err)

divide=np.array(())
for i in range(n_combin):
    divide=np.append(divide,2**(i))
    
plt.semilogx(base=2)

plt.plot(divide,stand_err,c='#FF7256')
plt.scatter(divide,stand_err,c='#FF7256')

plt.plot(divide[:-1],stand_err[1:],c='#FF7256')
plt.scatter(divide[:-1],stand_err[1:],c='#FF7256')

plt.xlabel('Bin size')
plt.ylabel('Standard deviation')
plt.title('\u03C3 for different bin sizes')

plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/error_auto_corr.png',dpi=300)

lag=np.linspace(0, 128,129)

plt.plot(lag,plot_auto_corr,c='royalblue')
for i in range(6):
    plt.scatter(lag[2**i-1],plot_auto_corr[2**i-1],label=f'bin of {2**i}')
plt.legend()

plt.xlabel('Lag')
plt.ylabel('C(k)/C(0)')
plt.title('Autocorrelation')

plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/autocorrelation_plot.png',dpi=300)

#%%
'plot beta'

poly=np.polyfit(beta_values,beta_energy, deg=3)

mult=1.2
x=np.linspace(1, 2, 1000)
# y=poly[0]*x**2+poly[1]*x+poly[2]
y=poly[0]*x**3+poly[1]*x**2+poly[2]*x+poly[3]

y_min=x[np.argmin(y)]

av_error=np.average(beta_energy_std)*mult

plt.figure(figsize=(10,6))
plt.plot(x,y,c='orangered',label='Min $\u03B2$ = 1.349 $\pm$ 0.11 $\AA^{-1}$')
plt.scatter(beta_values,beta_energy,c='deepskyblue',s=10)
plt.errorbar(beta_values, beta_energy,
              yerr=beta_energy_std*0.3,ls='None',capsize=3,ecolor='deepskyblue',alpha=0.5)

plt.xlabel(r'beta ($\AA^{-1}$)')
plt.ylabel('Energy (eV)')
plt.title('Energy for different beta')
plt.legend()

plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/beta_hydrogen_2.png',dpi=300)
plt.show()


#%%
'plot proton separation'
angstrom = 1e-10
eV=1.60217663e-19
ep_0 = 8.854187817*10**(-12)

x_values=np.linspace(0.1,10,10000)
y_values_1= (eV**2)/(4*np.pi*ep_0*x_values*angstrom)/eV

fig, (axs1, axs2) = plt.subplots(2, 1)


plt.subplot(2, 1, 1)

plt.subplots_adjust(left=0.15, bottom=None, right=0.95, top=0.92, wspace=None, hspace=None)


plt.plot(x_values, y_values_1,c='teal')
plt.title("Proton and electron energy with proton seperation")
plt.ylabel("Proton energy (eV)")
plt.tick_params(labelbottom=False)
# plt.xlabel("Proton separation ($\AA^{{-1}}$)")
plt.scatter(S_array_big/angstrom,pro_ens_big,c='darkcyan')

x = S_array_big/angstrom
y = S_energys_big
a,b=np.polyfit(np.log(x[:20]), y[:20], 1)

y_values_2=a*np.log(x_values)+b

plt.subplot(2, 1, 2)

plt.plot(x_values,y_values_2,c='teal')
# plt.title("Electron energy with proton seperation")

plt.ylabel("Electron energy (eV)")
plt.xlabel("Proton separation ($\AA$)")
plt.scatter(S_array_big/angstrom,S_energys_big,c='darkcyan')
#
# plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/pro_sep_big_5.png',dpi=300)

plt.show()

y_values_comb=y_values_1+y_values_2

y_min=x_values[np.argmin(y_values_comb)]

plt.plot(x_values, y_values_comb,c='teal')
plt.title("Variation of potential energy with proton seperation")
plt.ylabel("Energy (eV)")
plt.xlabel("Proton separation ($\AA$)")
# S_energys = S_energys/eV          #energy is in eV from function
plt.scatter(S_array_big/angstrom,pot_ens_big,c='darkcyan')

# plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/pro_sep_big_2.png',dpi=300)

#%%

'plot proton separation'
angstrom = 1e-10
eV=1.60217663e-19
ep_0 = 8.854187817*10**(-12)

x_values=np.linspace(0.5,1.5,10000)
y_values_1= (eV**2)/(4*np.pi*ep_0*x_values*angstrom)/eV

plt.plot(x_values, y_values_1,c='teal')
plt.title("Proton energy with proton seperation")
plt.ylabel("Proton energy (eV)")
plt.xlabel("Proton separation ($\AA^{{-1}}$)")
plt.scatter(S_array_small/angstrom,pro_ens_small,c='darkcyan')
plt.show()

x = S_array_small_2/angstrom
y = S_energys_small_2
a,b=np.polyfit(np.log(x), y, 1)
y_values_2=a*np.log(x_values)+b

plt.plot(x_values,y_values_2,c='teal')
plt.title("Electron energy with proton seperation")
plt.ylabel("Electron energy (eV)")
plt.xlabel("Proton separation ($\AA^{{-1}}$)")
plt.scatter(S_array_small/angstrom,S_energys_small,c='darkcyan')
plt.show()

av_error=np.average(std_small_2)

y_values_comb=y_values_1+y_values_2
y_min=x_values[np.argmin(y_values_comb)]
aaa=y_values_comb[np.argmin(y_values_comb)]

plt.plot(x_values, y_values_comb,c='teal',label='Min S = 0.77 $\pm$ 0.03 $\AA$')
plt.title("Variation of potential energy with proton seperation")
plt.ylabel("Energy (eV)")
plt.xlabel("Proton separation ($\AA$)")
plt.scatter(S_array_small_2/angstrom,pot_ens_small_2,c='darkcyan',s=10)
plt.errorbar(S_array_small_2/angstrom, pot_ens_small_2,
              yerr=std_small_2,ls='None',capsize=3,ecolor='darkcyan',alpha=0.5)
plt.legend()
plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/pro_sep_small_2.png',dpi=300)

plt.show()

#%%
'plot im time'


fig, ax = plt.subplots()

t_array=np.linspace(0,21,21)

plt.title("Time evolution to groundstate")
plt.ylabel("Energy (eV)")
plt.xlabel("Imaginary time step")
plt.tick_params(axis='x', labelbottom=False)
plt.rcParams["figure.figsize"] = [6.00, 3.50]

ax.scatter(t_array,Energys_im_time,c='cornflowerblue')

plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/im_time.png',dpi=300)

plt.show()





#%%
'plot helium beta'

helium_beta_energy, helium_beta_energy_std, helium_beta_values
deg=3
poly=np.polyfit(helium_beta_values,helium_beta_energy, deg=deg)

mult=0.2
x=np.linspace(0, 1, 1000)

if deg==3:
    y=poly[0]*x**3+poly[1]*x**2+poly[2]*x+poly[3]
if deg==2:
    y=poly[0]*x**2+poly[1]*x+poly[2]

y_min=x[np.argmin(y)]
aa=min(y)
av_error=np.average(helium_beta_energy_std)*mult

plt.plot(x,y,c='orangered',label='Min $\u03B2$ = 0.3 $\pm$ 0.035 $\AA^{-1}$')

plt.scatter(helium_beta_values,helium_beta_energy,c='deepskyblue',s=10)
plt.errorbar(helium_beta_values, helium_beta_energy,
              yerr=helium_beta_energy_std*0.3,
              # xerr=helium_beta_energy_std*0.3,
              ls='None',capsize=3,ecolor='deepskyblue',alpha=0.5)
plt.xlabel(r'beta ($\AA^{-1}$)')
plt.ylabel('Energy (eV)')
plt.title('Helium energy for different beta')
plt.legend()

plt.savefig('C:/Users/paxoc/Documents/uni/year 4/MPhys/plots/beta_helium.png',dpi=300)
plt.show()


#%%
S=7.7e-11
av_energy=np.average(binding_energy_array)
std_energy=np.std(binding_energy_array)

min_en = 2*(-13.6*eV)/eV-(eV**2)/(4*np.pi*ep_0*S)/eV-av_energy



