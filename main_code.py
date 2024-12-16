from math import exp
import numpy as np
import matplotlib.pyplot as plt
import random
import time

'relevant functions'

start_time = time.time()

h=6.62607015e-34
h_bar=h/(2*np.pi)
e_mass=9.1093837e-31 
e_charge=1.60217663e-19
eV=e_charge
ep_0 = 8.854187817*10**(-12)
a0 = 5.291772e-11  # m  Bohr radius hbar/(alpha m c)
angstrom = 1e-10  # m

def timee():
    'To calculate time'
    end_time = time.time()
    tot_time=end_time - start_time
    tot_time = "{:.2f}".format(tot_time)
    print(f"Calculation time in seconds is {tot_time}")

def func_to_get_a(S,reit=30,initial_guess_a=0.2):
    al=[]
    a_f=initial_guess_a
    for _ in range(reit): # iterate to a solution ...
        a_f=a0/(1.0+exp(-S/a_f))
        al.append(a_f)
    return(al)
 
def tot_fun(position_e_1,position_e_2,position_p_L,position_p_R,beta,a):
    r_12=position_e_2-position_e_1
    r_1L=position_p_L-position_e_1
    r_1R=position_p_R-position_e_1
    r_2L=position_p_L-position_e_2
    r_2R=position_p_R-position_e_2

    r_12_m=(r_12[0]**2+r_12[1]**2+r_12[2]**2)**0.5
    r_1L_m=(r_1L[0]**2+r_1L[1]**2+r_1L[2]**2)**0.5
    r_1R_m=(r_1R[0]**2+r_1R[1]**2+r_1R[2]**2)**0.5
    r_2L_m=(r_2L[0]**2+r_2L[1]**2+r_2L[2]**2)**0.5
    r_2R_m=(r_2R[0]**2+r_2R[1]**2+r_2R[2]**2)**0.5

    return phi(r_1L_m,r_1R_m,a)*phi(r_2L_m,r_2R_m,a)*f_func(r_12_m,beta)

def phi(r_iL,r_iR,a):
    return np.exp(-r_iL/a)+np.exp(-r_iR/a)

def f_func(r_12,beta):
    return np.exp(r_12/(2*a0*(1+beta*r_12)))   ##a0 in units of angstrom

def f_r_1st(r_12,beta):                           
    return f_func(r_12,beta)/(2*a0*(1+beta*r_12)**2) 

def trial_pos_function(pos_1,pos_2,delta):
    trial_1=np.array(())
    trial_2=np.array(())

    for i in range(3):
        trial_1=np.append(trial_1,pos_1[i]+2*delta*angstrom*(random.uniform(0,1)-0.5))
        trial_2=np.append(trial_2,pos_2[i]+2*delta*angstrom*(random.uniform(0,1)-0.5))
    return trial_1,trial_2

def get_points(n,L,R,a,beta,delta=0.45):
    positions_1=np.array(([[angstrom*(random.uniform(0,1)-0.5)],\
                            [angstrom*(random.uniform(0,1)-0.5)],\
                            [angstrom*(random.uniform(0,1)-0.5)]]))

    positions_2=np.array(([[angstrom*(random.uniform(0,1)-0.5)],\
                           [angstrom*(random.uniform(0,1)-0.5)],\
                           [angstrom*(random.uniform(0,1)-0.5)]]))

    count=0
    for _ in range(n):
        trial_pos_1,trial_pos_2=trial_pos_function(positions_1[:,-1],positions_2[:,-1], delta)
        if (tot_fun(trial_pos_1,trial_pos_2,L,R,beta,a))**2 >\
                (tot_fun(positions_1[:,-1],positions_2[:,-1],L,R,beta,a))**2:
            positions_1=np.c_[positions_1,[trial_pos_1[0],trial_pos_1[1],trial_pos_1[2]]]
            positions_2=np.c_[positions_2,[trial_pos_2[0],trial_pos_2[1],trial_pos_2[2]]]

        else:
            if random.uniform(0,1) <= (tot_fun(trial_pos_1,trial_pos_2,L,R,beta,a))**2/\
            (tot_fun(positions_1[:,-1],positions_2[:,-1],L,R,beta,a))**2:
                positions_1=np.c_[positions_1,[trial_pos_1[0],trial_pos_1[1],trial_pos_1[2]]]
                positions_2=np.c_[positions_2,[trial_pos_2[0],trial_pos_2[1],trial_pos_2[2]]]

            else:
                positions_1=np.c_[positions_1,[positions_1[0,-1],positions_1[1,-1],positions_1[2,-1]]]
                positions_2=np.c_[positions_2,[positions_2[0,-1],positions_2[1,-1],positions_2[2,-1]]]
                count+=1

    perc="{:.2f}".format(100-count/n*100)
    # print(perc)
    return positions_1,positions_2,perc

def grad_part(a,r_L,r_R):
    return 1/a**2-2/a*((1/r_L+1/r_R*np.exp((r_L-r_R)/a))/(1+np.exp((r_L-r_R)/a)))

def g_part(a,r_L,r_L_v,r_R,r_R_v,r_12_v):
    return -1/a*((np.dot(r_L_v,r_12_v))/r_L+(np.dot(r_R_v,r_12_v))/r_R*np.exp((r_L-r_R)/a))/(1+np.exp((r_L-r_R)/a))

def f_prime_by_f(r_12,beta):
    return 1/(2*a0*r_12*(1+beta*r_12)**2)   ##a0 in units of angstrom

def grad_f_part(r_12,beta):
    return 2/(2*a0*r_12*(1+beta*r_12)**3)+1/(4*a0**2*(1+beta*r_12)**4)

def g(r,a):                             #for part of kineic en calc
    return np.exp(-r/a)

def g_1st(r,a):
    return -(1/a)*g(r,a)
def energy(a,position_e_1,position_e_2,position_p_L,position_p_R,beta):
    h=6.62607015e-34
    h_bar=h/(2*np.pi)
    e_mass=9.1093837e-31
    kin_const=h_bar**2/e_mass

    e_charge=1.60217663e-19
    eps0= 8.854187817e-12
    pot_const=e_charge**2/(4*np.pi*eps0)

    r_12 = position_e_2-position_e_1
    r_1L = position_p_L-position_e_1
    r_1R = position_p_R-position_e_1
    r_2L = position_p_L-position_e_2
    r_2R = position_p_R-position_e_2

    r_12_m = (r_12[0]**2+r_12[1]**2+r_12[2]**2)**0.5
    r_1L_m = (r_1L[0]**2+r_1L[1]**2+r_1L[2]**2)**0.5
    r_1R_m = (r_1R[0]**2+r_1R[1]**2+r_1R[2]**2)**0.5
    r_2L_m = (r_2L[0]**2+r_2L[1]**2+r_2L[2]**2)**0.5
    r_2R_m = (r_2R[0]**2+r_2R[1]**2+r_2R[2]**2)**0.5

    K=-kin_const/2*(grad_part(a, r_1L_m, r_1R_m)+grad_part(a, r_2L_m, r_2R_m)\
    +2*(g_part(a,r_1L_m,r_1L,r_1R_m,r_1R,r_12)-g_part(a,r_2L_m,r_2L,r_2R_m,r_2R,r_12))*f_prime_by_f(r_12_m, beta)\
    +2*grad_f_part(r_12_m, beta))
    P=-pot_const*(1/r_1L_m+1/r_2L_m+1/r_1R_m+1/r_2R_m-1/r_12_m)
    return (P+K)/eV

def integral_cal(position_e_1, position_e_2, position_p_L, position_p_R, beta, n_combin, a, test_binsize=False, ind_resu=False):
    divide=np.array(())
    for i in range(n_combin):
        divide=np.append(divide,2**(i))
    print(divide)
    std=np.array(())
    result=np.array(())
    energy_array=np.array(())
    for j in range(len(position_e_1[0])):
        energy_step=energy(a, position_e_1[:,j], position_e_2[:,j], position_p_L, position_p_R, beta)
        energy_array=np.append(energy_array,energy_step)
    
    if test_binsize==True:
        for i in divide:
            counter=0
            temp_bin=np.zeros((3,0))
            fin_bin=np.zeros((3,0))
            for j in range(len(position_e_1[0])):
    
                temp_bin=np.append(temp_bin,energy_array[j])
                counter+=1
    
                if counter==i:
                    fin_bin=np.append(fin_bin,np.mean(temp_bin))
                    temp_bin=np.array(())
                    counter=0
            integral=1/(n/i)*sum(fin_bin)
            result=np.append(result,integral)
    
            'std calculation'
            std_sqrd=1/((n/i)-1)*(sum(fin_bin**2)/(n/i)-integral**2)
            std=np.append(std,std_sqrd**(1/2))
            print(i)
        return std,result
    else:
        divide=2**n_combin
        counter=0
        temp_bin=np.zeros((3,0))
        fin_bin=np.zeros((3,0))
        for j in range(len(position_e_1[0])):

            temp_bin=np.append(temp_bin,energy_array[j])
            counter+=1

            if counter==divide:
                fin_bin=np.append(fin_bin,np.mean(temp_bin))
                temp_bin=np.array(())
                counter=0
        if ind_resu==True:
            return fin_bin
        integral=1/(n/divide)*sum(fin_bin)
        result=np.append(result,integral)

        'std calculation'
        std_sqrd=1/((n/divide)-1)*(sum(fin_bin**2)/(n/divide)-integral**2)
        std=np.append(std,std_sqrd**(1/2))
    return std,integral

def auto_corr(func,
    position_e_1,
    position_e_2,
    position_p_L,
    position_p_R,
    beta,
    n_combin,
    a,
    lag,
    n):
    'find autocorrelation of a given function (func) in a lag (lag)'
    '1/(n-k)*sum((f_i-f_av)(f_i+k-f_av))'
    
    inte=func(position_e_1, position_e_2, position_p_L, position_p_R, beta, n_combin, a)
    av_f=inte[1]
    
    auto_corr=np.array(())
    for k in (lag):
        dif_f=func(position_e_1[:,:int(n-k)], position_e_2[:,:int(n-k)], position_p_L, position_p_R, beta, n_combin, a, ind_resu=True)-av_f   
        dif_k=func(position_e_1[:,int(k)+1:], position_e_2[:,int(k)+1:], position_p_L, position_p_R, beta, n_combin, a, ind_resu=True)-av_f

        auto_c_k=1/(n-k)*sum(dif_f*dif_k)
        auto_corr=np.append(auto_corr,auto_c_k)
    norm_auto=auto_corr/auto_corr[0]

    return norm_auto,auto_corr

'drift function'
dc = h_bar/e_mass  #drift constant (start of each drift component)

def coordinates(r_1,r_2,L_pro,R_pro):
    r_12= ((r_2[0]-r_1[0])**2 + (r_2[1]-r_1[1])**2 + (r_2[2]-r_1[2])**2)**(1/2)
    r_1L= ((L_pro[0]-r_1[0])**2+(L_pro[1]-r_1[1])**2+(L_pro[2]-r_1[2])**2)**(1/2)
    r_1R= ((R_pro[0]-r_1[0])**2+(R_pro[1]-r_1[1])**2+(R_pro[2]-r_1[2])**2)**(1/2)
    r_2L= ((L_pro[0]-r_2[0])**2+(L_pro[1]-r_2[1])**2+(L_pro[2]-r_2[2])**2)**(1/2)
    r_2R= ((R_pro[0]-r_2[0])**2+(R_pro[1]-r_2[1])**2+(R_pro[2]-r_2[2])**2)**(1/2)
    return r_12, r_1L, r_1R, r_2L, r_2R

def FA(r, r_iL, r_iR, a):
    return (1/phi(r_iL,r_iR,a))*(g_1st(r_iL,a)/r_iL + g_1st(r_iR,a)/r_iR)

def FB(r, r_iL, r_iR, a):
    return (1/phi(r_iL,r_iR,a))*(g_1st(r_iL,a)/r_iL - g_1st(r_iR,a)/r_iR)

def FE(r_12):
    # print("FE part 1:",f_r_1st(r_12))
    # print("FE part 2:",(f_r(r_12)*r_12))
    return f_r_1st(r_12,beta)/(f_func(r_12,beta)*r_12)
    
def D_12(r_12,r1,r_1L,r_1R,a,xy_1,xy_2):
    return dc*(FA(r1, r_1L, r_1R, a)*xy_1 + FE(r_12)*(xy_1-xy_2))

def D_3(r_12,r1,r_1L,r_1R,a,S,z_1,z_2):
    return dc*(FA(r1,r_1L,r_1R,a)*z_1+FE(r_12)*(z_1-z_2)+FB(r1,r_1L,r_1R,a)*(S/2))

def D_45(r_12,r2,r_2L,r_2R,a,xy_1,xy_2):
    return dc*(FA(r2, r_2L, r_2R, a)*xy_1 - FE(r_12)*(xy_1-xy_2))

def D_6(r_12,r2,r_2L,r_2R,a,S,z_1,z_2):
    # print("part 1:",dc*(FA(r2,r_2L,r_2R,a)*xy_1))
    # print("part 2:",FE(xy_1-xy_2))
    # print("part 3:",FB(r2,r_2L,r_2R,a)*(S/2))
    return dc*(FA(r2,r_2L,r_2R,a)*z_1-FE(r_12)*(z_1-z_2)+FB(r2,r_2L,r_2R,a)*(S/2))

def drift(r1,r2,L_pro,R_pro,a,S):
    r_12, r_1L, r_1R, r_2L, r_2R = coordinates(r1,r2,L_pro,R_pro)
    D = np.array([D_12(r_12,r1,r_1L,r_1R,a,r1[0],r2[0]),D_12(r_12,r1,r_1L,r_2R,a,r1[1],r2[1]),\
                     D_3(r_12,r1,r_1L,r_1R,a,S,r1[2],r2[2]),D_45(r_12,r2,r_2L,r_2R,a,r1[0],r2[0]),\
                     D_45(r_12,r2,r_2L,r_2R,a,r1[1],r2[1]),D_6(r_12,r2,r_2L,r_2R,a,S,r1[2],r2[2])])
    return D



def new_r(old_r1, old_r2, L_pro,R_pro,a,S, drift, t_step, n):
    chi = np.random.normal(0,np.sqrt(h_bar*t_step/e_mass),(6,n+1)) #*angstrom
    old_r = np.concatenate((old_r1,old_r2),axis=0)  #adding r1 + r2 column-wise
    # print('chi',chi[0])
    # print('drift',drift[0]*t_step)
    # print('old_r',old_r[0])

    return old_r + drift*t_step + chi 

#%%
'section 1'
'set parameters'

S=0.5*angstrom
a_values=func_to_get_a(S)
a_last=a_values[-1]
beta=1/angstrom
n = 2**16
n_combin=10   ##this tests all the way up to 2^10 bin size 

position_pro_L=np.array([-S/2,0,0])
position_pro_R=np.array([S/2,0,0])

'autocorrelation for beta=1A'
start_time = time.time()

position_ele_1,position_ele_2,perc=get_points(n,
                                              position_pro_L,
                                              position_pro_R,
                                              a_last,
                                              beta)
lag_points=2**7
lag=np.linspace(0, lag_points,lag_points+1)
stand_err,integral=integral_cal(position_ele_1,
                                position_ele_2,
                                position_pro_L,
                                position_pro_R,
                                beta,n_combin,
                                a_last,
                                test_binsize=True)
n_combin=1 ##for lag we dont want to see binning

plot_auto_corr, auto_correlation=auto_corr(integral_cal,
                                          position_ele_1,
                                          position_ele_2,
                                          position_pro_L,
                                          position_pro_R,
                                          beta,n_combin,
                                          a_last,lag,
                                          n)


# save plot_auto_corr, auto_correlation
np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/error_auto.csv',
            (stand_err, integral), delimiter=',')

#%%
for i in range(6):
    plt.scatter(lag[2**i-1],plot_auto_corr[2**i-1],label=f'bin of {2**i}')
plt.plot(lag,plot_auto_corr)
plt.legend()

timee()

#%%
'find lowest energy beta'
'section 2'
start_time = time.time()

'adjust before begin'
n_combin=5   ##from now on bin down to 2**5
beta_amount=30
betaa=np.linspace(1/angstrom, 2/angstrom,beta_amount)

at=100

resu=np.array(())
std_resu=np.array(())
beta_progress=0
for i in betaa:
    print('beta progress', beta_progress/beta_amount)
    beta_progress+=1
   
    ate=np.array(())
    for j in range(at):
        position_ele_1,position_ele_2,perc=get_points(n,
                                                      position_pro_L,
                                                      position_pro_R,
                                                      a_last,
                                                      i)
       
        stand_err,integral=integral_cal(position_ele_1,
                                        position_ele_2,
                                        position_pro_L,
                                        position_pro_R,
                                        i,
                                        n_combin,
                                        a_last)
        ate=np.append(ate, integral)

    av_ate=np.average(ate)
    st_dev=(sum((ate-av_ate)**2)/at)**0.5

    resu=np.append(resu,av_ate)
    std_resu=np.append(std_resu,st_dev)

betaa_ang = betaa * angstrom
best_beta=betaa_ang[np.argmin(resu)]


# save resu, st_dev and betaa_ang
np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/beta_data.csv',
           (resu, std_resu, betaa_ang), delimiter=',')

#%%
'section 3'
'find proton separation'
'range 0.1-10'
start_time = time.time()

best_beta=1.3103448275862069
n_combin=5
n = 2**16

n_proton_sep=100
S_array = np.linspace(0.1*angstrom,10*angstrom,n_proton_sep)

counter = 0

S_energys = np.array(())
pro_ens = np.array(())
pot_ens = np.array(())
for j in S_array:
    
    L_pro_S = np.array([-j/2,0,0]) 
    R_pro_S = np.array([j/2,0,0])
    a_S = func_to_get_a(j)
    a_last = a_S[-1]

    r1_S, r2_S, perc = get_points(n,L_pro_S,R_pro_S,a_last,2/angstrom)
    
    std,S_energy = integral_cal(r1_S,r2_S,L_pro_S,R_pro_S, best_beta, n_combin, a_last)
   
    pro_ens = np.append(pro_ens,((eV**2)/(4*np.pi*ep_0*j))/eV)      #proton potential in eV
    pot_ens = np.append(pot_ens,((eV**2)/(4*np.pi*ep_0*j))/eV + S_energy)   #overall potential
    
    S_energys=np.append(S_energys, S_energy)
    counter += 1 
    print(round(counter/len(S_array)*100,2),"%")

best_pro_sep=S_array[np.argmin(pot_ens)]/angstrom

##save S_energys, pro_ens, pot_ens, S_array
np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/proton_sep_data_big_10.csv',
            (S_energys, pro_ens, pot_ens, S_array), delimiter=',')
timee()

#%%
'section 4'
'find proton separation'
'range 0.5-1.5'
start_time = time.time()

n_proton_sep=100
S_array = np.linspace(0.5*angstrom,1.5*angstrom,n_proton_sep)
best_beta=1.3103448275862069
n_combin=5
n = 2**16

counter = 0

S_energys = np.array(())
std_array = np.array(())
pro_ens = np.array(())
pot_ens = np.array(())
for j in S_array:
    
    L_pro_S = np.array([-j/2,0,0]) 
    R_pro_S = np.array([j/2,0,0])
    a_S = func_to_get_a(j)
    a_last = a_S[-1]

    r1_S, r2_S, perc = get_points(n,L_pro_S,R_pro_S,a_last,2/angstrom)
    
    std,S_energy = integral_cal(r1_S,r2_S,L_pro_S,R_pro_S, best_beta, n_combin, a_last)
    
    std_array=np.append(std_array, std)
    pro_ens = np.append(pro_ens,((eV**2)/(4*np.pi*ep_0*j))/eV)      #proton potential in eV
    pot_ens = np.append(pot_ens,((eV**2)/(4*np.pi*ep_0*j))/eV + S_energy)   #overall potential
   
    S_energys=np.append(S_energys, S_energy)
    counter += 1 
    print(round(counter/len(S_array)*100,2),"%")

best_pro_sep=S_array[np.argmin(pot_ens)]/angstrom

#save S_energys, pro_ens, pot_ens, S_array
np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/proton_sep_data_small_2.csv',
            (S_energys, pro_ens, pot_ens, S_array, std_array), delimiter=',')

timee()

#%%

best_beta=1.3103448275862069
best_pro_sep=S_array[27]
best_pro_sep=7.727272727272728e-11
#%%
'section 5'
'imaginary time loop'
start_time = time.time()

n=2**16
def time_ev_f(best_beta,best_pro_sep):
    sigma=0.034*angstrom
    dtau=sigma**2*e_mass/h_bar
    time_step=dtau
   
    # time_step = 10**(-19)
   
    S=best_pro_sep
   
    a_values=func_to_get_a(S)
    a_last=a_values[-1]
    beta=best_beta/angstrom
    n_combin=5
   
    position_pro_L=np.array([-S/2,0,0])
    position_pro_R=np.array([S/2,0,0])
   
    num_steps = 30
    e_step=0.1
    error=np.array(e_step)
    for i in range(num_steps-1):
        e_step=e_step*0.6
        error =np.append(error,e_step)
       
    t_steps = np.array(time_step)
    for i in range(num_steps-1):
        time_step = time_step*(0.5)
        t_steps = np.append(t_steps,time_step)
   
    #original poisitions (G=PHI^2) 
    position_ele_1,position_ele_2,perc=get_points(n,
                                                  position_pro_L,
                                                  position_pro_R,
                                                  a_last,
                                                  beta)
   
    w = np.full(n, 1)   #setting original weights
   
    progress_counter = 0
    Energys = np.array(())
    E_n_array=np.array(())
   
    std, integral=integral_cal(position_ele_1,
                               position_ele_2,
                               position_pro_L,
                               position_pro_R,
                               beta,
                               n_combin,
                               a_last)
 
    Energys=np.append(Energys,integral*eV)
   
    E_attempt = np.array(())
    for t_step in t_steps:
        if t_step==10**(-19) or t_step==dtau:
            for i in range(2):
                # w = np.full(n, 1)   #setting original weights
                drift_v = drift(position_ele_1,
                                position_ele_2,
                                position_pro_L,
                                position_pro_R,
                                a_last,
                                S)            #drift for new posisitions
                new_pos = new_r(position_ele_1,
                                position_ele_2,
                                position_pro_L,
                                position_pro_R,
                                a_last,
                                S,
                                drift_v,
                                t_step,
                                n)
               
                #Finding E_n to maintain total weights = 1
                std,step_en = integral_cal(new_pos[:3,:],new_pos[3:,:],position_pro_L,position_pro_R,beta,n_combin,a_last)
                step_en=step_en*eV
               
                position_ele_1=new_pos[:3,:]
                position_ele_2=new_pos[3:,:]
               
                E_n = (h_bar/t_step)*np.log(n/(sum(w*np.exp(-step_en*t_step/h_bar)))) #using big T
                w = w*np.exp(-(step_en-E_n)*t_step/h_bar) #using big T
              
                E_attempt = np.append(E_attempt,(1/n)*sum(w*step_en))            #3rd step
                E_n_array=np.append(E_n_array,E_n)
    
        # print(abs(E_attempt[-1]-E_attempt[-2])/eV)
        while(abs(E_attempt[-1]-E_attempt[-2])>error[progress_counter]*eV):

            drift_v = drift(position_ele_1,
                            position_ele_2,
                            position_pro_L,
                            position_pro_R
                            ,a_last,
                            S)            #drift for new posisitions
            new_pos = new_r(position_ele_1,
                            position_ele_2,
                            position_pro_L,
                            position_pro_R,
                            a_last,
                            S,
                            drift_v,
                            t_step
                            ,n)
            
            #Finding E_n to maintain total weights = 1
            std,step_en = integral_cal(new_pos[:3,:],new_pos[3:,:],position_pro_L,position_pro_R,beta,n_combin,a_last)
            step_en=step_en*eV
            
            position_ele_1=new_pos[:3,:]
            position_ele_2=new_pos[3:,:]
            
            E_n = (h_bar/t_step)*np.log(n/(sum(w*np.exp(-step_en*t_step/h_bar)))) #using big T
            w = w*np.exp(-(step_en-E_n)*t_step/h_bar) #using big T
            
            E_attempt = np.append(E_attempt,(1/n)*sum(w*step_en))            #3rd step
            E_n_array=np.append(E_n_array,E_n)
            
        Energys = np.append(Energys, E_attempt[-1])

        progress_counter += 1
        # print(progress_counter/num_steps*100)
   
    Energys=Energys/eV
    E_n_array=E_n_array/eV
    
    min_en = min(Energys)          #extrapolated minimum energy

    E_b = 2*(-13.6*eV)-(eV**2)/(4*np.pi*ep_0*S)-min_en*eV
    E_b = E_b/eV           #bindng en (-4.75)
    print("Bindinng energy:",E_b)
    return E_b

bin_array=np.array(())
for i in range(100):
    bin_en=time_ev_f(best_beta, best_pro_sep)
    bin_array=np.append(bin_array,bin_en)
    print(i)

print(np.average(bin_array))

np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/im_time_data_mult.csv',
            (bin_array), delimiter=',')

#%%
'section 6'
'imaginary time 1 run'
best_beta=1.3103448275862069
best_pro_sep=7.727272727272728e-11
whileloop=0
while(whileloop==0):
    sigma=0.034*angstrom
    dtau=sigma**2*e_mass/h_bar
    time_step=dtau
    n=2**16
    n_combin=5
    # time_step = 10**(-19)
    
    S=best_pro_sep
    
    a_values=func_to_get_a(S)
    a_last=a_values[-1]
    beta=best_beta/angstrom

    n_combin=5
    
    position_pro_L=np.array([-S/2,0,0])
    position_pro_R=np.array([S/2,0,0])
    
    num_steps = 20

    e_step=0.1
    error=np.array(e_step)
    for i in range(num_steps-1):
        e_step=e_step*0.6
        error =np.append(error,e_step)
        
    t_steps = np.array(time_step)
    for i in range(num_steps-1):
        time_step = time_step*(0.5)
        t_steps = np.append(t_steps,time_step)
    break
    #original poisitions (G=PHI^2) 
    position_ele_1,position_ele_2,perc=get_points(n,
                                                  position_pro_L,
                                                  position_pro_R,
                                                  a_last,
                                                  beta)
    
    w = np.full(n, 1)   #setting original weights
    
    progress_counter = 0
    Energys = np.array(())
    E_n_array=np.array(())
    
    std, integral=integral_cal(position_ele_1,position_ele_2,position_pro_L,position_pro_R,beta,n_combin,a_last)
    
    Energys=np.append(Energys,integral*eV)
    E_attempt = np.array(())
    for t_step in t_steps:
        # if t_step==10**(-19) or t_step==dtau:
        for i in range(2):
            # w = np.full(n, 1)   #setting original weights
            drift_v = drift(position_ele_1,
                            position_ele_2,
                            position_pro_L,
                            position_pro_R,
                            a_last,
                            S)            #drift for new posisitions
            new_pos = new_r(position_ele_1,
                            position_ele_2,
                            position_pro_L,
                            position_pro_R,
                            a_last,
                            S,
                            drift_v,
                            t_step,
                            n)
            
            #Finding E_n to maintain total weights = 1
            std,step_en = integral_cal(new_pos[:3,:],new_pos[3:,:],position_pro_L,position_pro_R,beta,n_combin,a_last)
            step_en=step_en*eV
            
            position_ele_1=new_pos[:3,:]
            position_ele_2=new_pos[3:,:]
            
            E_n = (h_bar/t_step)*np.log(n/(sum(w*np.exp(-step_en*t_step/h_bar)))) #using big T
            w = w*np.exp(-(step_en-E_n)*t_step/h_bar) #using big T
            
            E_attempt = np.append(E_attempt,(1/n)*sum(w*step_en))            #3rd step
            E_n_array=np.append(E_n_array,E_n)
    
        # print(abs(E_attempt[-1]-E_attempt[-2])/eV)
        while(abs(E_attempt[-1]-E_attempt[-2])>error[progress_counter]*eV):
            print('a')
            drift_v = drift(position_ele_1,
                            position_ele_2,
                            position_pro_L,
                            position_pro_R,
                            a_last,
                            S)            #drift for new posisitions
            new_pos = new_r(position_ele_1,
                            position_ele_2,
                            position_pro_L,
                            position_pro_R,
                            a_last,
                            S,
                            drift_v,
                            t_step,
                            n)
            
            #Finding E_n to maintain total weights = 1
            std,step_en = integral_cal(new_pos[:3,:],new_pos[3:,:],position_pro_L,position_pro_R,beta,n_combin,a_last)
            step_en=step_en*eV
            
            position_ele_1=new_pos[:3,:]
            position_ele_2=new_pos[3:,:]
            
            E_n = (h_bar/t_step)*np.log(n/(sum(w*np.exp(-step_en*t_step/h_bar)))) #using big T
            w = w*np.exp(-(step_en-E_n)*t_step/h_bar) #using big T
            
            E_attempt = np.append(E_attempt,(1/n)*sum(w*step_en))            #3rd step
            E_n_array=np.append(E_n_array,E_n)
            
        Energys = np.append(Energys, E_attempt[-1])
        progress_counter += 1
    # print(progress_counter/num_steps*100)
    
    Energys=Energys/eV
    E_n_array=E_n_array/eV
    
    min_en = min(Energys)          #extrapolated minimum energy
    E_b = 2*(-13.6*eV)-(eV**2)/(4*np.pi*ep_0*S)-min_en*eV
    E_b = E_b/eV           #bindng en (-4.75)
    print("Bindinng energy:",E_b)
    
    if E_b<5 and E_b>4:
        whileloop=1
#%%

x=np.linspace(0, len(Energys), len(Energys))


plt.scatter(x,Energys)



save energys, t_steps

np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/im_time_data_1_run.csv',
            (Energys), delimiter=',')

timee()
#%%
'section 7'
'beta for helium'
start_time = time.time()

#%%

'adjust before begin'
n_combin=5   ##from now on bin down to 2**5
beta_amount=30
betaa=np.linspace(0.1/angstrom, 1/angstrom,beta_amount)

S=0*angstrom
position_pro_L=np.array([-S/2,0,0])
position_pro_R=np.array([S/2,0,0])
n=2**16

at=50

a_values=func_to_get_a(S)
a_last=a_values[-1]

resu=np.array(())
std_resu=np.array(())
beta_progress=0
for i in betaa:
    print('beta progress', beta_progress/30)
    beta_progress+=1
  
    ate=np.array(())
    for j in range(at):
        position_ele_1,position_ele_2,perc=get_points(n,
                                                      position_pro_L,
                                                      position_pro_R,
                                                      a_last,
                                                      i)
      
        stand_err,integral=integral_cal(position_ele_1,
                                        position_ele_2,
                                        position_pro_L,
                                        position_pro_R,
                                        i,
                                        n_combin,
                                        a_last)
        ate=np.append(ate, integral)

    av_ate=np.average(ate)
    st_dev=(sum((ate-av_ate)**2)/at)**0.5

    resu=np.append(resu,av_ate)
    std_resu=np.append(std_resu,st_dev)

betaa_ang = betaa * angstrom
best_beta=betaa_ang[np.argmin(resu)]
#%%
# save resu, st_dev and betaa_ang
np.savetxt('C:/Users/paxoc/Documents/uni/year 4/MPhys/Coding/beta_data_helium.csv',
           (resu, std_resu, betaa_ang), delimiter=',')
