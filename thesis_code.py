# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 11:29:07 2021

@author: Beth
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


Gas_type = 1 #1 Neon, 2 Nitrogen
N_simulation = 3000000 #define the number of particles in each simulation


if Gas_type == 1 : # Neon
        gamma = 5/3 # heat compacity ratio, 1.67
        C1 = 3.232 # table 2.2
        C2 = -0.7563
        C3 = 0.3937
        C4 = -0.0729
        A = 0.778 # 0.527; table2.3
        B = 0.495 # 0.545;
        C6_k = 0.758e-43*1e-12 #K*m^6, table 2.5
        sigma = 2.75e-10 # meter
        col_sigma = 5e-10 #metre square
        radius = 154e-12
        M = 20.1797e-3  # kg/mol  mole mass
        kappa = 1.98 #or 2.0, table 2.4
        n_s = 1.022e21 
        a = 0.806  #for Zref calculation
        d_Rn = 0  # for virtual source location;
elif Gas_type == 2: #Nitrogen
        gamma = 7/5 #heat compacity ratio, 1.4
        C1 = 3.606 #table 2.2
        C2 = -1.742
        C3 = 0.9226
        C4 = -0.2069
        A = 0.783
        B = 0.353
        C6_k = 6.2e-43*1e-12 #K*m^6, table 2.5
        radius = 182e-12 #kinetic radius/ m
        sigma = 3.85e-10 #meter particle cross section
        col_sigma = 8.8e-10 #meter square collision diameter
        M = 28.0134e-3 #kg/mol  mole mass
        kappa = 1.38; #or 1.47, table 2.4
        n_s = 5.92e20
        a = 0.591 #for Zref calculation
        d_Rn = 0.85 #for virtual source location
        

#%%physical constants
R = 8.31 #gas constant
k_B = 1.38e-23
NA = 6.022e23
pi = np.pi
        
#%%setup known parameters of the experiment
P0 = 5e5 #bar inlet pressure(stagnation pressure)
T0 = 300 #K stagnation temperature
n0 = P0/k_B/T0 #stagnation density
v0 = np.sqrt((3*R*T0)/M) #root mean square velocity
a0 = np.sqrt((gamma*R*T0)/M) #speed of sound

#initial pressures in each chamber
if Gas_type == 1:
    Pb1 = 2.21e-1 #Pa
    Pb2 = 4.87e-4*3.5  #Pa
    Pb3 = 6.37e-5*3.5 #Pa
elif Gas_type == 2:
    Pb1 = 6.17e-1 #Pa 3.8e-3 mbar
    Pb2 = 1.35e-3 #Pa
    Pb3 = 1.03e-4 #Pa
    
#initial number density in each chamber
nb1 = Pb1/k_B/T0
nb2 = Pb2/k_B/T0
nb3 = Pb3/k_B/T0

#skimmer geometries
d_nozzle = 30e-6 #m nozzle diameter
d_sk1 = 180e-6 #m diameter of skimmer 1
d_sk2 = 400e-6 #m diameter of skimmer 2 
d_det = 0.25e-3 #m diameter of detector pinhole
#the skimmers also have a length and outer radius as they are conical shaped
sk1_len = 6.58e-3 #m skimmer 1 length
sk1_end_r = 4.15e-3 #m end of skimmer 1 radius
sk2_len = 6.09e-3 #m skimmer 2 length
sk2_end_r = 4.15e-3 #m end of skimmer 2 radius

noz_sk1 = 0.004 #m nozzle to skimmer 1 distance
noz_sk2 = 0.025 #m nozzle to skimmer 2 distance
noz_det = 287e-3 #m nozzle to detector distance
sk1_det = noz_det - noz_sk1

#%% parameters at the nozzle
mass_flow = 0.25* 4 * n0*a0*(d_nozzle**2)*((2/(gamma+1))**((gamma+1)/(2*gamma - 2))) #particles per second
nozzle_flux = mass_flow/ (NA*M)
I_0 = (kappa*nozzle_flux)/pi  #centre-line intensity 
x_M = np.sqrt(P0/Pb1)*0.67*d_nozzle #position of Mach disk /m
v_inf = np.sqrt(((2*R*T0)/M)*(gamma/(gamma-1))) #terminal velocity 
S_parallel_final  = A*(np.sqrt(2)*n0*d_nozzle*(53/T0*C6_k)**(1/3))**B #parallel speed ratio  
M_parallel_final = np.sqrt(2/gamma) * S_parallel_final #final Mach number at quitting surface
x_quit = d_nozzle*(M_parallel_final/C1)**(1/(gamma-1)); #quitting surface location


#up to mach disk the particles travel in a continuos flow region
#effective Mach number at skimmer 1 for a 

# if noz_sk1 <= x_M:
#     if Gas_type == 2:
#         Mach = 3.65*(((noz_sk1/d_nozzle)-0.4)**0.4)-0.82*(((noz_sk1/d_nozzle)-0.4))**(-0.4) #diatomic gas
#     else:
#         Mach = 3.26*(((noz_sk1/d_nozzle)-0.075)**0.67) - 0.61*(((noz_sk1/d_nozzle)-0.075))**(-0.67) #monoatomic gas
#     #   Another way to calculate Mach may be
# else:    
#   Mach = 3.65*((noz_sk1/d_nozzle)**(gamma-1)) - ((gamma+1)/(2*(gamma-1)))*(3.65*((noz_sk1/d_nozzle)**(gamma-1)))**(-1)

x_q = x_quit #define the quitting surface


#Centre line Mach number
xoverd = x_q/d_nozzle
Mach = (xoverd**(gamma-1))*(C1 + (C2/xoverd) + (C3/xoverd**2) + (C4/xoverd**3))
    

#free ideal gas expansion
# laws of isentropic expansion
def isentropicT(T0,Mach): #isentropic temperature
    T = T0/(1 + (0.5*(gamma-1)*(Mach**2)))
    return T

def isentropicn(n0,Mach): #isentropic density
    n = n0/((1+((gamma-1)/2)*Mach**2)**(1/(gamma-1)))
    return n

def beamvelocity(Mach): #isentropic velocity
    v = Mach * np.sqrt((gamma*R*T0)/M) * (1+(0.5*(gamma-1))*(Mach**2))**(-0.5)
    return v

def isentropicP(Mach): #isentropic pressure
    P = (1+(0.5*(gamma-1))*(Mach**2))**(-gamma/(gamma-1))
    return P 
                         
#collision rate
z_jet = np.sqrt(2)* n0 * sigma*v0**(1+0.5*(gamma-1)*Mach**2)**-(0.5)*((gamma+1)/(gamma-1))

def meanfreepath(n):
    lamda = 1/sigma*np.sqrt(2)*n
    return lamda


T_q = isentropicT(T0,Mach)
#density at the quitting surface 
n_q = isentropicn(n0,Mach)
#speed of sound at the quitting surface 
c_q = np.sqrt((gamma*R*T_q)/M) #speed of sound
 

v_mean = beamvelocity(Mach) #mean velocity 
v_a = np.sqrt(R*T_q/M) #distribution of velocities

lamda = meanfreepath(n0)
lamda_sk1 = meanfreepath(n_q)
Kn = lamda/(d_nozzle) #knudsen number
Kn_back = Kn*(((2/5) * v_a**2))**(-2/12) #with backscattering
Mach_final = np.sqrt(2/(gamma-1))*np.sqrt(T0/T_q)

#%% collision term
#density at skimmer 2
nc1 = (n0*(Mach**2)*(d_sk1**2)*gamma)/(8*((sk1_det)**2))
nc21 = ((1+(gamma - 1)*(Mach**2))/2)**(1/(gamma-1))
#collisional loss terms
nc2a = (pi*((2*sigma)**2)*d_sk1*n0)/np.sqrt(2)
nc2b = (1 - (np.sqrt(gamma/(2*pi))*((Mach*d_sk1)/(2*(sk1_det)))))
n_det = nc1/ (nc21 + (nc2a*nc2b))
#magnitude of the collisional term
coll = (nc1/(nc21 + (nc2a*nc2b)))/(nc1/nc21)

        
#%% background attenuation    
Q = pi*sigma**2
att_eff = 1-np.exp(-S_parallel_final**2*(d_sk1/2/x_q)**2 * ((noz_sk1+sk1_det)/sk1_det)**2) *np.exp(-nb1*Q*noz_sk1-nb2*Q*0.025-nb3*Q*((noz_sk1+sk1_det)-0.025))
#this is the magnitude of the attenuation effect which will be
#applied to the final density
        
#%% single tracing       
##from quitting surface, partical ray tracing

  
z1 = noz_sk1 

def cosd(x):
    return np.cos(x * pi / 180)

def single_tracing(z,X0,Y0,Vx0,Vy0,z0,Vz0):
    T = (z-z0)/Vz0
    X = X0+Vx0*T
    Y = Y0+Vy0*T
    Z = z0+Vz0*T
    return X,Y,T,Z

#%% apply nozzle angle
#rotate all particles by angle:
a = 0  #degrees 
for angle in a:
    
    angle = angle*(np.pi/180) #radians
   
    XX1 = []
    YY1 = []
    Vx1 = []
    Vy1 = []
    Vz1 = []
    
    N_tot = 0
    theta0 = 180
        
    for n in range(10000): #deine the number of simulations to run
        # print(n)
        
        #for sperical quitting surface
        theta = np.arccos(np.random.rand(N_simulation))        
        phi = np.random.uniform(0,2*pi,N_simulation)
        
        #positions on quitting surface
        X0 = (x_q)*np.sin(theta)*np.cos(phi) #meter
        Y0 = (x_q)*np.sin(theta)*np.sin(phi) #meter
        Z0 = (x_q)*np.cos(theta) #meter
        
        position = [X0,Y0,Z0]
        
        # gaussian velocities
        VR0 = v_mean + np.random.normal(0,v_a,N_simulation)
        v_phi0 = np.random.normal(0,v_a,N_simulation)
        v_theta0 = np.random.normal(0,v_a,N_simulation)

        
        #spherical to cartesian        
        Vx0 = VR0*np.sin(theta)*np.cos(phi)  +v_theta0*np.cos(theta)*np.cos(phi)   -v_phi0*np.sin(phi)
        Vy0 = VR0*np.sin(theta)*np.sin(phi)  +v_theta0*np.cos(theta)*np.sin(phi)   +v_phi0*np.cos(phi)
        Vz0 = VR0*np.cos(theta)              -v_theta0*np.sin(theta);    
        velocities = [Vx0, Vy0, Vz0]
    
       
    #%% add a rotation
    #treat nozzle at an angle as a rotation of all particles about an axis
        
        #about x
        # rotation = np.array([[1, 0, 0],[0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])      
        # X0 = position[0]*rotation[0,0] + position[1]*rotation[0,1] + position[2]*rotation[0,2]
        # Y0 = position[0]*rotation[1,0] + position[1]*rotation[1,1] + position[2]*rotation[1,2]
        # Z0 = position[0]*rotation[2,0] + position[1]*rotation[2,1] + position[2]*rotation[2,2]
      
        #about y
        # rotation = np.array([[np.cos(angle), 0, np.sin(angle)],[0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
        # X0 = position[0]*rotation[0,0] + position[1]*rotation[0,1] + position[2]*rotation[0,2]
        # Y0 = position[0]*rotation[1,0] + position[1]*rotation[1,1] + position[2]*rotation[1,2]
        # Z0 = position[0]*rotation[2,0] + position[1]*rotation[2,1] + position[2]*rotation[2,2]
        
        #rotatedposition = [X0,Y0,Z0]

#%%        
        #new theta + phi
        
        # phi1 = np.arctan2(Y0,X0)
        
           
        # for i in range(len(X0)):
        #     if phi1[i]<0:
        #         phi1[i] = 2*pi + phi1[i]
                
        #     else:
        #         phi1[i] = phi1[i]
    
        
        
        # theta1 = np.arccos((np.sqrt(X0**2 + Y0**2))/(x_q+(np.zeros(N_simulation))))
   
        # #new velocities
        # Vx0 = VR0*np.sin(theta1)*np.cos(phi1)  +v_theta0*np.cos(theta1)*np.cos(phi1)   -v_phi0*np.sin(phi1)
        # Vy0 = VR0*np.sin(theta1)*np.sin(phi1)  +v_theta0*np.cos(theta1)*np.sin(phi1)   +v_phi0*np.cos(phi1)
        # Vz0 = VR0*np.cos(theta1)               -v_theta0*np.sin(theta1);  
        
        
        # #add rotation
        # Vx0 = velocities[0]*rotation[0,0] + velocities[1]*rotation[0,1] + velocities[2]*rotation[0,2]
        # Vy0 = velocities[0]*rotation[1,0] + velocities[1]*rotation[1,1] + velocities[2]*rotation[1,2]
        # Vz0 = velocities[0]*rotation[2,0] + velocities[1]*rotation[2,1] + velocities[2]*rotation[2,2]
#%%        
        
        z0 = Z0
        
        F = single_tracing(z1,X0,Y0,Vx0,Vy0,z0,Vz0)
        M1 =np.argwhere(((F[0])**2+(F[1])**2)**0.5 <= d_sk1/2)
        
        N_pass = len(M1)
        N_tot  = N_tot + N_pass
        
        
        #save the positions and locations of the particles that pass through the first skimmer orifice
        XX = F[0][M1]
        XX1 = np.append(XX1,XX)
        YY = F[1][M1]
        YY1 = np.append(YY1,YY)
        Vx = Vx0[M1]
        Vx1 = np.append(Vx1,Vx)
        Vy = Vy0[M1]
        Vy1 = np.append(Vy1,Vy)
        Vz = Vz0[M1]
        Vz1 = np.append(Vz1,Vz)


    X1 = XX1
    Y1 = YY1
    Vx1 = Vx1
    Vy1 = Vy1
    Vz1 = Vz1

    #partices through the end of the first skimmer       
    z2 = np.float64(z1 + sk1_len)
    F1 = single_tracing(z2,X1,Y1,Vx1,Vy1,z1,Vz1)
    M2 = np.array(np.argwhere(((F1[0])**2+(F1[1])**2)**0.5 <= 4.15e-3))
    X2 = F1[0][M2]
    Y2 = F1[1][M2]
    Vx2 = Vx1[M2]
    Vy2 = Vy1[M2]
    Vz2 = Vz1[M2]
    
    
    #%% second skimmer

    z3 = np.float64(z1 + 25e-3)
    F2 = single_tracing(z3,X2,Y2,Vx2,Vy2,z2,Vz2)
    M3 = np.array(np.argwhere(((F2[0])**2+(F2[1])**2)**0.5 <= 1e-3))
    
    X3 = F2[0][M3[:,0]]
    Y3 = F2[1][M3[:,0]]
    Vx3 = Vx2[M3[:,0]]
    Vy3 = Vy2[M3[:,0]]
    Vz3 = Vz2[M3[:,0]]
  
    #3rd skimmer
    z4 =np.float64(z1+300e-3)
    F3 = single_tracing(z4,X3,Y3,Vx3,Vy3,z4,Vz3)
    M4 = np.argwhere((F3[0]<= 1.8e-3) & (F3[0]>= -1.8e-3) & (F3[1]<= 0.5e-3) & (F3[1]>= -0.5e-3))

    X4 = F3[0][M4[:,0]]
    Y4 = F3[1][M4[:,0]]
    Vx4 = Vx3[M4[:,0]]
    Vy4 = Vy3[M4[:,0]]
    Vz4 = Vz3[M4[:,0]]

    #interaction point
    z5 =np.float64(z1+538e-3)
    F4 = single_tracing(z4,X4,Y4,Vx4,Vy4,z3,Vz4)
    M5 = np.array(np.argwhere(((F4[0])**2+(F4[1])**2)**0.5 <= 50e-3))
    
    X5 = F4[0][M5[:,0]]
    Y5 = F4[1][M5[:,0]]
    Vx5 = Vx4[M5[:,0]]
    Vy5 = Vy4[M5[:,0]]
    Vz5 = Vz4[M5[:,0]]

   
 #%% no skimmer 2 
    #looking at before the second skimmer
    # z6 =np.float64(z1+16e-3)
    # F5 = single_tracing(z6,X2,Y2,Vx2,Vy2,z2,Vz2)
    # M6 = np.argwhere((F5[0]**2+F5[1]**2)**0.5 <= 50e-3)
    # X6 = F5[0][M6[:,0]]
    # Y6 = F5[1][M6[:,0]]
    # Vx6 = Vx2[M6[:,0]]
    # Vy6 = Vy2[M6[:,0]]

    
    
    #%%density plot
    
plt.figure(dpi = 1200)
plt.clf()
x_min = np.min(X5) 
x_max = np.max(X5)  
y_min = np.min(Y5) 
y_max = np.max(Y5)   
x_bins = np.linspace(x_min, x_max, 25) 
y_bins = np.linspace(y_min, y_max, 25) 
plt.xlabel('X- position/m')  
plt.ylabel('Y- position/m') 
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
hist = plt.hist2d(X5[:,0],Y5[:,0], bins=(x_bins, y_bins),cmap = 'jet')

# compute the area of each bin
dx = np.diff(x_bins)
dy = np.diff(y_bins)
area = dx[0]*dy[0]

#hits = hist[0]/ (n*N_simulation)

#%% final density
fin_den_det = att_eff*n_q/area/n/N_simulation*2*x_q**2

# #%% add misalignment lines for inspection
(mux, sigmax) = norm.fit(X5)
(muy, sigmay) = norm.fit(Y5)

# mis = np.where(hist[0] == np.amax(hist[0])) 
misalignmentx = mux
misalignmenty = muy

x = np.linspace(x_min,x_max,100)
y = np.linspace(y_min,y_max,100)
plt.plot(x,np.zeros(100)+(misalignmenty),'b-')
plt.plot(np.zeros(100)+(misalignmentx),y,'b-')
plt.plot(x,np.zeros(100),'r-')
plt.plot(np.zeros(100),y,'r-')

#%% add density colour scale
plt.title('{} degrees'.format(angle/(np.pi/180)))
den = hist[0]*fin_den_det*2
plt.imshow(den,cmap  ='jet')
plt.colorbar(label = 'density m^-3')
plt.show()


#%% line scan at central bin
plt.figure(dpi = 1200)
plt.clf()
plt.xlabel('bin')
plt.ylabel('Intensity')
plt.plot(hist[0][12]/np.max(hist[0][25]))
plt.plot(hist[0][:,12]/np.max(hist[0][:,12]))


