import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Changing the default size
#fig_size = plt.rcParams["figure.figsize"]
#fig_size[0] = 20
#fig_size[1] = 16
#plt.rcParams["figure.figsize"] = fig_size

imax = 1001
imax = int(input("Enter imax: "))
dx = 10.0/(imax-1)
u = np.ndarray((imax),dtype=np.float64)
un = np.ndarray((imax),dtype=np.float64)
un1 = np.ndarray((imax),dtype=np.float64)

ud1 = np.zeros_like(u)
ud1n = np.zeros_like(u)
phi = np.zeros_like(u)

ud2 = np.zeros_like(u)
ud2n = np.zeros_like(u)

x = np.ndarray((imax),dtype=np.float64)

for i in range(imax):
    x[i] = i*dx
    u[i] = 0.0
    un[i] =0.0
    if x[i] > 4.0 and x[i] < 6.0:
        u[i] = 1.0
        un[i]=1.0

un1[:] = u[:]

#Initiate derivatives value
for i in range( 1, imax-1 ):
    ud1[i] = 0.5*(u[i+1] - u[i-1])/dx
    #ud1[i] = (u[i] - u[i-1])/dx

for i in range( 1, imax-1 ):
    ud2[i] = 0.5*(ud1[i+1] - ud1[i-1])/dx
    
    
dt = np.float64(input("Enter dt, dx=%s\n  "%dx ))
itermax = int( 2.0/dt ) 
c = 1.0
c = float(input("Enter c, +1.0 or -1.0 "))

switch = float( input("Artificial viscosity switch, on/off (1=on, 0=off) ") )

alpha = c*dt/dx
eps = 1.0e-8

xx = -c*dt
steps = itermax/10
print(itermax)

#calculating exact solution
uexact = np.zeros_like(u)
for i in range(imax):
    r1 = itermax*c*dt + 4.0
    r2 = r1 + 2.0 #did this on purpose, a reminder
    if x[i] >=r1 and x[i] <= r2:
        uexact[i] = 1.0

for iter in range(itermax):

    for i in range(1,imax-1):
        up = -np.sign(c)
        iup = i + int(up)
        xx = -c*dt
        
        udif = ( u[iup] - u[i] )/dx*up
        a2 = ( ud1[i] + ud1[iup] - 2.0*udif ) / (dx*dx)
        a3 = ( 3.0*udif - 2.0*ud1[i] - ud1[iup] ) / dx*up
        a4 = ud1[i]
                        
        un[i]  = u[i] + a4*xx + a3*xx**2 + a2*xx**3 
        ud1n[i] = ud1[i] + 2.0*a3*xx + 3.0*a2*xx**2 

    #artificial viscosity
    #Finding! Threshold at CFL number yields the best result; why?
    cfl = abs(c*dt/dx)
    threshold = cfl
    
    for i in range(2, imax-2):
        delta_plus_half = un[i+1] - un[i]
        delta_minus_half = un[i] - un[i-1]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i = (1.0 - cfl )*abs( abs(delta_plus_half) - abs(delta_minus_half)  )/(sensor + eps)
        else:
            seta_i = 0.0

        delta_plus_half = un[i+2] - un[i+1]
        delta_minus_half = un[i+1] - un[i]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i_plus_1 = (1.0 - cfl)*abs( abs(delta_plus_half) - abs(delta_minus_half)  )/(sensor + eps)
        else:
            seta_i_plus_1 = 0.0

        seta_plus_half = max( seta_i, seta_i_plus_1)

        delta_plus_half = un[i] - un[i-1]
        delta_minus_half = un[i-1] - un[i-2]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i_minus_1 = (1.0 - cfl)*abs( abs(delta_plus_half) - abs(delta_minus_half)  )/(sensor + eps)
        else:
            seta_i_minus_1 = 0.0

        seta_minus_half = max( seta_i, seta_i_minus_1)

        un1[i] = un[i] + switch*0.125*( seta_plus_half*(un[i+1]-un[i]) - seta_minus_half*(un[i]-un[i-1]))
        
    #update
    #u[:] = un[:]
    u[:] = un1[:]
    ud1[:] = ud1n[:]
    


    if iter%steps == 0:
        num = str(iter)
        filename = "./data1D/cip-ori-artvis-"+num.zfill(5)+".csv"
        df = pd.DataFrame({"x": x, "f": u})
        df.to_csv(filename, index=False)




    current = iter*dt + dt
    display = "t = %.4f"%(current)
    plt.axis([0.0, 10.0, -0.5, 1.5 ] )
    plt.title(display)
    plt.ylabel("U")
    plt.xlabel("x")
    plt.plot(x,u,'bo-')
    plt.pause(0.001)
    plt.clf() #clear drawing
    

filename = "final.png"
plt.axis([0.0, 10.0, -0.5, 1.5 ] )
plt.plot(x,u, 'bo-', x, uexact, 'kv-')
plt.title(display)
plt.ylabel("U")
plt.xlabel("x")
plt.savefig(filename)
plt.show()


filename = "final-cip-ori-artvis.csv"
df = pd.DataFrame({"x": x, "f": u})
df.to_csv(filename, index=False)
