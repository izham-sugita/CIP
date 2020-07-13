import numpy as np
import matplotlib.pyplot as plt

#Changing the default size
#fig_size = plt.rcParams["figure.figsize"]
#fig_size[0] = 20
#fig_size[1] = 16
#plt.rcParams["figure.figsize"] = fig_size

imax = 1001
dx = 10.0/(imax-1)
u = np.ndarray((imax),dtype=np.float64)
un = np.ndarray((imax),dtype=np.float64)

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

alpha = c*dt/dx
eps = 1.0e-6
xx = -c*dt
steps = 10
alpha = 1.0

#calculating exact solution
uexact = np.zeros_like(u)
#calculating exact solution
for i in range(imax):
    r1 = itermax*dt + 4.0
    r2 = r1 + (6.0 - 4.0) #did this on purpose, a reminder
    if x[i] >=r1 and x[i] <= r2:
        uexact[i] = 1.0

for iter in range(itermax):

    cfl = abs(c*dt/dx)
    threshold = 0.0

    for i in range(2,imax-2):
        delta_plus_half = u[i+1] - u[i]
        delta_minus_half = u[i] - u[i-1]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i = abs( abs(delta_plus_half) - abs(delta_minus_half)  )/(sensor + eps)
        else:
            seta_i = 0.0

        delta_plus_half = u[i+2] - u[i+1]
        delta_minus_half = u[i+1] - u[i]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i_plus_1 = abs( abs(delta_plus_half) - abs(delta_minus_half)  )/(sensor + eps)
        else:
            seta_i_plus_1 = 0.0

        seta_plus_half = max( seta_i, seta_i_plus_1)
        
        delta_plus_half = u[i] - u[i-1]
        delta_minus_half = u[i-1] - u[i-2]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i_minus_1 = abs( abs(delta_plus_half) - abs(delta_minus_half)  )/(sensor + eps)
        else:
            seta_i_minus_1 = 0.0
        
        seta_minus_half = max( seta_i, seta_i_minus_1)

        #allows for automatic alpha definition
        #the question is, what order of diffusion? 2nd order? 4th order?
        alpha = max(seta_plus_half, seta_minus_half) 

        #alpha = 0.1
        
        up = -np.sign(c)
        iup = i + int(up)
        xx = -c*dt
        
        udif = ( u[iup] - u[i] )/dx*up

        Beta =  ( abs( ( udif - ud1[i] ) / (ud1[iup] - udif + eps ) ) - 1.0 )/(dx*up);
        
        a2 = ( ud1[i] - udif + ( ud1[iup] - udif )*( 1.0 + alpha*Beta*dx*up ) ) /(dx*dx)

        a3 = alpha*Beta*udif + (udif - ud1[i] )/(dx*up) -a2*(dx*up)  

        a4 = ud1[i] + alpha*Beta*u[i]
                        
        un[i]  = ( u[i] + a4*xx + a3*xx**2 + a2*xx**3 ) / (1.0 + alpha*Beta*xx)
        
        ud1n[i] = ( (1.0 + alpha*Beta*xx)*(a4 + 2.0*a3*xx + 3.0*a2*xx**2) \
                    - ( u[i] + a4*xx + a3*xx**2 + a2*xx**3 ) *(alpha*Beta) ) \
                    / ( (1.0 + alpha*Beta*xx)*(1.0 + alpha*Beta*xx))             
        
    #update
    u[:] = un[:]
    ud1[:] = ud1n[:]


    #if iter%steps == 0:
    #    num = str(iter)
    #    filename = "./data1D/f"+num.zfill(5)+".csv"
    #    fp = open(filename, "w")
    #    fp.write("x, phi, 1minusphi\n")
    #    for i in range(imax):
    #        xp = round( x[i], 4 )
    #        rphi = round( phi[i], 4)
    #        mrphi = round ( (1.0 - phi[i]) , 4 )
    #        str1 = str(xp)
    #        str2 = str(rphi)
    #        str3 = str(mrphi)
    #        comma = ","
    #        nextline = "\n"
    #        strall = str1+comma+str2+comma+str3+nextline
    #        fp.write(strall)

    #    fp.close()

    current = iter*dt + dt
    display = "t = %.4f"%(current)
    
    
    #current = iter*dt + dt
    #display = "t = %.4f"%(current)
    plt.axis([0.0, 10.0, -0.5, 1.5 ] )
    plt.title(display)
    plt.ylabel("U")
    plt.xlabel("x")
    plt.plot(x,u,'bo-')
    plt.pause(0.001)
    plt.clf() #clear drawing
    

filename = "final.png"
plt.axis([0.0, 10.0, -0.5, 1.5 ] )
plt.plot(x,u, 'bo-')
plt.plot(x,uexact, 'kv-')
plt.title(display)
plt.ylabel("U")
plt.xlabel("x")
plt.savefig(filename)
plt.show()
#plt.show(block=False)
#plt.pause(3)
#plt.close()

filename = "final-cip-with-limiter.csv"
fp = open(filename, "w")
fp.write("x, u\n")
for i in range(imax):
    xp = round( x[i], 8 )
    rphi = round( u[i], 8)
    str1 = str(xp)
    str2 = str(rphi)
    comma = ","
    nextline = "\n"
    strall = str1+comma+str2+nextline
    fp.write(strall)

fp.close()
