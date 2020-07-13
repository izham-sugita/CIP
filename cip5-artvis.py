import numpy as np
import matplotlib.pyplot as plt

#Changing the default size
#fig_size = plt.rcParams["figure.figsize"]
#fig_size[0] = 20
#fig_size[1] = 16
#plt.rcParams["figure.figsize"] = fig_size

imax = 2001
imax = int( input("Enter imax ") )
dx = 10.0/(imax-1)
u = np.ndarray((imax),dtype=np.float64)
un = np.ndarray((imax),dtype=np.float64)
un1 = np.ndarray((imax),dtype=np.float64)

ud1 = np.zeros_like(u)
ud1n = np.zeros_like(u)

ud2 = np.zeros_like(u)
ud2n = np.zeros_like(u)

x = np.ndarray((imax),dtype=np.float64)

for i in range(imax):
    x[i] = i*dx
    u[i] = 0.0
    un[i] =0.0
    if x[i] >= 4.0 and x[i] <= 6.0:
        u[i] = 1.0
        un[i]=1.0

un1[:] = u[:]
        
#Initiate derivatives value
for i in range( 1, imax-1 ):
    ud1[i] = 0.5*(u[i+1] - u[i-1])/dx

for i in range( 1, imax-1 ):
    ud2[i] = 0.5*(ud1[i+1] - ud1[i-1])/dx
    
    
dt = np.float64(input("Enter dt, dx=%s\n  "%dx ))
itermax = int( 2.0/dt )
print("Maximum iteration: ", itermax)
c = 1.0
c = float(input("Enter c, +1.0 or -1.0 "))

alpha = c*dt/dx
eps = 1.0e-6


uexact = np.zeros_like(u)
#calculating exact solution
for i in range(imax):
    r1 = itermax*c*dt + 4.0
    r2 = r1 + 2.0 #did this on purpose, a reminder
    if x[i] >=r1 and x[i] <= r2:
        uexact[i] = 1.0
    


#matrix A
up = -np.sign(c)
A = np.array( [ [ (up*dx)**5, (up*dx)**4, (up*dx)**3],
                [5.0*(up*dx)**4, 4.0*(up*dx)**3, 3.0*(up*dx)**2],
                [20.0*(up*dx)**3, 12.0*(up*dx)**2, 6.0*up*dx] ] )

coef = np.array( [0.0, 0.0, 0.0] )
b = np.array( [0.0, 0.0, 0.0] )

xx = -c*dt
steps = 1
eps = 1.0e-8

phi = np.zeros_like(u)

for iter in range(itermax):

    for i in range(1,imax-1):
        up = -np.sign(c)
        iup = i + int(up)
        xx = -c*dt

        b[0] = ( u[iup] - u[i] ) -0.5*ud2[i]*dx*dx - ud1[i]*up*dx
        b[1] = ( ud1[iup] - ud1[i] ) - ud2[i]*up*dx
        b[2] = ud2[iup] - ud2[i]
        coef = np.linalg.solve(A, b)

        a0 = coef[0]
        a1 = coef[1]
        a2 = coef[2]
        a3 = ud2[i]*0.5
        a4 = ud1[i]

        #limiter
        udif = ( u[iup] - u[i] )/dx*up
        
        un[i]  = a0*xx**5 + a1*xx**4 + a2*xx**3 + a3*xx**2 + a4*xx + u[i]
        ud1n[i] = 5.0*a0*xx**4 + 4.0*a1*xx**3 + 3.0*a2*xx**2 + 2.0*a3*xx \
                  + ud1[i]  
        #the second derivative is not affected
        ud2n[i] = 20.0*a0*xx**3 + 12.0*a1*xx**2 + 6.0*a2*xx + ud2[i]


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

        un1[i] = un[i] + 0.125*( seta_plus_half*(un[i+1]-un[i]) - seta_minus_half*(un[i]-un[i-1]))
            

    #update
    #u[:] = un[:]
    u[:] = un1[:]
    ud1[:] = ud1n[:]
    ud2[:] = ud2n[:]
    

    if iter%steps == 0:
        print(iter)
        #num = str(iter)
        #filename = "./data1D/f"+num.zfill(5)+".csv"
        #fp = open(filename, "w")
        #fp.write("x, u\n")
        #for i in range(imax):
        #    str1 = str(x[i])
        #    str2 = str(u[i])
        #    comma = ","
        #    nextline = "\n"
        #    strall = str1+comma+str2+nextline
        #    fp.write(strall)

        #fp.close()

    current = iter*dt + dt
    display = "t = %.4f"%(current)
    
    
    
    #plt.axis([0.0, 10.0, -0.5, 1.5 ] )
    #plt.title(display)
    #plt.ylabel("U")
    #plt.xlabel("x")
    #plt.plot(x,u,'bo-')
    #plt.pause(0.001)
    #plt.clf() #clear drawing
    

filename = "final.png"
plt.axis([0.0, 10.0, -0.5, 1.5 ] )
plt.plot(x,u, 'bo-', x, uexact,'kv-')
plt.title(display)
plt.ylabel("U")
plt.xlabel("x")
plt.savefig(filename)
plt.show()
#plt.show(block=False)
#plt.pause(3)
#plt.close()

filename = "cip5-final.csv"
fp = open(filename, "w")
fp.write("x, u\n")
for i in range(imax):
    str1 = str(x[i])
    str2 = str(u[i])
    comma = ","
    nextline = "\n"
    strall = str1+comma+str2+nextline
    fp.write(strall)

fp.close()

