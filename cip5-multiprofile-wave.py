import numpy as np
import matplotlib.pyplot as plt

#Changing the default size
#fig_size = plt.rcParams["figure.figsize"]
#fig_size[0] = 20
#fig_size[1] = 16
#plt.rcParams["figure.figsize"] = fig_size

imax = 2001
imax = int( input("Enter imax ") )
length = 2.0 #-1<=x<=1
dx = length/(imax-1)
u = np.ndarray((imax),dtype=np.float64)
un = np.ndarray((imax),dtype=np.float64)

ud1 = np.zeros_like(u)
ud1n = np.zeros_like(u)

ud2 = np.zeros_like(u)
ud2n = np.zeros_like(u)

x = np.ndarray((imax),dtype=np.float64)

'''
for i in range(imax):
    x[i] = i*dx
    u[i] = 0.0
    un[i] =0.0
    if x[i] >= 4.0 and x[i] <= 6.0:
        u[i] = 1.0
        un[i]=1.0
'''

u[:] = 0.0
un[:] = 0.0
#multiple wave profile
for i in range(imax):
    x[i] = -1.0 + i*dx

    if x[i] >=-0.8 and x[i] <=-0.6:
        u[i] = np.exp( -np.log(2.0)*(x[i]+0.7)**2 / 0.0009  )
        un[i] = u[i]
    elif x[i] >=-0.5 and x[i] <=-0.2:
        u[i] = 1.0
        un[i] = u[i]
    elif x[i] >=0.0 and x[i] <=0.2:
        u[i] = 1.0 - abs(10.0*x[i] - 1.0)
        un[i] = u[i]
    elif x[i] >=0.4 and x[i] <=0.6:
        u[i] = np.sqrt( 1.0 - 100.0*(x[i] - 0.5)**2  )
        un[i] = u[i]


#Initiate derivatives value
for i in range( 1, imax-1 ):
    ud1[i] = 0.5*(u[i+1] - u[i-1])/dx

for i in range( 1, imax-1 ):
    ud2[i] = 0.5*(ud1[i+1] - ud1[i-1])/dx
    
    
dt = np.float64(input("Enter dt, dx=%s\n  "%dx ))
elapsed = 10.0
itermax = int( elapsed/dt )-int(elapsed/2.0) #adjusted timestep; don't know why
print("Maximum iteration: ", itermax)
c = 1.0
c = float(input("Enter c, +1.0 or -1.0 "))

alpha = c*dt/dx
eps = 1.0e-6

uexact = np.zeros_like(u)

'''
#calculating exact solution
for i in range(imax):
    r1 = itermax*dt + 4.0
    r2 = r1 + (6.0 - 4.0) #did this on purpose, a reminder
    if x[i] >=r1 and x[i] <= r2:
        uexact[i] = 1.0
'''
uexact[:] = u[:]


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

        #minmod limiter
        ratio = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps)
        phi0 = min(10.0*dx, ratio) #default is 1.0
        phi[iup] = max(0.0, phi0)
        #phi[iup] = 0.0

        #van Leer (continuous function) #very diffusive
        #ratio = (u[i] - u[i-1]) / (u[i+1] - u[i] + eps)
        #phi[iup] = (ratio + abs(ratio)) / (1.0 + ratio)

        
        #un[i]  = a0*xx**5 + a1*xx**4 + a2*xx**3 + a3*xx**2 + a4*xx + u[i]

        un[i]  = u[i] + (1.0-phi[iup])*(a4*xx + a3*xx**2 + a2*xx**3 + a1*xx**4  + a0*xx**5) \
                 + phi[iup]*(udif*xx) 

        ud1n[i] = (1.0 - phi[iup])*( 5.0*a0*xx**4 + 4.0*a1*xx**3 + 3.0*a2*xx**2 + 2.0*a3*xx \
                  + ud1[i] ) + phi[iup]*udif
        
        # weight 0.98, 0.01 is the least diffusive
        #putting weight only on the first derivative
        #un[i]  = u[i] + (1.0 - phi[iup])*(a4*xx) + a3*xx**2 + a2*xx**3 + a1*xx**4  + a0*xx**5 \
        #         + phi[iup]*(udif*xx) 
        #ud1n[i] =  5.0*a0*xx**4 + 4.0*a1*xx**3 + 3.0*a2*xx**2 + 2.0*a3*xx \
        #          + (1.0 - phi[iup])*ud1[i] + phi[iup]*udif

        #the second derivative is not affected
        ud2n[i] = 20.0*a0*xx**3 + 12.0*a1*xx**2 + 6.0*a2*xx + ud2[i]
        

    #update periodic BC
    u[0] = un[imax-2]
    ud1[0] = ud1n[imax-2]
    ud2[0] = ud2n[imax-2]
    
    u[imax-1] = un[imax-2]
    ud1[imax-1] = ud1n[imax-2]
    ud2[imax-1] = ud2n[imax-2]


    for i in range(1, imax-1):
        u[i] = un[i]
        ud1[i] = ud1n[i]
        ud2[i] = ud2n[i]

    #update periodic BC
    #u[imax-1] = un[imax-2]
    #ud1[imax-1] = ud1n[imax-2]
    #ud2[imax-1] = ud2n[imax-2]

    #u[0] = un[imax-2]
    #ud1[0] = ud1n[imax-2]
    #ud2[0] = ud2n[imax-2]
        
    '''    
    #update
    u[:] = un[:]
    ud1[:] = ud1n[:]
    ud2[:] = ud2n[:]
    '''

    #if iter%steps == 0:
    #    num = str(iter)
    #    filename = "./data1D/f"+num.zfill(5)+".csv"
    #    fp = open(filename, "w")
    #    fp.write("x, u\n")
    #    for i in range(imax):
    #        str1 = str(x[i])
    #        str2 = str(u[i])
    #        comma = ","
    #        nextline = "\n"
    #        strall = str1+comma+str2+nextline
    #        fp.write(strall)
    #    fp.close()

    current = iter*dt + dt
    display = "t = %.4f"%(current)
    phi[:] = 0.0
    
    current = iter*dt + dt
    display = "t = %.4f"%(current)
    #plt.axis([0.0, 10.0, -0.5, 1.5 ] )
    
    plt.axis([-2.0, 2.0, -0.5, 1.5 ] )
    plt.title(display)
    plt.ylabel("U")
    plt.xlabel("x")
    plt.plot(x,u,'bo-')
    plt.pause(0.001)
    plt.clf() #clear drawing
    

filename = "final.png"
#plt.axis([0.0, 10.0, -0.5, 1.5 ] )
plt.axis([-2.0, 2.0, -0.5, 1.5 ] )
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

