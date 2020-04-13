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

for i in range( 1, imax-1 ):
    ud2[i] = 0.5*(ud1[i+1] - ud1[i-1])/dx
    
    
dt = np.float64(input("Enter dt, dx=%s\n  "%dx ))
itermax = int( 1.0/dt ) 
c = 1.0
c = float(input("Enter c, +1.0 or -1.0 "))

alpha = c*dt/dx
eps = 1.0e-6

#matrix A
up = -np.sign(c)
A = np.array( [ [ (up*dx)**5, (up*dx)**4, (up*dx)**3],
                [5.0*(up*dx)**4, 4.0*(up*dx)**3, 3.0*(up*dx)**2],
                [20.0*(up*dx)**3, 12.0*(up*dx)**2, 6.0*up*dx] ] )

coef = np.array( [0.0, 0.0, 0.0] )
b = np.array( [0.0, 0.0, 0.0] )

xx = -c*dt


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

        un[i]  = a0*xx**5 + a1*xx**4 + a2*xx**3 + a3*xx**2 + a4*xx + u[i]
        ud1n[i] = 5.0*a0*xx**4 + 4.0*a1*xx**3 + 3.0*a2*xx**2 + 2.0*a3*xx + ud1[i]
        ud2n[i] = 20.0*a0*xx**3 + 12.0*a1*xx**2 + 6.0*a2*xx + ud2[i]

    #update
    u[:] = un[:]
    ud1[:] = ud1n[:]
    ud2[:] = ud2n[:]


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
plt.plot(x,u, 'bo-')
plt.title(display)
plt.ylabel("U")
plt.xlabel("x")
plt.savefig(filename)
plt.show()
#plt.show(block=False)
#plt.pause(3)
#plt.close()


