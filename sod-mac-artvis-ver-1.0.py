import numpy as np
import matplotlib.pyplot as plt

imax = 1001
L = 10.0 #to match canonical example
dx = L /(imax-1)

tend = 1.0
dt = float(input("Enter dt, dx=%.8f "%(dx) ))
print("dt = %.8f"%(dt))

itermax = tend/(dt)
print("Maximum iteration: %d "%(itermax) )

#Physical quantities vector
Q   = np.ndarray((imax, 3))
Qn = np.ndarray((imax, 3))
Qnn = np.ndarray((imax, 3))

#Flux vector
F  = np.ndarray((imax, 3))

#velocity
u   = np.ndarray((imax))
un = np.ndarray((imax))

#density
rou = np.ndarray((imax))
roun = np.ndarray((imax))

#internal energy
e =  np.ndarray((imax))
en =  np.ndarray((imax))

#pressure
p =  np.ndarray((imax))
pn =  np.ndarray((imax))

#specific heat ratio
gamma = 1.4

#Initial condition
xg =  np.ndarray((imax))
xstd =  np.ndarray((imax))

#Initial condition for physical quantities
for i in range(imax):
    xg[i] = i*dx
    if xg[i] <= 0.5*L:
        rou[i] = 1.0
        u[i] = 0.0
        p[i] = 1.0
    else:
        rou[i] = 0.125
        u[i] = 0.0
        p[i] = 0.1

#stdize grid
# making 0<=xstd<=1.0
x0 = xg[0]
for i in range(imax):
    xstd[i] = ( xg[i] - x0 ) / L

#Initial value for un, roun, pn
un[:] = u[:]
roun[:] = rou[:]
pn[:] = p[:]
        
#Initial condition for vectors
for i in range(imax):
    Q[i][0] = rou[i]
    Q[i][1] = rou[i]*u[i]
    Q[i][2] = p[i]/(gamma - 1.0) + 0.5*rou[i]*u[i]*u[i]

    F[i][0] = rou[i]*u[i]
    F[i][1] = rou[i]*u[i]*u[i] + p[i]
    F[i][2] = ( (p[i]/(gamma - 1.0) + 0.5*rou[i]*u[i]*u[i]) + p[i] ) * u[i]

#Initial value for Qn
Qn[:][:] = Q[:][:]
Qnn[:][:] = Q[:][:]

cfl = dt/dx
print("dt/dx: %.4f"%(cfl) )

ib0 = 0
ib1 = imax-1
eps = 1.0e-8
threshold = eps
iter = 0
while iter < itermax:
    #implementing MacCormack scheme
    #Predictor step
    for i in range(2, imax-2):
        #Qn[i][0] = Q[i][0]  - cfl*( F[i][0] - F[i-1][0] )
        #Qn[i][1] = Q[i][1]  - cfl*( F[i][1] - F[i-1][1] )
        #Qn[i][2] = Q[i][2]  - cfl*( F[i][2] - F[i-1][2] )

        #Better to use i,i-1 first
        Qn[i][0] = Q[i][0]  - cfl*( F[i+1][0] - F[i][0] )
        Qn[i][1] = Q[i][1]  - cfl*( F[i+1][1] - F[i][1] )
        Qn[i][2] = Q[i][2]  - cfl*( F[i+1][2] - F[i][2] )

    #update intermediate flux
    for i in range(2, imax-2):
        F[i][0] = Qn[i][1]
        press = (gamma-1.0)*( Qn[i][2] - 0.5*( (Qn[i][1]*Qn[i][1])/Qn[i][0] ) )
        F[i][1] = (Qn[i][1]*Qn[i][1])/Qn[i][0] + press 
        velocity = Qn[i][1]/Qn[i][0]
        F[i][2] = ( Qn[i][2] + press )*velocity 

    #Corrector step
    for i in range(2, imax-2):
        delta_plus_half = Q[i+1][0] - Q[i][0]
        delta_minus_half = Q[i][0] - Q[i-1][0]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i = abs( (abs(delta_plus_half) - abs(delta_minus_half) )/ (sensor+eps) )
        else:
            seta_i = eps

        delta_plus_half = Q[i+2][0] - Q[i+1][0]
        delta_minus_half = Q[i+1][0] - Q[i][0]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i_plus_1 = abs( ( abs(delta_plus_half) - abs(delta_minus_half) )/ (sensor+eps) )
        else:
            seta_i_plus_1 = eps

        seta_plus_half = max(seta_i, seta_i_plus_1 )

        delta_plus_half = Q[i][0] - Q[i-1][0]
        delta_minus_half = Q[i-1][0] - Q[i-2][0]
        sensor = abs(delta_plus_half) + abs(delta_minus_half)
        if sensor > threshold:
            seta_i_minus_1 = abs( (abs(delta_plus_half) - abs(delta_minus_half) )/(sensor+eps) )
        else:
            seta_i_minus_1 = eps

        seta_minus_half = max( seta_i, seta_i_minus_1)

                
        #Qnn[i][0] = 0.5*(Q[i][0] + Qn[i][0])  - 0.5*cfl*( F[i+1][0] - F[i][0] ) \
        #            + 0.125*( seta_plus_half*(Q[i+1][0]-Q[i][0]) - seta_minus_half*(Q[i][0]-Q[i-1][0])  )
        #Qnn[i][1] = 0.5*(Q[i][1] + Qn[i][1])  - 0.5*cfl*( F[i+1][1] - F[i][1] ) \
        #            + 0.125*( seta_plus_half*(Q[i+1][1]-Q[i][1]) - seta_minus_half*(Q[i][1]-Q[i-1][1])  )
        #Qnn[i][2] = 0.5*(Q[i][2] + Qn[i][2])  - 0.5*cfl*( F[i+1][2] - F[i][2] ) \
        #            + 0.125*( seta_plus_half*(Q[i+1][2]-Q[i][2]) - seta_minus_half*(Q[i][2]-Q[i-1][2])  )

        Qnn[i][0] = 0.5*(Q[i][0] + Qn[i][0])  - 0.5*cfl*( F[i][0] - F[i-1][0] ) \
                    + 0.125*( seta_plus_half*(Q[i+1][0]-Q[i][0]) - seta_minus_half*(Q[i][0]-Q[i-1][0])  )
        Qnn[i][1] = 0.5*(Q[i][1] + Qn[i][1])  - 0.5*cfl*( F[i][1] - F[i-1][1] ) \
                    + 0.125*( seta_plus_half*(Q[i+1][1]-Q[i][1]) - seta_minus_half*(Q[i][1]-Q[i-1][1])  )
        Qnn[i][2] = 0.5*(Q[i][2] + Qn[i][2])  - 0.5*cfl*( F[i][2] - F[i-1][2] ) \
                    + 0.125*( seta_plus_half*(Q[i+1][2]-Q[i][2]) - seta_minus_half*(Q[i][2]-Q[i-1][2])  )


    #update
    for i in range(2, imax-2):
        Q[i][0] = Qnn[i][0]
        Q[i][1] = Qnn[i][1]
        Q[i][2] = Qnn[i][2]

    #update new flux
    for i in range(2, imax-2):
        F[i][0] = Q[i][1]
        press = (gamma-1.0)*( Q[i][2] - 0.5*( (Q[i][1]*Q[i][1])/Q[i][0] ) )
        F[i][1] = (Q[i][1]*Q[i][1])/Q[i][0] + press 
        velocity = Q[i][1]/Q[i][0]
        F[i][2] = ( Q[i][2] + press )*velocity


    for i in range(imax):
        press = (gamma-1.0)*( Q[i][2] - 0.5*( (Q[i][1]*Q[i][1])/Q[i][0] ) )
        velocity = Q[i][1]/Q[i][0]
        roun[i] = press
        
        '''
    if iter%1 == 0:
        num = str(iter)
        filename = "./data/f"+num.zfill(5)+".csv"
        f = open(filename, "w")
        f.write("x, q0, q1, q2, velocity, pressure\n")
        for i in range(imax):
            velocity = Q[i][1]/Q[i][0]
            press = (gamma-1.0)*( Q[i][2] - 0.5*( (Q[i][1]*Q[i][1])/Q[i][0] ) )
            strx = str(xg[i])
            strq0 = str(Q[i][0])
            strq1 = str(Q[i][1])
            strq2 = str(Q[i][2])
            strv = str(velocity)
            strp = str(press)
            comma = ","
            nline ="\n"
            strall = strx+comma+strq0+comma+strq1+comma+strq2+comma+strv+comma+strp+nline
            f.write(strall)

        f.close()
        '''
    

    '''
    #plot
    current = iter*dt + dt
    display = "t = %.4f"%(current)
    plt.axis([0.0, 1.0+0.2, 0.0, 1.2 ] )
    plt.title(display)
    plt.ylabel("P")
    plt.xlabel("x")
    plt.plot(xstd, roun, 'bo-')
    plt.pause(0.001)
    plt.clf() #clear drawing
    '''

    if iter%100 == 0:
        print("Iteration : %d"%(iter))

    iter +=1


current = iter*dt + dt
display = "t = %.4f"%(current)
plt.axis([0.0, 1.0+0.2, 0.0, 1.2 ] )
plt.title(display)
plt.ylabel("P")
plt.xlabel("x")
plt.plot(xstd, roun, 'bo-')
plt.show()

