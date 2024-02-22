'''Modelling Elliptical Orbits by Solving ODEs using 4th-order Runge Kutta'''

import matplotlib.pyplot as plt
import math
from math import factorial, sin
from pylab import plot,xlabel,show
import numpy as np
from numpy import linspace,sin,cos,ones,copy,tan,pi,array,arange, power, argmin
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings("ignore")

# Using the fourth-order Runge-Kutta method with an adaptive step size to solve two second order differential equations,
# derived from Newton's Second law, that describe the dynamics of a comet in eliptical orbit around the Sun (in the xy plane)

G=66377.002 #defining universal gravitational constant, in m^3/ kg*yr^2
M=1.989*10**(30)  #defining mass of Sun
x0=4*10**12 #defining initial x-value
y0=0 #defining initial y-value
v_x0=0  #defining initial v_x value given in problem (initial velocity in x-direction)
v_y0=15778476000  ##defining initial v_y given in problem (initial velocity in y-direction), in m/yr

def f(r, t): #creating function that helps calculate partial derivatives
    x=r[0] #indexing the zero-th term in whatever r1 array the user hands the function, which will be the x position, assigning it to variable x
    v_x=r[1] #indexing the 1st term in whatever r1 array the user hands the function, which will be x component of velocity, assigning it to v_x
    dx_dt=v_x #setting v1=dx_dt
    y=r[2]  #indexing the zero-th term in whatever r2 array the user hands the function, which will be the y position, assigning it to variable y
    v_y=r[3]  #indexing the 1st term in whatever r2 array the user hands the function, which will be x component of velocity, assigning it to v_y
    dy_dt=v_y  #setting v2=dx_dt
    rad_=np.sqrt(x**2+y**2)
    dvx_dt=-G*M*(x/rad_**3) #calculating equation of motion for x direction
    dvy_dt=-G*M*(y/rad_**3)  #calculating equation of motion for y direction
    return(array([dx_dt, dvx_dt, dy_dt, dvy_dt])) #returns two arrays containing the derivatives, one with velocity and position in x direction, one with velocity and position in y-dir

x_values, vx_values=[], [] #creating empty arrays for x values and x components of velocity to be stored in
y_values, vy_values=[], [] #creating empty array for y values and y components of velocity to be stored in

r=array([x0, v_x0, y0, v_y0], float) #r vector with the first component giving the initial x coordinate, second giving the initial v_x velocity value, third giving initial y coordinate, fourth giving initial v_y velocity value
#T=2*np.pi*np.sqrt(x_0**3/(G*M))
#h=T/2000000000 #calculating step size between data points
T=165
N=500000
h=T/N #calculating timestep
timepoints=arange(0, T, h) #creating an array of time values within the range defined by limits of integration above

for t in timepoints: #for loop to perform calculations of fourth-order Runge-Kutta method
    x_values.append(r[0]) #adds current x value (corresponding to index 0/the 0th term in the x array) to empty array defined above
    vx_values.append(r[1]) #adds current v1 value (corresponding to index 1/the 1st term in the r1 array) to empty array defined above
    y_values.append(r[2]) #adds current y value (corresponding to index 0/the 0th term in the y array) to empty array defined above
    vy_values.append(r[3]) #adds current v2 value (corresponding to index 0/the 0th term in the y array) to empty array defined above
    #Runge-Kutta fourth-order calculations for f()
    k1=h*f(r, t) #calculating k1 in Runge-Kutta fourth-order method
    k2=h*f(r+0.5*k1,t+0.5*h) #calculating k2 in Runge-Kutta fourth-order method
    k3=h*f(r+0.5*k2,t+0.5*h) #calculating k3 in Runge-Kutta fourth-order method
    k4=h*f(r+k3,t+h) #calculating k4 in Runge-Kutta fourth-order method
    r+=(k1+2*k2+2*k3+k4)/6 #calculating r using k1,k2,k3,k4 and adding it to the current r value

plt.figure(1) #creating first figure
plt.plot(x_values, y_values, color='purple') #plotting x and y values generated with fixed step size
plt.title('Dynamics of Comet in Orbit Around Sun') #naming figure
plt.xlabel('x (km)') #labelling axes
plt.ylabel('y (km)')
plt.savefig('orbitaldynamics.png') #saving figure with given name
plt.show()


# Part c) modifying program in part b to use adaptive step size
delta = 1000  # Given target accuracy in m/yr

def step(x, t, h):
    def rk4(x, t, h):
        k1 = h * f(x, t) #calculating k1 in Runge-Kutta fourth-order method
        k2 = h * f(x + 0.5 * k1, t + 0.5 * h) #calculating k2 in Runge-Kutta fourth-order method
        k3 = h * f(x + 0.5 * k2, t + 0.5 * h) #calculating k3 in Runge-Kutta fourth-order method
        k4 = h * f(x + k3, t + h) #calculating k4 in Runge-Kutta fourth-order method
        return((k1 + 2 * k2 + 2 * k3 + k4)/6)  #returning finsl result of Runge-Kutta calculations

    step_h1 = rk4(x, t, h) #using adaptive step in rk4 function that calculates fourth-order Runge-Kutta method for h in 2 steps
    step_h2 = rk4(x + step_h1, t + h, h)

    step_h1_h2 = step_h1 + step_h2 #using adaptive step in rk4 function that calculates fourth-order Runge-Kutta method for 2h in 1 step
    step_2h = rk4(x, t, 2 * h)

    delta_x1 = step_h1_h2[0]
    delta_x2 = step_2h[0]
    delta_y1 = step_h1_h2[2]
    delta_y2 = step_2h[2]
    error = np.sqrt((delta_x1 - delta_x2) ** 2 + (delta_y1 - delta_y2) ** 2) / 30 #calculating/estimating the error for step sizes
    rho = h * delta / error #calculating rho value using target accuracy and error
    factor = np.power(rho, 1 / 4)  #calculating factor for multiplication of step size using rho


    if rho >= 1: #creating the condition for an adapted step size for the given target accuracy
        t = t + 2 * h
        if factor > 2: #limiting the step size
            h *= 2
        else:
            h *= factor
        return step_h1_h2, h, t
    else:  #program will recalculate with smaller step size of the given factor fails to statisfy the given target accuracy,
        return step(r, t, factor * h)

h = T/250000  #defining initial step size
tpoints = [] #creating empty array of time values
x_values2 = [] #creating empty array of x values
y_values2 = [] #creating empty array of y values
r = np.array([x0, v_x0, y0, v_y0], float)  #creating an array which stores initial conditions
t = 0 #initializing t

while t < T:
    tpoints.append(t)
    x_values2.append(r[0])
    y_values2.append(r[2])
    delta_r, h, t = step(r, t, h)
    r += delta_r

plt.figure(2) #creating second figure
plt.plot(array(x_values2), array(y_values2)) #plotting x and y values generated with adaptive step size
plt.title("Dynamics of Comet in Orbit with Target Accuracy of 1 km/yr") #naming plot
plt.xlabel("x (km)") #labelling axes
plt.ylabel("y (km)")
plt.savefig('orbitaldynamicsdelta.png') #saving figure with given title
plt.show()

# d) plotting the trajectory of the comet with dots showing the position at each Runge-Kutta step around a single orbit
plt.figure(3) #creating third figure
plt.plot(array(x_values), array(y_values), color="pink") #plotting x and y values generated with adaptive step size
plt.plot(array(x_values2)[::20], array(y_values2[::20]), "ro") #plotting dots
plt.title("Dynamics of Comet in Orbit with Dots at Each Runge-Kutta Step") #naming plot
plt.xlabel("x (km)") #labelling axes
plt.ylabel("y (km)")
plt.savefig('orbitaldynamicsdots.png') #saving figure with given title
plt.show()