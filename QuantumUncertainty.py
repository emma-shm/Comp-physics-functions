import matplotlib.pyplot as plt
import math
from math import factorial
import numpy as np
from numpy import linspace,sin,cos,ones,copy,tan,pi
import warnings
warnings.filterwarnings("ignore")

'''Quantum uncertainty in the harmonic oscillator'''


# Plotting harmonic oscillator wavefunction
def H(n,x): #defining function for the Hermite polynomial, then using that to write a function for psi --- the wavefunction of the nth energy level of the one-dimensional quantum harmonic oscillator
    if n==0:
        return(1)
    elif n==1:
        return(2*x)
    else:
        return(((2*x)*(H(n-1,x)))-(2*(n-1))*(H(n-2,x)))

def psi(n,x): #defining wavefunction of the nth energy level of the one-dimensional quantum harmonic oscillator using function for Hermite polynomial defined above
    psi=(1/(np.sqrt((2**n)*(math.factorial(n))*(np.sqrt(np.pi)))))*((np.e)**(-(x**2)/2))*(H(n,x))
    return(psi)

x_range_=linspace(-4,4,1000) #defining x-range to run wavefunction over using linspace(); will create array for domain of 1000 x-values between -4 and 4
psi_0=psi(0,x_range_) #calling the wavefunction for n=0, x=range defined above using linspace(), storing resulting array of wavefunctions with n=0 calculated using each of the 1000 x-values in the range
psi_1=psi(1,x_range_) #calling the wavefunction for n=1, x=range defined above using linspace(), storing resulting array of wavefunctions with n=1 calculated using each of the 1000 x-values in the range
psi_2=psi(2,x_range_) #calling the wavefunction for n=2, x=range defined above using linspace(), storing resulting array of wavefunctions with n=2 calculated using each of the 1000 x-values in the range
psi_3=psi(3,x_range_) #calling the wavefunction for n=3, x=range defined above using linspace(), storing resulting array of wavefunctions with n=3 calculated using each of the 1000 x-values in the range

plt.figure(2) #plotting figure 2, which will include the wavefunctions calculated for n=0,1,2,3
#plt.plot(x_range_,psi_0,color="red")
plt.plot(x_range_,psi_1,color="orange")
plt.plot(x_range_,psi_2,color="pink")
plt.plot(x_range_,psi_3,color="purple")
plt.title('Harmonic Oscillator Wavefunctions') #gives plot a title
plt.savefig("wavefunctions.png") #saves plot as a png image to whatever directory the user is currently in
plt.show() #shows the plot/figure and then clears it, so that data doesn't all get stacked

#Part b) Making separate plot of the wavefunction for n = 30 from x = −10 to x = 10.
x_range__=linspace(-10,10,1000)
psi_30=psi(30,x_range_) #calling the wavefunction for n=3, x=range defined above using linspace(), storing resulting array of wavefunctions with n=3 calculated using each of the 1000 x-values in the range

plt.figure(3) #plotting figure 2, which will include the wavefunctions calculated for n=0,1,2,3
plt.plot(x_range_,psi_30,color="blue")
plt.title('Harmonic Oscillator Wavefunction for n=30') #gives plot a title
plt.savefig("wavefunctions_30.png") #saves plot as a png image to whatever directory the user is currently in
plt.show() #shows the plot/figure and then clears it, so that data doesn't all get stacked



# Writing a program that evaluates the integral in the root-mean-square using Gaussian quadrature, in order to determine quantum uncertainty in particle position
from gaussxw import gaussxw

n=5 #defining n value you want to find uncertainty on
def f(x): #defining a function that calculates the integrand for psi
    int_=(x**2)*(abs(psi(n,x)))**2
    return(int_)

N=10000 #number of sample points
a=-1000 #choosing very large, finite limits of integration to avoid getting output of nan
b=1000

x,w = gaussxw(N) #calculating sample points and weights
xp = 0.5*(b-a)*x + 0.5*(b+a) #mapping sample points to required domain of integration
wp = 0.5*(b-a)*w #mapping weights to required domain of integration

s = 0.0 #initializing s at zero
for k in range(N): #Performing the integration
    s += wp[k]*f(xp[k]) #each time through for loop, adds to the previous term the newly calculated value of function at the given sample point, plus the corresponding weight
print(np.sqrt(s))