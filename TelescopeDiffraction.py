import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from numpy import linspace,sin,cos
import matplotlib.patches as mpatches
from pylab import imshow,show

'''Modeling the diffraction pattern of light in a telescope'''

#function that returns the Bessel function Jm(x) using Simpson's rule
#then making a plot of the Bessel functions for m=0, 1, and 2
def J(m,x): #defining function
    def f(theta): #defining nested function
        f=np.cos(m*theta-x*theta) #Defining the integrand of the bessel function as a nested function, to be used for Simpson's rule
        return(f) #returns the value of the integrand when f is called in Simpson's rule below
    a=0 #storing lower limit of integration to be used in Simpson's rule
    b=np.pi #storing upper limit of integration to be used in Simpson's rule
    N=1000
    h=(b-a)/1000 #defining the width of slices to be a 1000th of the total range of integration, making the slices smaller and thus increasing the accuracy of the estimate
    odd_terms=[] #defining empty array for the odd terms to be added to each time through the for loop
    for k in range(1,N,2): #using for loop to store, in the variable "odd", the term in Simpson's rule with a sum over odd values
        odd=f(a+k*h) #calls the function f(theta) that stores the integrand, using it in the sum
        odd_terms.append(odd) #adds most recently calculated value of f(x) to the empty array of odd terms
    odd_sum=sum(odd_terms) #calculating the sum of the now-populated array of odd terms calculated using the for loop above
    even_terms = []  #defining empty array for the even terms to be added to each time through the for loop
    for k in range(2,N,2):  #using for loop to store, in the variable "even", the term in Simpson's rule with a sum over even values
        even=f(a+k*h)  #calls the function f(theta) that stores the integrand, using it in the sum
        even_terms.append(even)  #adds most recently calculated value of f(x) to the empty array of even terms
    even_sum=sum(even_terms)  #calculating the sum of the now-populated array of even terms calculated using the for loop above
    J_s=(1/3)*h*(f(a)+f(b)+4*odd_sum+2*even_sum) #calculating Simpson's rule where we integrate from x=a to x=b for given function f(x) (f(x)=J_m),
    #where h=width of slices we're fitting quadratic to and integrating over.
    return((1/np.pi)*J_s)

x_range=linspace(0,20,1000) #defining x-range to run Bessel functions over using linspace(); will create array for domain of 1000 x-values between 0 and 2
J0=J(0,x_range) #calling the Bessel function for m=0, x=range defined above using linspace()
J1=J(1,x_range) #calling the Bessel function for m=1, x=range defined above using linspace()
J2=J(2,x_range) #calling the Bessel function for m=2, x=range defined above using linspace()

plt.figure(1) #plots the Bessel function as the first figure
plt.plot(x_range,J0,color="red")
plt.plot(x_range,J1,color="orange")
plt.plot(x_range,J2,color="pink")
plt.title('Bessel Function for m=0,1,2') #gives plot a title
plt.savefig("bessel.png") #saves plot as a png image to whatever directory the user is currently in
plt.show() #shows the plot/figure and then clears it, so that data doesn't all get stacked

#writing program that makes a density plot of the intensity of the circular diffraction pattern of a point light source with λ = 500 nm, in a square region of the focal plane
#first, defining constants given
lambda_=500/(1*10**9) #defining lambda, and converting from nm to m for calculations
x_=linspace(-(1/(1*10**6)),(1/(1*10**6)),100) #creating arrays of x and y coordinates using linspace to represent the plane in which diffraction pattern is observed
y_=linspace(-(1/(1*10**6)),(1/(1*10**6)),100)
X,Y=np.meshgrid(x_,y_)
r_=np.sqrt((X**2)+(Y**2)) #defining range of r values using linspace() and converting upper limit from micrometers to m
#r_= np.delete(r_, np.where(r_==0)) # excluding the value 0 from the range of r values calculated using linspace, because "r_" might include the value 0, which would cause a divide by zero error when calculating "k*r_".
k=(2*np.pi)/lambda_ #calculating k using formula given
I=((J(1,k*(r_)))/(k*r_))**2 #calculating intensity of light in the diffraction pattern using formula given, which will prodce an array of I values dimensionally identical to the array of r values previously defined
imshow(I)
show()