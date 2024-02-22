import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from numpy import linspace,sin,cos
import matplotlib.patches as mpatches
from pylab import imshow,show

'''Heat capacity of a solid'''

#Function that calculates the heat capacity of a solid at a given temperature T, using Debye's theory of solids
def Cv(T): #defining function that calculates Cv
    V=100/100 #defining the volume given, converting from cm to m
    rho_=6.022*10**28 #defining the number density of aluminum
    theta_D=429 #defining the Debye temperature
    k_B=1.380649*10**-23 #defining Boltzmann's constant
    def f(x):  # defining nested function
        if np.isclose(np.exp(x), 1).all(): #checking if the exponential function of x is close to 1 using the np.isclose function; returns 0 if it's true for all elements in the array, meaning np.exp(x) is very close to 1
            return 0
        else:
            f=((x**4)*((np.e)**x))/((((np.e)**x)-1)**2) #If np.exp(x) is not close to 1 for all elements, function calculates the value of the integrand using the input x and returns the result
            return (f)  #returns the value of the integrand when f is called in Simpson's rule below
    a=0 #storing lower limit of integration to be used in Simpson's rule
    b=theta_D/T #storing upper limit of integration to be used in Simpson's rule
    N=50 #defining N-value for the number of sample points as its given in the question
    h=(b-a)/50 #defining the width of slices to be a 50th of the total range of integration, making the slices smaller and thus increasing the accuracy of the estimate
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
    Cv_integral=(1/3)*h*(f(a)+f(b)+4*odd_sum+2*even_sum) #calculating Simpson's rule where we integrate from x=a to x=b for given function f(x) (f(x)=J_m),
    Cv_final=((9*V*rho_*k_B*((T/theta_D)**3)))*Cv_integral #calculating Cv by multiplying the integral in the heat capacity formula with the coefficient in front of it (coefficient
                                                            # #is calculated using constants defined in first 4 lines of Cv(T) function)
    return(Cv_final)

#Part b)
T_range=np.arange(5,500,5) #using np.arange to define range of temperatures T to run the function over, storing the range in an array
h_c=Cv(T_range) #calling the function for temperatures over the range defined abpve, calculating the heat capacity for each and storing the results in an array of heat capacity
# #values that each correspond to the entries in T_range
plt.figure(3)
plt.plot(T_range,h_c,'.') #plotting the range of temperatures defined above using np.arange VS the array of heat capacity values calculated using each of these temperatures
plt.title('Heat Capacity of Solid Aluminum as a Function of Temperature') #gives plot a title
plt.savefig("heatcapacity.png") #saves plot as a png image to whatever directory the user is currently in
plt.show() #shows the plot/figure and then clears it, so that data doesn't all get stacked