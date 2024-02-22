'''Altitude of a satellite'''

G=6.67*10**(-11) # defining gravitational constant
M=5.97*10**24 # defining mass of earth
R=6371*1000 # defining radius of earth and converting it from km to m
def satellite_altitude(T): #defining function that takes orbital period T, in minutes, and calculates altitude h
    h=((G*M*(T**2))/(4*(np.pi**2)))**(1/3)-R # formula for altitude of satellite orbiting earth
    return(h) #returns h

T=float(input("What is the orbital period (T) of the satellite, in minutes?"))*60 # prompts user with string to input the orbital period in minutes, converts the input from minutes to seconds for the calculations, and stores that number in the variable T
a=satellite_altitude(T) # uses the variable T defined by the user's input to call the function to calculate the value of h for that T; stores the result in the variable a
print("The altitude of the satellite, in km, is", a/1000) # prints results, dividing the output of the function by 1000 to convert from m to km


'''Planetary orbits'''
# Program that asks the user to enter the distance to the Sun and velocity
# at perihelion, then calculates and prints certain quantities
G=6.67*10**(-11) #defining gravitational constant
M=1.9891*10**30 # defining mass of orbiting planet

def planetary_orbit(v1, l1): # defining a function that takes the velocity, v1, and the distance to
                            # the Sun, l1, of a planet at perihelion
    v2=((2*G*M)/(l1*v1))-v1 # calculation of v2 using the root of the quadratic equation given in part a of the problem
    l2=(l1*v1)/v2 # calculation of distance from Sun to planet at aphelion
    a=(1/2)*(l1+l2) # calculation of semi major axis
    b=math.sqrt(l1*l2) # calculation of semi minor axis
    T=(2*(math.pi)*a*b)/(l1*v1) # calculation of orbital period using l1 and v1
    e=(l2-l1)/(l2+l1) # calculation of orbital eccentricity
    return(l2,v2,T,e) #returns tuple with values for l2,v2,T,e

v=float(input("What is the velocity of the planet at perihelion?"))  # prompts user with string to input velocity at perihelion v1, and stores that number in the variable v
l=float(input("What is the distance from the planet to the Sun at perihelion?")) # prompts user with string to input distance to sun at perihelion l1, and stores that number in the variable l
r=planetary_orbit(v,l) # uses the variables v and l defined by the user's input to call the function to calculate the values of l2,v2,T,e; stores the results in the a tuple assigned to the variable r
apheliondist, aphelionvel, orbitalperiod, eccentricity = r #breaks apart tuple that will be returned by the
                                                    # above line when the function is called
print("The distance from the planet to the Sun at aphelion is %d; the velocity of the planet at aphelion is %d; the orbital period is %d; and the orbital eccentricity is %d" %(apheliondist, aphelionvel, orbitalperiod, eccentricity)) # prints results