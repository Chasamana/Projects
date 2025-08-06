from scipy import integrate

import numpy as np

import matplotlib.pyplot as plt

def Euler(f, y0, t, h):

    """This function applies the Euler method to solve differential equations dy/dx = f(y, x).

    The inputs include the equation's right-hand side, the initial value of y, an array specifying solution points, and the step size h."""

   

    # Compute h as the difference between consecutive points in t

    h = t[1] - t[0]

   

    # Define array to store the solution

    y = np.zeros((len(t), len(y0)))

   

    # Set the initial condition

    y[0] = y0

   

    # Loop over elements of t up to the second last

    for i in range(len(t) - 1):



        # Compute each element based on the formula for the Euler method

        y[i + 1] = y[i] + h * f(y[i], t[i])

   

    # Return the transpose of y for compatibility with plotting

    return y.T



def Heun(f, y0, t, h):

    """This function applies the Heun's method to solve differential equations dy/dx = f(y, x).

    The inputs include the equation's right-hand side, the initial value of y, an array specifying solution points, and the step size h."""

   

   

    # Compute h as the difference between consecutive points in t

    h = t[1] - t[0]

   

    # Define array to store the solution

    y = np.zeros((len(t), len(y0)))

   

    # Set the initial condition

    y[0] = y0

   

    # Loop over elements of t up to the second last

    for i in range(len(t) - 1):

       

        k1 = f(y[i], t[i])

       

        k2 = f(y[i]+k1*h, t[i]+h)

       

        # Compute each element based on the formula for Heun's method

        y[i + 1] = y[i] + h * (0.5*(k1+k2))

   

    # Return the transpose of y for compatibility with plotting

    return y.T

#Initialise the step size

h = 0.1



# External force function (a unit impulse)

F = lambda t: t <= 0.1



# System dynamics as a lambda function for two masses

massSpringDamperTwoMasses = lambda y, t: np.array([y[1], y[3] - y[1] - 4*y[0] -2, y[3], y[1] - y[3] + F(t)])

 

# Initial conditions (initial displacements and velocities) for mass 1 and mass 2

y0 = np.array([0, 0, 0, 0])  



# Initialise Heun and Euler error values to infinity

errorHeun = float('inf')

errorEuler = float('inf')



# Initialize maximum step size to infinty

maxStepSize = float('inf')



# Loop containing Heun, Euler and Odeint solutions. Loop continues until both errors are less than 0.01

while errorHeun >= 0.01 or errorEuler >= 0.01:

   

    # Define array with points where the solution should be obtained - this line of code defines the end time of the simulation

    t = np.linspace(0, 17.5, int(17.5/ h))

   

    # Solution using Euler method with lambda function

    yEuler = Euler(massSpringDamperTwoMasses, y0, t, h)



    # Solution using Heun's method with lambda function

    yHeun = Heun(massSpringDamperTwoMasses , y0, t, h)

   

    #Solution using Odeint with lambda function

    yOdeint = integrate.odeint(massSpringDamperTwoMasses , y0, t)

   

    # Accessing the last values for all state variables

    lastValuesEuler = yEuler[:, -1]

    lastValuesHeun = yHeun[:, -1]

    lastValuesOdeint = yOdeint[-1]

   

    # Accessing displacement of mass 2 and calculating errors. These values have to be absolute or else the loop won't work properly.

    errorHeun = np.abs((lastValuesOdeint[2] - lastValuesHeun[2]) / lastValuesOdeint[2])

    errorEuler = np.abs((lastValuesOdeint[2] - lastValuesEuler[2]) / lastValuesOdeint[2])

   

    # This variable ensures that the step size that satisfies the condition of the loop is stored before it is decreased in the next step

    maxStepSize = h

   

    # Step size is decreased in case error is not small enough

    h = h * 0.99

   

# Turning error into percentage errors

errorHeunPercentage = 100*errorHeun

errorEulerPercentage = 100*errorEuler


print('Final Value of Odeint function:', lastValuesOdeint[2])

print('Final Value of Euler Method:', lastValuesEuler[2])

print('Final Value of Heun Method:', lastValuesHeun[2])

print("Final Step Size", maxStepSize)

print("Error Heun:", errorHeun)

print("Error Euler:", errorEuler)

print("Error Heun Percentage:", errorHeunPercentage)

print("Error Euler Percentage:", errorEulerPercentage)

# Plot the results

plt.plot(t, yOdeint[:, 2], label = 'ODEINT - Displacement Mass 2')

plt.plot(t, yEuler[2], label = 'Euler - Displacement Mass 2', linestyle = '--')

plt.plot(t, yHeun[2], label = 'Heun - Displacement Mass 2', linestyle = '--')

plt.xlabel('Time')

plt.ylabel('Displacement')

plt.legend()

plt.show()

