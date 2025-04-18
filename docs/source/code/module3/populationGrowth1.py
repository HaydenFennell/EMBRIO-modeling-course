import math
import numpy as np
import matplotlib.pyplot as plt

k = 0.693147
P = 100
dt = 0.5
t = 0
totalTime = 5 #hours

time = [t]
population = [P]

while t < totalTime:
    dP = k*(P)*dt #calculate our change in P over the timestep dt
    P += dP #add calculated change to running value
    population.append(int(P)) #add new population value to storage list
    t += dt #increment time by time step
    time.append(t) #add new time value to storage list

#print final value of storage arrays to check results
print("Estimated Value:", population[-1])
print("Final Time:", time[-1])

#Plot the data
plt.plot(time, population)
plt.title("dt = " + str(dt))
plt.xlabel("time (hours)")
plt.ylabel("population")
plt.show()