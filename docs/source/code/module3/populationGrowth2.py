import numpy as np
import matplotlib.pyplot as plt

k = 0.693147 #growth rate (found in example from Module 2)
P = 100
dt = 0.5
t = 0
totalTime = 5 #hours

time = [t]
population = [P]
populationSolved = [P] #adding a third storage list for the solved equation values

while t < totalTime:
    dP = k*(P)*dt
    P += dP
    population.append(int(P))
    Pt = 100*(2**(t+dt)) #adding in a line for calculating the exact result at each step
    populationSolved.append(int(Pt))
    t += dt
    time.append(t)

print("Estimated Value:", population[-1])
print("Solved Value:", populationSolved[-1])
print("Final Time:", time[-1])


plt.plot(time, population)
plt.plot(time, populationSolved)
plt.title("dt = " + str(dt))
plt.xlabel("time (hours)")
plt.ylabel("population")
plt.legend(["estimate", "solved"], loc="upper left") #adding a legend to the plot to distinguish the two lines
plt.show()