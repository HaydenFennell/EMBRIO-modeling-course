import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import scipy.sparse.linalg
from scipy import sparse
from IPython.display import Image
import math as ma

# example for increasing number of Taylor series terms to show how it affects accuracy
t = np.linspace(-4*np.pi,4*np.pi,1000)
f1 = t 
f2 = t - (t**3)/(ma.factorial(3))
f3 = t - (t**3)/(ma.factorial(3)) + (t**5)/(ma.factorial(5))
f4 = t - (t**3)/(ma.factorial(3)) + (t**5)/(ma.factorial(5)) - (t**7)/(ma.factorial(7))
f5 = t - (t**3)/(ma.factorial(3)) + (t**5)/(ma.factorial(5)) - (t**7)/(ma.factorial(7)) + (t**9)/(ma.factorial(9))
f6 = t - (t**3)/(ma.factorial(3)) + (t**5)/(ma.factorial(5)) - (t**7)/(ma.factorial(7)) + (t**9)/(ma.factorial(9)) - (t**11)/(ma.factorial(11))

plt.plot(t,np.sin(t),'k')
plt.plot(t,f1,'r--')
plt.plot(t,f2,color='darkorange',linestyle='--')
plt.plot(t,f3,'y--')
plt.plot(t,f4,'g--')
plt.plot(t,f5,'b--')
plt.plot(t,f6,color='violet',linestyle='--')
plt.plot(0,0,'ok')
plt.ylim(-1.5,1.5)