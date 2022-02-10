#!/usr/bin/env python
# coding: utf-8

# ### Harmonic Oscillator Energy By Monte Carlo Approach

# <span style="color:blue"> **Monte Carlo** </span> is an Algorithm that divides a function into two pieces, $f$, and $P$, as a probability function. In each step, select a random variable, such as $x$, and insert it into $P$, then assess its relationship divided by the probability of the preceding variable. If the result is more than one, this variable is added to the data list. In addition, if it is greater than another random integer, $x$ is added to the list. We now have a sequence of data in the $x$ list, which is fed into the main function as the energy function,<span style="color:blue"> which computes the sum of outputs divided by data numbers</span>. As a result, we have the system's energy.

# In[3]:


from random import random
import numpy as np
import matplotlib.pyplot as plt

x_list = []
x = random()
delta = 1
l = np.arange(0.1,1,0.05)
e_list = []
v_list = []
def p(x,l):
    return np.exp(-2*l*(x**2))

N = 10000
for j in l:
    x_list = []
    for i in range(N):
        xtrial = x+(delta*(2*random()-1)) #test sample
        w = p(xtrial,j)/p(x,j)
        
        if w > 1:
            x = xtrial
            x_list.append(x)
        else:
            r = random()
            if r <= w:
                x = xtrial 
                x_list.append(x)
            else:
                x_list.append(x)
                
    E,E1=0,0
    for y in x_list:
        E2 = j + ((y**2)*((1/2)-2*(j**2)))
        E = E + E2
        E1 = E1 + (E2**2)
        
    e_list.append(E/N)
    v_list.append(((E1/N)-((E/N)**2))) 
a = e_list.index(min(e_list)) # Minimum value of Lambda
print("lambda:",l[a])


# In[4]:


print(min(e_list))
plt.plot(l,e_list)
plt.ylabel("E")
plt.xlabel("lambda")
plt.show()
plt.plot(l,v_list)
plt.xlabel("lambda")
plt.ylabel("sinma^2")
plt.show()


# In[ ]:




