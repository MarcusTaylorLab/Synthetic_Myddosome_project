def oligomer (n, rlc, p1, p2, p3, p4):
    
    """ Stochastic simulation of binding to a Receptor-Ligand-Complex, with n Receptor-Ligand-Complexes."""  
    # rlc is an array of sites, of length n.
    # p1 and p2 are the probabilities of binding and dissociation.
    # Receptor can be bound (1) or free (0).
    
    # select random Receptor:
    rs = np.random.randint(n)  # generate a random integer between 0 and n
    s = rlc [rs] # select that position in the rlc array
    
    # generate random numbers r1 and r2
    
    r1 = np.random.rand() # generates a random number between 0 and 1; probability of association
    r2 = np.random.rand() # generates a random number between 0 and 1; probability of dissociation
    
     # modify the site 's' according to p1 and p2
  
    if s > 0 and s <= 5:           # If receptor 's' is bound, and the number of monomers in oligomer is between 1 and 5
        if r1 <= p1:
            rlc [rs] +=1        # increase bound oligomer by one monomer with probability p1
        
        if r2 <= p2:
            rlc [rs] -=1        # decrease bound oligomer by one monomer with probability p2
    
    elif s == 0:         # If receptor 's' is unbound, 
        if r1 <= p1:
            rlc [rs] +=1        # increase bound oligomer by one monomer with probability p1
    
    elif s >= 6 and s < 20:           # If receptor 's' is bound, and the number of monomers in oligomer is between 6 and 20
        if r1 <= p3:
            rlc [rs] +=1        # increase bound oligomer by one monomer with probability p3
        
        if r2 <= p4:
            rlc [rs] -=1        # decrease bound oligomer by one monomer with probability p4
    
    elif s == 20:           # If receptor 's' is bound by the maximal oligomer size
        if r2 <= p4:
            rlc [rs] -=1       # decrease bound oligomer by one monomer with probability p4

   # evaluate the Receptor array: count ones and zeros. Gives a single number for each state. 
    free = 0
    bound = 0
    total_bound = 0
   
    for site in rlc:
        if site == 0:
            free +=1            #sum up all the free sites
        
        elif site > 0:          #if site is occupied; sum up bound sites
            bound +=1    
            
    for site in range(len(rlc)):
        total_bound +=rlc[site]     #sum up the number of monomers in complex
        
    return free, bound, total_bound, rlc  # "rlc" is the new input for the next step. 


def oligomer_time (step, sec, n, rlc_0, p1, p2, p3, p4):
    
    """ Stochastic simulation of time series for oligomer. """
    
    t_list = np.linspace (0,sec,step)   # define the time vector, sec is the end point to evaluate and step the step size
    free_list = []                      # set up an empty list. Here you will collect the number of free receptors after each iteration. 
    bound_list = []                     # set up an empty list. Here you will collect the number of bound receptors after each iteration.
    total_bound_list = []               # set up an empty list. Here you will collect the number of monomers in complex after each iteration. 


    for i in range (len(t_list)):
        
        free, bound, total_bound, rlc_1 = oligomer (n, rlc_0,p1, p2, p3, p4)  # call the function "oligomer"
        rlc_0 = rlc_1                                      # take the array "rlc_1" from "oligomer" and use it as input to the next iteration
        free_list.append(free)                             # add the number of free receptors to the list "free_list"
        bound_list.append(bound)                           # add the number of bound receptors to the list "bound_list"
        total_bound_list.append(total_bound)               # add the number of of monomers in complex to the list "total_bound"
 
    
    return t_list, free_list, bound_list, total_bound_list


def oligomer_time_repeat (runs,step,sec,n, rlc_0, p1, p2, p3, p4):
    
    """Runs oligomer_time several times."""
    
    free_all  = []
    bound_all = []
    total_bound_all = []
    
    for j in range (runs):
        rlc_0 = [0 for _ in range(n)]  # alternative: rlc = np.zeros (n) # generates an array of n receptors, all initially unbound.
        t_list, free_list, bound_list, total_bound_list = oligomer_time(step,sec,n, rlc_0, p1, p2, p3, p4)
        free_all.append(free_list)
        bound_all.append(bound_list)
        total_bound_all.append(total_bound_list)
    return bound_all, free_all, total_bound_all, t_list
    

#main program:
#------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# set binding sites.
n = 20
# set time
sec = 15
sub = 1000
step = sec * sub

# set p1, p2 
p1 = 0.2 # probability of binding.
p2 = 0.15 # probability of dissociation.

p3 = p1
p4 = p2

# Run the stochastic model:
#------------------------------------------------------------------------
runs= 100
rlc_0 = [0 for _ in range(n)]   # alternative: rlc = np.zeros (n) # generates a oligomer of n binding sites, all initially unbound.

bound_all, free_all, total_bound_all, t_list = oligomer_time_repeat (runs,step,sec,n, rlc_0, p1, p2, p3, p4)  # call the function "oligomer_time_repeat"

#calculate mean size of MyD88 complexes
mean_complex_size = total_bound_all
for sublist in mean_complex_size:
    for i in range(len(sublist)):
        sublist[i] = sublist[i]/ n

mean_bound = np.mean(bound_all, axis=0)  # calculates the mean of all runs (bound)
mean_free = np.mean(free_all, axis=0)    # calculates the mean of all runs (free)
mean_total_bound = np.mean(total_bound_all, axis=0)    # calculates the mean of all runs (total_bound)
mean_complex_size_allruns = np.mean(mean_complex_size, axis=0)

# Now plot the data:
#------------------------------------------------------------------------

plt.subplot(3,1,1)
plt.plot(t_list, bound_all[0],'k',label='bound receptors')
plt.plot(t_list, free_all[0],'r',label='free receptors')
plt.plot(t_list, mean_complex_size[0],'g',label='avg complex size')
plt.title('Stochastic model: single cell')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(0,22)

plt.subplot(3,1,2)
plt.plot(t_list, mean_bound,'k', label='bound receptors')
plt.plot(t_list, mean_free,'r', label='free receptors')
plt.plot(t_list, mean_complex_size_allruns,'g', label='avg complex size')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Stochastic model: mean of {} cells'.format(runs))
plt.ylim(0,22)

plt.subplot(3,1,3)
for i in range(runs-2):
    plt.plot(t_list, mean_complex_size[i])
plt.plot(t_list, mean_complex_size[runs-1], label='avg complex size')
plt.title('Stochastic model: {} cells: avg complex size'.format(runs))
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(0,22)

plt.subplots_adjust(right=0.8, hspace=0.8)
#plt.show()

# save the plot as a PNG file
filename = '//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/python_scripts/Modeling_Mauriz/Figures/MyD88_growth.png'  # replace with your desired file path
plt.savefig(filename, dpi = 300, transparent = True, bbox_inches = 'tight')
