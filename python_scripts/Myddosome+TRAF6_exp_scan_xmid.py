def oligomer (n, s_max, rlc, traf6, p1, p2, p3, p4, p5, p6):
    
    """ Stochastic simulation of binding to a Receptor-Ligand-Complex, with n Receptor-Ligand-Complexes."""  
    # rlc is an array of sites, of length n.
    # p1 and p2 are the probabilities of binding and dissociation.
    # Receptor can be bound (1) or free (0).
    
    # select random Receptor:
    rs = np.random.randint(n)  # generate a random integer between 0 and n
    
    #Define the equations and probabilities for MyD88 growth:
        
    s = rlc [rs] # select that position in the rlc array
    
    # generate random numbers r1 and r2
    r1 = np.random.rand() # generates a random number between 0 and 1; probability of association
    r2 = np.random.rand() # generates a random number between 0 and 1; probability of dissociation
    
    p1_mod = 0.3+s*p1
    
     # modify the site 's' according to p1 and p2
  
    if s > 0 and s <= 5:           # If receptor 's' is bound, and the number of monomers in oligomer is between 1 and 5
        if r1 <= p1_mod:
            rlc [rs] +=1        # increase bound oligomer by one monomer with probability p1
        
        if r2 <= p2:
            rlc [rs] -=1        # decrease bound oligomer by one monomer with probability p2
    
    elif s == 0:         # If receptor 's' is unbound, 
        if r1 <= p1_mod:
            rlc [rs] +=1        # increase bound oligomer by one monomer with probability p1
    
    elif s >= 6 and s < s_max:           # If receptor 's' is bound, and the number of monomers in oligomer is between 6 and 20
        if r1 <= p3:
            rlc [rs] +=1        # increase bound oligomer by one monomer with probability p3
        
        if r2 <= p4:
            rlc [rs] -=1        # decrease bound oligomer by one monomer with probability p4
    
    elif s == s_max:           # If receptor 's' is bound by the maximal oligomer size
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
    
    # select random Receptor:
    rs = np.random.randint(n)  # generate a random integer between 0 and n
    
    #Define the equations and probabilities for TRAF6 growth:
    s = rlc [rs] # select that position in the rlc array
    t = traf6 [rs] #we look at the receptor s and will calculate the probability of TRAF6 association
    
    # generate random numbers r3 and r4
    r3 = np.random.rand() # generates a random number between 0 and 1; probability of association
    r4 = np.random.rand() # generates a random number between 0 and 1; probability of dissociation
    
    if s > 0:
        p5_mod = p5[0]/(1+np.exp((p5[1]-np.log(s))/p5[2])) #logistic growth model: p = Asym/(1+exp((xmid-log(input))/scal)) #-4.156397e-02 + 4.306951e-02*s + -1.327248e-03*s**2 + 5.393454e-06*s**3 + 3.639276e-07*s**4 + -4.279642e-09*s**5
        p6_mod = p6[0]/(1+np.exp((p6[1]-np.log(s))/p6[2])) #-1.755390e-02 + 2.102242e-02*s + 1.072152e-04*s**2 + -3.888753e-05*s**3 + 1.133250e-06*s**4 + -1.010123e-08*s**5
        
    elif s == 0:
        p5_mod = 0
        p6_mod = 0

    #we look at the same receptor as before and calculate the probability of TRAF6 association
    if s >= 6 and t > 0 and t < 10:  # If receptor 's' is bound by more than 6 monomers and the number of TRAF6 associated is greater than 0
        if r3 <= p5_mod:
            traf6 [rs] +=1        # increase bound TRAF6 by one TRAF6-monomer with probability p5
        
        if r4 <= p6_mod:
            traf6 [rs] -=1        # decrease bound TRAF6 by one TRAF6-monomer with probability p6
            
    elif s >= 6 and t == 0:       # If receptor 's' is bound by more than 6 monomers and one TRAF6 is associated
        if r3 <= p5_mod:
            traf6 [rs] +=1        # increase bound TRAF6 by one TRAF6-monomer with probability p5
            
    elif s == 0:                  # If the complex dissappeared, then TRAF6 will also dissociate
        traf6 [rs] = 0
            
    elif t == 10:                 # if receptor is bound by the maximal number of TRAF6
       if r4 <= p6_mod:
           traf6 [rs] -=1         # decrease bound TRAF6 by one TRAF6-monomer with probability p6
    
   # evaluate the Receptor array: count ones and zeros. Gives a single number for each state. 
    TRAF6_neg = 0
    TRAF6_pos = 0
    total_bound_TRAF6 = 0
   
    for site in traf6:
        if site == 0:
            TRAF6_neg +=1            #sum up all the free sites
        
        elif site > 0:          #if site is occupied; sum up bound sites
            TRAF6_pos +=1    
            
    for site in range(len(traf6)):
        total_bound_TRAF6 +=traf6[site]     #sum up the number of monomers in complex
        
    return free, bound, total_bound, rlc, TRAF6_neg, TRAF6_pos, total_bound_TRAF6, traf6  # "rlc" and "traf6" are the new input for the next step. 


def oligomer_time (step, sec, n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, p5, p6):
    
    """ Stochastic simulation of time series for oligomer. """
    
    t_list = np.linspace (0,sec,step)   # define the time vector, sec is the end point to evaluate and step the step size
    free_list = []                      # set up an empty list. Here you will collect the number of free receptors after each iteration. 
    bound_list = []                     # set up an empty list. Here you will collect the number of bound receptors after each iteration.
    total_bound_list = []               # set up an empty list. Here you will collect the number of monomers in complex after each iteration. 

    TRAF6_neg_list = []                 # set up an empty list. Here you will collect the number of complexes without TRAF6 after each iteration. 
    TRAF6_pos_list = []                 # set up an empty list. Here you will collect the number of complexes with TRAF6 after each iteration.
    total_bound_TRAF6_list = []         # set up an empty list. Here you will collect the number of total bound TRAF6 after each iteration. 
    
    size_list = [0 for _ in range(int(s_max))]  #set up a list to sum up the total number (cummulative number at the end of time) of complexes of certain size
    TRAF6_pos_size_list = [0 for _ in range(int(s_max))]



    for i in range (len(t_list)):
        
        free, bound, total_bound, rlc_1, TRAF6_neg, TRAF6_pos, total_bound_TRAF6, traf6_1 = oligomer (n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, p5, p6)  # call the function "oligomer"
        
        rlc_0 = rlc_1                                      # take the array "rlc_1" from "oligomer" and use it as input to the next iteration
        traf6_0 = traf6_1                                  # take the array "traf6_1" from "oligomer" and use it as input to the next iteration
        
        free_list.append(free)                             # add the number of free receptors to the list "free_list"
        bound_list.append(bound)                           # add the number of bound receptors to the list "bound_list"
        total_bound_list.append(total_bound)               # add the number of monomers in complex to the list "total_bound"
        
        TRAF6_neg_list.append(TRAF6_neg)                   # add the number of complexes without TRAF6 to the list "TRAF6_neg_list"
        TRAF6_pos_list.append(TRAF6_pos)                   # add the number of complexes with TRAF6 to the list "TRAF6_pos_list"
        total_bound_TRAF6_list.append(total_bound_TRAF6)   # add the total number of bound TRAF6 to the list "total_bound_TRAF6_list"
        
        for size in range(int(s_max)):
            for site in rlc_0:
                if site == size:
                    size_list[size] +=1                    # sum up the number of complexes with a distinct number of monomers
        
        for size in range(int(s_max)):
            for site in range(len(rlc_0)):
                if rlc_0[site] == size and traf6_0[site] > 1: # traf6_0[site] > 1 we count colocalization when two or more TRAF6 are associated
                    TRAF6_pos_size_list[size] +=1
                        

    return t_list, free_list, bound_list, total_bound_list, TRAF6_neg_list, TRAF6_pos_list, total_bound_TRAF6_list, size_list, TRAF6_pos_size_list


def oligomer_time_repeat_mean (runs, step, sec, n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, p5, p6):
    
    """Runs oligomer_time several times."""
    
    free_all  = []
    bound_all = []
    total_bound_all = []
    
    TRAF6_neg_all = []
    TRAF6_pos_all = []
    total_bound_TRAF6_all = []
    
    pct_traf6_all = []
    
    for j in range (runs):
        rlc_0 = [0 for _ in range(n)]  # alternative: rlc = np.zeros (n) # generates an array of n receptors, all initially unbound.
        traf6_0 = [0 for _ in range(n)] #so we start every run with a new array of 0
        
        t_list, free_list, bound_list, total_bound_list, TRAF6_neg_list, TRAF6_pos_list, total_bound_TRAF6_list, size_list, TRAF6_pos_size_list = oligomer_time(step, sec, n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, p5, p6) #call the function "oligomer_time"
        
        free_all.append(free_list)
        bound_all.append(bound_list)
        total_bound_all.append(total_bound_list)
        
        TRAF6_neg_all.append(TRAF6_neg_list)
        TRAF6_pos_all.append(TRAF6_pos_list)
        total_bound_TRAF6_all.append(total_bound_TRAF6_list)
        
        pct_traf6_list = []                                      # create an empty list to calculate the percentage of complexes recruiting TRAF6 for each run
        for i in range(len(size_list)):
            if size_list[i] == 0:
                quotient = np.nan                                # append NaN value if denominator is zero
            else:
                quotient = TRAF6_pos_size_list[i] / size_list[i] # divide the integers at the same position
            pct_traf6_list.append(quotient)                      # add the quotient to the result list
        
        pct_traf6_all.append(pct_traf6_list)                     # append the result for each run to the pct_traf6_all list
    
    mean_bound = np.mean(bound_all, axis=0)  # calculates the mean of all runs (bound)
    mean_TRAF6_pos = np.mean(TRAF6_pos_all, axis=0)
    mean_pct_traf6 = np.mean(pct_traf6_all, axis=0)

    return mean_bound, mean_TRAF6_pos, mean_pct_traf6, t_list

def oligomer_xmid_repeat_scan (runs, step, sec, n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, xmid_list, p6):
    
    """Runs oligomer_time_repeat with different time values"""
    
    pct_traf6_xmid = [] # pct of complexes with traf6 based on complex size after each simulation run with different p3
    TRAF6_pos_time = [] # complexes with TRAF6 over time
    bound_all_time = [] # bound receptors over time 
    
    for xmid in xmid_list:
        p5 = [0.4087897, xmid, 0.7862437] # p5 values for growth over 3 Frames #[0.4096933, 1.9628704, 0.4627068] #[Asym, xmid, scal] #p = Asym/(1+exp((xmid-log(input))/scal))
        mean_bound, mean_TRAF6_pos, mean_pct_traf6, t_list = oligomer_time_repeat_mean (runs, step, sec, n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, p5, p6)  # call the function "oligomer_time_repeat"
        TRAF6_pos_time.append(mean_TRAF6_pos)
        bound_all_time.append(mean_bound)
        pct_traf6_xmid.append(mean_pct_traf6)
        
    return pct_traf6_xmid, TRAF6_pos_time, bound_all_time, t_list

#main program:
#------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# set binding sites.
n = 50

# set time
sec = 30 #200
sub = 1000
step = int(sec*sub)

#set maximal oligomer size
s_max = 50
s_list = np.linspace (0, s_max, s_max)   # define the list for all sizes of oligomers

#set p3, p4
p3 = 0.6 
p4 = 0.4 #0.05

# set p1, p2
p1 = (p3-0.3)/5
p2 = 0.4 # probability of dissociation. #0.2

#Probabilities of TRAF6 association and dissociation
xmid_list = [np.log(_) for _ in np.linspace(2,50,13)]
# p5 = [0.4087897, 2.6175046, 0.7862437] # p5 values for growth over 3 Frames #[0.4096933, 1.9628704, 0.4627068] #[Asym, xmid, scal] #p = Asym/(1+exp((xmid-log(input))/scal))
p6 = [0.5326890, 3.4659471, 1.0717059] # p6 values for growth over 3 Frames #[0.4364759, 2.5783074, 0.6698108]

# Run the parameter model:
#------------------------------------------------------------------------
runs= 100
rlc_0 = [0 for _ in range(n)]   # alternative: rlc = np.zeros (n) # generates an array of n binding sites, all initially unbound.
traf6_0 = [0 for _ in range(n)]

#run the parameter scan
pct_traf6_xmid, TRAF6_pos_time, bound_all_time, t_list = oligomer_xmid_repeat_scan (runs, step, sec, n, s_max, rlc_0, traf6_0, p1, p2, p3, p4, xmid_list, p6)

# #export the percentage of puncta recruiting TRAF6 at distinct sizes as a csv file
# # Create an empty DataFrame with the desired column names
# df = pd.DataFrame(columns=['PCT_RECRUITMENT', 'SIZE', 'XMID'])

# # Loop through the sublists and add their data to the DataFrame        
# dataframes = []
# for i, sublist in enumerate(pct_traf6_xmid):
#     XMID = xmid_list[i]  # The index of the sublist
#     rows = [{'PCT_RECRUITMENT': value, 'SIZE': j, 'XMID': XMID} for j, value in enumerate(sublist)]
#     dataframes.append(pd.DataFrame(rows))
# df = pd.concat(dataframes, ignore_index=True)
        
# # # Save the DataFrame as a compressed CSV file
# file_path = '//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/python_scripts/Modeling_Mauriz/Figures/df_size2_3frames_scan_xmid.csv.gz'
# df.to_csv(file_path, compression='gzip', index=False)


plt.subplot(2,1,1)
for i in range(len(xmid_list)):
    plt.plot(s_list, pct_traf6_xmid[i])
#plt.plot(s_list, pct_traf6_time[i-1], label='pct of complexes with TRAF6')
plt.title('Parameter scan xmid: pct of complexes with TRAF6')
plt.xlabel('size of MyD88 complex')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(0,1.1)

plt.subplot(2,1,2)
for i in range(len(xmid_list)):
    plt.plot(t_list, TRAF6_pos_time[i])
plt.title('xmid; complexes with TRAF6')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(0,52)

plt.subplots_adjust(right=0.8, hspace=0.5)
plt.show()

