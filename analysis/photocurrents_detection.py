'''
Photocurrent analysis from ABF files. Designed to run in a Jupyter Notebook.

Written by A. Prominski
'''

import pyabf
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid

def get_data(filename):
    '''
    Fetches data from ABF files from all available sweeps
    '''
    abf = pyabf.ABF(filename)
    x = []
    I = []
    V = []
    analog =[]
    for sweep_number in abf.sweepList:   
        abf.setSweep(sweep_number, channel=0)
        x.append(abf.sweepX)
        I.append(abf.sweepY)
        abf.setSweep(sweep_number, channel=1)
        V.append(abf.sweepY)
        abf.setSweep(sweep_number, channel=2)
        analog.append(abf.sweepY)
    if len (x) == 1:
        return x[0], I[0], V[0], analog[0]   
    
    return x, I, V, analog



def find_all_pulses(analog,V):
    '''
    Finds indices of analog/TTL and voltage pulses in a clamp setup:
    TTL_up - analog/TTL pulse up
    TTL_down - analog/TTL pulse down
    TTL_level - maximum pulse level
    V_up - voltage pulse up
    V_down - voltage pulse down
    '''

    TTL_up = []
    TTL_down = []
    TTL_level = []
    V_up = []
    V_down = []

    ### FLAGS do not change
    tmp = 0
    UP = 1
    DOWN = 0
    direction = UP
    ###
    
    for n,_ in enumerate(analog):     
        if direction == UP and analog[n] - analog[tmp] > 0.1:
            TTL_up.append(n)
            TTL_level.append(analog[n])

            tmp = n
            direction = DOWN
            
        elif direction == DOWN and analog[n] - analog[tmp] > 0.1:
            tmp = n
            TTL_level.pop()
            TTL_level.append(analog[n])
        
        elif direction == DOWN and analog[n] - analog[tmp] < -0.1:
            TTL_down.append(n)

            tmp = n
            direction = UP
            
        elif direction == UP and analog[n] - analog[tmp] < -0.1:
            tmp = n
            
    ### FLAGS do not change
    tmp = 0
    UP = 1
    DOWN = 0
    direction = DOWN
    ### 
    
    for n,_ in enumerate(V[:-10]):
        tmpval = np.mean(V[tmp:tmp+10])
        currentval = np.mean(V[n:n+10])
        if direction == UP and currentval - tmpval > 0.4:
            V_up.append(n)

            tmp = n
            direction = DOWN

        elif direction == DOWN and currentval - tmpval > 0.4:
            tmp = n

        elif direction == DOWN and currentval - tmpval < -0.4:
            V_down.append(n)
            tmp = n
            direction = UP

        elif direction == UP and currentval - tmpval < -0.4:
            tmp = n
    
    ### Error handling
    if len(TTL_up) > 1:
        print('TTL_up too long')
        plt.plot(np.arange(0,len(analog),1), analog)
        plt.show()
    if len(V_up) > 1:
        print('V_up too long')
        plt.plot(np.arange(0,len(V),1), V)
        plt.show()
    
    if len(TTL_up)<1:
        return -1, -1, -1, V_up[0], V_down[0]
    
    if len(TTL_up) == 1:
        return TTL_up[0], TTL_down[0], TTL_level[0], V_up[0], V_down[0]
    
    return TTL_up, TTL_down, TTL_level, V_up, V_down


def getR(I,V,V_down):
    '''
    Function used for calculating pipette resistance
    '''

    ### raise error if voltage down step is far from expected position
    if abs(V_down - 10626) > 10:
        print('R calc error, V down step not where expected. Inspect the data.')
        print(V_down)
        plt.plot(np.arange(0,len(V),1),V)
        plt.show()     
       
    start = V_down - 1000
    end = V_down + 1000
    N=500

    I1 = np.sum(I[start:start+N])
    I2 = np.sum(I[end:end+N])
    V1 = np.sum(V[start:start+N])
    V2 = np.sum(V[end:end+N])
        
    Ip = abs(I2 - I1)*1e-12 ## mV
    Vp = abs(V2 - V1)*1e-3 ## pA
     
    R = Vp/Ip/1e6 ## MOhm

    return R


def photoresponse(t, If, A, k):
    '''
    Function fitting currents:
    If - faradaic current, i.e. exponential decay y-axis offset
    A - pre-exponential constant 
    k - rate constant
    t - time variable
    '''
    return A*np.exp(k*t)+If


def fit_currents(x, I, TTL_up, TTL_down):
    '''
    Model calculating the charges injected into the solution. 
    '''

    t_d = int(10/0.01) ### time before TTL for avg init calculation - 10 ms
    avg_init = np.mean(I[TTL_up-2*t_d:TTL_up-t_d])
    
    current_shifted = I[TTL_up:TTL_down] - avg_init ### holding current offset subtraction
    current_shifted = current_shifted/1000 ### pA to nA
    
    time = x[TTL_up:TTL_down] ### time axis for the pulse position
    time = time - time[0] ### first timepoint subtracted  
  
    n = int(len(current_shifted)*4/5) ### use only last 20% of the transient for calculations                       
    y_data = current_shifted[n:]
    x_data = time[n:]

    ### Intial guesses, magic numbers, but work well in practice
    If_guess = 0
    A_guess = np.amax(y_data)
    k_guess = -0.693/0.01
    p0 = [If_guess, A_guess, k_guess]

    b = ([0, 0, -np.inf],[np.inf, np.inf, np.inf] )
            
    ### Try fitting function
    try:
        popt, pcov, *rest = curve_fit(photoresponse, x_data, y_data, p0 = p0, method = 'trf', jac='3-point', bounds=b)
    except RuntimeWarning:
        print('No solution.')
        popt = [0,0,0]
  
    I_fit = popt[0]
            
    ### Plot data, initial guess and model for visual inspection of the fitting
    y_init = [photoresponse(t, *p0) for t in x_data]
    y_model = [photoresponse(t, *popt) for t in x_data]
    plt.scatter(x_data, y_data, s = 2)
    plt.plot(x_data, y_model, c='red')
    plt.plot(x_data, y_init, c = 'black')
    plt.show()    
    
    ### Calculate charges, through integration and subtraction
    If_final = I_fit
        
    faradaic = np.full_like(time, If_final)
    
    Q_tot = trapezoid(current_shifted, time)
    Qf_final = trapezoid(faradaic, time)
    Qc_final = Q_tot - Qf_final
        
    return Q_tot, Qc_final, Qf_final


def fitting(abf):
    '''
    Analyzes photocurrents in the file
    '''
    print(f'Analyzing: {abf}')
    file = abf
    x, I, V, analog = get_data(file)
    TTL_up, TTL_down, TTL_level, V_up, V_down = find_all_pulses(analog, V)
    tot, c, f = fit_currents(x, I, TTL_up, TTL_down)
    R = getR(I,V,V_down)
    
    return tot*R, c*R, f*R, R

'''
This script is designed to search all subfolders in a current folder for ABF 
files and analyze photocurrent transients in all found files. It outputs the 
total, capacitive, and faradaic charge for each measurement in the CSV files
grouped by the folder.
'''

if __name__ == '__main__':
    ### Find directory names
    dirs = sorted(glob.glob('*/'))
    ### Get abfs filenames in each directory
    abfs = [sorted(glob.glob(d + '*abf')) for d in dirs]
    ### generate sample names from directory names and print
    names = [x[:-1] for x in dirs]
    print(names)

    ### Calculate delivered charges for each measurement
    results = [[fitting(file) for file in a] for a in abfs]

    ### Reshape the array for better handling
    data = [np.array([[a[0],a[1],a[2],a[3]] for a in b]) for b in results]

    ### Export to csv, no header 
    for d,name in zip(data,names):
        np.savetxt(f'analysis_out_{name}-.csv', d, delimiter = ',')