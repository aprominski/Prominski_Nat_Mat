'''
Single channel EMG analysis

Originally run in Jupyter Notebook.

Written by A. Prominski
'''

import pyabf
import glob
import numpy as np
from matplotlib import pyplot as plt


def get_data(filename):
    '''
    get data, single EMG sweep
    '''
    abf = pyabf.ABF(filename)
    
    abf.setSweep(sweepNumber=0, channel=0)
    x = abf.sweepX
    TTL = abf.sweepY   
    abf.setSweep(sweepNumber=0, channel=1)
    EMG = abf.sweepY
    
   
    data = x,EMG,TTL
        
    return data

def find_all_pulses(a):
    '''
    detect all light pulses in the file
    '''
    pulse_up = []
    pulse_down = []
    tmp = 0
    UP = 1
    DOWN = 0
    direction = UP

    for n in range(1, len(a)):
        
        if direction == UP and a[n] - a[tmp] > 0.1:
            pulse_up.append(n)
            tmp = n
            direction = DOWN
            #print('switch down')
        elif direction == DOWN and a[n] - a[tmp] > 0.1:
            pass
        elif direction == DOWN and a[n] - a[tmp] < -0.1:
            pulse_down.append(n)
            tmp = n
            direction = UP
            #print('switch up')
        elif direction == UP and a[n] - a[tmp] < -0.1:
            pass    
    
    return pulse_up,pulse_down

def extract_pulses(abf,n):  
    '''
    find and analyze pulses
    '''
    print('{}'.format(n))
    
    ## get data
    x, EMG, TTL = get_data(abf)
    
    ## find where ligth pulses are
    pulse_up, pulse_down = find_all_pulses(TTL)
    
    
    ## extract pulses into the separate arrays
    time_base = 0.025 # 10 us
    length = pulse_down[0] - pulse_up[0]
    max_len = length*time_base
    
    t_extent =100 # datapoints ms
    extent = int(t_extent/time_base)
    
    x = np.arange(0,max_len+t_extent,time_base)
    
    data = None
    for up,down in zip(pulse_up,pulse_down):
        new_data = EMG[up:down+extent]
        data = np.vstack((data, new_data)) if data is not None else new_data
    subsets = data
    
    
    # calculate the and standard deviation of all the pulses
    mean = np.mean(subsets, axis=0)
    stdev = np.std(subsets, axis=0)
    
    ## plot the pulse +/- standard deviation
    plt.figure(figsize=(4,4))
    plt.fill_between(x, mean+stdev, mean-stdev, facecolor='blue', alpha=0.25)
    plt.plot(x,mean, c = 'red')
    plt.ylim(-0.3,0.7)
    plt.savefig('{}.svg'.format(n), dpi=300)
    plt.show()
    

    ## amplitude is defined as a miximum difference in the recorded signal
    amplitudes = []
    for x in subsets:
        amplitudes.append(np.amax(x) - np.amin(x))
    amplitudes = np.array(amplitudes)
    
    return amplitudes

if __name__ == 'main':
    abfs = glob.glob('*abf')

    names = [x[:-4] for x in abfs]

    for e,n, in enumerate(names):
        print(f'{e}: {n}')

    ### collect amplitudes from all recordings in the folder
    amplitudes = []
    for abf,n in zip(abfs,names):
        amplitudes.append(extract_pulses(abf,n))

    np.savetxt('analysis_out.csv', amplitudes, delimiter = ',')