'''
Isochronal map preparation

Originally run in Jupyter Notebook.

Written by A. Prominski
'''

# Load the background Python functions that allow for loading RHD files
from rhd_utilities import * # this is poor practice, should be fixed

from glob2 import glob
import matplotlib.pyplot as plt
import pyabf

import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

import pybaselines
from scipy.signal import find_peaks
from matplotlib import cm

def to_secs(s):
    '''
    changes 24hr time string to seconds
    '''

    a = int(s[0:2])
    b = int(s[2:4])
    c = int(s[4:6])
    return a*3600 + b*60 + c

if __name__ == '__main__':
    ### get filenames
    abfs = sorted(glob('03*/*.abf'))
    rhds = glob('RHD/**/*.rhd')

    ## read ABFs times
    abfs_time = []
    for a in abfs:
        f = pyabf.ABF(a)
        head = f.headerText
        a = re.search('abfDateTime.*$', head, re.MULTILINE)
        b = re.search('..:..:..', a[0])
        c = b[0].replace(':','')
        abfs_time.append(c)

    ## read RHDs times 
    rhds_time = []
    for r in rhds:
        rhds_time.append(r[-10:-4])

    matched_rhd = []
    for a in abfs_time:
        a_s = to_secs(a)
        for e,n in enumerate(rhds_time):
            flag = False
            n_s = to_secs(n)
            if abs(n_s-a_s) < 10:
                print('Found!')
                flag = True
                matched_rhd.append(e)
                break
        if flag == False:
            matched_rhd.append(None)
            print(f'Match not found for {a}')