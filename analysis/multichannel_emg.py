'''
Multichannel EMG analysis

Originally run in Jupyter Notebook.

Written by A. Prominski
'''

# Load the background Python functions that allow for data loading and plotting
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

plt.rcParams.update({'font.size': 5})
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['lines.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

def to_secs(s):
    '''
    changes 24hr time string to seconds
    '''

    a = int(s[0:2])
    b = int(s[2:4])
    c = int(s[4:6])
    return a*3600 + b*60 + c


def plot_sigs(n, xlim = None):

    '''
    plot signals for inspection
    '''
    result, data_present = load_file(rhds[match[n]], verbose = False)
    
    fig,axs = plt.subplots(3,1,figsize=(11,8.5))
    
    a002,a003,a004,ain = [get_tV(c, result) for c in chs]
    
    print('Raw data')
    for d,ax in zip((a002,a003,a004),axs):
        ax.plot(*d)
        #ax.set_xlim(-0.1,2)
    plt.show()
    
    
    baselines = [pybaselines.smooth.ipsa]

    for base in baselines:
        print(base.__name__)
        fig,axs = plt.subplots(3,1,figsize=(11,8.5))
        for d,ax in zip((a002,a003,a004),axs):
            x,y = d
            bkg_1 = base(y)[0]
            ax.plot(x, y-bkg_1)
            #ax.set_xlim(-0.1,2)
         
        plt.show()
    

def get_peaks(n):
    '''
    find and plot all the peaks from the file
    '''
    fig, axs = plt.subplots(11,3)
    
    result, data_present = load_file(rhds[match[n]], verbose = False)
    a002,a003,a004,ain = [get_tV(c, result) for c in chs]
    start = np.argwhere(a002[0] == 0)[0][0]
    
    for d in (a002,a003,a004):
        for n in range(0,12):
            x = (d[0][start+10000*n:start+250+10000*n]-n, d[1][start+10000*n:start+250+10000*n]-d[1][start+10000*n])
            
            plt.plot(*x)
        plt.show()

def plot_scan(l : list):
    '''
    Plots the entire series of recordings. Used for SI figure.
    '''
    N = len(l)
    
    length = 2500
    fig, axs = plt.subplots(3,N, figsize=(6,3))

    limits = [[0,0],[0,0],[0,0]]
    #fig.suptitle(f'{abfs[l[0]]}', fontsize=10)
    
    
    output_array = []
    for e,n in enumerate(l):
        scan_array = []
        
        result, data_present = load_file(rhds[match[n]], verbose = False)
        a002,a003,a004,ain = [get_tV(c, result) for c in chs]
        start = np.argwhere(a002[0] == 0)[0][0]
        
        x,y = a004
        bkg_1 = pybaselines.smooth.ipsa(y)[0]
        a004 = (x, y-bkg_1)
        
        window = 10000 ### in us, 25 ms
        
        for d,ax,lim in zip((a002,a003,a004),axs[:,e], limits):
            channel_array = []
            color = cm.Reds(np.linspace(0, 1, 12))
            for n,c in zip(range(0,12),color):
                x = (d[0][start+window*n:start+length+window*n]-n, d[1][start+window*n:start+length+window*n]-d[1][start+window*n])
                
                max_v = np.amax(x[1])
                min_v = np.amin(x[1])
                channel_array.append(max_v)
                
                if lim[0] > min_v:
                    lim[0] = min_v
                if lim[1] < max_v:
                    lim[1] = max_v
                    
                ax.plot(*x, c=c, linewidth=1)
            scan_array.append(channel_array)
        output_array.append(scan_array)
                
    print(limits)
    for ax,lim in zip(axs,limits):
        ylim = lim
        for sax in ax:
            sax.set_ylim(ylim[0]*1.1,ylim[1]*1.1)
            
            for axis in ['top','bottom','left','right']:
                sax.spines[axis].set_linewidth(0)
    

        
    for ax in axs:
        for sax in ax[1:]:
            
            sax.set_xticks([])
            sax.set_yticks([])
            
    
    
    axs[2][0].set_xlabel('Time (s)')
    
    axs[0][0].set_ylabel('EMG (μV)')
    axs[1][0].set_ylabel('EMG (μV)')
    axs[2][0].set_ylabel('EMG (μV)')
    plt.tight_layout()
    plt.savefig(f'EMG{str(l)}.svg', dpi= 600)
    plt.show()
    return output_array

## Fetches data from ABF file
def get_data(filename):
    abf = pyabf.ABF(filename)

    abf.setSweep(0, channel=0)
    t_vector = abf.sweepX
    TTL = abf.sweepY

    abf.setSweep(0, channel=1)
    TTL2 = abf.sweepY
    
    abf.setSweep(0, channel=2)
    LVP = abf.sweepY
    
    abf.setSweep(0, channel=3)
    ECG = abf.sweepY
    
    abf.setSweep(0, channel=4)
    PP = abf.sweepY

    return t_vector, PP, ECG, LVP, TTL, TTL2

def plot_full(abf_name):
    t_vector, PP, ECG, LVP, TTL, TTL2 = get_data(abf_name)
    x_data = t_vector
    y_data = [PP, ECG, LVP, TTL, TTL2]
    
    colors = ['black', 'blue', 'orange', 'red', 'green']
    
    fig, axs = plt.subplots(5,1,figsize=(11,8.5))
    
    for y_d,ax,col in zip(y_data,axs, colors):
        ax.plot(x_data, y_d, c = col)
        ax.set_xlim(4.5,6)
    plt.show()
    plt.close()


def isochrone_map(r, beats, offset = 0, name = 'map'):
    '''
    Main function for plotting isochorne maps
    '''
    result, data_present = load_file(r)
    
    ### electrode connection on MEA, was confirmed using laser irradiation
    
    ch_n = ['13','22','10','23','08','09','14','19','12', '17','18','15','20', '21', '16', '11']
    chnnls = [f'A-0{c}' for c in ch_n]


    T = []
    V = []
    fig, axs = plt.subplots(16,1)
    fig.set_size_inches(11, 8.5)
    for ch,ax in zip(chnnls,axs):

  
        t, v = get_tV(ch, result)
        T.append(t)
        V.append(v)
        ax.plot(t,v)
        ax.set_xlim(-2,4)

      
    plt.savefig(f'{r[-21:-4]}-signals.png', dpi=600)
    plt.show()  
    plt.close()
    Tstack = np.stack(T)
    print(f'tstack shape: {Tstack.shape}')
    Vstack = np.stack(V)
    
    reT = np.reshape(Tstack, (4,4,Tstack.shape[-1]))
    reV = np.reshape(Vstack, (4,4,Vstack.shape[-1]))
    print(f'reT shape: {reT.shape}')

    ### max plot
    zero = np.argwhere(reT[0,0,:] == 0)[0][0]

    mean = np.zeros((4,4))

    max_N = beats
    spacing = 10*10000//beats
    
    ### ignores first peak due to possible 
    for N in range(1,max_N):
        print(f'Signal: {N}')
        #print(f'map {N}')
        Vs = reV[:,:,zero+offset+spacing*N:zero+offset+spacing+spacing*N]
        ts = reT[:,:,zero+offset+spacing*N:zero+offset+spacing+spacing*N]*1000 ## to ms
        #plt.plot(ts[0,0,:],Vs[0,0,:])
        #plt.show()
        #print(np.amax(Vs, axis=-1))
        maxarg = np.argmax(Vs, axis = -1)
       
        out = np.empty((4,4))
        for n in range(4):
            for k in range(4):
                #tharg = np.argwhere(Vs[n,k,:] > 4000)[0]
                #out[k,n] = ts[k,n,tharg]
                out[k,n] = ts[k,n,maxarg[k,n]]

        img = out-ts[:,:,0]
        #print(reT[:,:,zero+offset+spacing*N])
        mean += img
        
        #img = img.round(4)

        fig, ax = plt.subplots()
        im = ax.imshow(img, interpolation='bicubic')

        for i in range(4):
            for j in range(4):
                text = ax.text(j, i, img[i, j], ha="center", va="center", color="w")

        plt.show()
        plt.close()
        
        

    plt.rcParams.update({'font.size': 5})
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    print('Mean map') 
    mean = mean/(max_N-1)
    img = mean-np.amin(mean)
    img = img.round(1)

    fig, ax = plt.subplots(figsize=(1.75,1.75))
    im = ax.imshow(img, cmap='coolwarm', interpolation = 'gaussian', vmin=0, vmax=4)

#     for i in range(4):
#         for j in range(4):
#             text = ax.text(j, i, img[i, j], ha="center", va="center", color="w")
    ax.axis('off')
    #plt.colorbar(im, fraction = 0.2)
    plt.savefig(f'{name}-map.png', dpi = 600)
    plt.show()

if __name__ == '__main__':
    abfs = sorted(glob('*.abf'))
    rhds = sorted(glob('*.rhd'))
    for e,a in enumerate(abfs):
        print(f'{e}: {a}')

    abfs_time = []
    for a in abfs:
        f = pyabf.ABF(a)
        head = f.headerText
        a = re.search('abfDateTime.*$', head, re.MULTILINE)
        b = re.search('..:..:..', a[0])
        c = b[0].replace(':','')
        abfs_time.append(c)

    rhds_time = []
    for r in rhds:
        rhds_time.append(r[-10:-4])

    match = []
    for a in abfs_time:
        a_s = to_secs(a)
        for e,n in enumerate(rhds_time):
            flag = False
            n_s = to_secs(n)
            if abs(n_s-a_s) < 5:
                print('Found!')
                flag = True
                match.append(e)
                break
        if flag == False:
            print(f'Match not found for {a}')

    chs = ['A-002','A-003','A-004','ANALOG-IN-07']

    l = range(0,10)
    output_array = plot_scan(l)

    # Normalizing signals' amplitudes
    EMG_amp = np.asarray(output_array)
    print(EMG_amp.shape)
    max_amp = np.amax(EMG_amp,axis = (0,2))
    print(max_amp)
    norm_amp = EMG_amp/max_amp[None,:,None]
    print(norm_amp.shape)

    ###selectivity index
    norm_amp_sum = np.sum(norm_amp, axis=1)

    selectivity_idx = norm_amp/norm_amp_sum[:,None,:]

    GM = [selectivity_idx[n,0,:] for n in range(8)]
    TA = [selectivity_idx[n,1,:] for n in range(8)]
    PI = [selectivity_idx[n,2,:] for n in range(8)]
    fig, axs = plt.subplots(1,3,figsize=(10,4))
    axs[0].boxplot(GM)

    axs[1].boxplot(TA)

    axs[2].boxplot(PI)

    for ax in axs:
        ax.set_ylim(0,1)
    plt.show()

    # preparing and plotting the final figure
    x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    gm = np.mean(selectivity_idx[:,0,:], axis=-1)
    ta = np.mean(selectivity_idx[:,1,:], axis=-1)
    pi = np.mean(selectivity_idx[:,2,:], axis=-1)

    gm[0] = 0
    pi[0] = 0

    gm[-1] = 0
    pi[-1] = 0

    fig,ax = plt.subplots(figsize=(1.5,1))

    ax.fill_between(x, gm[1:-1]*100,0, facecolor='blue', alpha=0.25)
    ax2=ax.twinx()

    ax2.fill_between(x, (pi[1:-1])*100, 0, facecolor='red', alpha=0.25)
    ax2.invert_yaxis()
    ax.set_ylabel('Muscle selectivity (%)')
    ax.set_ylim(0,100)
    ax2.set_ylim(100,0)

    ax.set_xlabel('Laser position (mm)')
    plt.tick_params(right = False, labelright = False)
    plt.xlim(0,0.9)
    plt.savefig('power7activation.png', dpi=600)
    plt.savefig('power7activation.svg', dpi=600, transparent=True)
    plt.show()

    name='LVonly'
    n = 0
    r = rhds[matched_rhd[n]]
    print(r[-21:-4])
    offset = 0
    beats = 10
    isochrone_map(r, beats, offset, name)