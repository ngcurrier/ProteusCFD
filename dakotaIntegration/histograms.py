#!/usr/bin/env python

import collections
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.stats
import os

# Reads a dakota tabular output file
def ReadDakotaTabFile(filename):
    dict = ReadTabularFile(filename)
    for key in list(dict.keys()):
        # eliminate the columns which contain no useful data
        if key == '%eval_id' or key == 'interface':
            del dict[key]
        else:
            tmp = []
            for val in dict[key]:
                tmp.append(float(val))
            dict[key] = tmp
    return dict

# Reads a tabular file and loads results into an ordered dict
def ReadTabularFile(filename):
    dict = collections.OrderedDict()
    f = open(filename,'r')
    headers = f.readline()
    headers = headers.split()
    for head in headers:
        dict[head] = []
    for line in f:
        contents = line.split()
        for value,key in zip(contents,dict):
            dict[key].append(value)
    f.close()
    return dict

# Plots histogram of data contained in dictionary under 'key'
# Will dump plots for ./pic_histograms folder
# showPlot = True will display plots one at a time, if False only write to disk
def plotHistogram(dict, key, showPlot=True):
    # check for directory histograms
    filename = './pic_histograms/'+key+'.png'
    dir = os.path.dirname(filename)
    try:
        os.stat(dir)
    except:
        os.mkdir(dir)
    
    arr = np.array(dict[key])
    sigma = np.std(arr)
    mu = np.mean(arr)
    
    # the histogram of the data
    n, bins, patches = plt.hist(arr, 50, density=1, facecolor='green', alpha=0.75)

    # add a 'best fit' line
    y = scipy.stats.norm.pdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)

    num = 1
    fig = plt.figure(num, facecolor='white', figsize=(15,10))
    ax = fig.add_subplot(111)
    
    plt.xlabel(key)
    plt.ylabel('Probability')
    plt.title(r'$\mathrm{Histogram\ of\ ' + key +':}\ \mu=' + str(mu) + ',\ \sigma=' + str(sigma) +'$')
    plt.grid(True)

    plt.savefig(filename, bbox_inches='tight')
    if showPlot:
        plt.show()

if __name__ == "__main__":
    nargs = len(sys.argv)
    if nargs != 2:
        print('USAGE: ' + sys.argv[0] + ' <filename>')
        exit()
    dict = ReadDakotaTabFile(sys.argv[1])
    for key in list(dict.keys()):
        plotHistogram(dict, key, True)
