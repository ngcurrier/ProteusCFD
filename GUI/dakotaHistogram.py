#!/usr/bin/python

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#reads a space delimited file with a header and returns a dictionary
#attempts to cast dictionary entries into floats, if it fails, leaves as is
def readSpaceDelimitedFile(filename):
    f = open(filename, 'r')
    headers = f.readline().split()
    dict = {}
    for header in headers:
        dict[header] = []
        

    for line in f:
        items = line.split()
        i = 0
        for header in headers:
            try:
                dict[header].append(float(items[i]))
            except:
                dict[header].append(items[i])
            i = i + 1
        
    f.close()
    return dict

#plots a histogram of data, computes basic stats, and labels chart
def plotHistogram(data, seriesName):
    # the histogram of the data
    n, bins, patches = plt.hist(data, 50, normed=1, facecolor='green', alpha=0.75)

    mu = np.mean(data)
    sigma = np.std(data)
    
    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)
    
    plt.xlabel(seriesName)
    plt.ylabel('Probability')
    plt.title(r'$\mathrm{Histogram\ of\ ' + seriesName + ':}\ \mu=' + str(mu) +',\ \sigma=' + str(sigma) +'$')
    plt.grid(True)
    
    plt.show()

if __name__ == '__main__':
    data = readSpaceDelimitedFile('dakota_tabular.dat')
    print data
    for idata in data:
        if idata != 'interface' and idata != '%eval_id':
            plotHistogram(data[idata], idata)
    
