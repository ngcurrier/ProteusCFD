#!/usr/bin/python
 
import matplotlib.pyplot as plt
import numpy as np
import math as math
from matplotlib.colors import BoundaryNorm
 
def readCorrelationFile(filename):
    # assumes file is structured as dakota output
    # first row is all variables in list
    # each successive row starts with variable name
    # and contains at first 1 value, then 2 on the second row and so on
    f = open(filename, 'r')
 
    dict = {}
    dict['labels'] = f.readline().split()
    max_values = len(dict['labels'])
 
    # create numpy matrix
    dict['data'] = np.zeros((max_values,max_values))
 
    values = 1
    for line in f:
        entries = line.split()
        for j in range(1,values+1):
            dict['data'][values-1, j-1] = float(entries[j])
        values = values + 1
   
    return dict
 
def plotHeatmap(column_labels, row_labels, data):
    # assume data is square
    mat_size = math.sqrt(data.size)
    plt.plot()
    cmap = plt.get_cmap('RdBu')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    bounds = np.linspace(-1,1,100)
    norm = BoundaryNorm(boundaries=bounds, ncolors=256)
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    heatmap = plt.pcolor(data, norm = norm, cmap=cmap)
 
    # put the major ticks at the cell middle
    plt.xticks(np.arange(0.5, mat_size+0.5, 1.0))
    plt.yticks(np.arange(0.5, mat_size+0.5, 1.0))
 
    # create a more table-like display
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.xaxis.set_ticks_position('top')
   
    ax.set_xticklabels(row_labels, minor=False)
    ax.set_yticklabels(column_labels, minor=False)

    plt.colorbar(format='%.2f')
    plt.xticks(rotation=45)
    plt.show()
 
if __name__ == "__main__":
 
    dict = readCorrelationFile('correlationTest.txt')
    plotHeatmap(dict['labels'], dict['labels'], dict['data'])
