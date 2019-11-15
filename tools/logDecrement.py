#!/usr/bin/env python

import sys
import math

def main():

    
    if len(sys.argv) != 3:
        print 'USAGE: ' + sys.argv[0] + ' filename ' + ' column of interest'
        return

    filename = sys.argv[1]
    col = sys.argv[2]
    col = int(col)

    try:
        f = open(filename, 'r')
    except:
        print 'File does not exist: ' + filename
        return

    # read in the lines
    lines = f.readlines()
    # remove ":" character from lines
    for i in range(0,len(lines)):
        if(lines[i].find(":")):
            lines[i] = lines[i].replace(":","")
    # do a double nested list comprehension to get the rest of the data
    data = [[float(val) for val in line.split()] for line in lines[0:]]
    
    # allocate a list
    column = []
    
    print '%d steps read from file' % len(data)

    # build out the list from the column in the data matrix
    for i in range(1, len(data)-1):
        column.append(data[i][col])

    # find a high spot and the next one, and calculate the log decrement
    # start at the end of the list, we assume the behavior is periodic there
    j = 0
    maxj = []
    maxi = []

    for i in range(len(column)-2 , -1, -1):
        if((column[i+1] < column[i]) and (column[i-1] < column[i])):
            # this is a debugging line, will print out line number and the three
            # points for which the max peak was detected, useful for weird data
            # print '%d %.8f %.8f %.8f)'% (i, column[i-1], column[i], column[i+1])
            maxj.append(column[i])
            maxi.append(i)
            j = j+1

    print '%d peaks found:' % j
    print '------------------------------'
    for i in range(0,j):
        print 'step: %d %.6f' % (maxi[i], maxj[i])
    if(j < 2):
        print 'Not enough data... patience is a virtue'
        return
        
    # look at the first two values, caculate log decrement
    delta = math.log(maxj[1]/maxj[0])
    pi = math.acos(-1.0)
    print 'delta: %.6f '% delta
    print 'damping ratio: %.6f ' % (1.0/(math.sqrt(1.0 + ((2.0*pi/delta) ** 2.0))) )
                        

if __name__ == '__main__':
    main()
