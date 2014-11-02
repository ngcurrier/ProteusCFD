#!/usr/bin/python

import sys
import math

def main():

    
    if len(sys.argv) != 3:
        print 'USAGE: ' + sys.argv[0] + ' filename ' + ' boundary id of interest'
        return

    targetString = 'Lift coefficient for body[' + sys.argv[2]
    targetTimestep = ': \n'

    filename = sys.argv[1]
    try:
        f = open(filename, 'r')
    except:
        print 'File does not exist: ' + filename
        return


    for line in f:
        if targetTimestep in line:
            'replace the colon and newline with a comma so we get a CSV file'
            line = line.replace (targetTimestep, ', ')
            print line,
            continue
        if targetString in line:
            'find the last colon and get the string after that which is the numeric value'
            'of the lift coefficient'
            pos = line.rfind(':')
            size = len(line)
            print line[pos+1:size],
            continue

if __name__ == '__main__':
    main()
