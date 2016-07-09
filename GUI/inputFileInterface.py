#!/usr/bin/python
import collections

# this function parses a proteus inputfile and returns a
# dictionary with keys as the keywords
def parseInputFile(filename):
    f = open(filename, 'r')
    dict = collections.OrderedDict()
    for line in f:
        line = line.strip()
        # skip empty lines
        if len(line) == 0:
            continue
        # skip comment lines and commands
        if line[0] != '#' and line[0] != '<':
            # skip any lines that are not definitions, identified by '='
            if '=' in line:
                # end of line comments are allowed, strip them off
                line = line.split('#')[0]
                # split the line on the = sign
                if not line.find('='):
                    continue
                parts = line.split('=')
                key = parts[0].strip()
                val = parts[1].strip()
                dict[key] = val
    f.close()
    return dict

# this function takes a dictionary generated from parseInputFile
# which may have modified values and writes the changes back to
# this inputfile for solution
def writeInputFile(inputFile, outputFile, modifiedDict):
    fin = open(inputFile, 'r')
    fout = open(outputFile, 'w')
    
    for line in fin:
        for key in modifiedDict:
            # if the key is an exact match, modify the line
            if key == line.split('=')[0].strip():
                line = key + ' = ' + str(modifiedDict[key]) +'\n'
        fout.write(line)
                
    fin.close()
    fout.close()
    
if __name__ == "__main__":

    filename = 'flatplate.param'
    outputFile = 'written.param'
    
    dict = parseInputFile(filename)
    for item in dict:
        print item + ' = ' + dict[item]

    #change a key to test
    dict['refLength'] = 20000
        
    writeInputFile(filename, outputFile, dict)
