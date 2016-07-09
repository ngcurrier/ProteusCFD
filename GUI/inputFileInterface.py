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
            # skip any lines that are not definitions
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
    return dict
            
if __name__ == "__main__":
    filename = 'flatplate.param'
    dict = parseInputFile(filename)
    for item in dict:
        print item + ' = ' + dict[item]
