#!/usr/bin/python
import xml.etree.ElementTree

def readParameterDefinitions(filename):
    
    tree = xml.etree.ElementTree.parse(filename)
    root = tree.getroot()

    for node in root:
        print node.tag
        try:
            print node.find('description').text
        except:
            print 'no description'
            
if __name__ == "__main__":
    filename = '../parameters.xml'
    readParameterDefinitions(filename)    
