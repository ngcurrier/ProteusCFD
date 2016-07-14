#!/usr/bin/python
import xml.etree.ElementTree

def readParameterDefinitions(filename):
    
    tree = xml.etree.ElementTree.parse(filename)
    root = tree.getroot()

    for node in root:
        print node.tag
        try:
            description = node.find('description').text
            keyword = node.find('keyword').text
            default = node.find('default').text
            min = node.find('minimum').text
            max = node.find('maximum').text
            default = node.find('default').text
            print keyword + '=' + default
        except:
            print 'no description'
            
if __name__ == "__main__":
    filename = '../parameters.xml'
    readParameterDefinitions(filename)    
