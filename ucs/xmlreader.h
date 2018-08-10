#ifndef XML_READER_H__
#define XML_READER_H__

#include "general.h"
#include <string>
#include <vector>
#include <map>
#include "tinyxml.h"

class XMLReaderCallbackInterface
{
public:
    // The prefix "cbi" is to prevent naming clashes.
    virtual void cbiCallbackFunction(int) = 0;
};

class CallbackCalleeExample : public XMLReaderCallbackInterface
{
public:
    // The callback function that Caller will call.
    void cbiCallbackFunction(int i)  
    { 
      printf("  Callee::cbiCallbackFunction() inside callback\n");
    }
};

class XMLReader
{
 public:
  XMLReader(std::string filename);
  ~XMLReader();

  void RegisterCallback(std::string path, XMLReaderCallbackInterface* callback);
  void ProcessXML();
  
 protected:
  
 private:
  XMLReader(){}; //no use
  TiXmlDocument* doc;
  std::map<std::string, XMLReaderCallbackInterface*> callbacks;
};

// a utility function defining a very simple method to indent a line of text
const char * getIndent( unsigned int numIndents );
// dumps out an xml strcture to stdout
void dump_to_stdout( TiXmlNode * pParent, unsigned int indent = 0 );

#endif
