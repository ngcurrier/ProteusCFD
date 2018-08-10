#ifndef XML_READER_H__
#define XML_READER_H__

#include "general.h"
#include <string>
#include <vector>
#include <map>
#include "tinyxml.h"
#include "exceptions.h"

class XMLReaderCallbackInterface
{
public:
    virtual void callbackFunction(TiXmlNode* node) = 0;
};

class CallbackCalleeDouble : public XMLReaderCallbackInterface
{
 public:
 CallbackCalleeDouble(double& storagePtr, std::string key):
  storagePtr(storagePtr), key(key) { };
  void callbackFunction(TiXmlNode* node)
  {
    int t = node->Type();
    if(t == TiXmlNode::TINYXML_TEXT){
      TiXmlPrinter printer;
      node->Accept(&printer);
      std::stringstream ss(printer.CStr());
      ss >> storagePtr;
    }
    else{
      Abort << "CallbackCalleeDouble::calbackFunction found wrong type";
    }
  };
 private:
  double& storagePtr;
  std::string key;
};

class CallbackCalleeExample : public XMLReaderCallbackInterface
{
public:
    // The callback function that Caller will call.
    void callbackFunction(TiXmlNode* node)  
    { 
      printf("  Callee::callbackFunction() inside callback\n");
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
