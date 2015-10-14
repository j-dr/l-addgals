#ifndef SD_H
#define SD_H

#include "Logger.h"

using std::string;

#define MAXLENGTH 200

class StringNode
{
private:
	const string parameterName;
	string parameterValue;
	StringNode* const nextStringNode;
public:
	StringNode(StringNode* nextStringNode, string parameterName, string parameterValue=NULL): parameterName(parameterName), nextStringNode(nextStringNode), parameterValue(parameterValue){}
	const string getParameterName();
	string getParameterValue();
	void setParameterValue(string ParameterValue);
	StringNode* getNextStringNode();
};

class StringDatabase
{

private:
	StringNode* rootStringNode;

public:
	StringDatabase(char* fileName);
	void print();
	int fillTable(StringNode* samRootStringNode);
	string findParameterValue(string parameterName);
	~StringDatabase();
};

#endif

