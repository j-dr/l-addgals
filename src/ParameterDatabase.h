#ifndef PD_H
#define PD_H
#include <string>
#include <simulation.h>
#include "Logger.h"

#include "List.h"

//using std::string;
using namespace std;

//#define MAXLENGTH 200 //defined in List.h

class ParameterNode
{
private:
	const string parameterName;
	double parameterValue;
	ParameterNode* const nextParameterNode;
public:
	ParameterNode(ParameterNode* nextParameterNode, string parameterName, double parameterValue=0): parameterName(parameterName), nextParameterNode(nextParameterNode), parameterValue(parameterValue){}
	const string getParameterName();
	double getParameterValue();
	void setParameterValue(double ParameterValue);
	ParameterNode* getNextParameterNode();
};

class ParameterDatabase
{

private:
	ParameterNode* rootParameterNode;

public:
	ParameterDatabase(char* fileName);
	void print();
	int fillTable(ParameterNode* samRootParameterNode);
	double findParameterValue(string parameterName);
	~ParameterDatabase();
};

#endif

