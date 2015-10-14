#include "StringDatabase.h"


const string StringNode::getParameterName()
{
	return parameterName;

}

string StringNode::getParameterValue()
{
	return parameterValue;
}

void StringNode::setParameterValue(string parameterValue)
{
	this->parameterValue = parameterValue; 
}

StringNode* StringNode::getNextStringNode()
{
	return nextStringNode;
}


StringDatabase::StringDatabase(char* fileName)
{
	FILE *inputFile = fopen(fileName, "r");
	rootStringNode = NULL;

	char parameterName[MAXLENGTH];
	char parameterValue[MAXLENGTH];
	
	if(inputFile == NULL)
	{
		Logger::print(1, "Could not open parameter input file");
	}
	else
	{
		while (fscanf(inputFile, "%s", parameterName) == 1) 
		{
			fscanf(inputFile, " %s", parameterValue);
			rootStringNode = new StringNode(rootStringNode, parameterName, parameterValue);
		}
	}
	fclose(inputFile);
	//parameterTable->printList();
}


int StringDatabase::fillTable(StringNode* samStringNode)
{
	int unfilled = 0;
	while( samStringNode != NULL)
	{
		string parameterValue = findParameterValue(samStringNode->getParameterName());
		if(parameterValue != "")
		{
			samStringNode->setParameterValue(parameterValue);
		}
		else
		{
			unfilled++;
		}
		samStringNode = samStringNode->getNextStringNode();
	}
	return unfilled;
}

string StringDatabase::findParameterValue(const string parameterName)
{
	StringNode* ptr = rootStringNode;
	
	while( ptr != NULL)
	{
		if( Utility::caselessStringCompare(ptr->getParameterName(), parameterName) == 0)
		{
			return ptr->getParameterValue();
		}
		else
		{
			ptr = ptr->getNextStringNode();
		}
	}
	Logger::print(2, "StringDatabase::findParameterValue(string parameterName) parameterName=%s not found\n", parameterName.c_str());
	//return NULL;
	return "not_found";
}

void StringDatabase::print()
{
	StringNode* ptr = rootStringNode;
	
	while( ptr != NULL)
	{
		Logger::print(1, " %s %s", ptr->getParameterName().c_str(), ptr->getParameterValue().c_str() );
		ptr = ptr->getNextStringNode();
	}
}

StringDatabase::~StringDatabase()
{
	StringNode* ptr;
	
	while( rootStringNode != NULL)
	{
		ptr = rootStringNode;
		rootStringNode = rootStringNode->getNextStringNode();
		delete ptr;
	}
}
