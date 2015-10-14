#include "ParameterDatabase.h"


const string ParameterNode::getParameterName()
{
	return parameterName;

}

double ParameterNode::getParameterValue()
{
	return parameterValue;
}

void ParameterNode::setParameterValue(double parameterValue)
{
	this->parameterValue = parameterValue; 
}

ParameterNode* ParameterNode::getNextParameterNode()
{
	return nextParameterNode;
}


ParameterDatabase::ParameterDatabase(char* fileName)
{
	FILE *inputFile = fopen(fileName, "r");
	rootParameterNode = NULL;

	char parameterName[MAXLENGTH];
	char parameterValue[MAXLENGTH];
	double parameterDoubleValue;

	if(inputFile == NULL)
	{
		Logger::print(1, "Could not open parameter input file");
	}
	else
	{
		while (fscanf(inputFile, "%s", parameterName) == 1) 
		{
			fscanf(inputFile, " %s", parameterValue);
			parameterDoubleValue = atof(parameterValue);
			rootParameterNode = new ParameterNode(rootParameterNode, parameterName, parameterDoubleValue);
		}
	}
	fclose(inputFile);
	//parameterTable->printList();
}


int ParameterDatabase::fillTable(ParameterNode* samParameterNode)
{
	int unfilled = 0;
	while( samParameterNode != NULL)
	{
		double parameterValue = findParameterValue(samParameterNode->getParameterName());
		if(parameterValue!=-1)
		{
			samParameterNode->setParameterValue(parameterValue);
		}
		else
		{
			unfilled++;
		}
		samParameterNode = samParameterNode->getNextParameterNode();
	}
	return unfilled;
}

double ParameterDatabase::findParameterValue(const string parameterName)
{
	ParameterNode* ptr = rootParameterNode;
	
	while( ptr != NULL)
	{
		if( Utility::caselessStringCompare(ptr->getParameterName(), parameterName) == 0)
		{
			return ptr->getParameterValue();
		}
		else
		{
			ptr = ptr->getNextParameterNode();
		}
	}
	Logger::print(2, "ParameterDatabase::findParameterValue(string parameterName) parameterName=%s not found", parameterName.c_str());
	return -1;
}

void ParameterDatabase::print()
{
	ParameterNode* ptr = rootParameterNode;
	
	while( ptr != NULL)
	{
		Logger::print(1, " %s %f", ptr->getParameterName().c_str(), ptr->getParameterValue() );
		ptr = ptr->getNextParameterNode();
	}
}

ParameterDatabase::~ParameterDatabase()
{
	ParameterNode* ptr;
	
	while( rootParameterNode != NULL)
	{
		ptr = rootParameterNode;
		rootParameterNode = rootParameterNode->getNextParameterNode();
		delete ptr;
	}
}
