#include "Utility.h"

int Utility::caselessStringCompare(const char* string1, const char* string2)
{
	int length1 = strlen(string1);
	int length2 = strlen(string2);
	if(length1 != length2)
		return length1 - length2;
	for(int i=0;i<length1;i++)
	{
		if(lowerChar(string1[i]) != lowerChar(string2[i]))
		{
			return (string1[i] - string2[i]);
		}
	}
	return 0;
}

int Utility::caselessStringCompare(const string string1, const string& string2)
{
	int length1 = string1.length();
	int length2 = string2.length();
	if(length1 != length2)
		return length1 - length2;
	for(int i=0;i<length1;i++)
	{
		if(lowerChar(string1[i]) != lowerChar(string2[i]))
		{
			return (string1[i] - string2[i]);
		}
	}
	return 0;
}

char Utility::lowerChar(char ch)
{
	if ( ( ch >= 'A' && ch <= 'Z' ) )
	{
		return ch + 'a' - 'A';
	}
	else
		return ch;
}

