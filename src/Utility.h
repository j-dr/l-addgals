#ifndef UTILITY
#define UTILITY


#include "Logger.h"

class Utility
{
public:
	static int caselessStringCompare(const char* string1, const char* string2);
	static int caselessStringCompare(const std::string string1, const std::string& string2);
	static char lowerChar(char ch);
};


#endif

