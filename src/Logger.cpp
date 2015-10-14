#include "Logger.h"


void Logger::print(int level, char* msg, ...)
{
	if(level>0)
	{
		va_list argptr;
		va_start(argptr,msg);
		fprintf(stderr, "\n");
		vfprintf(stderr, msg,argptr);
		va_end(argptr);
	}
}

