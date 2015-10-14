#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "Utility.h"

#ifndef NULL
#define NULL 0
#endif

using namespace std;


class Logger
{
public:
	static void print(int level, char* msg, ...);
};

#endif
