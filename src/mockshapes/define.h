/*
 				define.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	mapmap
*
*	Author:		J. P. Dietrich, University of Michigan
*
*	Contents:	global definitions.
*
*	Last modify:	2009-12-31
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------ what, who, when and where ------------------------*/

#define		BANNER		"mockshapes"
#define		COPYRIGHT	"J.P. Dietrich <jorgd@umich.edu>"
#define		INSTITUTE	"University of Michigan"


/*----------------------------- Internal constants --------------------------*/

#define	MAXCHAR			256		/* max. number of characters */
#define	BIG			1e+30		/* a huge number */


/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)


/*------------------------------- Other Macros -----------------------------*/

#define MALLOC(ptr, size, msg)                  \
    do {                                        \
        (ptr) = malloc(size*sizeof(*ptr));      \
        if ((ptr) == NULL){                     \
            perror(msg);                        \
            exit(EXIT_FAILURE);                 \
        }                                       \
    } while(0)

#define CALLOC(ptr, size, msg)                  \
    do {                                        \
        (ptr) = calloc(size, sizeof(*ptr));     \
        if ((ptr) == NULL){                     \
            perror(msg);                        \
            exit(EXIT_FAILURE);                 \
        }                                       \
    } while(0)



#define FITS_ERR_LEN (128)

#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif

#ifndef SQRT2l
# define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
#endif

#ifndef PATH_MAX
# define PATH_MAX       4096
#endif

#define SYNTAX \
"mockshapes catalog [-c <configuration_file>] \\\n" \
"                       [-<keyword> <value>]\n" \
"                   or, to dump a default configuration file:\n" \
"                   mockshapes -d \n"

