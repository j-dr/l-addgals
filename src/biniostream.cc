// biniostream.cc
//
// The very few non-inlined parts of the biniostream classes
//
// Patrik Jonsson
//
// $id$
// $Log: biniostream.cc,v $
// Revision 1.1.1.1  2004/04/18 23:31:57  risa
// First sansa import
//
// Revision 1.1.1.1  2002/05/17 21:20:10  risa
// Imported sources
//
// Revision 1.2  1998/10/10 08:20:51  patrik
// Fixed bugs in the new version of the biniostream classes.
//
// Revision 1.1  1998/10/09 05:00:18  patrik
// Biniostream is rewritten as a base class binfstream : public fstream,
// and binifstream and binofstream are derived from binfstream. They are opened
// with default modes like ofstream and ifstream, and binfstream can do I/O
// both ways. They also have open() and default constructors.
//

#include "biniostream.h"

// Define the byte order of the machine
#ifdef BINIO_BIG_ENDIAN
const int binfstream::machine_big_endian=1;
#else
const int binfstream::machine_big_endian=0;
#endif
