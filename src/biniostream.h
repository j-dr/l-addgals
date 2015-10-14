// biniostream.h
//
// Declaration and inlines of binary iofilestreams
//
// Patrik Jonsson
//
// 020328: Risa Wechsler modification for gcc 3.0.4
//
// 981001: Added default constructors to anable array allocation, and
//         members for i/o of null-terminated strings
// $Id: biniostream.h,v 1.2 2004/07/04 22:23:43 risa Exp $
// $Log: biniostream.h,v $
// Revision 1.2  2004/07/04 22:23:43  risa
// removed long double from biniostream class
// added simple cosmology definition class
// added simple simulation class, which depends on cosmology
// added simple point class, which depends on simulation class
// for periodic distance measurements.
// added physical constants definition file.  you should take cspeed and G
// out of various programs.
//
// Revision 1.1.1.1  2004/04/18 23:31:57  risa
// First sansa import
//
// Revision 1.1.1.1  2002/05/17 21:20:10  risa
// Imported sources
//
//
// Revision 1.4  1998/10/09 05:00:18  patrik
// Biniostream is rewritten as a base class binfstream : public fstream,
// and binifstream and binofstream are derived from binfstream. They are opened
// with default modes like ofstream and ifstream, and binfstream can do I/O
// both ways. They also have open() and default constructors.
//
// Revision 1.3  1998/10/02 19:11:31  patrik
// Changed the biniostream classes to correctly handle endianness.
// The endianness of the files can now be set, and the define BINIO_BIG_ENDIAN
// must be defined on big endian machines to correctly handle this.
//

#ifndef __biniostream__
#define __biniostream__

#include <iostream>
#include <fstream>

using namespace std;

// Uncomment to read files written on a system that has
// different byte ordering (ie Sun-DEC)!
//#define BINIO_BIG_ENDIAN


// **********************
// *                    *
// * Binary file stream *
// *                    *
// **********************

// This is the base class for binary file streams and can be
// used for both reading and writing

class binfstream : public fstream {
private:
  // The byte order of the machine
  static const int machine_big_endian;

  // Byte order of the file
  int file_native;

public:
  binfstream();
  //  binfstream(const char*,int=0,int=0664);
  binfstream(const char*);
  ~binfstream();

  // Set file type - defaults to native
  void bigendian();
  void littleendian();
  void native();

  // These output operators replace the ofstream ones
  binfstream& operator<<(const char);
  binfstream& operator<<(const char*);
  binfstream& operator<<(const short int);
  binfstream& operator<<(const int);
  binfstream& operator<<(const long int);
  binfstream& operator<<(const float);
  binfstream& operator<<(const double);

  // These input operators replace the ifstream ones
  binfstream& operator>>(char&);
  binfstream& operator>>(char*);
  binfstream& operator>>(int&);
  binfstream& operator>>(short int&);
  binfstream& operator>>(long int&);
  binfstream& operator>>(float&);
  binfstream& operator>>(double&);
  binfstream& skip(int);
};


// Binofstream is derived from binfstream

class binofstream : public binfstream {
public:
  binofstream();
  binofstream(const char*);
  ~binofstream();

  void open(const char*);
};


// Binifstream is derived from binfstream

class binifstream : public binfstream {
public:
  binifstream();
  //  binifstream(const char*,int=ios::in);
  binifstream(const char*);
  ~binifstream();

  void open(const char*);
};


// *** Inline definitions ***

// *** binfstream ***


// Default constructor
inline binfstream::binfstream() : fstream()
{
  file_native=1;
}

// Constructor to open an ofstream to filename and send output there
//inline binfstream::binfstream(const char* filename,int mode,int prot) :
// fstream(filename,mode,prot)
inline binfstream::binfstream(const char* filename) :
  fstream(filename)
{
  file_native=1;
}

// Destructor 
inline binfstream::~binfstream()
{
}

// Byte order stuff
inline void binfstream::bigendian()
{
  file_native = machine_big_endian;
}
  
inline void binfstream::littleendian()
{
  file_native = !machine_big_endian;
}

inline void binfstream::native()
{
  file_native = 1;
}


// Outputs a char to the binary ofstream
inline binfstream& binfstream::operator<<(const char c)
{
  put(c);
  return *this;  
}

// Outputs a null-terminated string of char to the binary ofstream
inline binfstream& binfstream::operator<<(const char* c)
{
  while(*c)
    put(*c++);

  return *this;  
}

// Outputs a short int to the binary ofstream
inline binfstream& binfstream::operator<<(const short int d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  

  // Put bytes on output stream
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      put(*p++);
    else
      put(*p--);
  
  // Return binostream for continuing I/O 
  return *this;  
}

// Outputs an int to the binary ofstream
inline binfstream& binfstream::operator<<(const int d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;


  // Put bytes on output stream
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      put(*p++);
    else
      put(*p--);


  // Return binostream for continuing I/O 
  return *this;  
}

// Outputs a long int to the binary ofstream
inline binfstream& binfstream::operator<<(const long int d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;


  // Put bytes on output stream
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      put(*p++);
    else
      put(*p--);
  

  // Return binostream for continuing I/O 
  return *this;  
}

// Outputs a float to the binary ofstream
inline binfstream& binfstream::operator<<(const float d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  
  
  // Put bytes on output stream
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      put(*p++);
    else
      put(*p--);
  
  
  // Return binostream for continuing I/O 
  return *this;  
}

// Outputs a double to the binary ofstream
inline binfstream& binfstream::operator<<(const double d)
{
  char* p;

  // Get a pointer to the double
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  
  
  // Put bytes on output stream
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      put(*p++);
    else
      put(*p--);
  
  
  // Return binostream for continuing I/O 
  return *this;  
}

// Inputs a char from the binary ifstream
inline binfstream& binfstream::operator>>(char& c)
{
  get(c);
  return *this;
}

// Inputs a null-terminated string of char from the binary ifstream
inline binfstream& binfstream::operator>>(char* c)
{
  get(*c);
  while(*c++)
    get(*c);
  return *this;
}

// Inputs a short from the binary ifstream
inline binfstream& binfstream::operator>>(short int& d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  

  // Now read size number of bytes and copy unsigned into d's address
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      get(*p++);
    else
      get(*p--);
  

  // Return binistream for continuing I/O
  return *this;
}

// Inputs an int from the binary ifstream
inline binfstream& binfstream::operator>>(int& d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  

  // Now read size number of bytes and copy into d's address
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      get(*p++);
    else
      get(*p--);
  

  // Return binistream for continuing I/O
  return *this;
}

// Inputs a long int from the binary ifstream
inline binfstream& binfstream::operator>>(long int& d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  

  // Now read size number of bytes and copy into d's address
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      get(*p++);
    else
      get(*p--);
  
  
  // Return binistream for continuing I/O
  return *this;
}

// Inputs a float from the binary ifstream
inline binfstream& binfstream::operator>>(float& d)
{
  char* p;

  // Get a pointer
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;


  // Now read size number of bytes and copy into d's address
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      get(*p++);
    else
      get(*p--);


  // Return binistream for continuing I/O
  return *this;
}

// Inputs a double from the binary ifstream
inline binfstream& binfstream::operator>>(double& d)
{
  char* p;

  // Get a pointer to the double
  if(file_native)
    p = (char*)&d;
  else
    p = (char*)&d+sizeof(d)-1;
  

  // Now read size number of bytes and copy into d's address
  for(unsigned int i=0;i<sizeof(d);i++)
    if(file_native)
      get(*p++);
    else
      get(*p--);
  

  // Return binistream for continuing I/O
  return *this;
}

// Skips n bytes in stream
inline binfstream& binfstream::skip(int n)
{
  char c;
  
  while(n--)
    get(c);
  
  return *this;
}

// *** binofstream ***
inline binofstream::binofstream() : binfstream()
{}

inline binofstream::binofstream(const char* name) :
  binfstream(name)
  //  binfstream(name,ios::out)
{}

// Destructor 
inline binofstream::~binofstream()
{
}

inline void binofstream::open(const char* name)
{
  binfstream::open(name,ios::out);
}

// *** binifstream ***
inline binifstream::binifstream() : binfstream()
{}

//inline binifstream::binifstream(const char* name, int mode) :
//  binfstream(name,mode)
inline binifstream::binifstream(const char* name) :
  binfstream(name)
{}

// Destructor 
inline binifstream::~binifstream()
{
}

inline void binifstream::open(const char* name)
{
  binfstream::open(name,ios::in);
}

#endif

