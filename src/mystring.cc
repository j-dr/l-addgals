
#include <string>
#include <strstream>

using namespace std;

//just like string, but truncates or adds zeros to the string
//so that it's of length size.

string MakeString(int i, int size)
{
  static const int StringBufferSize = (1 << 20);
  static char stringBuffer[StringBufferSize];
  ostrstream p(stringBuffer, size);
  p << i;
  string q = "";
  for(int i = 0; i< size-p.pcount();i++){
    q = "0"+q;
  }
  stringBuffer[p.pcount()] = 0;
  return q+string(stringBuffer);
}


