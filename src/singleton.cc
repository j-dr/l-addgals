#include "singleton.h"

Zmax* Zmax::_instance = 0;

Zmax* Zmax::Instance(){
  if (_instance == 0) {
    _instance = new Zmax;
  }
  return _instance;
}

Zmax::Zmax()
{
  _zmax = 0;
}
