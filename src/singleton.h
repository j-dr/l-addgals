#ifndef singleton_h
#define singleton_h



class Zmax{
 public:
  static Zmax* Instance();
  float GetVal(){return _zmax;};
  void SetVal(float val){_zmax=val;};
 protected:
  Zmax();
 private:
  static Zmax* _instance;
  float _zmax; 
};


static Zmax *zmax = Zmax::Instance();
/*
#define NP 7
class ChainCov{
 public:
  static NRMat <double> rotation_matrix(NP,NP);  
  static NRVec <double> eigenvect(NP);
  const double StepScaleFactor = 2.0/(sqrt(NP));  
}
*/

#endif
