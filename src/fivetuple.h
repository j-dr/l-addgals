#include <iostream>

using std::endl;
using std::cout;

class FiveTuple{
public:
  FiveTuple(){};
  FiveTuple(float tcm, float tfm, float tfs, float tp):
    cm(tcm),fm(tfm),fs(tfs),p(tp){ cs = 0.35; alpha = 0.0; };
    //cm(tcm),fm(tfm),fs(tfs),p(tp){cs = 0.6;};
  FiveTuple(float tcm, float tcs, float tfm, float tfs, float tp):
    cm(tcm),cs(tcs),fm(tfm),fs(tfs),p(tp){alpha = 0.0;};
  FiveTuple(float tcm, float tcs, float tfm, float tfs, float tp, float talpha):
    cm(tcm),cs(tcs),fm(tfm),fs(tfs),p(tp),alpha(talpha){};
  void Print(){cout<<cm<<" "<<cs<<" "<<fm<<" "<<fs<<" "<<p<<endl;};
  float operator[](int index) const{
    float return_me=0;
    assert( (index >= 0) && (index <6) ); // bounds checking
    if(index==0) return_me = cm;
    else if(index==1) return_me = cs;
    else if(index==2) return_me = fm;
    else if(index==3) return_me = fs;
    else if(index==4) return_me = p;
    else if(index==5) return_me = alpha;
    return return_me;
  }
  float LocalDens() const;
private:
  float cm;
  float cs;
  float fm;
  float fs;
  float p;
  float alpha;
};
