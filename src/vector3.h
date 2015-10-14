#ifndef vector3_h
#define vector3_h

#include <vector>
#include <iostream>
#include <fstream>


using namespace std;

template <class R> class vector3;

template <class T>  
class vector3 {
 public:
  vector< vector < vector <T> > >v;
  vector<vector <T> > & operator[] (int i){return v[i];};
  vector<vector <T> > const & operator[] (int i)const{return v[i];};
  size_t Xdim(){return xdim;};
  size_t Ydim(){return ydim;};
  size_t Zdim(){return zdim;};

 private:
  size_t xdim;
  size_t ydim;
  size_t zdim;
 public:
  template<class R> operator vector3<R>();
  vector3(){
    xdim=0;
    ydim=0;
    zdim=0;
  };
  vector3(int i){ Init(i);};
  vector3(int i, int j, int k){Init(i,j,k);};
  void Init(int i){Init(i,i,i);};
  void Init(int i, int j, int k){
    xdim=i; ydim=j; zdim=k;
    v.resize(i); 
    for(unsigned int i=0; i<v.size();i++){
      v[i].resize(j);
      for(unsigned int j=0; j<v.size();j++){
	v[i][j].resize(k);
      }
    } 
  };
  T Sum()const{
    T sum = 0;
    for(unsigned int i=0; i<v.size();i++)
      for(unsigned int j=0; j<v[i].size();j++)
	for(unsigned int k=0; k<v[i][j].size();k++)
	  sum += v[i][j][k];
    return sum;
  };
  double Ave()const{ return 1.0*Sum()/Size();};
  double AveFilled()const{
    return 1.0*Sum()/FilledElements();
  };
  int FilledElements()const{
    int nonzero = 0;
    for(unsigned int i=0; i<v.size();i++)
      for(unsigned int j=0; j<v[i].size();j++)
	for(unsigned int k=0; k<v[i][j].size();k++){
	  if(v[i][j][k]>0)
	    nonzero++;
	}
    return nonzero;
  };
  int Size()const{ 
    return xdim*ydim*zdim;
  };

  void Write(ofstream &file)const{
    for(unsigned int i=0; i<xdim;i++)
      for(unsigned int j=0; j<ydim;j++)
	for(unsigned int k=0; k<zdim;k++)
	  file<<v[i][j][k]<<"\n";
  };

  vector3 <T> KeepCells(int i1, int i2, int j1, int j2, int k1, int k2){
    vector3 <T> newvec3(i2-i1, j2-j1, k2-k1);
    //    cout<<newvec3.size()<<" "<<newvec3[newvec3.size()-1].size()<<" "<<newvec3[newvec3.size()-1][newvec3.size()-1].size()<<endl;
    unsigned int newi=0;
    for(unsigned int i=i1; i<i2;i++){
      unsigned int newj = 0;
      for(unsigned int j=j1; j<j2;j++){
	unsigned int newk = 0;
	for(unsigned int k=k1; k<k2;k++){
	  //	cout<<i<<" "<<j<<" "<<k<<" "<<newi<<" "<<newj<<" "<<newk<<endl;
	  newvec3[newi][newj][newk] = v[i][j][k];
	  newk++;
	}
	newj++;
      }
      newi++;
    }
    return newvec3;
  };



 // vector3 scalar operators (element-wise operations)
  vector3 operator+(T) const;
  //const vector3 operator-(T) const;
  //const vector3 operator*(T) const;
  //const vector3 operator/(T) const;

  //friend const vector3 operator+(T, const vector3&);
  //friend const vector3 operator-(T, const vector3&);
  //friend const vector3 operator*(T, const vector3&);
  //friend const vector3 operator/(T, const vector3&);

  vector3& operator=(T);
  vector3& operator+=(T);
  vector3& operator-=(T);
  vector3& operator*=(T);
  vector3& operator/=(T);
};
  // Implicit conversion operator between different vector3<>'s
template<class T> template<class R> 
inline vector3<T>::operator vector3<R>()
    {
      vector3<R> w(xdim,ydim,zdim);
      for(unsigned int i=0; i<xdim;i++)
	for(unsigned int j=0; j<ydim;j++)
	  for(unsigned int k=0; k<zdim;k++)
	    w[i][j][k] = static_cast<R>(v[i][j][k]);
      
      return w;
    }
template<class T> inline vector3<T> vector3<T>::operator+(T d)const
{
  vector3<T> w;
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	w[i][j][k] = d+v[i][j][k];

  return w;
}



template<class T> inline vector3<T>& vector3<T>::operator=(T d)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	v[i][j][k] = d;

  return *this;
}

template<class T> inline vector3<T>& vector3<T>::operator+=(T d)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	v[i][j][k] += d;

  return *this;
}

template<class T> inline vector3<T>& vector3<T>::operator-=(T d)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	v[i][j][k] -= d;

  return *this;
}

template<class T> inline vector3<T>& vector3<T>::operator*=(T d)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	v[i][j][k] *= d;

  return *this;
}
template<class T> inline vector3<T>& vector3<T>::operator/=(T d)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	v[i][j][k] /= d;

  return *this;
}


// Adds vec3 m to this
/* 

template<class T> inline vector3<T>& vector3<T>::operator++(vector3<T>& m)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	m[i][j][k]++;

  return m;
}


vector3<T>& vector3<T>::operator+=(double a){
    for(unsigned int i=0; i<xdim;i++)
      for(unsigned int j=0; j<ydim;j++)
	for(unsigned int k=0; k<zdim;k++)
	  v[i][j][k] += a;
    
    return *this;
  }

// Assigns constant value to vec3
template<class T> inline vec3<T>& vec3<T>::operator=(T d)
{
  c[0]=d;
  c[1]=d;
  c[2]=d;

  return *this;
}

// Adds d to all elements of vec3
template<class T> inline vec3<T>& vec3<T>::operator+=(T d)
{
  c[0]+=d;
  c[1]+=d;
  c[2]+=d;

  return *this;
}

template<class T> inline vector3<T>& vector3<T>::operator+=(int a)
{
  for(unsigned int i=0; i<xdim;i++)
    for(unsigned int j=0; j<ydim;j++)
      for(unsigned int k=0; k<zdim;k++)
	v[i][j][k] += a;

  return *this;
}

*/


#endif
