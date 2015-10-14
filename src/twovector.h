#include <vector>
using std::vector;

template <class T>
class twovector
{
 public:
  twovector():dimRow(0), dimCol(0){;}
  twovector(const twovector& rhs){
    dimRow = rhs.Rows();
    dimCol = rhs.Cols();
    V.resize(dimRow);
    for(int i=0;i < dimRow; i++){
      V[i].resize(dimCol);
    }
    for(int i=0;i<dimRow;i++){
      for(int j=0;j<dimRow;j++){
	V[i][j] = rhs[i][j];
      }
    }
  }
  twovector(int rows, int columns){
    dimRow = rows;
    dimCol = columns;
    V.resize(rows);
    for(int i=0;i < rows; i++){
      V[i].resize(columns);
    }
  }
  int Rows() const{
    return dimRow;
  }
  int Cols() const{
    return dimCol;
  }
  void Resize(int rows, int columns){
   dimRow = rows;
   dimCol = columns;
    V.resize(rows);
    for(int i=0;i < rows; i++){
      V[i].resize(columns);
    }
  }
  void GrowRow(int newSize) {
    if (newSize <= dimRow)
      return;
    dimRow = newSize;
    for(int i = 0 ; i < newSize - dimCol; i++) {
      vector<T> x(dimRow);
      V.push_back(x);
    }
  } 

  void GrowCol(int newSize) {
    if(newSize <= dimCol)
      return;
    dimCol = newSize;
    for (int i=0; i <dimRow; i++)
      V[i].resize(newSize);
  }

  vector<T>& operator[](int x) {
    return V[x];
  }
 private:
  vector <vector<T> > V;
  unsigned int dimRow;
  unsigned int dimCol;
};
