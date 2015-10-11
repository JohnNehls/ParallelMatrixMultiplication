#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <complex>

using std::vector;
using std::complex;

// simple wrapper of a vector
template <class T>
class Matrix {
public:
  // create an empty matrix
  Matrix(int numrows, int numcols)
    :Nrow(numrows), Ncol(numcols), elements(Nrow*Ncol) {}

  // construct it from existing data
  Matrix(int numrows, int numcols, T* data)
    :Nrow(numrows), Ncol(numcols), elements(data, data+numrows*numcols) {}

  int rows() {return Nrow;}
  int cols() {return Ncol;}

  // access to elements, col is the fast axis
  T operator() (int row, int col) const {return elements[Ncol*row + col];}
  T& operator() (int row, int col) {return elements[Ncol*row + col];}

  
  // get raw pointer to elements[0]
  T* data() {return elements.data();}
  const vector<T>& elem() {return elements;}

private:
  int Nrow, Ncol;
  vector<T> elements;
};


#endif // MATRIX_H_
