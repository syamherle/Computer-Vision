#include <stdio.h>
#include <iostream>
#include "CImg.h"
#include<algorithm>
using namespace cimg_library;
using namespace std;
int main(){
  
  //CImg<float> temp_mat(1,2);
  CImg<float> A_mat(2,2);
  //CImg<float> x_mat(1,2);

  A_mat(0,0) = 1;
  A_mat(1,0) = 2;
  A_mat(0,1) = 3;
  A_mat(1,1) = 4;
  A_mat.invert();
	cout << "0,0:" << A_mat(0, 0) << endl;
	cout << "0,1:" << A_mat(0, 1) << endl;
	cout << "1,0:" << A_mat(1, 0) << endl;
	cout << "1,1:" << A_mat(1, 1) << endl;
  /*
  x_mat(0,0) = 1;
  x_mat(0,1) = 2;
  

  temp_mat = A_mat * x_mat;
  
  cout << "0,0" << temp_mat(0,0) << endl;
  cout << "0,1" << temp_mat(0,1) << endl;
*/
 /*
  int a[] = {3, 6, 8, 33};
  int x = 8;
  size_t myArraySize = sizeof(a) / sizeof(int);
  int *end = a + myArraySize;
// find the value 0:
  //int *result = std::find(a, end, 0);
  bool exists =  std::find(a, end, x) != end;
  cout << exists;
*/
  return 0;
}
 
