#include<iostream>
#include<cmath>

using namespace std;

/*type-safe version*/
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main()
{
  float a;
  cout<<"Enter a"<<endl;
  cin>>a;
  
  float c = 1.0;
  int iup = (a > 0) - (a < 0); //C++ style for sign() in Fortran
  cout<<iup<<endl;

  iup = sgn(a);
  cout<<iup<<endl;
  
  
}
