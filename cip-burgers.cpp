#include<iostream>
#include<fstream>   //file I/O
#include<vector>    //vector container
#include<algorithm> //swap()
#include<iomanip>   //setw(), setprecision()
#include<sstream> //stringstream
#include<cmath>

#define pi 4.0*atan(1.0)

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void output(int steps, vector<float> x, vector<float> fx)
{
  ofstream fp;
  stringstream buf;
  string filenumber;

  buf<<setfill('0');
  filenumber = to_string(steps);
  buf<<setw(5)<<filenumber;

  string filename="f"+buf.str()+".csv";
  fp.open(filename, ios::out);
  fp<<"x, fx\n";
  for(int i=0; i<x.size(); ++i){
    fp<<x[i]<<", "
      <<fx[i]<<"\n";
  }
  fp.close();
  buf.str(string()); //clear buffer
}

int main()
{
  float dx, dt;
  float cx=1.0;
  vector<float> x, fx, gx, fxn, gxn;

  
  int i, imax;
  int iter, itermax;
  int steps;

  cout<<"CIP method for Burger's equation."<<endl;
  cout<<"Enter imax: "<<endl;
  cin>>imax;

  dx = (2.0*pi)/(imax-1);
  cout<<"Enter dt (dx= "<<dx<<" ) "<<endl;
  cin>>dt;

  cout<<"Enter maximum iteration: "<<endl;
  cin>>itermax;
  cout<<"Enter steps for record: "<<endl;
  cin>>steps;

  //  cout<<sgn(steps)<<endl;
  
  x.resize(imax);
  fx.resize(imax);
  gx.resize(imax);
  fxn.resize(imax);
  gxn.resize(imax);

  for(i=0; i<imax; ++i){
    x[i] = i*dx; 
  }

  
  for(i=0; i<imax; ++i){
    fx[i] = sin(x[i]);
    fxn[i] = fx[i];  //fx[i];
    gxn[i] = 0.0;
  }  
  

  
  for(i=1; i<imax-1; ++i){
    gx[i] = (fx[i]-fx[i-1])/dx;
  }
  
  iter = 0;
  int iup;
  float fdif,xam1,xbm1,xx;
  float up;
  float alpha =1.0;
  cout<<"Enter alpha(0< , <=1.0), control parameter:"<<endl;
  cin>>alpha;
  float B;
  float a,b,c,d;
  do{

    for(i=1; i<imax-1; ++i){
      up = -sgn(fx[i]); // df/dt+fdf/dx=0; 
      iup = i + (int)up;
      xx = -fx[i]*dt;
      fdif = (fx[iup]-fx[i])/dx*up;

      B = ( abs( (fdif-gx[i])/(gx[iup]-fdif +1.0e-6) ) - 1.0 )/(dx*up);

      a = ( gx[i]-fdif + (gx[iup]-fdif)*(1.0 + alpha*B*(dx*up) ) )/(dx*dx);

      b = alpha*B*fdif + (fdif-gx[i])/(dx*up) -a*(dx*up);

      c = gx[i] + alpha*B*fx[i];
      

      fxn[i] = ( ( (a*xx + b)*xx + c )*xx + fx[i] )/(1.0 + alpha*B*xx);

      gxn[i] =( (1.0 + alpha*B*xx)*((3.0*a*xx + 2.0*b)*xx + c)
		- ( ((a*xx + b)*xx + c )*xx + fx[i] )*(alpha*B) )
	/( (1.0 + alpha*B*xx)*(1.0 + alpha*B*xx) );

      //update non-advective term
      gxn[i] = gxn[i] - dt*(gxn[i]*gxn[i]);
      
    }

    if(iter%steps == 0) output(iter,x,fxn);
    iter +=1;
    swap(fx,fxn);
    swap(gx,gxn);
    
  }while(iter<itermax);
  
  
}
