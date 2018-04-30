//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
// Syam Sundar Herle
// Neelam Tikone
// Gulshan Madhwani
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>
#include<cmath>
#include<vector>
#include<time.h>
using namespace std;

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
  
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
  
};

static double calculate_avg(vector<double> &v){
        double return_value = 0.0;
        int n = v.size();
       
        for ( int i=0; i < n; i++)
        {
            
            return_value += v[i];
        }
        
        return ( return_value / n);
};


static double calculate_var(vector<double> &v,double mean){
        double sum = 0.0;
        double temp =0.0;
        double var =0.0;
       
        for ( int j =0; j <= v.size()-1; j++)
        {
            
            temp = pow((v[j] - mean),2);
            sum += temp;
        }
        
        return var = sum;   
}

// Write this in Part 1.1
 
 //
 //going pixel by pixel get magnitude value by using this formulae
 //E(i,j)= log(√(R(i,j)^2 )+I(i,j)^2 ) 
 //and save it in a image buffer and return the image buffer used to save magnitude
 //
 static SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag){
     int i,j;
     SDoublePlane output_real=fft_real;
     for(i=0;i<=fft_real.rows()-1;i++){
         for(j=0;j<=fft_real.cols()-1;j++){
             output_real[i][j] = log(sqrt(pow(fft_real[i][j],2)+pow(fft_imag[i][j],2)));
             
         }
     }
    return  output_real;
 }


// Write this in Part 1.2
//
// To remove interfernce get input image, use FFT to get real and imaginay part
//in real and imaginary part leaving the central pixel value, check each pixel of real and imaginary part for this condition
//√(R^2 )+I^2>1 if this condition is satisfied. then this is a noise having high frequency, so change this value with neighbouring value
//
//

SDoublePlane remove_interference(const SDoublePlane &input){
      SDoublePlane input_img=input;
      SDoublePlane fft_real,fft_imag,output_real;
      fft(input_img, fft_real, fft_imag);
      
      for(int i=0;i<=((fft_real.rows()-1)/2)-20;i++){
          for(int j=0;j<=(fft_real.cols()-1)/2;j++){
              
              if(sqrt(pow(fft_real[i][j],2)+pow(fft_imag[i][j],2))>1.00){
                  fft_real[i][j] = fft_real[i-5][j-5];
                  fft_imag[i][j] = fft_imag[i-5][j-5];
              }
          }
      }
      
      
      for(int i=((fft_real.rows()-1)/2)+20;i<=(fft_real.rows()-1);i++){
          for(int j=(fft_real.cols()-1)/2;j<=(fft_real.rows()-1);j++){
              
              if(sqrt(pow(fft_real[i][j],2)+pow(fft_imag[i][j],2))>=1.00){
                  fft_real[i][j] = fft_real[i-10][j-10]; 
                  fft_imag[i][j] = fft_imag[i-10][j-10];
              }
          }
      }
      
      ifft(fft_real, fft_imag, output_real);
      
     return output_real;
     
}   

// Write this in Part 1.3 -- add watermark N to image

//
// For marking a image get binary vector for integer N of length magic_l and access magic_l number of distinctive point in real
//change the real part pixel value and simillarly access the symmetry of real part pixel value using teta+pi. To change the pixel value use
//the following formulae
//R(i,j)=R(i,j)+ alpha|R(i,j)| v_i
//
//Once done get idft and get watermarked image
//
//Here we have a circular radius of r(magic)=128 length l = 16 and alph = 5
//

SDoublePlane mark_image(const SDoublePlane &input, int N){
    

    SDoublePlane fft_real,fft_imag,output_real;

    //  magic parameters  
    int magic_l=16;
    
    int alpha=5;
    //    vector generator need to change it to random
    vector<unsigned>  bin;
    
    srand(N);
    for(int i=0;i<=magic_l-1;i++){
       
       bin.push_back(rand()%2);
    }
    
    
//    Call fft to have real and Imaginary part
    fft(input,fft_real,fft_imag);
    
//    center of the real matrix
    int c_x = fft_real.rows()/2;
    int c_y = fft_real.cols()/2;
    int r=128;
    for(int teta=0+(180/magic_l),i=0;teta < 180,i<=magic_l-1; teta+= (180/magic_l),i++){
        int x = c_x + (r*cos(teta));
        int y = c_y + (r*sin(teta));
        int opp_x= c_x + (r*cos(teta+M_PI));
        int opp_y = c_y + (r*sin(teta+M_PI));
        fft_real[x][y] = fft_real[x][y]+( alpha *(fft_real[x][y])*bin[i]);
        fft_real[opp_x][opp_y] = fft_real[opp_x][opp_y]+( alpha *(fft_real[opp_x][opp_y])*bin[i]);
        
    }
    
    ifft(fft_real, fft_imag, output_real);
    //SDoublePlane mag_output_real=fft_magnitude(fft_real,fft_imag);
    
    return output_real;
    
    
};

//
// For checking a image get binary vector for integer N of length magic_l and access magic_l number of distinctive point in real
//and extract them as vector makedBin. Find correlation coefficient using the following formulae
//cr= Σ (c_i-μ_c)(v_i-μ_v)/√var(c)var(v)
// if the correlation is greater than t(0.3)(magic) declare it as watermarked. and print correlation value 
//and as well as True or False denotingwater marked or not
//the following formulae

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N){
    
    SDoublePlane fft_real,fft_imag,output_real;

    //  magic parameters  
    int magic_l=16;
    
    int alpha=5;
    
//    vector generator need to change it to random
    vector<double>  bin;
    
    srand(N);
    for(int i=0;i<=magic_l-1;i++){
       
       bin.push_back(rand()%2);
    }
    
    
    vector<double> markedBin;
    
    
//    Call fft to have real and Imaginary part
    fft(input,fft_real,fft_imag);
    
//    center of the real matrix
    int c_x = fft_real.rows()/2;
    int c_y = fft_real.cols()/2;
    int r=128;
    for(int teta=0+(180/magic_l);teta < 180; teta+= (180/magic_l)){
        int x = c_x + (r*cos(teta));
        int y = c_y + (r*sin(teta));
		
        markedBin.push_back(abs(fft_real[x][y]));        
    }
    double normalise_l=r*M_PI;


    //cout<<markedBin.size()<<endl;
    double marked_avg = calculate_avg(markedBin);
    double bin_avg = calculate_avg(bin);
    
    double marked_var =calculate_var(markedBin,marked_avg);
    double bin_var =calculate_var(bin,bin_avg);
    double covariance_val=0.0;
   
    for(int i=0;i<=magic_l-1;i++){
        covariance_val += (markedBin[i]-marked_avg)*(bin[i]-bin_avg);
    }

    double correlation_coefficient = covariance_val/(sqrt(marked_var)*sqrt(bin_var));
	cout<< correlation_coefficient<<endl;
    if(correlation_coefficient > 0.30){cout<<"True"<<endl;}else{cout<<"False"<<endl;}
    
	return fft_real;
    
};




int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }
    
    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
    
    SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    
    if(part == "1.1")
      {
	// do something here!
        SDoublePlane input_fft_real,input_fft_imag;
        
        fft(input_image,input_fft_real, input_fft_imag);
        
        SDoublePlane output_real=fft_magnitude(input_fft_real,input_fft_imag);
        SImageIO::write_png_file(outputFile.c_str(),output_real,output_real,output_real);    
        
      }
    else if(part == "1.2")
      {
        // do something here!
        SDoublePlane output_real = remove_interference(input_image);
        SImageIO::write_png_file(outputFile.c_str(),output_real,output_real,output_real);    
      }
    else if(part == "1.3")
      {
        
	if(argc < 6)
	  {
	    cout << "Need 6 parameters for watermark part:" << endl;
	    cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    return -1;
	  }
	string op(argv[4]);
        int N = atoi(argv[5]);
        if(op == "add")
	  {
	    
            SDoublePlane output=mark_image(input_image,N);
            
            SImageIO::write_png_file(outputFile.c_str(),output,output,output);    
            // add watermark
            
           
        }
	else if(op == "check")
	  {
	    // check watermark
			SDoublePlane output = check_image(input_image, N);
            
	  }
	else
	  throw string("Bad operation!");
       
	
      }
    else
      throw string("Bad part!");

  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








