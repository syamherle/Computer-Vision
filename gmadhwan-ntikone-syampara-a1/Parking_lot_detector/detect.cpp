//
// detect.cpp : Detect cars in satellite images.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include<math.h>
#include <cmath>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <ctime>

int orientationFlag = 0;
double orientationVal;

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//


void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &cars)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &cars, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.
//
/*
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output(input.rows(), input.cols());

  // Convolution code here
  
  return output;
}
*/



//to handle border pixels.
int border(int x, int m)
{
    if(x < 0)
    {
        return -x - 1;
    }
    if(x >= m)
    {
        return 2*m - x - 1;
    }
    return x;
    
}

//Part 2.2
// Convolve an image with a separable convolution kernel
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
    SDoublePlane output(input.rows(), input.cols());
    SDoublePlane intermediate(input.rows(), input.cols());
    output=input;
    intermediate = input;		//to store the intermediate values for separable kernel
    
    int m = input.rows();
    int n = input.cols();
    int v = row_filter.rows()/2;
    int h = col_filter.cols()/2;
    int i = 0,j = 0,f = 0;
    double sum=0.0;
    /*Horizontal multiplication*/
    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            sum=0;
            for(f=-h; f<=h; f++)
            {
                int newIndex = border(j-f,n);
                sum+=col_filter[0][f+h] * input[i][newIndex];
            }
            intermediate[i][j] = std::abs(sum);		//store the value in intermediate
        }
    }
    
    /*Vertical Multiplication*/
    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            sum=0;
            for(f=-v; f<=v; f++)
            {
                int newIndex = border(i-f,m);
                sum+=row_filter[f+v][0] * intermediate[newIndex][j];
            }
            output[i][j] = std::abs(sum);
        }
    }
    
    // Convolution code here
    
    return output;
}

// Convolve an image with a  convolution kernel
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  // Convolution code here
   
    SDoublePlane output(input.rows(), input.cols());
    output = input;
    int M = input.rows();
    int N = input.cols();
    int h = filter.rows()/2;
    int l = filter.cols()/2;
    
    int sum;
    int x,y,j,k;
    int i;
    
    for( x = 0; x < M; x++)
    {
        for( y = 0; y < N; y++)
        {
            double sum = 0;
            for( j = -h; j<= h; j++)
            {
                for( k = -l; k <= l; k++)
                {
                    int x1 = border(x - j, M);
                    int y1 = border(y - k, N);
                    sum = sum + filter[j+h][k+l] * input[x1][y1];
                }
            }
            output[x][y] = std::abs(sum);
        }
    }
  return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
    int i,j;

    SDoublePlane output(input.rows(), input.cols());
    double sobelX[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
    double sobelY[3][3] = {{1,2,1}, {0,0,0}, {-1,-2,-1}};
    SDoublePlane sobelx(3,3);
    SDoublePlane sobely(3,3);
    if(_gx == 1)
    {
        for(i = 0; i<3; i++)
        {
            for(j = 0;j<3; j++)
            {
                sobelx[i][j]= sobelX[i][j]/8;
            }
        }
        
        output = convolve_general(input, sobelx);


    }
    else{
        for(i = 0; i<3; i++)
        {
            for(j = 0;j<3; j++)
            {
                sobely[i][j]= sobelY[i][j]/8;
            }
        }
        output = convolve_general(input, sobely);
    }
    

  // Implement a sobel gradient estimation filter with 1-d filters

  return output;
}
// Function to add two images after taking derivative in x and y direction.
SDoublePlane addTwo(SDoublePlane &input1, SDoublePlane &input2)
{
    SDoublePlane resulting(input2.rows(), input2.cols());
    
    int M = input1.rows();
    int N = input1.cols();
    int i,j;
    for(i = 0;i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            double diff = (input1[i][j] + input2[i][j]);
            resulting[i][j]  = diff;
        }
    }
    return resulting;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(SDoublePlane &input_gX, SDoublePlane &input_gY,  double thresh=0)
{
    int M = input_gX.rows();
    int N = input_gX.cols();
    int i,j,k, rowStrt, rowend, colStrt, colEnd;
    double pi = 3.14159265;
    
    SDoublePlane output(M, N);
    SDoublePlane gradientEdgesImg = addTwo(input_gX,input_gY);

    
    double magnitude[M][N];
    int orientation[M][N];
    float orient;
    float maxm = 75;
    int angle;
    float min = -180;
    int ival,jval;
    int maxcol;
    int maxrow;
    float maxMag = 0;
    int count0 = 0;
    int count90 = 0;
    
    for( i = 0; i<M; i++)
    {
        for(j = 0;j<N; j++)
        {
        magnitude[i][j] = sqrt(input_gX[i][j]*input_gX[i][j] + input_gY[i][j]*input_gY[i][j]);
        if(maxm > magnitude[i][j])
        {
            maxm = magnitude[i][j];
        }
        orient = atan2(input_gY[i][j],input_gX[i][j]) * 180.0/pi;
            
            if((orient > 45 && orient <= 135) || (orient > 225 && orient <= 315))
                angle = 90;
            if(orient <= 45 || orient > 315 || (orient > 135 && orient <= 225))
                angle = 0;
            orientation[i][j] = angle;
        }
    }
    for( i = 0; i<M; i++)
    {
        for(j = 0;j<N; j++)
        {
            if(orientation[i][j] ==0)
                count0++;
            else
                count90++;
        }
    }
    
    
    for(i = 0; i<M ; i++)
    {
        for(j = 0; j<N; j++)
        {
            if(magnitude[i][j] > 21) // if magnitude greater than 21 then consider neighboring pixels
            {
                if(orientation[i][j] == 0)
                {
                    maxcol = j;
                    if(i >= 2 && i <= M -2)
                    {
                        maxrow = i-2;
                        rowStrt = i-2;
                        rowend = i +1;
                    }
                    else if( i < 2)
                    {
                        maxrow = i;
                        rowStrt = i;
                        rowend = i +1;
                    }
                    else if(i >= M - 2)
                    {
                        maxrow = i-2;
                        rowStrt = i-2;
                        rowend = i;
                    }
                    else
                    {
                        maxrow = i-2;
                        rowStrt = i;
                        rowend = i+2;
                    }
                    maxMag = magnitude[maxrow][maxcol];
                    for(k=rowStrt; k<=rowend; k++)
                    {
                        if(maxMag < magnitude[k][maxcol])
                        {
                            gradientEdgesImg[maxrow][maxcol] = 0;
                            maxMag = magnitude[k][maxcol];
                            maxrow = k;
                        }
                    }
                    gradientEdgesImg[maxrow][maxcol]= 255;

                    
                    
                }
                else if(orientation[i][j] ==90)
                {
                
                    maxrow = i;
                    if(j >= 2 && j <= N -2)
                    {
                        maxcol = j-2;
                        colStrt = j-2;
                        colEnd = j +1;
                    }
                    else if( j < 2)
                    {
                        maxcol = j;
                        colStrt = j;
                        colEnd = j +1;
                    }
                    else if(j >= N - 2)
                    {
                        maxcol = j-2;
                        colStrt = j-2;
                        colEnd = j;
                    }
                    else
                    {
                        maxcol = j-2;
                        colStrt = i;
                        colEnd = i+2;
                    }
                    maxMag = magnitude[maxrow][maxcol];
                    for(k=colStrt; k<=colEnd; k++)
                    {
                        if(maxMag < magnitude[maxrow][k])
                        {
                            gradientEdgesImg[maxrow][maxcol] = 0;
                            maxMag = magnitude[maxrow][k];
                            maxcol = k;
                        }
                    }
                    gradientEdgesImg[maxrow][maxcol]= 255;
                    
                
                }
                
            }
            
            else
            {
                gradientEdgesImg[i][j] = 0;// if magnitude less than 21 then update grayscale value to 0
            }
        }
    }
    
    
    

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
    
    if(orientationFlag ==0){
        orientationVal =(double) count0/count90;
        orientationFlag = 1;
    }
    
  return gradientEdgesImg;
}



//function to create a gaussian Kernel.

SDoublePlane createGaussKernel(double (&gaussKernel)[5][5], double sigma)
{
    double div = 2*sigma*sigma;
    double sum = 0.0;
    double pi = 22/7;
    double first, sec;
    SDoublePlane gaussK(5,5);
    
    int i,j;
    
    for(i= -2;i<=2;i++)
    {
        for(j = -2; j <= 2; j++)
        {
            first = 1/(pi* div);
            sec = exp(-((i*i)+(j*j))/(div));
            gaussKernel[i+2][j+2] = first * sec;
            sum+=gaussKernel[i+2][j+2];
        }
    }
     
    
    for( i= 0;i<5;i++)
    {
        for( j = 0; j <5; j++)
        {
            gaussKernel[i][j] /=sum;
        }
    }
    for( i= 0;i<5;i++)
    {
        for( j = 0; j <5; j++)
        {
            gaussK[i][j] = gaussKernel[i][j];
        }
    }

    
    return gaussK;
    
}

/*
 
 This function is used to match the Image to the Template. If the image matches the template to a certain threshold, Then it is considered as a match.
 
 */

vector<DetectedBox> match_template(SDoublePlane img, SDoublePlane temp, vector<DetectedBox> cars)
{
    DetectedBox s;
    int i,j,x,y,k,l;
    double max = 0;
    double sum = 0;
    int maxrow, maxcol;
    int N = temp.rows()*temp.rows()* temp.cols()*temp.cols();
    int N1 = img.rows()* img.cols();
    int N2 = temp.rows()* temp.cols();

    int flag = 0;
    int count = 0;
    int countGlobal = 0;
    
    for(x = 0; x< temp.rows();x++)
        {
            for(y = 0; y<temp.cols();y++)
            {
                    if(temp[x][y] == 255)
                    {
                        countGlobal++;
                    }
            }}
    
    if(orientationVal > 0.7) // if greater than 0.7, then consider Vertical car template
    {
    for(i = 0; i< (img.rows()- 40);i++)
    {
        for(j = 5; j< (img.cols()- temp.cols());j++)
        {
            sum = 0;
            for(x = 5; x< 40-5;x++)
            {
                for(y = 5; y<temp.cols();y++)
                {
                    if(img[i+x][j+y-5] == 255 && temp[x][y] == 255)
                        sum+=1;
                }
            }
            if(sum > 22)
            {
                maxrow = i;
                maxcol = j+15;
                s.row = maxrow-10;
                s.col = maxcol;
                s.width = 18;
                s.height = 35;
                s.confidence = sum;
                cars.push_back(s);
                j = j+22;
                flag = 1;
                count++;
            }
        }
        if(flag == 1)
        {
            if(count > 20)
            {
            i+=53;
            }
            else
            {
                i+=29;
            }
            flag = 0;
        }
        count = 0;
    
    }
    }
    if(orientationVal <=0.7) // if less than 0.7, then consider Horizontal car template
    {
    
        for(i = 0; i< (img.cols() - temp.cols());i++)
        {
            for(j = 0; j< (img.rows()- temp.rows());j++)
            {
                sum = 0;
                for(x = 0; x< temp.cols()-5;x++)
                {
                    for(y = 0; y<temp.rows()-5;y++)
                    {
                        if(img[j+y][i+x+5] == 255 && temp[y][x] == 255)
                            sum+=1;
                    }
                }
                if(sum > 7)
                {
                    maxrow = j;
                    maxcol = i;
                    s.row = maxrow+10;
                    s.col = maxcol+10;
                    s.width = 40;
                    s.height = 25;
                    s.confidence = sum;
                    cars.push_back(s);
                    j = j+23;
                    flag = 1;
                    count++;
                }
            }
            if(flag == 1)
            {
                if(count > 8)
                {
                    i+=52;
                    
                }
                else
                {
                    i+=41;
                }
                flag = 0;
            }
            count = 0;
        }
    }
    return cars;
}

SDoublePlane getColFilter(double input[3][3])
{
    SDoublePlane colKernel(1,3);
    //col_filter
    for(int i=0; i<3; i++)
    {
        colKernel[0][i] = input[2][i]/3;
    }
    
    return colKernel;
}

SDoublePlane getRowFilter(double input[3][3])
{
    SDoublePlane rowKernel(3,1);
    //row_filter
    for(int i=0; i<3; i++)
    {
        rowKernel[i][0] = input[i][2]/3;
    }
    
    return rowKernel;
}

//
// This main file just outputs a few test images. You'll want to change it to do
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
    string tempName = "temp.png";
    SDoublePlane temp = SImageIO::read_png_file(tempName.c_str());
    SDoublePlane sobelxImg1  ;
    SDoublePlane sobelyImg1 ;
    SDoublePlane edgeImg1 ;
    string tempName2 = "temp2.png";
    SDoublePlane temp2 = SImageIO::read_png_file(tempName2.c_str());
    double sobelX[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
    double sobelY[3][3] = {{-1,-2,-1}, {0,0,0}, {1,2,1}};
 
    
    //SDoublePlane output_image = convolve_general(input_image, mean_filter);
  
    double gauss[5][5];
    SDoublePlane gaussKernel(5,5);
    
    double sigma1= 1.0;
    gaussKernel = createGaussKernel(gauss, sigma1);
    
    SDoublePlane gaussImg = convolve_general(input_image, gaussKernel);
    SDoublePlane sobelxImg  = sobel_gradient_filter(gaussImg,1);
    SDoublePlane sobelyImg  = sobel_gradient_filter(gaussImg,0);
	
	//uncomment following two lines to convolve using separable kernel
    /*SDoublePlane sobelxImg = convolve_separable(gaussImg, getRowFilter(sobelX), getColFilter(sobelX));
    SDoublePlane sobelyImg = convolve_separable(gaussImg, getRowFilter(sobelY), getColFilter(sobelY));*/

    SDoublePlane edgeImg = find_edges(sobelxImg,sobelyImg);
    
    if(orientationVal > 0.7)
    {
     sobelxImg1  = sobel_gradient_filter(temp,1);
     sobelyImg1  = sobel_gradient_filter(temp,0);
     edgeImg1 = find_edges(sobelxImg1,sobelyImg1);
    }
    else
    {
        sobelxImg1  = sobel_gradient_filter(temp2,1);
        sobelyImg1  = sobel_gradient_filter(temp2,0);
        edgeImg1 = find_edges(sobelxImg1,sobelyImg1);
    
    }
    vector<DetectedBox> cars;
    vector<DetectedBox> carBox = match_template(edgeImg, edgeImg1, cars);
  
  
  write_detection_txt("detected.txt", carBox);
  write_detection_image("detected.png", carBox, input_image);
   
    
    
}
