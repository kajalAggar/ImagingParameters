// To estimate SNR values for IQP 3D reconstructed image

#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdint.h>
#include <cmath>
#include <algorithm>
#undef min
#undef max
#endif


using namespace std;


// Number of voxels (present in the uniform region) used to estimate SNR
const int Digi = 64480;
const int Dual = 64480;
const int IRIS = 64480;
const int Jack = 64480;
const int Axial = 64480;




float calculateSD(float data[])
{
    float sumSD = 0.0, meanSD, standardDeviation = 0.0;

    int i;

    for(i = 0; i < 64480; ++i)
    {
        sumSD += data[i];
    }

    meanSD = sumSD/64480;

    for(i = 0; i < 64480; ++i)
        standardDeviation += pow(data[i] - meanSD, 2);

    return sqrt(standardDeviation / 64480);
}

// Number of voxels in the 3D image matrix
// IRIS = 170; 170; 250
// Digi, dual = 100; 100; 275
// Axial, Jack = 100; 100; 270


void SNR(void)
{
  
 float Value;
 int z_pixels = 275;
 int y_pixels = 100;
 int x_pixels = 100;
 float VoxelSize = 0.4;
 
 float Rod1_radii = 2.5;
 float Rod2_radii = 2.0;
 float Rod3_radii = 1.5;
 float Rod4_radii = 1.0;
 float Rod5_radii = 0.5;
 
   
 
 
 std::string fileName("/biosave/kaggarwa/IRIS_12June/h_NORMALIZATION/DIGI/11nMore/11nMore_it1.img");
 



 // ------- Generating a 3D array as a function of z, y, x ----------------
   float*** array_img; 
   array_img = new float**[z_pixels];
   for (int z = 0; z < z_pixels; ++z)
   {
      array_img[z] = new float*[y_pixels];
      for(int y = 0; y < y_pixels; ++y)
      {
	  array_img[z][y] = new float[x_pixels];
	  for(int x = 0; x < x_pixels; ++x)
	  {
	      array_img[z][y][x] = 0.0;
	  }
      }
   }



// ------- Generating a 2D array as a function of y, x ----------------
   float array_sum[64480] = {}; 
      
   
      
   FILE* r = fopen(fileName.c_str(), "rb");

   for(int za=0; za<z_pixels; ++za){
    for(int ya=0; ya<y_pixels; ++ya){
      for(int xa=0; xa<x_pixels; ++xa){

        fread( &Value, sizeof(float), 1, r);
       
	array_img[za][ya][xa] = Value;
	
	} } }
	
	
   int count = 0;
   float sum = 0.0;
 
   // Identification of voxels lying in the ROI of the uniform region
   // voxels are considered on the basis of their centre position
   for(int zb=0; zb<z_pixels; ++zb)
   {
      for(int yb=0; yb<y_pixels; ++yb)
      {
	  for(int xb=0; xb<x_pixels; ++xb)
	  {
 
	      float z_voxel_cordinate = (zb * (VoxelSize/2.0)) + ((zb + 1) * (VoxelSize/2.0)) - ((z_pixels/2) * VoxelSize);
	      float y_voxel_cordinate = (yb * (VoxelSize/2.0)) + ((yb + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);
	      float x_voxel_cordinate = (xb * (VoxelSize/2.0)) + ((xb + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);
	      
	      float raddi_of_roi = 11.25;
	      float factor;
	      factor = (x_voxel_cordinate * x_voxel_cordinate) + (y_voxel_cordinate * y_voxel_cordinate) - (raddi_of_roi * raddi_of_roi);
 
	          
	      //if ((factor <= 0) && (zb <= 154) && (zb >= 130) )
	      if ((factor <= 0) && (z_voxel_cordinate <= 9.0) && (z_voxel_cordinate >= (-1.0)) )
	      {
		  //cout<<zb<<"	"<<yb<<"	"<<xb<<"	"<<array_img[zb][yb][xb]<<endl; 
		  sum += array_img[zb][yb][xb];
		  count++;
	          array_sum[count] += array_img[zb][yb][xb];
		
	      }
	

	
	
	  }
      }
  }
   
   
  
  float uniformity = sum/count;
  
  //--> cout<<"uniformity =		"<<uniformity<<endl;
  //--> cout<<"SNR = "<<uniformity/calculateSD(array_sum)<<endl;
// ----- Mean -------------- Sigma ------------------------- SNR -----
cout<<uniformity<<"	"<<calculateSD(array_sum)<<"	"<<uniformity/calculateSD(array_sum)<<endl;
} 
  
  
 
