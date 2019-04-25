// To estimate the RC for each rod of IQP



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


void RecoveryCoefficient(void)
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
 

 
 std::string fileName("/biosave/kaggarwa/IRIS_12June/e.DualLayer_simul/dual_imgs_100,100,275_moreit/dual_imgs_100,100,275_moreit_it39.img");
 
 
 



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
   float** array_sum; 
   array_sum = new float*[y_pixels];
   for (int y = 0; y < y_pixels; ++y)
   {
		array_sum[y] = new float[x_pixels];
		for(int x = 0; x < x_pixels; ++x)
		{
			array_sum[y][x] = 0.0;
		}
    }   
   


   float max;
   float max2, max3, max4, max5;
   
  
   
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
	      }
	

	
	
	  }
      }
  }
   
  cout<<sum<<"	"<<count<<endl; 
  float uniformity = sum/count;
  //cout<<"uniformity =		"<<uniformity<<endl;
  


 
// ----------------------------------------------------------- 
// ------- R C part begins here ------------------------------- 
// -----------------------------------------------------------	
  
FILE* r = fopen(fileName.c_str(), "rb");
	for(int za=0; za<z_pixels; ++za)
	{
		for(int ya=0; ya<y_pixels; ++ya)
		{
			for(int xa=0; xa<x_pixels; ++xa)
			{
				fread( &Value, sizeof(float), 1, r);
				array_img[za][ya][xa] = Value;
				float z_voxel_cordinate = (za * (VoxelSize/2.0)) + ((za + 1) * (VoxelSize/2.0)) - ((z_pixels/2) * VoxelSize);
				float y_voxel_cordinate = (ya * (VoxelSize/2.0)) + ((ya + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);
				float x_voxel_cordinate = (xa * (VoxelSize/2.0)) + ((xa + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);
				
				float Rod1; // for radii = 2.5 mm rod
			   	Rod1 = ((x_voxel_cordinate - 2.16) * (x_voxel_cordinate - 2.16)) + ((y_voxel_cordinate + 6.66)  * (y_voxel_cordinate + 6.66)) - (Rod1_radii * Rod1_radii);
				
				if ((Rod1 <= 0) && (z_voxel_cordinate <= (-8.5)) && (z_voxel_cordinate >= (-18.5)) ) // this za limits will remains constant
				{
					//cout<<za<<"	"<<ya<<"	"<<xa<<"	"<<array_img[za][ya][xa]<<endl;
					
					
					double temp = array_sum[ya][xa];
					temp = temp + Value;
					//cout<<"temp "<<temp<<endl;
					if ( temp > max )
					{
						max = temp;
					}
					array_sum[ya][xa] = temp;
				}
				
				float Rod2; // for radii = 2 mm rod
			   	Rod2 = ((x_voxel_cordinate + 5.66) * (x_voxel_cordinate + 5.66)) + ((y_voxel_cordinate + 4.11)  * (y_voxel_cordinate + 4.11)) - (Rod2_radii * Rod2_radii);
				
				if ((Rod2 <= 0) && (z_voxel_cordinate <= (-8.5)) && (z_voxel_cordinate >= (-18.5)) ) // this za limits will remains constant
				{
					//cout<<za<<"	"<<ya<<"	"<<xa<<"	"<<array_img[za][ya][xa]<<endl;
					
					
					double temp2 = array_sum[ya][xa];
					temp2 = temp2 + Value;
					//cout<<"temp "<<temp<<endl;
					if ( temp2 > max2 )
					{
						max2 = temp2;
					}
					array_sum[ya][xa] = temp2;
				}
				
				float Rod3; // for radii = 1.5 mm rod
			   	Rod3 = ((x_voxel_cordinate + 5.66) * (x_voxel_cordinate + 5.66)) + ((y_voxel_cordinate - 4.11)  * (y_voxel_cordinate - 4.11)) - (Rod3_radii * Rod3_radii);
				
				if ((Rod3 <= 0) && (z_voxel_cordinate <= (-8.5)) && (z_voxel_cordinate >= (-18.5)) ) // this za limits will remains constant
				{
					//cout<<za<<"	"<<ya<<"	"<<xa<<"	"<<array_img[za][ya][xa]<<endl;
					
					
					double temp3 = array_sum[ya][xa];
					temp3 = temp3 + Value;
					//cout<<"temp "<<temp<<endl;
					if ( temp3 > max3 )
					{
						max3 = temp3;
					}
					array_sum[ya][xa] = temp3;
				}
				
				float Rod4; // for radii = 1.0 mm rod
			   	Rod4 = ((x_voxel_cordinate - 2.16) * (x_voxel_cordinate - 2.16)) + ((y_voxel_cordinate - 6.66)  * (y_voxel_cordinate - 6.66)) - (Rod4_radii * Rod4_radii);
				
				if ((Rod4 <= 0) && (z_voxel_cordinate <= (-8.5)) && (z_voxel_cordinate >= (-18.5)) ) // this za limits will remains constant
				{
					//cout<<za<<"	"<<ya<<"	"<<xa<<"	"<<array_img[za][ya][xa]<<endl;
					
					
					double temp4 = array_sum[ya][xa];
					temp4 = temp4 + Value;
					//cout<<"temp "<<temp<<endl;
					if ( temp4 > max4 )
					{
						max4 = temp4;
					}
					array_sum[ya][xa] = temp4;
				}
				
				float Rod5; // for radii = 0.5 mm rod
			   	Rod5 = ((x_voxel_cordinate - 7.0) * (x_voxel_cordinate - 7.0)) + ((y_voxel_cordinate - 0.0)  * (y_voxel_cordinate - 0.0)) - (Rod5_radii * Rod5_radii);
				
				if ((Rod5 <= 0) && (z_voxel_cordinate <= (-8.5)) && (z_voxel_cordinate >= (-18.5)) ) // this za limits will remains constant
				{
					//cout<<za<<"	"<<ya<<"	"<<xa<<"	"<<array_img[za][ya][xa]<<endl;
					
					
					double temp5 = array_sum[ya][xa];
					temp5 = temp5 + Value;
					//cout<<"temp "<<temp<<endl;
					if ( temp5 > max5 )
					{
						max5 = temp5;
					}
					array_sum[ya][xa] = temp5;
				}
				
			} 
		} 
	}
	
	 //cout<<"max for rod1 =		"<< max <<endl;
	 //cout<<"max for rod2 =		"<< max2 <<endl;
	 //cout<<"max for rod3 =		"<< max3 <<endl;
	 //cout<<"max for rod4 =		"<< max4 <<endl;
	 //cout<<"max for rod5 =		"<< max5 <<endl;
	
	int maximum_y_voxel;
	int maximum_x_voxel;
	
	int maximum_y_voxel2, maximum_y_voxel3, maximum_y_voxel4, maximum_y_voxel5;
	int maximum_x_voxel2, maximum_x_voxel3, maximum_x_voxel4, maximum_x_voxel5;
	
// After summation, determining which x and y voxel summation has the max value
	for( int yb=0; yb < y_pixels; yb++)
	{
		for( int xb=0; xb < x_pixels; xb++)
		{
		  if (max == (array_sum[yb][xb]))
		  {
			//cout << "Sum id maximum for voxel y = " << yb << " and voxel x = " << xb << " : " << array_sum[yb][xb] << endl;
			maximum_y_voxel = yb;
			maximum_x_voxel = xb;
		  }  
		
		  if (max2 == (array_sum[yb][xb]))
		  {
			//cout << "Sum id maximum for voxel y = " << yb << " and voxel x = " << xb << " : " << array_sum[yb][xb] << endl;
			maximum_y_voxel2 = yb;
			maximum_x_voxel2 = xb;
		  } 
		  
		  if (max3 == (array_sum[yb][xb]))
		  {
			//cout << "Sum id maximum for voxel y = " << yb << " and voxel x = " << xb << " : " << array_sum[yb][xb] << endl;
			maximum_y_voxel3 = yb;
			maximum_x_voxel3 = xb;
		  } 
		  
		  if (max4 == (array_sum[yb][xb]))
		  {
			//cout << "Sum id maximum for voxel y = " << yb << " and voxel x = " << xb << " : " << array_sum[yb][xb] << endl;
			maximum_y_voxel4 = yb;
			maximum_x_voxel4 = xb;
		  } 
		  
		  if (max5 == (array_sum[yb][xb]))
		  {
			//cout << "Sum id maximum for voxel y = " << yb << " and voxel x = " << xb << " : " << array_sum[yb][xb] << endl;
			maximum_y_voxel5 = yb;
			maximum_x_voxel5 = xb;
		  } 
		}
	}
	
	


float x_5dia = (maximum_x_voxel * (VoxelSize/2.0)) + ((maximum_x_voxel + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);
float x_4dia = (maximum_x_voxel2 * (VoxelSize/2.0)) + ((maximum_x_voxel2 + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);
float x_3dia = (maximum_x_voxel3 * (VoxelSize/2.0)) + ((maximum_x_voxel3 + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);
float x_2dia = (maximum_x_voxel4 * (VoxelSize/2.0)) + ((maximum_x_voxel4 + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);
float x_1dia = (maximum_x_voxel5 * (VoxelSize/2.0)) + ((maximum_x_voxel5 + 1) * (VoxelSize/2.0)) - ((x_pixels/2) * VoxelSize);

float y_5dia = (maximum_y_voxel * (VoxelSize/2.0)) + ((maximum_y_voxel + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);
float y_4dia = (maximum_y_voxel2 * (VoxelSize/2.0)) + ((maximum_y_voxel2 + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);
float y_3dia = (maximum_y_voxel3 * (VoxelSize/2.0)) + ((maximum_y_voxel3 + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);
float y_2dia = (maximum_y_voxel4 * (VoxelSize/2.0)) + ((maximum_y_voxel4 + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);
float y_1dia = (maximum_y_voxel5 * (VoxelSize/2.0)) + ((maximum_y_voxel5 + 1) * (VoxelSize/2.0)) - ((y_pixels/2) * VoxelSize);


cout<<endl;

cout<<"			max voxels	 ; 	voxel centre"<<endl;
cout<<"			  x:y	 	 ; 	x	:	  y"<<endl;
cout<<"5 mm dia rod :		"<<maximum_x_voxel<<" : "<<maximum_y_voxel<<"		 ;  	"<< x_5dia<<"	:	"<<y_5dia<<endl;
cout<<"4 mm dia rod :		"<<maximum_x_voxel2<<" : "<<maximum_y_voxel2<<"		 ;	"<< x_4dia<<"	:	"<<y_4dia<<endl;
cout<<"3 mm dia rod :		"<<maximum_x_voxel3<<" : "<<maximum_y_voxel3<<"		 ;	"<< x_3dia<<"	:	"<<y_3dia<<endl;
cout<<"2 mm dia rod :		"<<maximum_x_voxel4<<" : "<<maximum_y_voxel4<<"	 ;  	"<< x_2dia<<"	:	"<<y_2dia<<endl;
cout<<"1 mm dia rod :		"<<maximum_x_voxel5<<" : "<<maximum_y_voxel5<<"	 ;	"<< x_1dia<<"	:	"<<y_1dia<<endl;

int count_RC = 0;
float sum_RC = 0.0;
  
int count_RC2 = 0;
float sum_RC2 = 0.0;

int count_RC3 = 0;
float sum_RC3 = 0.0;

int count_RC4 = 0;
float sum_RC4 = 0.0;

int count_RC5 = 0;
float sum_RC5 = 0.0;



cout<<((-15.4)/VoxelSize) - 0.5 + (z_pixels/2.0)<<endl;
cout<<((-5.4)/VoxelSize) - 0.5 + (z_pixels/2.0)<<endl;

int a = ((-15.4)/VoxelSize) - 0.5 + (z_pixels/2.0);
int b = ((-5.4)/VoxelSize) - 0.5 + (z_pixels/2.0);

 for(int zb=(a-1); zb<(b+1); ++zb)
  //for(int zb=(-15.4); zb<(-5.3); zb+=0.4)   
  {
      for(int yb=0; yb<y_pixels; ++yb)
      {
	  for(int xb=0; xb<x_pixels; ++xb)
	  {
	      // assume that x = 101 and y = 84 have the maxmimu value, then what?
	      if (zb >= a && zb <= b && yb == maximum_y_voxel && xb == maximum_x_voxel)
	      {
		      //cout<<zb<<"	"<<yb<<"	"<<xb<<"	"<<(array_img[zb][yb][xb])/uniformity<<endl;
		      sum_RC += ((array_img[zb][yb][xb])/uniformity);
		      count_RC++;
	      }
	     
	      if (zb >= a && zb <= b && yb == maximum_y_voxel2 && xb == maximum_x_voxel2)
	      {
		      //cout<<zb<<"	"<<yb<<"	"<<xb<<"	"<<(array_img[zb][yb][xb])/uniformity<<endl;
		      sum_RC2 += ((array_img[zb][yb][xb])/uniformity);
		      count_RC2++;
		
	      }	
	      
	      if (zb >= a && zb <= b && yb == maximum_y_voxel3 && xb == maximum_x_voxel3)
	      {
		      //cout<<zb<<"	"<<yb<<"	"<<xb<<"	"<<(array_img[zb][yb][xb])/uniformity<<endl;
		      sum_RC3 += ((array_img[zb][yb][xb])/uniformity);
		      count_RC3++;
		
	      }	
	      
	      if (zb >= a && zb <= b && yb == maximum_y_voxel4 && xb == maximum_x_voxel4)
	      {
		      //cout<<zb<<"	"<<yb<<"	"<<xb<<"	"<<(array_img[zb][yb][xb])/uniformity<<endl;
		      sum_RC4 += ((array_img[zb][yb][xb])/uniformity);
		      count_RC4++;
		
	      }	
	      
	      if (zb >= a && zb <= b && yb == maximum_y_voxel5 && xb == maximum_x_voxel5)
	      {
		      //cout<<zb<<"	"<<yb<<"	"<<xb<<"	"<<(array_img[zb][yb][xb])/uniformity<<endl;
		      sum_RC5 += ((array_img[zb][yb][xb])/uniformity);
		      count_RC5++;
		
	      }	
	      
	  }
      }

  }   
 
 cout<<endl;
 
  
 //cout<<sum_RC<<"	"<<count_RC<<endl; 
 //cout<<sum_RC2<<"	"<<count_RC2<<endl; 
 //cout<<sum_RC3<<"	"<<count_RC3<<endl; 
 //cout<<sum_RC4<<"	"<<count_RC4<<endl; 
 //cout<<sum_RC5<<"	"<<count_RC5<<endl; 


  cout<<"RC for 5 dia rod =	"<<sum_RC/count_RC<<endl; 
  cout<<"RC for 4 dia rod =	"<<sum_RC2/count_RC2<<endl; 
  cout<<"RC for 3 dia rod =	"<<sum_RC3/count_RC3<<endl; 
  cout<<"RC for 2 dia rod =	"<<sum_RC4/count_RC4<<endl;
  cout<<"RC for 1 dia rod =	"<<sum_RC5/count_RC5<<endl; 
  cout<<endl;
} 
  
  
 
