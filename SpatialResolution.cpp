// To estimate the spatial resolution (FWHM mm) for DigiPET and Dual layer PET following the NEMA standards


#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <stdint.h>
#include <algorithm>
#include <iterator>

#undef min
#undef max
#endif

using namespace std;

 
// Providing information about image space
// Number of 3D voxels in 3D space 
// Size of each side of voxel in mm.

const int z_pixels = 270;
const int y_pixels = 200;
const int x_pixels = 200;
const float VoxelSize = 0.2; 
const float zVoxelSize = 0.4;




std::string fileName("/biosave/kaggarwa/IRIS_12June/f.NEMA_sec3/DualLayer/FinalSimulations/1500_5/1500_5_it10.img");
cout<<"castor reconstructed -- 20"<<endl;




const int MAX2z = z_pixels;   
const int MAX2y = y_pixels;   
const int MAX2x = x_pixels;   


// Function to find the nearest value to a given value

float searchNearest1stz(float anArray[],float key, float t)
{
    float value = (key - anArray[0]); // input value - 1st element of the array
    float num1z = anArray[0]; 
    for(int x = 0;x < t; x++) 
    {
        if(value > abs(key - anArray[x]))
	{
	    value = (key - anArray[x]);
	    num1z = anArray[x];
	}
    }
    return num1z;
}


float searchNearest2ndz(float anArray[],float key, float t) 
{
    float value = (key - anArray[MAX2z-1]); // input value - 1st element of the array
    float num2z = anArray[MAX2z-1]; 
    for(int x = (MAX2z - 1); x > t; x--) // Going in opposite direction (final to mid)
    {
        if(value > abs(key - anArray[x]))
	{
	    value = (key - anArray[x]);
	    num2z = anArray[x];
	}
    }
    return num2z;
 }



float searchNearest1sty(float anArray[],float key, float t) 
{
    float value = (key - anArray[0]); 
    float num1y = anArray[0]; 
    for(int x = 0;x < t; x++) 
    {
        if(value > abs(key - anArray[x]))
	{
	    value = (key - anArray[x]);
	    num1y = anArray[x];
	}
    }
    return num1y;
}


 float searchNearest2ndy(float anArray[],float key, float t) 
 {
    float value = (key - anArray[MAX2y - 1]); 
    float num2y = anArray[MAX2y - 1]; 
    for(int x = (MAX2y -1); x > t; x--) 
    {
        if(value > abs(key - anArray[x]))
	{
	    value = (key - anArray[x]);
	    num2y = anArray[x];
	}
    }
    return num2y;
 }




float searchNearest1stx(float anArray[],float key, float t) 
{
    float value = (key - anArray[0]); 
    float num1x = anArray[0]; 
    for(int x = 0;x < t; x++) 
    {
        if(value > abs(key - anArray[x]))
	{
	    value = (key - anArray[x]);
	    num1x = anArray[x];
	}
    }
    return num1x;
}



float searchNearest2ndx(float anArray[],float key, float t) 
{
    float value = (key - anArray[MAX2x-1]); 
    float num2x = anArray[MAX2x-1]; 
    for(int x = (MAX2x -1); x > t; x--) 
    {
        if(value > abs(key - anArray[x]))
	{
	    value = (key - anArray[x]);
	    num2x = anArray[x];
	}
    }
    return num2x;
}



int SR_dual_2()
{
  
  TCanvas *c1 = new TCanvas("c1","1D Response function",900,700);
  gStyle->SetOptStat(0);
  c1->Divide(2,2,0,0);
  //c1->Divide(1,1,0,0);
  
  
  Double_t x_graph_Z[z_pixels], y_graph_Z[z_pixels];
  Double_t x_graph_Y[y_pixels], y_graph_Y[y_pixels];
  Double_t x_graph_X[x_pixels], y_graph_X[x_pixels];
  
  
  

  //----------------------------------------------------------------------------------
  //-------- Generating a 3D array as a function of z, y, x ----------------------------
  //----------------------------------------------------------------------------------
   
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



   
  // -----------------------------------------------------------------------------------
  // ------- Reading the image file and storing the value of each voxel ------------------  
  // -----------------------------------------------------------------------------------
  
    float VoxelContent;
    int x, y, z;
    FILE* r = fopen(fileName.c_str(), "rb");
    for(z=0; z<z_pixels; ++z)
    {
	for(y=0; y<y_pixels; ++y)
	{
	    for(x=0; x<x_pixels; ++x)
	    {
		fread( &VoxelContent, sizeof(float), 1, r);
		array_img[z][y][x] = VoxelContent;
	    } 
	}
    }
    
    


   // ------------------------------------------------------------------------------------
   // ------- Defining 1D z response function -----------------------------------------------
   // ------------------------------------------------------------------------------------
    float* array_zResponseFunction; 
    array_zResponseFunction = new float[z_pixels];
    for (int t = 0; t < z_pixels; ++t)
    {
	array_zResponseFunction[t] = 0.0;
    }
    
    float* array_yResponseFunction; 
    array_yResponseFunction = new float[y_pixels];
    for (int t = 0; t < y_pixels; ++t)
    {
	array_yResponseFunction[t] = 0.0;
    }
    
    float* array_xResponseFunction; 
    array_xResponseFunction = new float[x_pixels];
    for (int x = 0; x < x_pixels; ++x)
    {
	array_xResponseFunction[x] = 0.0;
    }
    
    
    
    // Add the bin content for all the voxels corresponding to same z, same y and same x in respective arrays
    for (int z = 0; z < z_pixels; ++z)
    {
	for (int y = 0; y < y_pixels; ++y)
	{
	    for (int x = 0; x < x_pixels; ++x)
	    {
	      array_zResponseFunction[z] = array_zResponseFunction[z] + array_img[z][y][x];
	      array_yResponseFunction[y] = array_yResponseFunction[y] + array_img[z][y][x];
	      array_xResponseFunction[x] = array_xResponseFunction[x] + array_img[z][y][x];
	    }
	}
    }
    
    
    
    for (int t = 0; t < z_pixels; ++t)
    {
        x_graph_Z[t] = t;
        y_graph_Z[t] = array_zResponseFunction[t];
    }
    
    TGraph* gr1 = new TGraph(z_pixels, x_graph_Z, y_graph_Z);
    c1->cd(1);
    gr1->SetTitle("z 1D response function");
    gr1->Draw("AL*");
   
   


   // ------------------------------------------------------------------------------------
   // ------- Defining 1D y response function -----------------------------------------------
   // ------------------------------------------------------------------------------------
    
    for (int t = 0; t < y_pixels; ++t)
    {
       //cout<<t<<"	"<<array_yResponseFunction[t]<<endl;
       x_graph_Y[t] = t;
       y_graph_Y[t] = array_yResponseFunction[t];
     }
   
    TGraph* gr2 = new TGraph(y_pixels, x_graph_Y, y_graph_Y);
    c1->cd(2);
    gr2->SetTitle("y 1D response function");
    gr2->Draw("AL*");
   
     

   // ------------------------------------------------------------------------------------
   // ------- Done defining 1D y response function ------------------------------------------
   // ------------------------------------------------------------------------------------
   
   



   // ------------------------------------------------------------------------------------
   // ------- Defining 1D x response function -----------------------------------------------
   // ------------------------------------------------------------------------------------
        
    for (int x = 0; x < x_pixels; ++x)
    {
       // cout<<x<<"	"<<array_xResponseFunction[x]<<endl;
       x_graph_X[x] = x;
       y_graph_X[x] = array_xResponseFunction[x];
     }
     
   TGraph* gr3 = new TGraph(x_pixels, x_graph_X, y_graph_X);
   c1->cd(3);
   gr3->SetTitle("x 1D response function");
   gr3->Draw("AL*");
   
   
   // ------------------------------------------------------------------------------------
   // ------- Done defining 1D x response function ------------------------------------------
   // ------------------------------------------------------------------------------------
   
   



   // ------------------------------------------------------------------------------------
   // Find the maximum value stored in the array_zResponseFunction[z]
   // ------------------------------------------------------------------------------------
   
    float maximum_z = *max_element(array_zResponseFunction, array_zResponseFunction + z_pixels); 
         
    int Max_z_voxel;
    for (int z = 0; z < z_pixels; ++z)
    {
	if (array_zResponseFunction[z] == maximum_z)
	{
	    Max_z_voxel = z; 
	    
	}
    }
    
    
 
    
   // ------------------------------------------------------------------------------------
   // Find the maximum value stored in the array_yResponseFunction[y]
   // ------------------------------------------------------------------------------------
   
    float maximum_y = *max_element(array_yResponseFunction, array_yResponseFunction + y_pixels); 
         
    int Max_y_voxel;
    for (int t = 0; t < y_pixels; ++t)
    {
	if (array_yResponseFunction[t] == maximum_y)
	{
	    Max_y_voxel = t; 
	    // cout<<"Maximum value is at = "<<Max_z_voxel<<endl;
	}
    }
    
   
    

    
   // ------------------------------------------------------------------------------------
   // Find the maximum value stored in the array_xResponseFunction[x]
   // ------------------------------------------------------------------------------------
   
    float maximum_x = *max_element(array_xResponseFunction, array_xResponseFunction + x_pixels); 
         
    int Max_x_voxel;
    for (int t = 0; t < x_pixels; ++t)
    {
	if (array_xResponseFunction[t] == maximum_x)
	{
	    Max_x_voxel = t; 
	    // cout<<"Maximum value is at = "<<Max_z_voxel<<endl;
	}
    }
    
    
    
    cout<<endl;
    cout<<" ***** Maximum value and its voxel ID for each response function ***** "<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"Maximum value:			"<<maximum_z<<"		|		"<<maximum_y<<"		|		"<<maximum_x<<endl;
    cout<<"Voxel n:			"<<Max_z_voxel<<"		|		"<<Max_y_voxel<<"		|		"<<Max_x_voxel<<endl;
    
    
    




    // -------------------------------- Parabola for z ------------------------------------
    // Neighbouring voxels of Max_z_voxel and their cordinates to begin with parabola thing
    // ------------------------------------------------------------------------------------
    
    int Neighbour_z1 = Max_z_voxel - 1;
    int Neighbour_z3 = Max_z_voxel + 1;
    
    float zx1, zx2, zx3;
    float zy1, zy2, zy3;
    
    zx1 = Neighbour_z1 - 0.5;	
    zx2 = Max_z_voxel - 0.5;	
    zx3 = Neighbour_z3 - 0.5;
        
    zy1 = array_zResponseFunction[Neighbour_z1]; 
    zy2 = maximum_z;
    zy3 = array_zResponseFunction[Neighbour_z3];
    
      
    
    // Need to fit parabola to determine the vertex or maximum value of Z
    double zdenom = (zx1 - zx2) * (zx1 - zx3) * (zx2 - zx3);
    double zA     = (zx3 * (zy2 - zy1) + zx2 * (zy1 - zy3) + zx1 * (zy3 - zy2)) / zdenom;
    double zB     = (zx3*zx3 * (zy1 - zy2) + zx2*zx2 * (zy3 - zy1) + zx1*zx1 * (zy2 - zy3)) / zdenom;
    double zC     = (zx2 * zx3 * (zx2 - zx3) * zy1 + zx3 * zx1 * (zx3 - zx1) * zy2 + zx1 * zx2 * (zx1 - zx2) * zy3) / zdenom;

    float zVertex_x = -zB / (2*zA);
    float zVertex_y = zC - zB*zB / (4*zA);
    
     


    
    // -------------------------------- Parabola for y ------------------------------------
    // Neighbouring voxels of Max_y_voxel and their cordinates to begin with parabola thing
    // ------------------------------------------------------------------------------------
    
    int Neighbour_y1 = Max_y_voxel - 1;
    int Neighbour_y3 = Max_y_voxel + 1;
    
    float yx1, yx2, yx3;
    float yy1, yy2, yy3;
    
    yx1 = Neighbour_y1 - 0.5;	
    yx2 = Max_y_voxel - 0.5;	
    yx3 = Neighbour_y3 - 0.5;
        
    yy1 = array_yResponseFunction[Neighbour_y1]; 
    yy2 = maximum_y;
    yy3 = array_yResponseFunction[Neighbour_y3];
    
    
    // Need to fit parabola to determine the vertex or maximum value of Y
    double ydenom = (yx1 - yx2) * (yx1 - yx3) * (yx2 - yx3);
    double yA     = (yx3 * (yy2 - yy1) + yx2 * (yy1 - yy3) + yx1 * (yy3 - yy2)) / ydenom;
    double yB     = (yx3*yx3 * (yy1 - yy2) + yx2*yx2 * (yy3 - yy1) + yx1*yx1 * (yy2 - yy3)) / ydenom;
    double yC     = (yx2 * yx3 * (yx2 - yx3) * yy1 + yx3 * yx1 * (yx3 - yx1) * yy2 + yx1 * yx2 * (yx1 - yx2) * yy3) / ydenom;

    float yVertex_x = -yB / (2*yA);
    float yVertex_y = yC - yB*yB / (4*yA);
    
   
    


    
    
    // -------------------------------- Parabola for x ------------------------------------
    // Neighbouring voxels of Max_y_voxel and their cordinates to begin with parabola thing
    // ------------------------------------------------------------------------------------
    
    int Neighbour_x1 = Max_x_voxel - 1;
    int Neighbour_x3 = Max_x_voxel + 1;
    
    float xx1, xx2, xx3;
    float xy1, xy2, xy3;
    
    xx1 = Neighbour_x1 - 0.5;	
    xx2 = Max_x_voxel - 0.5;	
    xx3 = Neighbour_x3 - 0.5;
        
    xy1 = array_xResponseFunction[Neighbour_x1]; 
    xy2 = maximum_x;
    xy3 = array_xResponseFunction[Neighbour_x3];
    
    
    
    // Need to fit parabola to determine the vertex or maximum value of X
    double xdenom = (xx1 - xx2) * (xx1 - xx3) * (xx2 - xx3);
    double xA     = (xx3 * (xy2 - xy1) + xx2 * (xy1 - xy3) + xx1 * (xy3 - xy2)) / xdenom;
    double xB     = (xx3*xx3 * (xy1 - xy2) + xx2*xx2 * (xy3 - xy1) + xx1*xx1 * (xy2 - xy3)) / xdenom;
    double xC     = (xx2 * xx3 * (xx2 - xx3) * xy1 + xx3 * xx1 * (xx3 - xx1) * xy2 + xx1 * xx2 * (xx1 - xx2) * xy3) / xdenom;

    float xVertex_x = -xB / (2*xA);
    float xVertex_y = xC - xB*xB / (4*xA);
    
     cout<<"	12	"<<endl;
    
    cout<<endl;
    cout<<endl;
    cout<<" ***** Coordinates of max point and two nearest neighbours for parabola *****"<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"				x  :  y		|		x  :  y		|		x  :  y		"<<endl;
    cout<<"1st point on pb: 	"<<zx1<<"  :  "<<zy1<<"	|	"<<yx1<<"  :  "<<yy1<<"	|	"<<xx1<<"  :  "<<xy1<<endl;
    cout<<"2nd point on pb: 	"<<zx2<<"  :  "<<zy2<<"	|	"<<yx2<<"  :  "<<yy2<<"	|	"<<xx2<<"  :  "<<xy2<<endl;
    cout<<"3rd point on pb: 	"<<zx3<<"  :  "<<zy3<<"	|	"<<yx3<<"  :  "<<yy3<<"	|	"<<xx3<<"  :  "<<xy3<<endl;
   
    
    cout<<endl;
    cout<<endl;
    cout<<" ***** Parabola Vertices estimated by fitting it on above three points ***** "<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"			x Vertex  :  y Vertex	|	x Vertex  :  y Vertex	|	x Vertex  :  y Vertex "<<endl;
    cout<<"Vertices:		   "<<zVertex_x<<"   :  "<<zVertex_y<<"	|	 "<<yVertex_x<<"   :  "<<yVertex_y<<"	|	 "<<xVertex_x<<"  :  "<<xVertex_y<<endl;
    
    
    cout<<endl;
    cout<<endl;
    cout<<" ***** Half the Maximum value or vertex of the parabola  ***** "<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"Half the maxima: 		"<<zVertex_y/2.0<<"		|		"<<yVertex_y/2.0<<"		|		"<<xVertex_y/2.0<<endl;
    
    
    
   
   
    
    // ------------------------------------------------------------------------------------
    // Finding the nearest value (in first half of the curve) to the half of the maximum value of z
    // ------------------------------------------------------------------------------------
      
    float znearest1;
    float znearest1_location1;
    znearest1 = searchNearest1stz(array_zResponseFunction, zVertex_y/2.0, zVertex_x);
    
    // Determining the location of 1st nearest value of half the maximum value
    for (int t = 0; t < (zVertex_x); ++t)
    {
	if (array_zResponseFunction[t] == znearest1)
	{
	    znearest1_location1 = t;
	    // cout<<"(greater) (1/2)maxima : voxel n\B0 "<<znearest1_location1<<" : "<<array_zResponseFunction[znearest1_location1]<<endl;
	}
    }
    
   
    
    // Smaller value for the first nearest point
    float znearest1_location2;
    
    if (znearest1 - (zVertex_y/2.0) <= 0.001)
    {
        znearest1_location2 = znearest1_location1 + 1;
    } else
    {
      znearest1_location2 = znearest1_location1 - 1;
    }
    
   
    if (array_zResponseFunction[znearest1_location1] <= zVertex_y/2.0 && array_zResponseFunction[znearest1_location2] <= zVertex_y/2.0)
    {
     znearest1_location1 = znearest1_location1 + 1;
     znearest1_location2 = znearest1_location2 + 1;  
    }
    
    if (array_zResponseFunction[znearest1_location1] <= zVertex_y/2.0 && array_zResponseFunction[znearest1_location2] <= zVertex_y/2.0)
    {
     znearest1_location1 = znearest1_location1 + 1;
     znearest1_location2 = znearest1_location2 + 1;  
    }
    
     if (array_zResponseFunction[znearest1_location1] >= zVertex_y/2.0 && array_zResponseFunction[znearest1_location2] >= zVertex_y/2.0)
    {
     znearest1_location1 = znearest1_location1 - 1;
     znearest1_location2 = znearest1_location2 - 1;  
    }
      
    // coordinates of neighbours of half the maximum value
    float zhalf1_x1, zhalf1_x2, zhalf1_y1, zhalf1_y2;
    zhalf1_x1 = znearest1_location1 - 0.5;
    zhalf1_x2 = znearest1_location2 - 0.5;
    zhalf1_y1 = array_zResponseFunction[znearest1_location1];
    zhalf1_y2 = array_zResponseFunction[znearest1_location2];
    
    
    
   
    
    
    
    
    
    // ------------------------------------------------------------------------------------
    // Finding the nearest value (in first half of the curve) to the half of the maximum value of y
    // ------------------------------------------------------------------------------------
      
    float ynearest1;
    float ynearest1_location1;
    ynearest1 = searchNearest1sty(array_yResponseFunction, yVertex_y/2.0, yVertex_x);
    
    // Determining the location of 1st nearest value of half the maximum value
    for (int t = 0; t < (yVertex_x); ++t)
    {
	if (array_yResponseFunction[t] == ynearest1)
	{
	    ynearest1_location1 = t;
	    // cout<<"(greater) (1/2)maxima : voxel n\B0 "<<znearest1_location1<<" : "<<array_zResponseFunction[znearest1_location1]<<endl;
	}
    }
    
    
    // Smaller value for the first nearest point
    float ynearest1_location2;
    
    if (ynearest1 - (yVertex_y/2.0) <= 0.001)
    {
        ynearest1_location2 = ynearest1_location1 + 1;
    } else
    {
      ynearest1_location2 = ynearest1_location1 - 1;
    }
	
	
    if (array_yResponseFunction[ynearest1_location1] <= yVertex_y/2.0 && array_yResponseFunction[ynearest1_location2] <= yVertex_y/2.0)
    {
     ynearest1_location1 = ynearest1_location1 + 1;
     ynearest1_location2 = ynearest1_location2 + 1;  
    }
    
    if (array_yResponseFunction[ynearest1_location1] <= yVertex_y/2.0 && array_yResponseFunction[ynearest1_location2] <= yVertex_y/2.0)
    {
     ynearest1_location1 = ynearest1_location1 + 1;
     ynearest1_location2 = ynearest1_location2 + 1;  
    }
	
    if (array_yResponseFunction[ynearest1_location1] >= yVertex_y/2.0 && array_yResponseFunction[ynearest1_location2] >= yVertex_y/2.0)
    {
     ynearest1_location1 = ynearest1_location1 - 1;
     ynearest1_location2 = ynearest1_location2 - 1;  
    }
    
    
    // coordinates of neighbours of half the maximum value
    float yhalf1_x1, yhalf1_x2, yhalf1_y1, yhalf1_y2;
    yhalf1_x1 = ynearest1_location1 - 0.5;
    yhalf1_x2 = ynearest1_location2 - 0.5;
    yhalf1_y1 = array_yResponseFunction[ynearest1_location1];
    yhalf1_y2 = array_yResponseFunction[ynearest1_location2];
    
    



    // ------------------------------------------------------------------------------------
    // Finding the nearest value (in first half of the curve) to the half of the maximum value of x
    // ------------------------------------------------------------------------------------
      
    float xnearest1;
    float xnearest1_location1;
    xnearest1 = searchNearest1stx(array_xResponseFunction, xVertex_y/2.0, xVertex_x);
    
    // Determining the location of 1st nearest value of half the maximum value
    for (int t = 0; t < (xVertex_x); ++t)
    {
	if (array_xResponseFunction[t] == xnearest1)
	{
	    xnearest1_location1 = t;
	    // cout<<"(greater) (1/2)maxima : voxel n\B0 "<<znearest1_location1<<" : "<<array_zResponseFunction[znearest1_location1]<<endl;
	}
    }
    
   
    
    // Smaller value for the first nearest point
    float xnearest1_location2;
    
    if (xnearest1 - (xVertex_y/2.0) <= 0.001)
    {
        xnearest1_location2 = xnearest1_location1 + 1;
    } else
    {
      xnearest1_location2 = xnearest1_location1 - 1;
    }
	
    if (array_xResponseFunction[xnearest1_location1] <= xVertex_y/2.0 && array_xResponseFunction[xnearest1_location2] <= xVertex_y/2.0)
    {
     xnearest1_location1 = xnearest1_location1 + 1;
     xnearest1_location2 = xnearest1_location2 + 1;  
    }
   
   
    if (array_xResponseFunction[xnearest1_location1] <= xVertex_y/2.0 && array_xResponseFunction[xnearest1_location2] <= xVertex_y/2.0)
    {
     xnearest1_location1 = xnearest1_location1 + 1;
     xnearest1_location2 = xnearest1_location2 + 1;  
    }
    
    
    if (array_xResponseFunction[xnearest1_location1] >= xVertex_y/2.0 && array_xResponseFunction[xnearest1_location2] >= xVertex_y/2.0)
    {
     xnearest1_location1 = xnearest1_location1 - 1;
     xnearest1_location2 = xnearest1_location2 - 1;  
    }
    
    
    
    
    // coordinates of neighbours of half the maximum value
    float xhalf1_x1, xhalf1_x2, xhalf1_y1, xhalf1_y2;
    xhalf1_x1 = xnearest1_location1 - 0.5;
    xhalf1_x2 = xnearest1_location2 - 0.5;
    xhalf1_y1 = array_xResponseFunction[xnearest1_location1];
    xhalf1_y2 = array_xResponseFunction[xnearest1_location2];
   
    
   


        
    
    // ------------------------------------------------------------------------------------
    //  Finding the point of intersection for Z by interpolating the above two found cordinates
    // ------------------------------------------------------------------------------------
          
    float zFWHM_x1;
    zFWHM_x1 = ((zVertex_y/2.0)- zhalf1_y1+ (zhalf1_y2-zhalf1_y1)/(zhalf1_x2-zhalf1_x1)*zhalf1_x1)/((zhalf1_y2 - zhalf1_y1)/(zhalf1_x2 - zhalf1_x1));
    
    


    // ------------------------------------------------------------------------------------
    //  Finding the point of intersection for Y by interpolating the above two found cordinates
    // ------------------------------------------------------------------------------------
          
    float yFWHM_x1;
    yFWHM_x1 = ((yVertex_y/2.0) - yhalf1_y1 + (yhalf1_y2 - yhalf1_y1)/(yhalf1_x2 - yhalf1_x1) * yhalf1_x1)/((yhalf1_y2 - yhalf1_y1)/(yhalf1_x2 - yhalf1_x1));
   


    
    // ------------------------------------------------------------------------------------
    //  Finding the point of intersection for X by interpolating the above two found cordinates
    // ------------------------------------------------------------------------------------
          
    float xFWHM_x1;
    xFWHM_x1 = ((xVertex_y/2.0) - xhalf1_y1 + (xhalf1_y2 - xhalf1_y1)/(xhalf1_x2 - xhalf1_x1) * xhalf1_x1)/((xhalf1_y2 - xhalf1_y1)/(xhalf1_x2 - xhalf1_x1));
   


    
    // ------------------------------------------------------------------------------------
    // --------------------- S U M M A R Y  on 1st FWHM point------------------------------------------ 
    // ------------------------------------------------------------------------------------
        
    cout<<endl;
    cout<<endl;
    cout<<" ***** Finding two nearest values to 1/2 maxima on LHS of the response function ***** "<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"				x  :  y		|		x  :  y		|		x  :  y"<<endl;
    cout<<"Greater point:		"<<znearest1_location1<<"  :  "<<array_zResponseFunction[znearest1_location1]<<"		|	"<<ynearest1_location1<<"  :  "<<array_yResponseFunction[ynearest1_location1]<<"		|	"<<xnearest1_location1<<"  :  "<<array_xResponseFunction[xnearest1_location1]<<endl;
    cout<<"Smaller point:		"<<znearest1_location2<<"  :  "<<array_zResponseFunction[znearest1_location2]<<"		|	"<<ynearest1_location2<<"  :  "<<array_yResponseFunction[ynearest1_location2]<<"		|	"<<xnearest1_location2<<"  :  "<<array_xResponseFunction[xnearest1_location2]<<endl;
    
    cout<<endl;
    cout<<endl;
    cout<<" ***** Taking the center cordinate of each voxel to do interpolation *****"<<endl;
    cout<<" ***** Coordinates for two points to do interpolation *****"<<endl;
    cout<<endl;
    cout<<"				x  :  y		|		x  :  y		|		x  :  y	"<<endl;
    cout<<"Greater coordinates:	"<<zhalf1_x1<<"  :  "<<zhalf1_y1<<"	|	"<<yhalf1_x1<<"  :  "<<yhalf1_y1<<"	|	"<<xhalf1_x1<<"  :  "<<xhalf1_y1<<endl;
    cout<<"Smaller coordinates:	"<<zhalf1_x2<<"  :  "<<zhalf1_y2<<"	|	"<<yhalf1_x2<<"  :  "<<yhalf1_y2<<"	|	"<<xhalf1_x2<<"  :  "<<xhalf1_y2<<endl;
    
    cout<<endl;
    cout<<endl;
    cout<<" ***** Coordinated calculated after interpolation *****"<<endl;
    cout<<endl;
    cout<<"FWHM 1st point:			"<<zFWHM_x1<<"		|		"<<yFWHM_x1<<"		|		"<<xFWHM_x1<<endl;
    cout<<endl;
   
    
    cout<<"**********************************************"<<endl;   
    
    // ------------------------------------------------------------------------------------
    // Finding the nearest value (in the second half of the curve) to the half of the maximum value
    // ------------------------------------------------------------------------------------
     
    // Z- Response function
    float znearest2;
    float znearest2_location1;
    znearest2 = searchNearest2ndz(array_zResponseFunction, zVertex_y/2.0, zVertex_x);
  
    // Determining the location of 2nd nearest value of half the maximum value
    for (int t = zVertex_x; t < z_pixels; ++t)
    {
	if (array_zResponseFunction[t] == znearest2)
	{
	    znearest2_location1 = t;
	    //cout<<"(greater) (1/2)maxima : voxel n\B0 "<<znearest2_location1<<" : "<<array_zResponseFunction[znearest2_location1]<<endl;
	}
    }
    
    
    // Smaller value for the second nearest point
    float znearest2_location2;
    
     if (znearest2 - (zVertex_y/2.0) <= 0.001)
    {
        znearest2_location2 = znearest2_location1 - 1;
    } else
    {
      znearest2_location2 = znearest2_location1 + 1;
    }
    
    
    
    
    if (array_zResponseFunction[znearest2_location1] <= zVertex_y/2.0 && array_zResponseFunction[znearest2_location2] <= zVertex_y/2.0)
    {
     znearest2_location1 = znearest2_location1 - 1;
     znearest2_location2 = znearest2_location2 - 1;  
    }
    
    if (array_zResponseFunction[znearest2_location1] <= zVertex_y/2.0 && array_zResponseFunction[znearest2_location2] <= zVertex_y/2.0)
    {
     znearest2_location1 = znearest2_location1 - 1;
     znearest2_location2 = znearest2_location2 - 1;  
    }
   
    if (array_zResponseFunction[znearest2_location1] >= zVertex_y/2.0 && array_zResponseFunction[znearest2_location2] >= zVertex_y/2.0)
    {
     znearest2_location1 = znearest2_location1 + 1;
     znearest2_location2 = znearest2_location2 + 1;  
    }
 
   
    // coordinates of neighbours of half the maximum value
    float zhalf2_x1, zhalf2_x2, zhalf2_y1, zhalf2_y2;
    zhalf2_x1 = znearest2_location1 - 0.5;
    zhalf2_x2 = znearest2_location2 - 0.5;
    zhalf2_y1 = array_zResponseFunction[znearest2_location1];
    zhalf2_y2 = array_zResponseFunction[znearest2_location2];
    
 
    
    
    
     // ----------------------
     // Y- Response function
    float ynearest2;
    float ynearest2_location1;
    ynearest2 = searchNearest2ndy(array_yResponseFunction, yVertex_y/2.0, yVertex_x);
 
    // Determining the location of 2nd nearest value of half the maximum value
    for (int t = yVertex_x; t < y_pixels; ++t)
    {
	if (array_yResponseFunction[t] == ynearest2)
	{
	    ynearest2_location1 = t;
	    //cout<<"(greater) (1/2)maxima : voxel n\B0 "<<znearest2_location1<<" : "<<array_zResponseFunction[znearest2_location1]<<endl;
	}
    }
  
  
    // Smaller value for the second nearest point
    float ynearest2_location2;
    
     if (ynearest2 - (yVertex_y/2.0) <= 0.001)
    {
        ynearest2_location2 = ynearest2_location1 - 1;
    } else
    {
      ynearest2_location2 = ynearest2_location1 + 1;
    }
    
    
    
    
     if (array_yResponseFunction[ynearest2_location1] <= yVertex_y/2.0 && array_yResponseFunction[ynearest2_location2] <= yVertex_y/2.0)
    {
     ynearest2_location1 = ynearest2_location1 - 1;
     ynearest2_location2 = ynearest2_location2 - 1;  
    }
    
    if (array_yResponseFunction[ynearest2_location1] <= yVertex_y/2.0 && array_yResponseFunction[ynearest2_location2] <= yVertex_y/2.0)
    {
     ynearest2_location1 = ynearest2_location1 - 1;
     ynearest2_location2 = ynearest2_location2 - 1;  
    }
    
    if (array_yResponseFunction[ynearest2_location1] >= yVertex_y/2.0 && array_yResponseFunction[ynearest2_location2] >= yVertex_y/2.0)
    {
     ynearest2_location1 = ynearest2_location1 + 1;
     ynearest2_location2 = ynearest2_location2 + 1;  
    }
    
    

    
    // coordinates of neighbours of half the maximum value
    float yhalf2_x1, yhalf2_x2, yhalf2_y1, yhalf2_y2;
    yhalf2_x1 = ynearest2_location1 - 0.5;
    yhalf2_x2 = ynearest2_location2 - 0.5;
    yhalf2_y1 = array_yResponseFunction[ynearest2_location1];
    yhalf2_y2 = array_yResponseFunction[ynearest2_location2];
    
 
    
    
     // ----------------------
     // X- Response function
    float xnearest2;
    float xnearest2_location1;
    xnearest2 = searchNearest2ndx(array_xResponseFunction, xVertex_y/2.0, xVertex_x);
 
    cout<<" ************************** "<<endl;
    cout<<xnearest2<<"	xVertex_x = "<<xVertex_x<<endl;
    cout<<" ************************** "<<endl;
    // Determining the location of 2nd nearest value of half the maximum value
    for (int t = xVertex_x; t < x_pixels; ++t)
    {
	if (array_xResponseFunction[t] == xnearest2)
	{
	    xnearest2_location1 = t;
	    //cout<<"(greater) (1/2)maxima : voxel n\B0 "<<znearest2_location1<<" : "<<array_zResponseFunction[znearest2_location1]<<endl;
	}
    }
    
    
    // Smaller value for the second nearest point
    float xnearest2_location2;
    
     if (xnearest2 - (xVertex_y/2.0) <= 0.001)
    {
        xnearest2_location2 = xnearest2_location1 - 1;
    } else
    {
      xnearest2_location2 = xnearest2_location1 + 1;
    }
    
  
  
    
    if (array_xResponseFunction[xnearest2_location1] <= xVertex_y/2.0 && array_xResponseFunction[xnearest2_location2] <= xVertex_y/2.0)
    {
     xnearest2_location1 = xnearest2_location1 - 1;
     xnearest2_location2 = xnearest2_location2 - 1;  
    }
   
   
    if (array_xResponseFunction[xnearest2_location1] <= xVertex_y/2.0 && array_xResponseFunction[xnearest2_location2] <= xVertex_y/2.0)
    {
     xnearest2_location1 = xnearest2_location1 - 1;
     xnearest2_location2 = xnearest2_location2 - 1;  
    }
    
    if (array_xResponseFunction[xnearest2_location1] >= xVertex_y/2.0 && array_xResponseFunction[xnearest2_location2] >= xVertex_y/2.0)
    {
     xnearest2_location1 = xnearest2_location1 + 1;
     xnearest2_location2 = xnearest2_location2 + 1;  
    }
    
    
  
  
  
  
    // coordinates of neighbours of half the maximum value
    float xhalf2_x1, xhalf2_x2, xhalf2_y1, xhalf2_y2;
    xhalf2_x1 = xnearest2_location1 - 0.5;
    xhalf2_x2 = xnearest2_location2 - 0.5;
    xhalf2_y1 = array_xResponseFunction[xnearest2_location1];
    xhalf2_y2 = array_xResponseFunction[xnearest2_location2];
    
    
    
    // ------------------------------------------------------------------------------------
    //  Finding the point of intersection for Z by interpolating the above two found cordinates
    // ------------------------------------------------------------------------------------
          
    float zFWHM_x2;
    zFWHM_x2 = ((zVertex_y/2.0)- zhalf2_y1+ (zhalf2_y2-zhalf2_y1)/(zhalf2_x2-zhalf2_x1)*zhalf2_x1)/((zhalf2_y2 - zhalf2_y1)/(zhalf2_x2 - zhalf2_x1));
   
    
    // ------------------------------------------------------------------------------------
    //  Finding the point of intersection for Y by interpolating the above two found cordinates
    // ------------------------------------------------------------------------------------
          
    float yFWHM_x2;
    yFWHM_x2 = ((yVertex_y/2.0)- yhalf2_y1+ (yhalf2_y2-yhalf2_y1)/(yhalf2_x2-yhalf2_x1)*yhalf2_x1)/((yhalf2_y2 - yhalf2_y1)/(yhalf2_x2 - yhalf2_x1));
    
    
    // ------------------------------------------------------------------------------------
    //  Finding the point of intersection for Y by interpolating the above two found cordinates
    // ------------------------------------------------------------------------------------
          
    float xFWHM_x2;
    xFWHM_x2 = ((xVertex_y/2.0)- xhalf2_y1+ (xhalf2_y2-xhalf2_y1)/(xhalf2_x2-xhalf2_x1)*xhalf2_x1)/((xhalf2_y2 - xhalf2_y1)/(xhalf2_x2 - xhalf2_x1));
   
   
    // ------------------------------------------------------------------------------------
    // --------------------- S U M M A R Y  on 2nd FWHM point------------------------------------------ 
    // ------------------------------------------------------------------------------------
        
    cout<<endl;
    cout<<endl;
    cout<<" ***** Finding two nearest values to 1/2 maxima on RHS on the response function ***** "<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"				x  :  y		|		x  :  y		|		x  :  y"<<endl;
    cout<<"Greater point:		"<<znearest2_location1<<"  :  "<<array_zResponseFunction[znearest2_location1]<<"		|	"<<ynearest2_location1<<"  :  "<<array_yResponseFunction[ynearest2_location1]<<"		|	"<<xnearest2_location1<<"  :  "<<array_xResponseFunction[xnearest2_location1]<<endl;
    cout<<"Smaller point:		"<<znearest2_location2<<"  :  "<<array_zResponseFunction[znearest2_location2]<<"		|	"<<ynearest2_location2<<"  :  "<<array_yResponseFunction[ynearest2_location2]<<"		|	"<<xnearest2_location2<<"  :  "<<array_xResponseFunction[xnearest2_location2]<<endl;
    cout<<endl;
    cout<<endl;
    cout<<" ***** Taking the center cordinate of each voxel to do interpolation *****"<<endl;
    cout<<" ***** Coordinates for two points to do interpolation *****"<<endl;
    cout<<endl;
    cout<<"				x  :  y		|		x  :  y		|		x  :  y	"<<endl;
    cout<<"Greater coordinates:	"<<zhalf2_x1<<"  :  "<<zhalf2_y1<<"	|		"<<yhalf2_x1<<"  :  "<<yhalf2_y1<<"	|		"<<xhalf2_x1<<"  :  "<<xhalf2_y1<<endl;
    cout<<"Smaller coordinates:	"<<zhalf2_x2<<"  :  "<<zhalf2_y2<<"	|		"<<yhalf2_x2<<"  :  "<<yhalf2_y2<<"	|		"<<xhalf2_x2<<"  :  "<<xhalf2_y2<<endl;
    cout<<endl;
    cout<<endl;
    cout<<" ***** Coordinated calculated after interpolation *****"<<endl;
    cout<<endl;
    cout<<"FWHM 2nd point:			"<<zFWHM_x2<<"		|		"<<yFWHM_x2<<"		|		"<<xFWHM_x2<<endl;
    
    
    cout<<"**********************************************"<<endl; 
    cout<<endl;
    cout<<endl;
    cout<<" ***** FWHM ***** "<<endl;
    cout<<endl;
    cout<<"			Z- Response function	|	Y- Response function	|	X- Response function"<<endl;
    cout<<"			"<<(zFWHM_x2 - zFWHM_x1) * zVoxelSize<<"		|		"<<(yFWHM_x2 - yFWHM_x1) * VoxelSize<<"		|		"<<(xFWHM_x2 - xFWHM_x1) * VoxelSize<<endl;
    
    
     
    return 0;
}




  
