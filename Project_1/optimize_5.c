#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#define im_rows         288  //image size 
#define im_collumns     352
#define padded_rows     292  // padded image size
#define padded_collumns 356  // added two zeros on each side (for the case where N=5)

#define N 5     // filter class
#define STD 0.88  // standard deviation  
#define K 3     // constant
#define gauss_const 501187233627271,5

#define filename      "image_in.yuv"  //input file name define
#define finalfilename "image_out5.yuv" //output file name define

// Global Variables
#pragma arm section zidata="ram"
double masktemp[N][N];
int tmp_img1[N][padded_collumns];
double tmp_img2[N][im_collumns];	
#pragma arm section
int image[im_rows ][im_collumns]; //the image we are filtering
int padded_image[padded_rows][padded_collumns]={0}; // padded image
double g_flipped[N][N]; // Gauss filter mask/kernell NXN size
double conv_out[im_rows ][im_collumns]; //output of convolution
//conv_out is the same size as the initial image 

// end initialization

//Image Read
void read()  //image read
{
  int i,j;  //indexes //error when initialized inside for 
  FILE *frame_c;
  if((frame_c=fopen(filename,"rb"))==NULL)
  {
    printf("current frame doesn't exist\n");
    exit(-1);
  }

  for(i=0;i<im_rows;i++)  // extracts lines
  {
    for(j=0;j<im_collumns;j++) //extracts collumns
    {
      image[i][j] =fgetc(frame_c);
    }
  }
  fclose(frame_c);
}
// end read

// Gauss Filter Kernell
void gauss(){
  int i,j;  //indexes //error when initialized inside for
  int x,y;      //needed for gauss meshgrid
  double c;        //gauss filter coefficient
  double g[N][N];  // Gauss filter mask/kernell NXN size

  x =-N/2;
    c =exp(N*N/(STD*STD));

    for(i=0;i<N;i++){
        y =-N/2;
        for(j=0;j<N;j++){ //we avoided using the pow function, in order to make it
            g[i][j] =c*exp(-(x*x+y*y)/(2*STD*STD)); // more efficient/faster
            g[i][j] = g[i][j]/gauss_const;
            y++;
        }
        x++;
    }
    for(i=0; i<N; i++){
      for( j=0; j<N; j++){
        g_flipped[i][j] =g[N-i-1][N-j-1];
      }
    }
   
}//end gauss


// convolution 2D
void convolution2D(){
  int i,j; //indexes //error when initialized inside for
    int center =N/2; //automata kovei dekadika dld floor
    int left =center-1;
    int top =center-1;
    int x,y;

    //padding
    for( i=top+1;i<=im_rows+top;i++){
        for( j=1+left;j<=im_collumns+left;j=j+8){
           padded_image[i][j]=     image[i-top-1][j-left-1];
            padded_image[i][j+1]= image[i-top-1][j+1-left-1];
            padded_image[i][j+2]= image[i-top-1][j+2-left-1];
            padded_image[i][j+3]= image[i-top-1][j+3-left-1];
            padded_image[i][j+4]= image[i-top-1][j+4-left-1];
            padded_image[i][j+5]= image[i-top-1][j+5-left-1];
            padded_image[i][j+6]= image[i-top-1][j+6-left-1];
            padded_image[i][j+7]= image[i-top-1][j+7-left-1]; 
        }  
    }
    //padding  x,y define the center pixel
    for(x=0;x<im_rows;x++){//more efficient to read per row
        for( y=0;y<im_collumns;y++){  //kernell is applied on each pixel
            for(i=0;i<N;i++){ //center of the kernel is centered on the pixel we are reading
                
                    conv_out[x][y]=conv_out[x][y]+( padded_image[i+x][y]   *g_flipped[i][0] );
                    conv_out[x][y]=conv_out[x][y]+( padded_image[i+x][y+1] *g_flipped[i][1]);
                    conv_out[x][y]=conv_out[x][y]+( padded_image[i+x][y+2] *g_flipped[i][2]);
                    conv_out[x][y]=conv_out[x][y]+( padded_image[i+x][y+3] *g_flipped[i][3]);
                    conv_out[x][y]=conv_out[x][y]+( padded_image[i+x][y+4] *g_flipped[i][4]);
                //mask is applied(multiplication) on the surrounding pixels of the centered one
            } // then the results of each mult. are summed together
        } 
    }      
}    
//end  convolution 2D

// final steps of algorithm

void final_steps(){
  int i,j;  //indexes //error when initialized inside for
// image values range from 0 to 255 thus we need clipping
  for(i=0;i<im_rows;i++){
        for(j=0;j<im_collumns;j=j+2){ // Imask=I-Is , Iunsharp=I+K*Imask
            image[i][j] =image[i][j] + K*(image[i][j]-conv_out[i][j]);
            if(image[i][j]>255){image[i][j] =255;}// maximum pixel intensity is 255
            else if(image[i][j]<0){image[i][j] =0;}// minimum pixel intensity is 0
            image[i][j+1] =image[i][j+1] + K*(image[i][j+1]-conv_out[i][j+1]);
            if(image[i][j+1]>255){image[i][j+1] =255;}// maximum pixel intensity is 255
            else if(image[i][j+1]<0){image[i][j+1] =0;}// minimum pixel intensity is 0
        }
    }
}
//end final_steps


//write file    
void write(){
  int i,j;  //indexes //error when initialized inside for
  FILE *frame_y;
  frame_y=fopen(finalfilename,"wb");

  for(i=0;i<im_rows ;i++)
  {
    for(j=0;j<im_collumns;j++)
    {
      fputc(image[i][j],frame_y);
    }
  }
  fclose(frame_y);
}

// end write    

void main(){
  read();
  gauss();
  convolution2D();
  final_steps();
  write();
} 
