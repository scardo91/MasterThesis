#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>


 
#define BLOCK_SIZE 16
#define MAX_IMG_PEL 255
#define min(a,b) ((a<b)?a:b)
#define max(a,b) ((a>b)?a:b)
#define absm(a) ((a)<(0) ? (-(a)):(a))
#define log2(a) (log((float) a)/log(2.0))

 int** mvvet;
				
void get_block_masks(unsigned char **pict_dvc, int x_pos, int y_pos, 
	int width, int height, int **block, int block_w, int block_h)
{

  int dx, dy;
  int x, y;
  int i, j;
  int maxold_x,maxold_y;
  int result;
  int pres_x;
  int pres_y; 
  int tmp_res[4][9];
  static const int COEF[6] = {    1, -5, 20, 20, -5, 1  };

  dx = x_pos&3;
  dy = y_pos&3;
  x_pos = (x_pos-dx)/4;
  y_pos = (y_pos-dy)/4;

  maxold_x = width-1;
  maxold_y = height-1;


  if (dx == 0 && dy == 0) 
  {  /* fullpel position */
    for (j = 0; j < block_h; j++)
      for (i = 0; i < block_w; i++)
        block[j][i] = pict_dvc[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
  }
  else 
  { /* other positions */

    if (dy == 0) 
    { /* No vertical interpolation */

      for (j = 0; j < block_h; j++) 
      {
        for (i = 0; i < block_w; i++) 
        {
          for (result = 0, x = -2; x < 4; x++)
            result += pict_dvc[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[j][i] = max(0, min(MAX_IMG_PEL, (result+16)/32));
        }
      }

      if ((dx&1) == 1) 
      {
        for (j = 0; j < block_h; j++)
          for (i = 0; i < block_w; i++)
            block[j][i] = (block[j][i] + pict_dvc[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+dx/2))] +1 )/2;
      }
    }
    else if (dx == 0) 
    {  /* No horizontal interpolation */

      for (j = 0; j < block_h; j++) 
      {
        for (i = 0; i < block_w; i++) 
        {
          for (result = 0, y = -2; y < 4; y++)
            result += pict_dvc[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
          block[j][i] = max(0, min(MAX_IMG_PEL, (result+16)/32));
        }
      }

      if ((dy&1) == 1) 
      {
        for (j = 0; j < block_h; j++)
          for (i = 0; i < block_w; i++)
           block[j][i] = (block[j][i] + pict_dvc[max(0,min(maxold_y,y_pos+j+dy/2))][max(0,min(maxold_x,x_pos+i))] +1 )/2;
      }
    }
    else if (dx == 2) 
    {  /* Vertical & horizontal interpolation */

      for (j = -2; j < (block_h+3); j++) 
      {
        for (i = 0; i < block_w; i++)
          for (tmp_res[j][i+2] = 0, x = -2; x < 4; x++)
            tmp_res[j][i+2] += pict_dvc[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
      }

      for (j = 0; j < block_h; j++) 
      {
        for (i = 0; i < block_w; i++) 
        {
          for (result = 0, y = -2; y < 4; y++)
            result += tmp_res[j][i+y+2]*COEF[y+2];
          block[j][i] = max(0, min(MAX_IMG_PEL, (result+512)/1024));
        } 
      }

      if ((dy&1) == 1)
      {
        for (j = 0; j < block_h; j++)
          for (i = 0; i < block_w; i++)
            block[j][i] = (block[j][i] + max(0, min(MAX_IMG_PEL, (tmp_res[j][i+2+dy/2]+16)/32)) +1 )/2;
      }
    }
    else if (dy == 2)
    {  /* Horizontal & vertical interpolation */

      for (j = 0; j < block_h; j++)
      {
        for (i = -2; i < (block_w+3); i++)
          for (tmp_res[j][i+2] = 0, y = -2; y < 4; y++)
            tmp_res[j][i+2] += pict_dvc[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
      }

      for (j = 0; j < block_h; j++)
      {
        for (i = 0; i < block_w; i++)
        {
          for (result = 0, x = -2; x < 4; x++)
            result += tmp_res[j][i+x+2]*COEF[x+2];
          block[j][i] = max(0, min(MAX_IMG_PEL, (result+512)/1024));
        }
      }

      if ((dx&1) == 1)
      {
        for (j = 0; j < block_h; j++)
          for (i = 0; i < block_w; i++)



            block[j][i] = (block[j][i] + max(0, min(MAX_IMG_PEL, (tmp_res[j][i+2+dx/2]+16)/32))+1)/2;
      }
    }
    else
    {  /* Diagonal interpolation */

      for (j = 0; j < block_h; j++)
      {
        for (i = 0; i < block_w; i++)
        {
          pres_y = dy == 1 ? y_pos+j : y_pos+j+1;
          pres_y = max(0,min(maxold_y,pres_y));
          for (result = 0, x = -2; x < 4; x++)
            result += pict_dvc[pres_y][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[j][i] = max(0, min(MAX_IMG_PEL, (result+16)/32));
        }
      }

      for (j = 0; j < block_h; j++)
      {
        for (i = 0; i < block_w; i++)
        {
          pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
          pres_x = max(0,min(maxold_x,pres_x));
          for (result = 0, y = -2; y < 4; y++)
            result += pict_dvc[max(0,min(maxold_y,y_pos+j+y))][pres_x]*COEF[y+2];
          block[j][i] = (block[j][i] + max(0, min(MAX_IMG_PEL, (result+16)/32)) +1 ) / 2;
        }
      }

    }
  }
}
 
 void find_mv(unsigned char** image, int **block, int** pred, int block_w, int block_h,
	int width, int height, int limdx, int limsx, int limup, int limdown, int mv[2])

{

int x, y, pointer, dx, dy, di, dj, best_dx, best_dy, c;
float dist_block,best_dist_block;
			
best_dist_block=1000000.0;
for (dy=(limup*4);dy<(limdown*4);dy+=4)
for (dx=(limsx*4);dx<(limdx*4);dx+=4)
	{	
	get_block_masks(image,dx,dy,width,height,pred,block_w,block_h);
				
	dist_block=0.0;
	c=0;
	for (di=0;di<block_w;di++)
	for (dj=0;dj<block_h;dj++)
			{
			dist_block=dist_block+
				absm(pred[dj][di]-block[dj][di]);
			c++;
			}
	
	dist_block=dist_block/c;
				
	/*
	dist_block=dist_block+5.0*((mv[0]-(dx/4))*(mv[0]-(dx/4)));
	dist_block=dist_block+5.0*((mv[1]-(dy/4))*(mv[1]-(dy/4)));
	*/
				
	if (dist_block<best_dist_block)
		{
		best_dx=dx;
		best_dy=dy;
		best_dist_block=dist_block;
		}
				
	}

get_block_masks(image,best_dx,best_dy,width,height,pred,block_w,block_h);

best_dx=best_dx/4;
best_dy=best_dy/4;


mv[0]=best_dx;
mv[1]=best_dy;
	
}

void compute_pred(unsigned char** image, unsigned char** image2, unsigned char** image3, unsigned short** masks, int width, int height, int Nelem)

{

int x, y, xb, yb, c, pointer, dx, dy, di, dj, best_dx, best_dy;
int limdx, limsx, limup, limdown, block_w, block_h, cnt;
int mv[2];
float dist_block,best_dist_block;
int **block, **bmasks, **pred;		

block_w=block_h=4;

block=(int**) malloc(sizeof(int*)*block_h);
for (y=0;y<block_h;y++)
	block[y]=(int*) malloc(sizeof(int)*block_w);
    
pred=(int**) malloc(sizeof(int*)*block_h);
for (y=0;y<block_h;y++)
	pred[y]=(int*) malloc(sizeof(int)*block_w);

cnt=0;
for (yb=0;yb<height;yb+=4)
for (xb=0;xb<width;xb+=4,cnt++)
if ((masks[yb][xb]<65535)|0)
{
    
    limsx=xb;
    limup=yb;
    
    mv[0]=limsx;
    mv[1]=limup;

    for (y=0;y<block_h;y++)
        for (x=0;x<block_w;x++)
            block[y][x]=image[limup+y][limsx+x];

    printf("Find best prediction for superpixel %d @ ( %d x %d )\n",c,limsx,limup);
    
    limdown=limup+32;
    limup=limup-32;
    limdx=limsx+32;
    limsx=limsx-32;
    
    limup=(limup<0)?0:limup;
    limdown=(limdown>(height-1))?(height-1):limdown;
    limsx=(limsx<0)?0:limsx;
    limdx=(limdx>(width-1))?(width-1):limdx;
    
    printf("Find best prediction for superpixel %d ( %d x %d ) in window [ %d , %d ] x [ %d , %d ]\n",
        c,block_h,block_w,limup,limdown,limsx,limdx);
    
    find_mv(image2,block,pred,block_w,block_h,width,height,limdx,limsx,limup,limdown,mv);

    printf("Best prediction for superpixel MV [ %d , %d ]\n",mv[0],mv[1]);
    
    mvvet[cnt][0]=mv[0]-xb;
    mvvet[cnt][1]=mv[1]-yb;

    for (y=0;y<block_h;y++)
        for (x=0;x<block_w;x++)
            image3[yb+y][xb+x]=pred[y][x];
    
}

   for (y=0;y<block_h;y++)
        free(block[y]);
    free(block);
    
    for (y=0;y<block_h;y++)
        free(pred[y]);
    free(pred);

	
}


/*****************************************************************/
/*!
 *  \fn int main (int argc, char** argv)
 *  \brief Main function for the whole program 
 *
 * \bug None (already known)
 *
 *  \author <a href="mailto: sim1mil@dei.unipd.it"> Simone Milani
 *                    (sim1mil@dei.unipd.it) </a>
 *  \date   2004-05-12
 *
 * \param argc : number of input parameter
 * \param argv : input parameter string
 * \return void
 *
******************************************************************/
 int main (int argc, char** argv){
 
 int i, j, width, height;
 int  N;
 unsigned char** image;
 unsigned char** image2;
 unsigned char** image3;
 unsigned short** masks;
 FILE *fp_cur, *fp_ref, *fp_masks, *fp_out;



 if (argc != 7) {
 	fprintf(stderr,"Wrong number of input parameters \n");
 	fprintf(stderr,"Usage: <width> <height> <cur> <ref> <masks> <out>\n");
 	exit(1);
 	}
 
 width=atoi(argv[1]);
 height=atoi(argv[2]);
 fp_cur=fopen(argv[3],"r");
 fp_ref=fopen(argv[4],"r");
 fp_masks=fopen(argv[5],"r");
 fp_out=fopen(argv[6],"w");
 N=(width*height)>>4;
 
  image=(unsigned char**) malloc(sizeof(unsigned char*)*height);
  for (j=0;j<height;j++)
      image[j]=(unsigned char*) malloc(sizeof(unsigned char)*width);
  
  image2=(unsigned char**) malloc(sizeof(unsigned char*)*height);
  for (j=0;j<height;j++)
      image2[j]=(unsigned char*) malloc(sizeof(unsigned char)*width);
  
  
  image3=(unsigned char**) malloc(sizeof(unsigned char*)*height);
  for (j=0;j<height;j++)
      image3[j]=(unsigned char*) malloc(sizeof(unsigned char)*width);
  
  mvvet=(int**) malloc(sizeof(int*)*N);
  for (j=0;j<N;j++)
      mvvet[j]=(int*) malloc(sizeof(int)*2);
  
  masks=(unsigned short**) malloc(sizeof(unsigned short*)*height);
  for (j=0;j<height;j++)
      masks[j]=(unsigned short*) malloc(sizeof(unsigned short)*width);
  
  
  for (j=0;j<height;j++)
  {
  fread(&(image2[j][0]),sizeof(char),width,fp_ref);
  fread(&(image[j][0]),sizeof(char),width,fp_cur);
  fread(&(masks[j][0]),sizeof(unsigned short),width,fp_masks);
  }
  
  for (j=0;j<height;j++)
 	for (i=0;i<width;i++)
		image3[j][i]=0;

  /*  call the C subroutine */
  
  printf("try to predict the image %d x %d\n",width,height);
  
  compute_pred(image,image2,image3,masks,width,height,N);
  
 
 for (j=0;j<height;j++)
  fwrite(&(image3[j][0]),sizeof(char),width,fp_out);

 
 for (j=0;j<height;j++)
    free(image[j]);
  free(image);
  
  for (j=0;j<height;j++)
    free(image2[j]);
  free(image2);
  
  for (j=0;j<height;j++)
    free(masks[j]);
  free(masks);
  
  for (j=0;j<height;j++)
    free(image3[j]);
  free(image3);
  
  for (j=0;j<N;j++)
    free(mvvet[j]);
  free(mvvet);
 
 
 fclose(fp_ref);
 fclose(fp_cur);
 fclose(fp_masks);
 fclose(fp_out);

}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    if(nrhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:full_src:nrhs",
                      "Three inputs required.");
        }


    if(nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:full_src:nlhs",
		"One output required.");
        }
    

    if( !mxIsUint8(prhs[0]) || mxIsComplex(prhs[0])) {
            mexErrMsgIdAndTxt("MyToolbox:full_src:notDouble",
                "Input matrix must be type uint8.");
            }

       /* check that number of rows in first input argument is 1 */
    if(mxGetM(prhs[0]) == 1 ) {
        mexErrMsgIdAndTxt("MyToolbox:full_src:notRowVector",
                      "Input must be a matrix.");
            }

    if( !mxIsUint8(prhs[1]) || mxIsComplex(prhs[1])) {
            mexErrMsgIdAndTxt("MyToolbox:full_src:notDouble",
                "Input matrix must be type uint8.");
            }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1]) == 1 ) {
        mexErrMsgIdAndTxt("MyToolbox:full_src:notRowVector",
                      "Input must be a matrix.");
            }
        
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
            mxGetNumberOfElements(prhs[2]) != 1 ) {
                mexErrMsgIdAndTxt("MyToolbox:full_src:notScalar",
                      "Input factor must be a scalar.");
                }


    /*DEFINE VARIABLES*/

    unsigned char *img_in;       /* MxN input matrix */
    unsigned char *img_ref;       /* MxN input matrix */
    int qp;
    int n_mv,x,y,xb,yb,w4,h4,ind,coeff_cost,i,j;
    int *coeff_array; 

    mwSize ndesc_row;           /* size of matrix */
    mwSize ndesc_col;           /* size of matrix */
    
    int *outCent;      /* output vet */
    int block[4][4];
    /*END OF DEFINITION*/

    /*MEMORY ALLOCATION*/
    /* get the value of the scalar input  */
    qp = (int) ((double) mxGetScalar(prhs[2]));
 
    /* create a pointer to the real data in the input matrix  */
    img_in = (unsigned char *) mxGetData(prhs[0]);
    img_ref = (unsigned char *) mxGetData(prhs[1]);

    /* get dimensions of the input matrix */
    ndesc_col = mxGetN(prhs[0]);
    ndesc_row = mxGetM(prhs[0]);
    
    n_mv=(ndesc_col*ndesc_row)/(16);
    w4=ndesc_col/4;
    h4=ndesc_row/4;
    
    /* create the output image */
    plhs[0] = mxCreateNumericMatrix(n_mv,16,mxINT32_CLASS,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outCent = (int *) mxGetData(plhs[0]);

    coeff_array=(int*) malloc(17*sizeof(int));
    
    /*Compute*/
    
    for (yb=0;yb<h4;yb++)
	for (xb=0;xb<w4;xb++)
        {
		ind=xb+w4*yb;
		for (y=0;y<4;y++)
			for (x=0;x<4;x++)
                {
                j=(yb<<2)+y;
                i=(xb<<2)+x;
				block[y][x]=img_in[j*ndesc_row+i]-img_ref[j*ndesc_row+i];
                }
				
		coeff_cost=0;
		dct_luma(block,&(coeff_array[0]),&coeff_cost,qp,ind);
        memcpy(&(outCent[ind*16]),&(coeff_array[0]),(16*sizeof(int)));
        }


    free(coeff_array);
}
