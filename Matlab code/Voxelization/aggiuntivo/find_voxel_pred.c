#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#define absm(a)  ((a<0)?-1*a:a)

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    if(nrhs != 5) {
        mexErrMsgIdAndTxt("MyToolbox:full_src:nrhs",
                      "Four inputs required.");
        }


    if(nlhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:full_src:nlhs",
		"Two output required.");
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
    
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
            mxGetNumberOfElements(prhs[3]) != 1 ) {
                mexErrMsgIdAndTxt("MyToolbox:full_src:notScalar",
                      "Input factor must be a scalar.");
                }


    /*DEFINE VARIABLES*/

    unsigned char *img_in;       /* MxN input matrix */
    unsigned char *img_ref;       /* MxN input matrix */
    int Nvx,blk_size,win_size,nvx_b,nvx_c;
    int limx1, limx2, limy1, limy2, limz1, limz2, i, j, dx, dy, dz;
    int n_mv,x,y,z,xb,yb,zb,w4,w4b,w4c,tmp_c,npx_blk,cnt_blk,npred;
    unsigned char *pos,*pos_ref, *pos1, *pos2, *best_pos, *pos_out;
    int *pos3;
    int *pos_blk;
    double sum_diff,best_diff;

    mwSize ndesc_row;           /* size of matrix */
    mwSize ndesc_col;           /* size of matrix */
    
    unsigned char *outCent;      /* output vet */
    int *outMV;      /* output vet */
    /*END OF DEFINITION*/

    /*MEMORY ALLOCATION*/
    /* get the value of the scalar input  */
    Nvx = (int) ((double) mxGetScalar(prhs[2]));
    blk_size = (int) ((double) mxGetScalar(prhs[3]));
    win_size = (int) ((double) mxGetScalar(prhs[4]));
 
    /* create a pointer to the real data in the input matrix  */
    img_in = (unsigned char *) mxGetData(prhs[0]);
    img_ref = (unsigned char *) mxGetData(prhs[1]);

    /* get dimensions of the input matrix */
    ndesc_col = mxGetN(prhs[0]);
    ndesc_row = mxGetM(prhs[0]);
    
    n_mv=(ndesc_row)/(blk_size*blk_size*blk_size);
    w4=Nvx/blk_size;
    
    /* create the output image */
    plhs[0] = mxCreateNumericMatrix(ndesc_row,1,mxINT8_CLASS,mxREAL);
    /*create mv array*/
    plhs[1] = mxCreateNumericMatrix(3,n_mv,mxINT32_CLASS,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outCent = (unsigned char *) mxGetData(plhs[0]);
    outMV = (int*) mxGetData(plhs[1]);
    
   
    w4b=Nvx*blk_size-Nvx;
    w4c=(Nvx*Nvx)*(blk_size-1);
    nvx_b=Nvx;
    nvx_c=Nvx*Nvx;
    npx_blk=blk_size*blk_size*blk_size;
    
    pos_blk=(int*) malloc(npx_blk*sizeof(int));
    
    
    /*define coordinates ofr block voxels*/
    i=0;
    tmp_c=0;
    for (z=0;z<blk_size;z++)
        {
        for (y=0;y<blk_size;y++)
            {
            for (x=0;x<blk_size;x++,i++)
                {
                pos_blk[i]=tmp_c;
                tmp_c++;
                }
            tmp_c+=(Nvx-blk_size);
            }
        tmp_c+=(nvx_c-blk_size*Nvx);
        }
    
    /*for (i=0;i<npx_blk;i++)
        printf("%d ",pos_blk[i]);
    printf("\n");*/

    /*Compute*/
    cnt_blk=0;
    pos=&(img_in[0]);
    pos_out=&(outCent[0]);
    for (zb=0;zb<Nvx;zb+=blk_size)
        {
        for (yb=0;yb<Nvx;yb+=blk_size)
            {
            for (xb=0;xb<Nvx;xb+=blk_size)
                {
                
                limx1=(xb-win_size);
                limx1=(limx1<0)?0:limx1;  
                limx2=(xb+win_size+1);
                limx2=(limx2>(Nvx-blk_size+1))?(Nvx-blk_size+1):limx2;
                
                limy1=(yb-win_size);
                limy1=(limy1<0)?0:limy1;  
                limy2=(yb+win_size+1);
                limy2=(limy2>(Nvx-blk_size+1))?(Nvx-blk_size+1):limy2;
                
                limz1=(zb-win_size);
                limz1=(limz1<0)?0:limz1;  
                limz2=(zb+win_size+1);
                limz2=(limz2>(Nvx-blk_size+1))?(Nvx-blk_size+1):limz2;
                
                /*printf("%d %d(%d,%d,%d) [%d,%d]x[%d,%d]x[%d,%d]\n",
                        cnt_blk,pos-&(img_in[0]),xb,yb,zb,limx1,limx2,limy1,limy2,limz1,limz2);*/

                /*find best pred for the current block*/
                best_diff=10000000;
                pos_ref=&(img_ref[0])+limz1*nvx_c+limy1*nvx_b+limx1;
                for (z=limz1;z<limz2;z++)
                    {
                    for (y=limy1;y<limy2;y++)
                        {
                        for (x=limx1;x<limx2;x++)
                            {
                            sum_diff=0;
                            pos3=&(pos_blk[0]);
                            for (i=0;i<npx_blk;i++,pos3++)
                                {
                                pos1=pos+(*pos3);
                                pos2=pos_ref+(*pos3);
                                sum_diff=sum_diff+((*pos1)^(*pos2));
                                }
                            
                            dx=x-xb;
                            dy=y-yb;
                            dz=z-zb;
                            
                            sum_diff=sum_diff+0.1*(absm(dx)+absm(dy)+absm(dz));
                            /*printf("%d %d(%d,%d,%d) %d\n",cnt_blk,pos_ref-&(img_ref[0]),
                                    x-xb,y-yb,z-zb,sum_diff);*/
                            
                            if (sum_diff<best_diff)
                                {
                                best_diff=sum_diff;
                                best_pos=pos_ref;
                                /*printf("%d (%d,%d,%d) %d\n",cnt_blk,x-xb,y-yb,z-zb,best_diff);*/
                                outMV[cnt_blk*3]=dx;
                                outMV[cnt_blk*3+1]=dy;
                                outMV[cnt_blk*3+2]=dz;
                                }
                            
                            pos_ref++;
                            }
                        pos_ref+=(nvx_b+limx1-limx2);
                        }
                    pos_ref+=(nvx_c+((limy1-limy2)*nvx_b));
                    }
                /*end of prediction for the current block*/
                pos3=&(pos_blk[0]);
                pos_out=&(outCent[0])+(pos-&(img_in[0]));
                for (i=0;i<npx_blk;i++,pos3++)
                	{
                    pos1=best_pos+(*pos3);
                    pos2=pos_out+(*pos3);
                    *pos2=*pos1;
                    }
                
                /*printf("\n");*/
                cnt_blk=cnt_blk+1;
                
                pos+=blk_size;
                }
            pos+=w4b;
            }
        pos+=w4c;
        }

    free(pos_blk);
}