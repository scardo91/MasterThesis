#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

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
    int Nvx,blk_size,win_size,nvx_b,nvx_c,sum_diff,best_diff;
    int limx1, limx2, limy1, limy2, limz1, limz2, i, j;
    int n_mv,x,y,z,xb,yb,zb,w4,w4b,w4c,tmp_c,npx_blk,cnt_blk,npred,tmp_c2;
    unsigned char *pos,*pos_ref, *pos_ref2, *pos1, *pos2;
    int *pos3;;
    int *pos_blk;
    int *pos_pred;

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
    npred=(2*win_size+1)*(2*win_size+1)*(2*win_size+1);
    
    pos_blk=(int*) malloc(npx_blk*sizeof(int));
    pos_pred=(int*) malloc(npred*sizeof(int));
    
    printf("%d \n",n_mv);
    
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
    
    tmp_c=-win_size*(Nvx*Nvx)-win_size*Nvx-win_size;
    i=0;
    for (z=0;z<(2*win_size+1);z++)
    	{
    	for (y=0;y<(2*win_size+1);y++)
            {
            for (x=0;x<(2*win_size+1);x++,i++)
                {
                pos_pred[i]=tmp_c;
                tmp_c++;
                }
            tmp_c+=(Nvx-(2*win_size+1));
            }
        tmp_c+=(nvx_c-(2*win_size+1)*Nvx);
        }
    
    tmp_c2=Nvx*Nvx*blk_size;
    
    for (i=0;i<npred;i++)
        printf("%d ",pos_pred[i]);
    printf("\n");
    
    /*Compute*/
    cnt_blk=0;
    pos=&(img_in[0]);
    for (zb=0;zb<Nvx;zb+=blk_size)
        {
        for (yb=0;yb<Nvx;yb+=blk_size)
            {
            for (xb=0;xb<Nvx;xb+=blk_size)
                {
                best_diff=10000000;
                pos_ref=&(img_ref[0])+(pos-&(img_in[0]));
                for (j=0;j<npred;j++)
                    {
                    pos_ref2=pos_ref+pos_pred[j];
                    tmp_c=pos_ref2-&(img_ref[0]);
                    if ((tmp_c>0)&((tmp_c+tmp_c2)<ndesc_row))
                        {
                        sum_diff=0;
                        pos3=&(pos_blk[0]);
                        for (i=0;i<npx_blk;i++,pos3++)
                            {
                            pos1=pos+(*pos3);
                            pos2=pos_ref2+(*pos3);
                            sum_diff=sum_diff+((*pos1)^(*pos2));
                            }
                        tmp_c=pos_pred[j];
                        x=tmp_c%Nvx;
                        tmp_c=tmp_c/Nvx;
                        y=tmp_c%Nvx;
                        tmp_c=tmp_c/Nvx;
                        z=tmp_c%Nvx;
                        printf("%d %d(%d,%d,%d) %d\n",cnt_blk,pos_ref2-&(img_ref[0]),
                                    x,y,z,sum_diff);
                            
                        if (sum_diff<best_diff)
                            {
                            best_diff=sum_diff;
                            /*printf("%d (%d,%d,%d) %d\n",cnt_blk,x-xb,y-yb,z-zb,best_diff);*/
                            outMV[cnt_blk*3]=x;
                            outMV[cnt_blk*3+1]=y;
                            outMV[cnt_blk*3+2]=z;
                            }
                        }
                    }

                printf("\n");
                cnt_blk=cnt_blk+1;
                
                pos+=blk_size;
                }
            pos+=w4b;
            }
        pos+=w4c;
        }

    free(pos_blk);
    free(pos_pred);
}