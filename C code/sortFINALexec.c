//! File for RAHT-based color coding
/*!Author: Enrico Polo
Date: Dec. 2017
RAHT-based coding routine for PCL coding (color)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>

#include "ac.h" /*!< header files for arith. coder functions */

//to store occupied voxels attributes
struct RGBvalue {
	double R;
	double G;
	double B;
};

//to navigate occupied voxels
struct tosort {
	int idx;
	int colidx;
};


//finds linear indexes of occupied bits (in raster scan order) from input and return their number
int find_occupied(unsigned char* voxel, int dim, struct tosort* occ)	
{								
	int i, j =0;
	#pragma acc kernels
	for (i = 0; i < dim*dim*dim; i++)
	{
		if (*voxel == 1)
		{
			occ[j].idx = i+1;
			occ[j].colidx = j;
			j++;
		}
	voxel++;
	}
	return j;
}

//une step of the RAHT transform + coordinates update (for the next iteration)
int colorTF(struct tosort *occ, int size, struct RGBvalue *colorRGB, unsigned char *weight, int *msx, int *mdx, int ptr, int k, int *actualdim)
{
	int i =0, j = 0,sx,dx,mod,div,tmp,tmp2,index;
	double a, b;
	struct RGBvalue csx, cdx;
	
	while (i+ptr < size)
		{
		if ((occ[i+ ptr].idx % 2 == 1) & (occ[i+ptr].idx == occ[i +ptr+ 1].idx - 1))	//adjacent voxels
			{
			sx = occ[ptr+i].colidx;
			dx = occ[ptr+i+1].colidx;
			a = (double)weight[sx] / (weight[sx] + weight[dx]);
			b = (double)weight[dx] / (weight[sx] + weight[dx]);
			a = sqrt(a);
			b = sqrt(b);
			csx = colorRGB[sx];
		        cdx = colorRGB[dx];
			//DC in sx
			colorRGB[sx].R = b*cdx.R + a*csx.R;
			colorRGB[sx].G = b*cdx.G + a*csx.G;
			colorRGB[sx].B = b*cdx.B + a*csx.B;
			//AC in dx 
			colorRGB[dx].R = a*cdx.R - b*csx.R;
			colorRGB[dx].G = a*cdx.G - b*csx.G;
			colorRGB[dx].B = a*cdx.B - b*csx.B;

			//update coordinates
			occ[ptr+i+1].idx = 0;
			occ[ptr+i].idx = ceil((double)occ[ptr+i].idx / 2);
			tmp = actualdim[(k - 1) % 3] * actualdim[k % 3];
			tmp2 = ((double)(occ[ptr+i].idx-1) / tmp);
			div = floor((double)tmp2);
			mod = (occ[ptr+i].idx-1) % tmp;
			occ[ptr+i].idx = div*actualdim[k % 3] + floor((double)mod / actualdim[(k - 1) % 3]) + mod%actualdim[(k - 1) % 3] * actualdim[k % 3] * actualdim[(k + 1) % 3] +1;

			//increment weights
			weight[sx] = weight[sx]+1;		
			weight[dx] = weight[dx] + 1;
			
			mdx[ptr + j] = dx;
			msx[ptr + j] = sx;
			j++;
			i++;
		}
		else {  //only one occupied voxel in a considered couple
		//update coordinates
			occ[ptr+i].idx = ceil((double)occ[ptr+i].idx / 2);	
			tmp = actualdim[(k - 1) % 3] * actualdim[k % 3];
			tmp2 = ((double)(occ[ptr+i].idx-1) / tmp);
			div = floor((double)tmp2);
			mod = (occ[ptr+i].idx-1) % tmp;
			occ[ptr+i].idx = div*actualdim[k % 3] + floor((double)mod / actualdim[(k - 1) % 3]) + mod%actualdim[(k - 1) % 3] * actualdim[k % 3] * actualdim[(k + 1) % 3] +1;
		}
		i++;
	}
	return j;
}


int tosort_compare(const void* a, const void* b){	
	return ((struct tosort*)a)->idx - ((struct tosort*)b)->idx;
}



void fullcolorTF(struct tosort *occ, int size, struct RGBvalue *colorRGB, unsigned char *weight, int iteration, int *actualdim, int *msx, int *mdx, int *ptr)
{
	int i = 0, j = 0;
	struct tosort * aux;
	for (i = 0; i < 3 * iteration; i++) {
		if (i == 0) {
			//init unit weight
			while (j < size) {		
				weight[j] = 1;	
				j++;
			}
			//half the resolution of the considered dimension (to update indexes)	
			actualdim[i % 3] = actualdim[i % 3] / 2;  
			//transform and update linear indexes (one iteration)	
			ptr[i] = colorTF(occ, size, colorRGB, weight,msx,mdx,0,i+1,actualdim);				
		}
		else {	
		clock_t tStart = clock();
		//sort survived voxels according to new linear indexes ()					
		if (i == 1){		
			clock_t tStart2 = clock();
			qsort(occ, size, sizeof(struct tosort), &tosort_compare);
		}
		else{
			aux = &(occ[ptr[i-2]]);
			qsort(aux, size-ptr[i-2], sizeof(struct tosort), &tosort_compare);			
		}
		//half the resolution of the considered dimension (to update indexes)
		actualdim[i % 3] = actualdim[i % 3] / 2;		
		//transform and update linear indexes (one iteration)		
		ptr[i] = colorTF(occ, size, colorRGB, weight, msx, mdx,ptr[i-1],i+1,actualdim); 
		ptr[i] = ptr[i-1]+ ptr[i];
		
		}
	}
}

void colorITF(struct RGBvalue * rgb, int * sxidx, int * dxidx, unsigned char* we, int gap, unsigned char* decomposed) {
	double a, b;
	int sx, dx, j=0;
	struct RGBvalue csx,cdx;
	//for all merged voxel at a certain level
	for (int i = 0; i < gap; i++,sxidx++,dxidx++) {
			sx = sxidx[0];
			dx = dxidx[0];
			we[dx] = we[dx] - 1;
			we[sx] = we[sx] - 1;

			a = (double)we[sx] / (we[sx] + we[dx]);
			b = (double)we[dx] / (we[sx] + we[dx]);
			a = sqrt(a);
			b = sqrt(b);

			csx = rgb[sx];
			cdx = rgb[dx];
			decomposed[sx] = 1;
			decomposed[dx] = 1;
				rgb[sx].R = -b*cdx.R + a*csx.R;
				rgb[sx].G = -b*cdx.G + a*csx.G;
				rgb[sx].B = -b*cdx.B + a*csx.B;
				
				rgb[dx].R = a*cdx.R + b*csx.R;
				rgb[dx].G = a*cdx.G + b*csx.G;
				rgb[dx].B = a*cdx.B + b*csx.B;
				j++;		
		}
}

void fullcolorITF(unsigned char *weight, struct RGBvalue * rgb, int* mergeddx, int* mergedsx, int *nonzero,int iter2,int iter, unsigned char* fl) {
	
	int sx, dx, aux, aux2, *sss, *ddd,ff;
	nonzero = &(nonzero[3*iter-1]);		//points at the end
	for (int i = 3 * iter -1; i>3 * iter -3 * iter2 -1; i--) {
		aux = (*nonzero);
		if (i != 0) {
			nonzero--;
			aux2 = (*nonzero);
		}
		else{
			aux2 = 0;
		}
		sss = &(mergedsx[aux2]);	//first sx merged voxel at a certain decomposition level
		ddd = &(mergeddx[aux2]);	//first dx merged voxel at a certain decomposition level
		//inverse transform (one step)
		colorITF(rgb, sss, ddd, weight, aux-aux2,fl);		
	}
}

void quantizeYUV(struct RGBvalue *color, struct RGBvalue *quant, int stp, int occ, int *mx, int *mn, int* st, int flag){
	int min[3];
	mn[0] = 0;
	mn[1] = 0;
	mn[2] = 0;
	int max[3];
	mx[0] = 0;
	mx[1] = 0;
	mx[2] = 0;
	int step[3];
	st[0] = stp;
	st[1] = stp;
	st[2] = stp;
	//if YUV
	if(flag){
		st[1] = st[0]*2;
		st[2] = st[0]*2;
	}
	for (int i = 0; i<occ;i++){
		quant[i].R = round((double) color[i].R/st[0])*st[0];
		if(quant[i].R<mn[0]){
			mn[0] = (int) quant[i].R;
		}
		else if(quant[i].R>mx[0]){
			mx[0] = (int) quant[i].R;
		}
		quant[i].G = round((double) color[i].G/st[1])*st[1];
		if(quant[i].G<mn[1]){
			mn[1] = (int) quant[i].G;
		}
		else if(quant[i].G>mx[1]){
			mx[1] = (int) quant[i].G;
		}
		quant[i].B = round((double) color[i].B/st[2])*st[2];
		if(quant[i].B<mn[2]){
			mn[2] = (int) quant[i].B;
		}	
		else if(quant[i].B>mx[2]){
			mx[2] = (int) quant[i].B;
		}	
	}		
}

double compMSE(struct RGBvalue *orig, struct RGBvalue *quant, int occ){
	double MSE = 0;
	for (int i=0;i<occ;i++){
		MSE = MSE + (orig[i].R-quant[i].R)*(orig[i].R-quant[i].R);
		MSE = MSE + (orig[i].G-quant[i].G)*(orig[i].G-quant[i].G);
		MSE = MSE + (orig[i].B-quant[i].B)*(orig[i].B-quant[i].B);
	}
	MSE = (double) MSE/(3*occ);
	return MSE;
}


int main(int argc, char **argv)
{
	struct RGBvalue *color, *quant_color, *quant_color2, *color_orig;
	struct tosort *occupied;
	int ccc, ddd,eee,dim, occsize, *actdim, iter, iter2,*mergedsx,*mergeddx,*nz,*nz2, symbols, i, *max, *offs, step, *st;
	unsigned char *input, *w, *inputcol, *decomp;
	double PSNR,MSE,tftime, enttime;
	FILE* fp, *fc, *fout;


	tftime=0;
	enttime=0;
	for(int x=0;x<30;x++){

	dim = (int)atoi(argv[2]);	//grid resolution
	step = (int)atoi(argv[4]);	//quantization step
	int dimm[] = { dim, dim, dim };  //actual resolution along {x,y,z}
	actdim = &dimm[0];

	input = (unsigned char*)malloc(sizeof(unsigned char)*dim*dim*dim);
	occupied = (struct tosort *)malloc((dim*dim*dim) * sizeof(struct tosort));	
	
	/*read the input voxel volume*/
	fp = fopen(argv[1], "r");
	ccc = fread(input, sizeof(unsigned char), dim*dim*dim, fp);
	fclose(fp);
		
	iter = ((double)log10((double)dim) / log10(2.0));	//number of decomposition levels
	
	clock_t tStart4 = clock();
	occsize = find_occupied(input, dim, occupied);
	printf("\n\nTime taken (only initial scan): %.10fs\n", (double)(clock() - tStart4) / CLOCKS_PER_SEC);
	inputcol = (unsigned char*)malloc(sizeof(unsigned char)*3*occsize);

	/*read the input color*/
	fc = fopen(argv[3], "rb");
	ddd = fread(inputcol, sizeof(unsigned char),3*occsize, fc);
	fclose(fc);

	//init color arrays and weight
	w = (unsigned char*)malloc(sizeof(unsigned char)*occsize);	
	color = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);	
	quant_color = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);	
        quant_color2 = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);	
	color_orig = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
	decomp = (unsigned char*)malloc(sizeof(unsigned char)*3*occsize);

	
	//fill color structs
	for (int i = 0; i < occsize; i++) {
		color[i].R = (double) inputcol[i];
		color[i].G = (double) inputcol[i+occsize];
		color[i].B = (double) inputcol[i + 2 * occsize];
		color_orig[i].R = (double) inputcol[i];
		color_orig[i].G = (double) inputcol[i+occsize];
		color_orig[i].B = (double) inputcol[i + 2 * occsize];
	}

	//RGB->YUV
	if ((int)atoi(argv[5]) == 1){
		for (int i = 0; i < occsize; i++) {
		color[i].R = (double) 0.299*inputcol[i] + 0.587*inputcol[i+occsize] + 0.114*inputcol[i + 2 * occsize];
		color[i].G = (double) 0.492*(inputcol[i + 2 * occsize]-color[i].R);
		color[i].B = (double) 0.877*(inputcol[i]-color[i].R);
		}
	}

	//init arrays in order to decode
	mergeddx = (int*)malloc(sizeof(int)*occsize);
	mergedsx = (int*)malloc(sizeof(int)*occsize);
	nz = (int*)malloc(sizeof(int) * 3 * iter);

	ac_encoder ace1, ace2, ace3;
	ac_model acm1, acm2, acm3;
	
	//transform
	clock_t tStart = clock();
	fullcolorTF(occupied, occsize, color, w, iter, actdim, mergedsx, mergeddx, nz);
	
	//quantization
	max = (int*)malloc(sizeof(int)*3);
	offs = (int*)malloc(sizeof(int)*3);
	st = (int*)malloc(sizeof(int)*3);
	quantizeYUV(color, quant_color, step, occsize, max, offs, st, (int)atoi(argv[5]));
	tftime = tftime+ (double)(clock() - tStart) / CLOCKS_PER_SEC;
	printf("Time taken (only tf): %.10fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	int symbols1 = (max[0] - offs[0])/st[0] +1;
	int symbols2 = (max[1] - offs[1])/st[1] +1;
	int symbols3 = (max[2] - offs[2])/st[2] +1;			
	
	//arith coding
	clock_t t2Start = clock();
	ac_encoder_init (&ace1, "foo1");
	ac_encoder_init (&ace2, "foo2");
	ac_encoder_init (&ace3, "foo3");
    	ac_model_init (&acm1, symbols1, NULL, 1);
	ac_model_init (&acm2, symbols2, NULL, 1);
	ac_model_init (&acm3, symbols3, NULL, 1);
	for (i=0; i<occsize; i++){
        	ac_encode_symbol (&ace1,&acm1, (int) (quant_color[i].R - offs[0])/st[0]);		//offs to adjust negative values 
		ac_encode_symbol (&ace2,&acm2, (int) (quant_color[i].G - offs[1])/st[1]);		//step to adjust real different values
		ac_encode_symbol (&ace3,&acm3, (int) (quant_color[i].B - offs[2])/st[2]);
	}
	enttime = enttime+ (double)(clock() - t2Start) / CLOCKS_PER_SEC;

	//coded bits
	double outlength1=(double) ac_encoder_bits (&ace1);
	double outlength2=(double) ac_encoder_bits (&ace2);
	double outlength3=(double) ac_encoder_bits (&ace3);
	ac_encoder_done (&ace1);
	ac_encoder_done (&ace2);
	ac_encoder_done (&ace3);
	ac_model_done (&acm1);
	ac_model_done (&acm2);
	ac_model_done (&acm3);

	printf("Time taken (entropy coding): %.10fs\n", (double)(clock() - t2Start) / CLOCKS_PER_SEC);
	printf("Time taken (full coding): %.10fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("Volume dimension: %d*%d*%d\n",dim,dim,dim);
	printf("Bits written (Y):%f\n",outlength1);
	printf("Arith. coder different symbols: %d\n",symbols1);
	printf("Bits written (U):%f\n",outlength2);
	printf("Arith. coder different symbols: %d\n",symbols2);
	printf("Bits written (V):%f\n",outlength3);
	printf("Arith. coder different symbols: %d\n",symbols3);
	printf("Bits written (total):%f\n",outlength3+outlength2+outlength1);
	printf("Occupied voxels: ");
	printf("%d\n", occsize);


	//decoding routine
	clock_t t3Start = clock();
	ac_decoder acd1,acd2,acd3;
	ac_decoder_init (&acd1, "foo1");
	ac_decoder_init (&acd2, "foo2");
	ac_decoder_init (&acd3, "foo3");
    	ac_model_init (&acm1, symbols1, NULL, 1);
	ac_model_init (&acm2, symbols2, NULL, 1);
	ac_model_init (&acm3, symbols3, NULL, 1);


	//entropy decoding
	for (i=0; i<occsize; i++){
	quant_color2[i].R = (double) ac_decode_symbol (&acd1, &acm1)*st[0] + offs[0];
		/*if (quant_color2[i].R!=quant_color[i].R)
				printf("Error!");*/
	quant_color2[i].G = (double) ac_decode_symbol (&acd2, &acm2)*st[1] + offs[1];
	/*if (quant_color2[i].G!=quant_color[i].G)
				printf("Error!");	*/
	quant_color2[i].B = (double) ac_decode_symbol (&acd3, &acm3)*st[2] + offs[2];
	/*if (quant_color2[i].B!=quant_color[i].B)
				printf("Error!");*/
	}

	ac_decoder_done (&acd1);
	ac_decoder_done (&acd2);
	ac_decoder_done (&acd3);
	ac_model_done (&acm1);
	ac_model_done (&acm2);
	ac_model_done (&acm3);
	printf("Time taken (entropy decoding): %.10fs\n", (double)(clock() - t3Start) / CLOCKS_PER_SEC);

	//inverse transform
	iter2 = (int) atoi(argv[6]);  //number of decoding levels (enables partial decomposition)
	fullcolorITF(w, quant_color2, mergeddx, mergedsx, nz, iter2, iter, decomp);
	printf("Time taken (entropydecoding+itf): %.10fs\n", (double)(clock() - t3Start) / CLOCKS_PER_SEC);



	//YUV --> RGB
	if ((int)atoi(argv[5]) == 1){
		for (int i = 0; i < occsize; i++) {
		double Y = quant_color2[i].R; 
		double U = quant_color2[i].G; 
		double V = quant_color2[i].B; 
		quant_color2[i].R = Y+1.14*V;
		quant_color2[i].G = Y-0.395*U-0.581*V;
		quant_color2[i].B = Y+2.032*U;
		}
	}
	
	//  to compare original and partially decoded colors (doesn't work, we have original color attributes only at full resolution) 
	int outsize = 0;
	for (int i = 0; i < occsize; i++) {
		if(decomp[i] == 1){
			/*quant_color2[outsize].R = (double) quant_color2[i].R;///((iter-iter2)*sqrt(2));
			quant_color2[outsize].G = (double) quant_color2[i].G;///((iter-iter2)*sqrt(2));
			quant_color2[outsize].B = (double) quant_color2[i].B;///((iter-iter2)*sqrt(2));
			color_orig[outsize].B = color_orig[i].B;
			color_orig[outsize].R = color_orig[i].R;
			color_orig[outsize].G = color_orig[i].G;*/
			outsize++;		
		}
	}

	
	//compute PSNR (works only if iter = iter2, otherwise computed PSNR is meaningless)
	MSE = compMSE(color_orig,quant_color2,outsize);
	PSNR = 20*log10((double) 255/sqrt(MSE));
	printf("PSNR: %f dB\n", PSNR);
	printf("Decoded voxel at level %d: %d\n",iter2, outsize);


	//write output
	fout = fopen("tf_color.col", "wb");
	eee = fwrite(quant_color2, sizeof(double), 3 * outsize, fout);
	fclose(fout);

	fout = fopen("orig_color.col", "wb");
	eee = fwrite(color_orig, sizeof(double), 3 * outsize, fout);
	fclose(fout);
	
	free(occupied);
	free(w);
	free(color);
	free(mergeddx);
	free(mergedsx);
	free(nz);
	free(quant_color);
	free(quant_color2);
	free(color_orig);
	free(input);
	//free(occupied);
	free(inputcol);
}
printf("\nmean transfrom time:%f\n",tftime/30 );
printf("mean entropy coder time:%f\n",enttime/30 );
	return 0;
}
