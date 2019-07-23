//! File for RAHT-based color coding
/*!Author: Enrico Polo
Date: Dec. 2017
RAHT-based coding routine for dynamic PCL coding (color)
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
	int i, j = 0;
	for (i = 0; i < dim*dim*dim; i++)
	{
		if (*voxel == 1)
		{
			occ[j].idx = i + 1;
			occ[j].colidx = j;
			j++;
		}
		voxel++;
	}
	return j;
}

int colorTF(struct tosort *occ, int size, struct RGBvalue *colorRGB,unsigned char * weight, int *msx, int *mdx, int ptr, int k, int *actualdim)
{
	int i = 0, j = 0, sx, dx, mod, div, tmp, tmp2, index;
	double a, b;
	struct RGBvalue csx, cdx;
	while (i + ptr < size)
	{
		if ((occ[i + ptr].idx % 2 == 1) & (occ[i + ptr].idx == occ[i + ptr + 1].idx - 1))  //adjacent voxels
		{
			sx = occ[ptr + i].colidx;
			dx = occ[ptr + i + 1].colidx;
			a = (double)weight[sx] / (weight[sx] + weight[dx]);
			b = (double)weight[dx] / (weight[sx] + weight[dx]);
			a = sqrt(a);
			b = sqrt(b);
			csx = colorRGB[sx];
			cdx = colorRGB[dx];
			//transform
			//DC in sx
			colorRGB[sx].R = b*cdx.R + a*csx.R;
			colorRGB[sx].G = b*cdx.G + a*csx.G;
			colorRGB[sx].B = b*cdx.B + a*csx.B;
			//AC in dx 
			colorRGB[dx].R = a*cdx.R - b*csx.R;
			colorRGB[dx].G = a*cdx.G - b*csx.G;
			colorRGB[dx].B = a*cdx.B - b*csx.B;
			/*uncomment to test approximated transform
			colorRGB[sx].R = (double)(cdx.R + csx.R) / 2;
			colorRGB[sx].G = (double)(cdx.G + csx.G) / 2;
			colorRGB[sx].B = (double)(cdx.B + csx.B) / 2;
			//AC in dx 
			colorRGB[dx].R = (double)(cdx.R - csx.R) / 2;
			colorRGB[dx].G = (double)(cdx.G - csx.G) / 2;
			colorRGB[dx].B = (double)(cdx.B - csx.B) / 2;
			*/
			//update coordinates
			occ[ptr + i + 1].idx = 0;
			occ[ptr + i].idx = ceil((double)occ[ptr + i].idx / 2);
			tmp = actualdim[(k - 1) % 3] * actualdim[k % 3];
			tmp2 = ((double)(occ[ptr + i].idx - 1) / tmp);
			div = floor((double)tmp2);
			mod = (occ[ptr + i].idx - 1) % tmp;
			occ[ptr + i].idx = div*actualdim[k % 3] + floor((double)mod / actualdim[(k - 1) % 3]) + mod%actualdim[(k - 1) % 3] * actualdim[k % 3] * actualdim[(k + 1) % 3] + 1;
			//increment weights
			weight[sx] = weight[sx]+1;		
			weight[dx] = weight[dx] + 1;

			mdx[ptr + j] = dx;
			msx[ptr + j] = sx;

			j++;
			i++;
		}
		else {	//only one occupied voxel in a considered couple
		//update coordinates
			occ[ptr + i].idx = ceil((double)occ[ptr + i].idx / 2);	
			tmp = actualdim[(k - 1) % 3] * actualdim[k % 3];
			tmp2 = ((double)(occ[ptr + i].idx - 1) / tmp);
			div = floor((double)tmp2);
			mod = (occ[ptr + i].idx - 1) % tmp;
			occ[ptr + i].idx = div*actualdim[k % 3] + floor((double)mod / actualdim[(k - 1) % 3]) + mod%actualdim[(k - 1) % 3] * actualdim[k % 3] * actualdim[(k + 1) % 3] + 1;
		}
		i++;
	}
	return j;
}

int tosort_compare(const void* a, const void* b) {
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


void colorITF(struct RGBvalue * rgb, unsigned char *we, int * sxidx, int * dxidx, int gap, unsigned char* decomposed) {
	double a, b;
	int sx, dx, j = 0;
	struct RGBvalue csx, cdx;
	//printf("itf\n");
	for (int i = 0; i < gap; i++, sxidx++, dxidx++) {
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
		
		/*uncomment to test approximated ITF
		rgb[sx].R = (double)(-cdx.R + csx.R);
		rgb[sx].G = (double)(-cdx.G + csx.G);
		rgb[sx].B = (double)(-cdx.B + csx.B);

		rgb[dx].R = (double)(cdx.R + csx.R);
		rgb[dx].G = (double)(cdx.G + csx.G);
		rgb[dx].B = (double)(cdx.B + csx.B);*/
		
		//ITF
		rgb[sx].R = -b*cdx.R + a*csx.R;
		rgb[sx].G = -b*cdx.G + a*csx.G;
		rgb[sx].B = -b*cdx.B + a*csx.B;

		rgb[dx].R = a*cdx.R + b*csx.R;
		rgb[dx].G = a*cdx.G + b*csx.G;
		rgb[dx].B = a*cdx.B + b*csx.B;
		j++;
	}
}

void fullcolorITF(struct RGBvalue * rgb,unsigned char *weight, int* mergeddx, int* mergedsx, int *nonzero, int iter2, int iter, unsigned char* fl) {

	int sx, dx, aux, aux2, *sss, *ddd, ff;
	nonzero = &(nonzero[3 * iter - 1]);		//points at the end
	for (int i = 3 * iter - 1; i>3 * iter - 3 * iter2 - 1; i--) {
		aux = (*nonzero);
		if (i != 0) {
			nonzero--;
			aux2 = (*nonzero);
		}
		else {
			aux2 = 0;
		}
		sss = &(mergedsx[aux2]); 	//first sx merged voxel at a certain decomposition level
		ddd = &(mergeddx[aux2]);	//first dx merged voxel at a certain decomposition level
		//inverse transform (one step)
		colorITF(rgb,weight, sss, ddd, aux - aux2, fl);
	}
}

void quantizeYUV(struct RGBvalue *color, struct RGBvalue *quant, int stp, int occ, int *mx, int *mn, int* st, int fl) {
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
	//RGB	
	st[1] = st[0];
	st[2] = st[0];
	//if YUV
	if (fl) {
		st[1] = st[0] * 2;
		st[2] = st[0] * 2;
	}
	for (int i = 0; i<occ; i++) {
		quant[i].R = round((double)color[i].R / st[0])*st[0];
		if (quant[i].R<mn[0]) {
			mn[0] = (int)quant[i].R;
		}
		else if (quant[i].R>mx[0]) {
			mx[0] = (int)quant[i].R;
		}
		quant[i].G = round((double)color[i].G / st[1])*st[1];
		if (quant[i].G<mn[1]) {
			mn[1] = (int)quant[i].G;
		}
		else if (quant[i].G>mx[1]) {
			mx[1] = (int)quant[i].G;
		}
		quant[i].B = round((double)color[i].B / st[2])*st[2];
		if (quant[i].B<mn[2]) {
			mn[2] = (int)quant[i].B;
		}
		else if (quant[i].B>mx[2]) {
			mx[2] = (int)quant[i].B;
		}
	}
}

double compMSE(struct RGBvalue *orig, struct RGBvalue *quant, int occ) {
	double MSE = 0;
	for (int i = 0; i<occ; i++) {
		MSE = MSE + (orig[i].R - quant[i].R)*(orig[i].R - quant[i].R);
		MSE = MSE + (orig[i].G - quant[i].G)*(orig[i].G - quant[i].G);
		MSE = MSE + (orig[i].B - quant[i].B)*(orig[i].B - quant[i].B);
	}
	MSE = (double)MSE / (3 * occ);
	return MSE;
}


struct RGBvalue* RAHTcoder(struct tosort * geom, int res, struct RGBvalue * icol, int qst, char * yuv, char * decstp, int occsize, double *mr, double *mt, double *me)
{

	struct RGBvalue *color, *quant_color, *quant_color2, *color_orig, *outcol;
	int ccc, ddd, eee, dim,  *actdim, iter, iter2, *mergedsx, *mergeddx, *nz, *nz2, symbols, i, *max, *offs, step, *st;
	int symbols1, symbols2, symbols3;
	unsigned char *input, *w,  *decomp;
	double PSNR, MSE;
	FILE* fp, *fc, *fout;


	int dimm[] = { res, res, res };	//actual resolution along {x,y,z}
	actdim = &dimm[0];
	iter = ceil(((double)log10((double)res) / log10(2.0))); 	//number of decomposition levels


	//init color and weight
	w = (unsigned char*)malloc(sizeof(unsigned char)*occsize);
	color = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
	quant_color = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
	quant_color2 = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
	decomp = (unsigned char*)malloc(sizeof(unsigned char) * 3 * occsize);


	//init array in order to decode (actually not coded)
	mergeddx = (int*)malloc(sizeof(int)*occsize);
	mergedsx = (int*)malloc(sizeof(int)*occsize);
	nz = (int*)malloc(sizeof(int) * 3 * iter);

	ac_encoder ace1, ace2, ace3;
	ac_model acm1, acm2, acm3;

	//transform
	clock_t tStart = clock();
	fullcolorTF(geom, occsize, icol,w, iter, actdim, mergedsx, mergeddx, nz);
	
	//quantization
	max = (int*)malloc(sizeof(int) * 3);
	offs = (int*)malloc(sizeof(int) * 3);
	st = (int*)malloc(sizeof(int) * 3);
	quantizeYUV(icol, quant_color, qst, occsize, max, offs, st, (int) atoi(yuv));

	double sc = (double)(clock() - tStart) / CLOCKS_PER_SEC;
	*mt = *mt + sc; //transform total time
	printf("Time taken (only tf): %.10fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	symbols1 = (max[0] - offs[0]) / st[0] + 1;
	symbols2 = (max[1] - offs[1]) / st[1] + 1;
	symbols3 = (max[2] - offs[2]) / st[2] + 1;

	//arith coding
	clock_t t2Start = clock();
	ac_encoder_init(&ace1, "foo1");
	ac_encoder_init(&ace2, "foo2");
	ac_encoder_init(&ace3, "foo3");
	
	ac_model_init(&acm1, symbols1, NULL, 1);
	ac_model_init(&acm2, symbols2, NULL, 1);
	ac_model_init(&acm3, symbols3, NULL, 1);

	for (i = 0; i < occsize; i++) {
		ac_encode_symbol(&ace1, &acm1, (int)(quant_color[i].R - offs[0]) / st[0]);		//offs to adjust negative values 
		ac_encode_symbol(&ace2, &acm2, (int)(quant_color[i].G - offs[1]) / st[1]);		//step to adjust real different values
		ac_encode_symbol(&ace3, &acm3, (int)(quant_color[i].B - offs[2]) / st[2]);
	}
	*me = *me + (double)(clock() - t2Start) / CLOCKS_PER_SEC; //entropy coding total time
	//coded bits
	double outlength1 = (double)ac_encoder_bits(&ace1);
	double outlength2 = (double)ac_encoder_bits(&ace2);
	double outlength3 = (double)ac_encoder_bits(&ace3);
	ac_encoder_done(&ace1);
	ac_encoder_done(&ace2);
	ac_encoder_done(&ace3);
	ac_model_done(&acm1);
	ac_model_done(&acm2);
	ac_model_done(&acm3);

	printf("Time taken (entropy coding): %.10fs\n", (double)(clock() - t2Start) / CLOCKS_PER_SEC);
	printf("Time taken (full coding): %.10fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("Volume dimension: %d*%d*%d\n", res, res, res);
	printf("Bits written (Y):%f\n", outlength1);
	printf("Arith. coder different symbols: %d\n", symbols1);
	printf("Bits written (U):%f\n", outlength2);
	printf("Arith. coder different symbols: %d\n", symbols2);
	printf("Bits written (V):%f\n", outlength3);
	printf("Arith. coder different symbols: %d\n", symbols3);
	printf("Bits written (total):%f\n", outlength3 + outlength2 + outlength1);
	printf("Occupied voxels: ");
	printf("%d\n", occsize);
	*mr =  *mr + outlength3 + outlength2 + outlength1;	//overall bit stream generated

	////decoding routine
	clock_t t3Start = clock();
	ac_decoder acd1, acd2, acd3;
	ac_decoder_init(&acd1, "foo1");
	ac_decoder_init(&acd2, "foo2");
	ac_decoder_init(&acd3, "foo3");
	
	ac_model_init(&acm1, symbols1, NULL, 1);
	ac_model_init(&acm2, symbols2, NULL, 1);
	ac_model_init(&acm3, symbols3, NULL, 1);
	//entropy decoding
	for (i = 0; i < occsize; i++) {
		quant_color2[i].R = (double)ac_decode_symbol(&acd1, &acm1)*st[0] + offs[0];
		/*if (quant_color2[i].R!=quant_color[i].R)
		printf("Error!");*/
		quant_color2[i].G = (double)ac_decode_symbol(&acd2, &acm2)*st[1] + offs[1];
		/*if (quant_color2[i].G!=quant_color[i].G)
		printf("Error!");*/
		quant_color2[i].B = (double)ac_decode_symbol(&acd3, &acm3)*st[2] + offs[2];
		/*if (quant_color2[i].B!=quant_color[i].B)
		printf("Error!");*/
	}

	ac_decoder_done(&acd1);
	ac_decoder_done(&acd2);
	ac_decoder_done(&acd3);
	ac_model_done(&acm1);
	ac_model_done(&acm2);
	ac_model_done(&acm3);
	printf("Time taken (entropy decoding): %.10fs\n", (double)(clock() - t3Start) / CLOCKS_PER_SEC);
	iter2 = (int) atoi(decstp); //number of decoding levels (enables partial decomposition)
	printf("%d,%d,\n", iter, iter2);
	fullcolorITF(quant_color2,w, mergeddx, mergedsx, nz, iter2, iter, decomp);	//inverse transform
	printf("Time taken (entropydecoding+itf): %.10fs\n", (double)(clock() - t3Start) / CLOCKS_PER_SEC);
	
	
	return quant_color2;		//decoded color
}


void predict(struct RGBvalue * act,struct tosort * actgeom,struct RGBvalue * prev,unsigned char * prevgeom, double min, struct RGBvalue * out, int occ, int dim){
	int cnt = 0, *cpoint, help;
	struct RGBvalue *mean;
	mean = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*1);
	mean[0].B = 0;
	mean[0].R = 0;
	mean[0].G = 0;
		for (int i = 0; i < occ; i++) {
			if (prevgeom[actgeom[i].idx - 1] == 0) {	//if voxel in prev. frame is empty
				for (int j2 = -min; j2 <= min; j2++) {
					for (int j1 = -min; j1 <= min; j1++) {
						for (int j = -min; j <= min; j++) {
								help = actgeom[i].idx - 1 - j - j1*dim - j2*dim*dim;
								if(help>-1 && help<dim*dim*dim){
									mean[0].R = mean[0].R + prev[help].R;
									mean[0].G = mean[0].G + prev[help].G;
									mean[0].B = mean[0].B + prev[help].B;
									cnt = cnt + prevgeom[help];
								}
						}
					}
				}
				out[i].R = act[i].R - mean[0].R / cnt;
				out[i].G = act[i].G- mean[0].G / cnt;
				out[i].B = act[i].B -mean[0].B / cnt;
				//if no occ.voxel found in the neighborhood, skip pred.
				if (cnt == 0) {		
					out[i].R = act[i].R;
					out[i].G = act[i].G;
					out[i].B = act[i].B;
				}
				//reset quantities
				cnt = 0;
				mean[0].R = 0;
				mean[0].G = 0;
				mean[0].B = 0;
			}
			else {	//if voxel in prev. frame is occupied
				out[i].R = act[i].R - prev[actgeom[i].idx - 1].R;
				out[i].G = act[i].G- prev[actgeom[i].idx - 1].G;
				out[i].B = act[i].B- prev[actgeom[i].idx - 1].B;
			}
		}
	}

void invpred(struct RGBvalue * decoded,struct tosort * actgeom,struct RGBvalue * prev, unsigned char * prevgeom, double min, struct RGBvalue * out, int occ, int dim){
	int cnt = 0, *cpoint,help;
	struct RGBvalue *mean;
	mean = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*1);
	mean[0].B = 0;
	mean[0].R = 0;
	mean[0].G = 0;
	
		for (int i = 0; i < occ; i++) {
			if (prevgeom[actgeom[i].idx - 1] == 0) {	//if voxel in prev. frame is empty
				for (int j2 = -min; j2 <= min; j2++) {
					for (int j1 = -min; j1 <= min; j1++) {
						for (int j = -min; j <= min; j++) {
								help = actgeom[i].idx - 1 - j - j1*dim - j2*dim*dim;
								if(help>-1 && help<dim*dim*dim){
									mean[0].R = mean[0].R + prev[help].R;
									mean[0].G = mean[0].G + prev[help].G;
									mean[0].B = mean[0].B + prev[help].B;
									cnt = cnt + prevgeom[help];
								}
						}
					}
				}
				out[i].R = decoded[i].R + mean[0].R / cnt;
				out[i].G = decoded[i].G+ mean[0].G / cnt;
				out[i].B = decoded[i].B+ mean[0].B / cnt;
				//if no occ.voxel found in the neighborhood, skip pred.
				if (cnt == 0) {
					out[i].R = decoded[i].R;
					out[i].G = decoded[i].G;
					out[i].B = decoded[i].B;
				}
				//reset quantities
				cnt = 0;
				mean[0].R = 0;
				mean[0].G = 0;
				mean[0].B = 0;
			}
			else {	//if voxel in prev. frame is occupied
				out[i].R = decoded[i].R + prev[actgeom[i].idx - 1].R;
				out[i].G = decoded[i].G+ prev[actgeom[i].idx - 1].G;
				out[i].B = decoded[i].B+prev[actgeom[i].idx - 1].B;
			}
		}
		
	}



int main(int argc, char **argv)
{
	struct RGBvalue* previous, *mean;
	struct RGBvalue *color,*color2, *color_orig, *pdiff, *estim, *diff, *toPSNR, *color3D, *color3D2;
	struct tosort *occupied, *occupied2, *tmp;
	int ccc, ddd, eee, dim, occsize, occsize_old, *actdim, iter, iter2, *nz, i, step, *offs, *max, *st, fl, tot_bit,reference, predstatfl;
	int symbols1, symbols2,symbols3, nframes;
	unsigned char *input, *input2, *inputcol, *inputcol2;
	char next[80], nextcol[80];
	char *aux, *aux2;
	double PSNR, MSE, meanr, *meanrate, PSNRtot, predtime, *tf_time, ttime, *ent_time, etime;
	FILE* fp, *fc, *fout;

	ttime =  0;
	etime = 0;
	predtime = 0;
	PSNRtot = 0;
	meanr = 0;
	meanrate = &meanr;
	tf_time = &ttime;
	ent_time = &etime;
	tot_bit = 0;
	dim = (int)atoi(argv[2]);		//grid resolution along one axis
	step = (int)atoi(argv[4]);	//quantization step	
	fl = (int)atoi(argv[5]);
	nframes = (int)atoi(argv[7]);	//total number of frames
	reference = (int)atoi(argv[9]);		//GOP value
	predstatfl = 1;

	occupied = (struct tosort *)malloc((dim*dim*dim) * sizeof(struct tosort));
	tmp = (struct tosort *)malloc((dim*dim*dim) * sizeof(struct tosort));
	mean = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*1);
	for (int w = 0; w < nframes; w++) {
		predstatfl = 1;
		if (w%reference != 0 && w != 0) {
			predstatfl = 0;	//do prediction
		}
		else if (w != 0 && predstatfl ==1) { //no prediction for this frame
			free(input);		//otherwise no-prediction coding crashes
			free(inputcol),
			free(color);
			free(color_orig);
			free(color3D);
			free(previous);
			free(toPSNR);
		}
		//load actual frame
		if (w < 10) {
			sprintf(next,"frame000%d.bin",w);
			sprintf(nextcol,"frame000%d.col",w);
		}
		else if (w<100) {
			sprintf(next,"frame00%d.bin",w);
			sprintf(nextcol,"frame00%d.col",w);			
		}
		else {
			sprintf(next,"frame0%d.bin",w);
			sprintf(nextcol,"frame0%d.col",w);
		}
		printf("\n\n\n");
		printf("Frame number: %d\n", w);
		
		if (predstatfl) {  //static encoding routine
			/*read the input voxel volume*/
			input = (unsigned char*)malloc(sizeof(unsigned char)*dim*dim*dim);
			fp = fopen(next, "r");
			ccc = fread(input, sizeof(unsigned char), dim*dim*dim, fp);
			fclose(fp);
		
			
			occsize = find_occupied(input, dim, occupied);
			occsize = find_occupied(input, dim, tmp);
			inputcol = (unsigned char*)malloc(sizeof(unsigned char) * 3 * occsize);
			/*read the input color*/
			fp = fopen(nextcol, "r");
			ccc = fread(inputcol, sizeof(unsigned char), 3 * occsize, fp);
			fclose(fp);
			
			tot_bit = tot_bit + occsize;
			color = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
			color_orig = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
			color3D = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*dim*dim*dim);
			
			//fill color structs				
			i=0;
			for (int j = 0; j < dim*dim*dim; j++) {
				if(input[j]==1){
				color[i].R = (double)inputcol[i];
				color[i].G = (double)inputcol[i + occsize];
				color[i].B = (double)inputcol[i + 2 * occsize];
				color_orig[i].R = (double)color[i].R;
				color_orig[i].G = (double)color[i].G;
				color_orig[i].B = (double)color[i].B;
				i++;
				}
			}

			//RGB --> YUV
			i=0;
			if ((int)atoi(argv[5]) == 1) {
				for (int j = 0; j < dim*dim*dim; j++) {
					if(input[j]==1){
					color[i].R = (double) 0.299*inputcol[i] + 0.587*inputcol[i + occsize] + 0.114*inputcol[i + 2 * occsize];
					color[i].G = (double) 0.492*(inputcol[i + 2 * occsize] - color[i].R);
					color[i].B = (double) 0.877*(inputcol[i] - color[i].R);
					i++;
					}
				}
			}
			
			//coding and decoding routine
			previous = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
			previous = RAHTcoder(tmp, dim, color, step, argv[5], argv[6], occsize, meanrate, tf_time, ent_time);

			//backup decoded colors
			i=0;
			for (int j = 0; j < dim*dim*dim; j++) {
				if(input[j]==1){
				color3D[j].R = (double)previous[i].R;
				color3D[j].G = (double)previous[i].G;
				color3D[j].B = (double)previous[i].B;
				i++;
				}
			}
			i=0;

		toPSNR = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
		//YUV-->RGB
		if ((int) atoi(argv[5]) == 1) {
			for (int j = 0; j < occsize; j++) {
			double Y = previous[j].R;
			double U = previous[j].G;
			double V = previous[j].B;
			toPSNR[j].R = Y + 1.14*V;
			toPSNR[j].G = Y - 0.395*U - 0.581*V;
			toPSNR[j].B = Y + 2.032*U;
			
			}
		}
		else{ //already RGB
			for (int j = 0; j < occsize; j++) {
			toPSNR[j].R = previous[j].R;
			toPSNR[j].G = previous[j].G;
			toPSNR[j].B = previous[j].B;
			}
		}
			//compute PSNR
			MSE = compMSE(color_orig, toPSNR, occsize);
			PSNR = 20 * log10((double)255 / sqrt(MSE));
			printf("PSNR: %f dB\n", PSNR);
			PSNRtot = PSNRtot + PSNR;	//sum all PSNRs
			
			/*append decoded frame to the ones already done*/
			fout = fopen("tf_color.col", "a");
			eee = fwrite(toPSNR, sizeof(double), 3 * occsize, fout);
			fclose(fout);
			/*append original frame to the ones already done*/
			fout = fopen("orig_color.col", "a");
			eee = fwrite(color_orig, sizeof(double), 3 * occsize, fout);
			fclose(fout);
		}
		else {	//predictive coding routine
			/*read the input voxel volume*/
			input2 = (unsigned char*)malloc(sizeof(unsigned char)*dim*dim*dim);
			fp = fopen(next, "r");
			ccc = fread(input2, sizeof(unsigned char), dim*dim*dim, fp);
			fclose(fp);

			occupied2 = (struct tosort *)malloc((dim*dim*dim) * sizeof(struct tosort));
			occsize = find_occupied(input2, dim, occupied2);
			occsize = find_occupied(input2, dim, tmp);
			tot_bit = tot_bit + occsize;

			/*read the input color*/
			inputcol2 = (unsigned char*)malloc(sizeof(unsigned char) * 3 * occsize);
			fp = fopen(nextcol, "r");
			ccc = fread(inputcol2, sizeof(unsigned char), 3*occsize, fp);
			fclose(fp);

			//fill color structs
			color2 = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
			i=0;
			for (int j = 0; j < dim*dim*dim; j++) {
				if(input2[j]==1){
				color2[i].R = (double)inputcol2[i];
				color2[i].G = (double)inputcol2[i + occsize];
				color2[i].B = (double)inputcol2[i + 2 * occsize];
				i++;
				}
			}
			//RGB-->YUV
			i=0;
			if ((int)atoi(argv[5]) == 1) {
				for (int j = 0; j < dim*dim*dim; j++) {
					if(input2[j]==1){
					color2[i].R = (double) 0.299*inputcol2[i] + 0.587*inputcol2[i + occsize] + 0.114*inputcol2[i + 2 * occsize];
					color2[i].G = (double) 0.492*(inputcol2[i + 2 * occsize] - color2[i].R);
					color2[i].B = (double) 0.877*(inputcol2[i] - color2[i].R);
					i++;
					}
				}
			}

			//make prediction
			diff = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
			clock_t t3Start = clock();
			double min = atol(argv[8]);			
			double coor = 0.0;
			int cnt=0;
			predict(color2,occupied2,color3D, input, min, diff, occsize, dim);
			predtime = predtime + (double)(clock() - t3Start) / CLOCKS_PER_SEC;
			printf("Time taken (pred): %.10fs\n", (double)(clock() - t3Start) / CLOCKS_PER_SEC);
			
			pdiff = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);			
			estim = (struct RGBvalue*)malloc(sizeof(struct RGBvalue)*occsize);
			
			//code and decode prediction residuals
			pdiff = RAHTcoder(tmp, dim, diff, step, argv[5], argv[6], occsize, meanrate, tf_time, ent_time);
			
			//inverse prediction
			invpred(pdiff,occupied2,color3D, input, min, estim, occsize,dim);

			//update old frame (backup decoded colors)
			i = 0;
			for (int j = 0; j < dim*dim*dim; j++) {
				input[j] = 0;
				color3D[j].R = 0;
				color3D[j].G = 0;
				color3D[j].B = 0;
				if (input2[j] == 1) {
					color3D[j].R = (double)estim[i].R;
					color3D[j].G = (double)estim[i].G;
					color3D[j].B = (double)estim[i].B;
					input[j] = 1;
					i++;
				}
			}
			//YUV-->RGB
			if ((int) atoi(argv[5]) == 1) {
				for (int i = 0; i < occsize; i++) {
				double Y = color2[i].R;
				double U = color2[i].G;
				double V = color2[i].B;
				color2[i].R = Y + 1.14*V;
				color2[i].G = Y - 0.395*U - 0.581*V;
				color2[i].B = Y + 2.032*U;
				}

				for (int i = 0; i < occsize; i++) {
				double Y = estim[i].R;
				double U = estim[i].G;
				double V = estim[i].B;
				estim[i].R = Y + 1.14*V;
				estim[i].G = Y - 0.395*U - 0.581*V;
				estim[i].B = Y + 2.032*U;
				}
			}
			//compute PSNR
			MSE = compMSE(color2, estim, occsize);
			PSNR = 20 * log10((double)255 / sqrt(MSE));
			printf("PSNR: %f dB\n", PSNR);
			PSNRtot = PSNRtot + PSNR;	//sum PSNRs


			/*append decoded frame to the ones already done*/
			fout = fopen("tf_color.col", "a");
			eee = fwrite(estim, sizeof(double), 3 * occsize, fout);
			fclose(fout);

			/*append original frame to the ones already done*/
			fout = fopen("orig_color.col", "a");
			eee = fwrite(color2, sizeof(double), 3 * occsize, fout);
			fclose(fout);


			free(estim);
			free(occupied2);
			free(color2);
			free(pdiff);
			free(diff);
			free(inputcol2);
			free(input2);
		}
		
	}
	printf("average prediction time:%f\n", (double)predtime / (nframes - floor(nframes / reference)));
	printf("average transform time:%f\n", (double)(*tf_time) / nframes);
	printf("average entropy time:%f\n", (double)(*ent_time) / nframes);
	printf("average rate:%f\n", (double) (*meanrate) / nframes);
	printf("average bit/voxel:%f\n", (double)(*meanrate) /tot_bit);
	printf("average PSNR:%f", PSNRtot /nframes);
	return 0;
}
