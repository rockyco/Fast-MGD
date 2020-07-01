#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include<ctime>
#include "AGD.h"

bool saveFile_band(const char* name, gg_O *result);
bool loadFile(const char* name1, attribute_16nt *img_data_in);
int main(void)
{
	attribute__t j,k;

	attribute_16nt * src_img=new attribute_16nt[SIZE*BAND];
	gg_O * out_img=new gg_O[SIZE*aver_BNAD];
	clock_t start,finish;
	double totaltime;
	loadFile("Segundo.txt",src_img);
	start = clock();
	parallel_AGD((data_g *)src_img,out_img);
	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	std::cout << "\n runing time is£º" << totaltime << "Ãë£¡" << std::endl;
	saveFile_band("hls_aver.txt", out_img);
    delete [] out_img;
	delete [] src_img;
	    return 0;
}

bool loadFile(const char* name1, attribute_16nt *img_data_in){
	FILE* fp1 = fopen(name1, "rb");
	int i = 0;
	int j = 0;
	double temp;
	if ((fp1 == NULL)) {
		std::cout << "Load Error!" << std::endl;
		return false;
	}
	for (i = 0; i<SIZE; i++) {
		for (j = 0; j<BAND; j++) {
			if(j<ori_BAND){
				fscanf(fp1, "%lf", &temp);
				img_data_in[i*BAND+j] = (attribute_16nt)temp;
			}
			else{
				img_data_in[i*BAND+j] = (attribute_16nt)0;
			}
		}
	}
	fflush(fp1);
	fclose(fp1);
	std::cout << "Load Success!" << std::endl;
	return true;
}
bool saveFile_band(const char* name, gg_O *result){
	FILE* fp = fopen(name, "wb");
	int j = 0;
	double res_temp;
	if (fp == NULL) {
		std::cout << "Save Error!" << std::endl;
		return false;
	}
	for (int i = 0; i<H; i++) {
		for (int j = 0; j<W; j++) {
				res_temp = (double)result[i*W+j];
				fprintf(fp, "%lf\n", res_temp);
		}
	}
	fflush(fp);
	fclose(fp);
	std::cout << "Save Success!" << std::endl;
	return true;
}

