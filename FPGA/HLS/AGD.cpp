
#include "AGD.h"
void Data_in(data_g *img_data, stream<data_g> ImgdataBWin[aver_BNAD])
{
	fusion_dix i;
	mode__1 j;
	mode__2 m;
	unsigned int startAddr = 0;
	data_g imSrcMem[BAND/PINNUM];
#pragma HLS ARRAY_PARTITION variable=imSrcMem complete dim=1
	data_g tmpDat;

	pix_loop :
	for (i=0; i<SIZE; i++) {
		memcpy(imSrcMem, &img_data[startAddr], BAND*sizeof(attribute_16nt));
#pragma HLS PIPELINE
		Data_in_label7:for(m=0;m<aver_BNAD;m++){
			Data_in_label8:for(j=0;j<BAND/(PINNUM*aver_BNAD);j++){
				 tmpDat = imSrcMem[m*(BAND/(PINNUM*aver_BNAD))+j];
				 ImgdataBWin[m]<<(tmpDat);
			}
		}
		startAddr += BAND/PINNUM;
	}
}

void Data_sum(stream<data_g> ImgdataBWin[aver_BNAD], stream<fusion_sum> A[aver_BNAD])
{

	fusion_dix i;
	mode__1 j;
	mode__1 k;
	mode__2 m;
	data_g tempat;
	attribute_16nt subPixDat = 0;
	data_16nt subPix_u16;
	attribute_16nt mean=0;
	fusion_sum sum=0;
	Data_AVER_label3:for (i=0; i<SIZE; i++) {
#pragma HLS PIPELINE
		Data_AVER_label5:for(m=0;m<aver_BNAD;m++){
			for(j=0;j<BAND/(PINNUM*aver_BNAD);j++){
				 ImgdataBWin[m]>>(tempat);
				 Data_AVER_label9:for (k=0; k<PINNUM; k++) {
						subPix_u16 = tempat(15, 0);
						subPixDat(15, 0) = subPix_u16(15, 0);
						sum=sum+subPixDat ;
						tempat >>= 16;
				 }
			}
			 A[m]<<(sum);
			 sum=0;
		}
	}
}


void Data_AVER(stream<fusion_sum> ImgdataBWin[aver_BNAD],stream<attribute_16nt> A[aver_BNAD]
		,stream<attribute_16nt> B[aver_BNAD],
		stream<attribute_16nt> C[aver_BNAD],stream<attribute_16nt> D[aver_BNAD])
{

	fusion_dix i;
	mode__2 m;
	fusion_sum tempat;

	attribute_16nt mean=0;
	fusion_var SS=1.0/aver;
	Data_AVER_label12:for (i=0; i<SIZE; i++) {
		Data_AVER_label11:for(m=0;m<aver_BNAD;m++){
#pragma HLS UNROLL
			ImgdataBWin[m]>>(tempat);
			 mean=tempat*SS;
			 A[m]<<(mean);
			 B[m]<<(mean);
			 C[m]<<(mean);
			 D[m]<<(mean);
		}
	}
}



void windowshift(attribute_16nt window[attributeR-1],attribute_16nt tmp,attribute__t index)
{
#pragma HLS INLINE
	attribute__t j;
	windowshift_label1:for (j = 1;j<attributeR-1; j++){
#pragma HLS UNROLL
		if(index==0){
			window[j - 1] = tmp;
		}
		else{
		window[j - 1] = window[j];
		}
	}
	window[attributeR - 2] = tmp;
}

attribute_16nt com_max(attribute_16nt window[attributeR-1],attribute_16nt in)
{
#pragma HLS INLINE
	unsigned char k;
	attribute_16nt max;
	max=in;
    com_max_label24:for(k=0;k<attributeR-1;k++){
  				    	if(max<window[k])
  				    		max=window[k];
  				    }
    return max;
}

attribute_16nt com_min(attribute_16nt window[attributeR-1],attribute_16nt in)
{
#pragma HLS INLINE
	unsigned char k;
	attribute_16nt min;
	min=in;
    			for(k=0;k<attributeR-1;k++){
  				    	if(min>window[k])
  				    		min=window[k];
  				    }
    return min;
}

void col_windowshift(attribute_16nt window[attributeR-1][W],attribute_16nt tmp,
					 attribute__t index,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	col_windowshift_label2:for (k = 1;k<attributeR-1; k++){
#pragma HLS UNROLL
		if(index==0){
			window[k-1][j] = tmp;
		}
		else{
		window[k-1][j] = window[k][j];
		}
	}
	window[attributeR - 2][j] = tmp;
}

attribute_16nt col_com_max(attribute_16nt window[attributeR-1][W],attribute_16nt in,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	attribute_16nt max;
	max=in;
    col_com_max_label26:for(k=0;k<attributeR-1;k++){
  				    	if(max<window[k][j])
  				    		max=window[k][j];
  				    }
    return max;
}

attribute_16nt col_com_min(attribute_16nt window[attributeR-1][W],attribute_16nt in,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	attribute_16nt min;
	min=in;
    for(k=0;k<attributeR-1;k++){
#pragma HLS PIPELINE
  				    	if(min>window[k][j])
  				    		min=window[k][j];
  				    }
    return min;
}

void dliate_max(stream<attribute_16nt> lineImgMem[aver_BNAD],stream<attribute_16nt> ori[aver_BNAD],
		stream<attribute_16nt> res_result[aver_BNAD],stream<attribute_16nt> ori_copy[aver_BNAD] )
{
	attribute_16nt slide_win[aver_BNAD][attributeR-1];
#pragma HLS ARRAY_PARTITION variable=slide_win complete dim=0
	attribute_16nt line_mem[aver_BNAD][attributeR-1][W];
#pragma HLS ARRAY_RESHAPE variable=line_mem complete dim=2
	attribute_16nt row_Max[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=row_Max complete dim=0
	attribute_16nt line_mem_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=line_mem_tmp complete dim=0
	attribute_16nt res_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=res_tmp complete dim=0

	stream<attribute_16nt> buf_yunatu[aver_BNAD];
#pragma HLS STREAM variable=buf_yunatu depth=1205
#pragma HLS ARRAY_PARTITION variable=buf_yunatu complete dim=0
	attribute_16nt ori_tmp_1[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=ori_tmp_1 complete dim=0
	attribute_16nt ori_tmp_2[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=ori_tmp_2 complete dim=0
	attribute__t i, j;
	mode__2 m;
	for (i = 0; i<H+attributer; i++){
		for (j = 0; j<W+attributer; j++){
#pragma HLS DEPENDENCE array inter false
#pragma HLS LOOP_FLATTEN off
#pragma HLS PIPELINE
			for(m=0;m<aver_BNAD;m++){
				if(j<W&&i<H){
					lineImgMem[m]>>(line_mem_tmp[m]);
					ori[m]>>(ori_tmp_1[m]);
					buf_yunatu[m]<<(ori_tmp_1[m]);
				}
				if(i<H){
					row_Max[m]=com_max(slide_win[m],line_mem_tmp[m]);
					windowshift(slide_win[m],line_mem_tmp[m],j);
				}
				else{
					row_Max[m]=(attribute_16nt)attributeMa;
				}
				if(j>=attributer){
					res_tmp[m]=col_com_max(line_mem[m],row_Max[m],j-attributer);
					col_windowshift(line_mem[m],row_Max[m],i,j-attributer);
					if(i>=attributer){
						res_result[m]<<(res_tmp[m]);
						buf_yunatu[m]>>(ori_tmp_2[m]);
					    ori_copy[m]<<(ori_tmp_2[m]);
					}
				}
			}
		}
	}
}


void erosion_min(stream<attribute_16nt> lineImgMem[aver_BNAD],stream<attribute_16nt> ori[aver_BNAD],
		stream<attribute_16nt> res_result[aver_BNAD],stream<attribute_16nt> ori_copy[aver_BNAD] )
{
	attribute_16nt slide_win[aver_BNAD][attributeR-1];
#pragma HLS ARRAY_PARTITION variable=slide_win complete dim=0
	attribute_16nt line_mem[aver_BNAD][attributeR-1][W];
#pragma HLS ARRAY_RESHAPE variable=line_mem complete dim=2
	attribute_16nt copy_RAM[aver_BNAD][attributeR-1][W];
#pragma HLS ARRAY_RESHAPE variable=line_mem complete dim=2
	attribute_16nt row_Min[aver_BNAD];
	attribute_16nt line_mem_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=line_mem_tmp complete dim=0
	attribute_16nt res_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=res_tmp complete dim=0
//	attribute_16nt ori_tmp[aver_BNAD];
//#pragma HLS ARRAY_PARTITION variable=ori_tmp complete dim=0

	attribute_16nt ori_tmp_1[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=ori_tmp_1 complete dim=0
	attribute_16nt ori_tmp_2[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=ori_tmp_2 complete dim=0

	stream<attribute_16nt> buf_yunatu[aver_BNAD];
#pragma HLS STREAM variable=buf_yunatu depth=1205  
#pragma HLS ARRAY_PARTITION variable=buf_yunatu complete dim=0

	attribute__t i, j;
	mode__2 m;
	for (i = 0; i<H+attributer; i++){
		for (j = 0; j<W+attributer; j++){
#pragma HLS DEPENDENCE array inter false
#pragma HLS LOOP_FLATTEN off
#pragma HLS PIPELINE
			for(m=0;m<aver_BNAD;m++){
				if(j<W&&i<H){
					lineImgMem[m]>>(line_mem_tmp[m]);
					ori[m]>>(ori_tmp_1[m]);
					buf_yunatu[m]<<(ori_tmp_1[m]);
				}
				if(i<H){
					row_Min[m]=com_min(slide_win[m],line_mem_tmp[m]);
					windowshift(slide_win[m],line_mem_tmp[m],j);
				}
				else{
					row_Min[m]=(attribute_16nt)attributeMi;
				}
				if(j>=attributer){
					res_tmp[m]=col_com_min(line_mem[m],row_Min[m],j-attributer);
					col_windowshift(line_mem[m],row_Min[m],i,j-attributer);
					if(i>=attributer){
					res_result[m]<<(res_tmp[m]);
					buf_yunatu[m]>>(ori_tmp_2[m]);
				    ori_copy[m]<<(ori_tmp_2[m]);
					}
				}
			}
		}
	}
}
void compare_min(stream<attribute_16nt> data_1[aver_BNAD],stream<attribute_16nt> yuan[aver_BNAD],
				 stream<attribute_16nt> result[aver_BNAD], stream<attribute_16nt> ori[aver_BNAD]){
	   attribute_16nt I_tmp;
	   attribute_16nt a_tmp;
	   attribute_16nt b_tmp;
	   attribute__t i,j;
	   mode__2 m;
	   for (i = 0; i<H; i++){
			 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
				 for(m=0;m<aver_BNAD;m++){
					 data_1[m]>>(I_tmp);
					 yuan[m]>>(a_tmp);

					 if(I_tmp<=a_tmp){
					  result[m]<<(I_tmp);
					  ori[m]<<a_tmp;
					 }
					 else{
					  result[m]<<(a_tmp);
					  ori[m]<<a_tmp;
					 }
				 }
			 }
	   }
}
void compare_min_last(stream<attribute_16nt> data_1[aver_BNAD],stream<attribute_16nt> yuan[aver_BNAD],
		         stream<attribute_16nt> result[aver_BNAD])
{
	   attribute_16nt I_tmp;
	   attribute_16nt a_tmp;
	   attribute_16nt b_tmp;
	   attribute__t i,j;
	   mode__2 m;
	   for (i = 0; i<H; i++){
			 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
				 for(m=0;m<aver_BNAD;m++){
					 data_1[m]>>(I_tmp);
					 yuan[m]>>(a_tmp);
					 if(I_tmp<=a_tmp){
					  result[m]<<(I_tmp);
					 }
					 else{
					  result[m]<<(a_tmp);
					 }
				 }
			 }
	   }
}

void compare_max(stream<attribute_16nt> data_1[aver_BNAD],stream<attribute_16nt> yuan[aver_BNAD],
				 stream<attribute_16nt> result[aver_BNAD], stream<attribute_16nt> ori[aver_BNAD]){
	   attribute_16nt I_tmp;
	   attribute_16nt a_tmp;
	   attribute_16nt b_tmp;
	   attribute__t i,j;
	   mode__2 m;
	   for (i = 0; i<H; i++){
			 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
				 for(m=0;m<aver_BNAD;m++){
					 data_1[m]>>(I_tmp);
					 yuan[m]>>(a_tmp);
					 if(I_tmp<=a_tmp){
					  result[m]<<(a_tmp);
					  ori[m]<<a_tmp;
					 }
					 else{
					  result[m]<<(I_tmp);
					  ori[m]<<a_tmp;
					 }
				 }
			 }
	   }
}
void compare_max_last(stream<attribute_16nt> data_1[aver_BNAD],stream<attribute_16nt> yuan[aver_BNAD],
		         stream<attribute_16nt> result[aver_BNAD])
{
	   attribute_16nt I_tmp;
	   attribute_16nt a_tmp;
	   attribute_16nt b_tmp;
	   attribute__t i,j;
	   mode__2 m;
	   for (i = 0; i<H; i++){
			 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
				 for(m=0;m<aver_BNAD;m++){
					 data_1[m]>>(I_tmp);
					 yuan[m]>>(a_tmp);
					 if(I_tmp<=a_tmp){
					  result[m]<<(a_tmp);
					 }
					 else{
					  result[m]<<(I_tmp);
					 }
				 }
			 }
	   }
}
void sub(stream<attribute_16nt> data_max[aver_BNAD],stream<attribute_16nt> data_min[aver_BNAD],
		         stream<attribute_16nt> result[aver_BNAD])
{
	   attribute_16nt I_tmp;
	   attribute_16nt a_tmp;
	   attribute_16nt b_tmp;
	   attribute__t i,j;
	   mode__2 m;
	   for (i = 0; i<H; i++){
			 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
				 for(m=0;m<aver_BNAD;m++){
					 data_max[m]>>(I_tmp);
					 data_min[m]>>(a_tmp);
					 b_tmp=I_tmp-a_tmp;
					 result[m]<<b_tmp;
				 }
			 }
	   }
}

void morph_open(stream<attribute_16nt> in[aver_BNAD],stream<attribute_16nt> yuan[aver_BNAD]
				,stream<attribute_16nt> out[aver_BNAD]){
#pragma HLS DATAFLOW
	stream<attribute_16nt> mid1[iteration][aver_BNAD];
#pragma HLS STREAM variable=mid1 depth=1 dim=2
	stream<attribute_16nt> ori_1[iteration][aver_BNAD];
#pragma HLS STREAM variable=ori_1 depth=1 dim=2

	stream<attribute_16nt> mid2[iteration][aver_BNAD];
#pragma HLS STREAM variable=mid2 depth=1 dim=2
	stream<attribute_16nt> ori_2[iteration][aver_BNAD];
#pragma HLS STREAM variable=ori_2 depth=1 dim=2

	stream<attribute_16nt> mid3[aver_BNAD];
#pragma HLS STREAM variable=mid3 depth=1

	stream<attribute_16nt> ori_3[aver_BNAD];
#pragma HLS STREAM variable=ori_3 depth=1


	erosion_min(in,yuan,mid1[0],ori_1[0]);
	for(unsigned char i=0;i<iteration-1;i++){
#pragma HLS UNROLL
		dliate_max(mid1[i],ori_1[i],mid2[i],ori_2[i]);
		compare_min(mid2[i],ori_2[i],mid1[i+1],ori_1[i+1]);
	}
	 dliate_max(mid1[iteration-1],ori_1[iteration-1],mid3,ori_3);
	 compare_min_last(mid3,ori_3,out);

}

void morph_close(stream<attribute_16nt> in[aver_BNAD],stream<attribute_16nt> yuan[aver_BNAD]
				,stream<attribute_16nt> out[aver_BNAD]){
#pragma HLS DATAFLOW
	stream<attribute_16nt> mid1[iteration][aver_BNAD];
#pragma HLS STREAM variable=mid1 depth=1 dim=2
	stream<attribute_16nt> ori_1[iteration][aver_BNAD];
#pragma HLS STREAM variable=ori_1 depth=1 dim=2

	stream<attribute_16nt> mid2[iteration][aver_BNAD];
#pragma HLS STREAM variable=mid2 depth=1 dim=2
	stream<attribute_16nt> ori_2[iteration][aver_BNAD];
#pragma HLS STREAM variable=ori_2 depth=1 dim=2

	stream<attribute_16nt> mid3[aver_BNAD];
#pragma HLS STREAM variable=mid3 depth=1

	stream<attribute_16nt> ori_3[aver_BNAD];
#pragma HLS STREAM variable=ori_3 depth=1


	dliate_max(in,yuan,mid1[0],ori_1[0]);
	for(unsigned char i=0;i<iteration-1;i++){
#pragma HLS UNROLL
		erosion_min(mid1[i],ori_1[i],mid2[i],ori_2[i]);
		compare_max(mid2[i],ori_2[i],mid1[i+1],ori_1[i+1]);
	}
	erosion_min(mid1[iteration-1],ori_1[iteration-1],mid3,ori_3);
	compare_max_last(mid3,ori_3,out);

}

void morph_recon(stream<attribute_16nt> A[aver_BNAD],stream<attribute_16nt> B[aver_BNAD],
				stream<attribute_16nt> C[aver_BNAD],stream<attribute_16nt> D[aver_BNAD]
				,stream<attribute_16nt> out[aver_BNAD]){

#pragma HLS DATAFLOW

	stream<attribute_16nt> k_ero_res[aver_BNAD];
	#pragma HLS STREAM variable=k_ero_res depth=1
	stream<attribute_16nt> d_ero_res[aver_BNAD];
	#pragma HLS STREAM variable=d_ero_res depth=1 dim=1
	morph_open(A,B,k_ero_res);
	morph_close(C,D,d_ero_res);
	sub(d_ero_res,k_ero_res,out);


}

g_row com_SUM(attribute_16nt window[gR-1],attribute_16nt in)
{
#pragma HLS INLINE
	unsigned char k;
	g_row SUM;
	SUM=in;
    com_max_label24:for(k=0;k<gR-1;k++){
  				    	SUM=SUM+window[k];
  				    }
    return SUM;
}

void g_windowshift(attribute_16nt window[gR-1],attribute_16nt tmp,attribute__t index)
{
#pragma HLS INLINE
	attribute__t j;
	for (j = 1;j<gR-1; j++){
#pragma HLS UNROLL
		if(index==0){
			window[j - 1] = (attribute_16nt)0;
		}
		else{
		window[j - 1] = window[j];
		}
	}
	window[gR - 2] = tmp;
}


g_col col_com_SUM(g_row window[gR-1][W],g_row in,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	g_col SUM;
	SUM=in;
    for(k=0;k<gR-1;k++){
    			SUM=SUM+window[k][j];
    }
    return SUM;
}


void g_col_windowshift(g_row window[gR-1][W],g_row tmp,
					 attribute__t index,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	for (k = 1;k<gR-1; k++){
#pragma HLS UNROLL
		if(index==0){
			window[k-1][j] = (g_row)0;
		}
		else{
		window[k-1][j] = window[k][j];
		}
	}
	window[gR - 2][j] = tmp;
}

void boxfilter_sum(stream<attribute_16nt> lineImgMem[aver_BNAD],
		stream<attribute_16nt> lineImgMem_copy[aver_BNAD],stream<g_2m> lineImgMem_II[aver_BNAD],
		stream<g_col> col_result[aver_BNAD]){

	attribute_16nt slide_win[aver_BNAD][gR-1];
#pragma HLS ARRAY_PARTITION variable=slide_win complete dim=0
	g_row line_mem[aver_BNAD][gR-1][W];
#pragma HLS ARRAY_RESHAPE variable=line_mem complete dim=2
	g_row row_Sum[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=row_Sum complete dim=0
	attribute_16nt line_mem_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=line_mem_tmp complete dim=0
	g_col res_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=res_tmp complete dim=0


		attribute_16nt ori_tmp_2[aver_BNAD];
	#pragma HLS ARRAY_PARTITION variable=ori_tmp_2 complete dim=0

		stream<attribute_16nt> buf_yunatu[aver_BNAD];
	#pragma HLS STREAM variable=buf_yunatu depth=302  
	#pragma HLS ARRAY_PARTITION variable=buf_yunatu complete dim=0
	attribute__t i, j,k;
	mode__2 m;
	for (i = 0; i<H+gr; i++){
		for (j = 0; j<W+gr; j++){
#pragma HLS DEPENDENCE array inter false
#pragma HLS LOOP_FLATTEN off
#pragma HLS PIPELINE
			for(m=0;m<aver_BNAD;m++){
				if(j<W&&i<H){
					lineImgMem[m]>>(line_mem_tmp[m]);
					buf_yunatu[m]<<(line_mem_tmp[m]);
				}
				if(i<H){
					row_Sum[m]=com_SUM(slide_win[m],line_mem_tmp[m]);
					g_windowshift(slide_win[m],line_mem_tmp[m],j);
				}
				else{
					row_Sum[m]=(g_row)0;
				}
				if(j>=gr){
					res_tmp[m]=col_com_SUM(line_mem[m],row_Sum[m],j-gr);
					g_col_windowshift(line_mem[m],row_Sum[m],i,j-gr);
					if(i>=gr){
						col_result[m]<<(res_tmp[m]);
						buf_yunatu[m]>>(ori_tmp_2[m]);
						lineImgMem_copy[m]<<(ori_tmp_2[m]);
						lineImgMem_II[m]<<g_2m(ori_tmp_2[m]*ori_tmp_2[m]);
					}
				}
			}
		}
	}
}

void box_divide(stream<attribute_16nt> g_in_copy[aver_BNAD],stream<g_2m> II_1[aver_BNAD],
		stream<g_col> img_col_result[aver_BNAD],stream<attribute_16nt> g_in1[aver_BNAD],
		stream<g_2m> g_II[aver_BNAD],stream<g_div> g_MI[aver_BNAD])
{
	g_col img_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=img_tmp complete dim=0
	attribute_16nt yunatu[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=yunatu complete dim=0
	g_2m II_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=II_tmp complete dim=0
	g_div divide_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=divide_tmp complete dim=0
   g_eps invers_tmp;
   attribute__t i,j;
   mode__2 m;
   invers_tmp = 1.0 / gR2;
   for (i=0; i<H; i++) {
	   for(j=0;j<W;j++){
#pragma HLS PIPELINE
		   for(m=0;m<aver_BNAD;m++){
			   img_col_result[m]>>(img_tmp[m]);
			   g_in_copy[m]>>(yunatu[m]);
			   II_1[m]>>(II_tmp[m]);
			   divide_tmp[m] = img_tmp[m]*invers_tmp;

			   g_MI[m]<<(divide_tmp[m]);
			   g_in1[m]<<(yunatu[m]);
			   g_II[m]<<(II_tmp[m]);
		   }
	 }
   }
}


gg_row com_SUM_II(g_2m window[gR-1],g_2m in)
{
#pragma HLS INLINE
	unsigned char k;
	gg_row SUM;
	SUM=in;
    com_max_label24:for(k=0;k<gR-1;k++){
  				    	SUM=SUM+window[k];
  				    }
    return SUM;
}

void g_windowshift_II(g_2m window[gR-1],g_2m tmp,attribute__t index)
{
#pragma HLS INLINE
	attribute__t j;
	for (j = 1;j<gR-1; j++){
#pragma HLS UNROLL
		if(index==0){
			window[j - 1] = (g_2m)0;
		}
		else{
		window[j - 1] = window[j];
		}
	}
	window[gR - 2] = tmp;
}


gg_col col_com_SUM_II(gg_row window[gR-1][W],gg_row in,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	gg_col SUM;
	SUM=in;
    for(k=0;k<gR-1;k++){
    			SUM=SUM+window[k][j];
    }
    return SUM;
}


void g_col_windowshift_II(gg_row window[gR-1][W],gg_row tmp,
					 attribute__t index,attribute__t j)
{
#pragma HLS INLINE
	attribute__t k;
	for (k = 1;k<gR-1; k++){
#pragma HLS UNROLL
		if(index==0){
			window[k-1][j] = (gg_row)0;
		}
		else{
		window[k-1][j] = window[k][j];
		}
	}
	window[gR - 2][j] = tmp;
}

void boxfilter_sum_II(stream<attribute_16nt> g_in1[aver_BNAD],
		stream<g_2m> g_II[aver_BNAD],stream<g_div> g_MI[aver_BNAD],
		stream<attribute_16nt> g_in_copy[aver_BNAD],
		stream<gg_col> col_result[aver_BNAD],stream<g_div> g_MI_copy[aver_BNAD]
	    ,stream<gg_MM> g_MM[aver_BNAD]

		){

	g_2m slide_win[aver_BNAD][gR-1];
#pragma HLS ARRAY_PARTITION variable=slide_win complete dim=0
	gg_row line_mem[aver_BNAD][gR-1][W];
#pragma HLS ARRAY_RESHAPE variable=line_mem complete dim=2
	gg_row row_Sum[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=row_Sum complete dim=0
	g_2m line_mem_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=line_mem_tmp complete dim=0
	gg_col res_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=res_tmp complete dim=0

	attribute_16nt ori_tmp_1[aver_BNAD];
	#pragma HLS ARRAY_PARTITION variable=ori_tmp_1 complete dim=0
		attribute_16nt ori_tmp_2[aver_BNAD];
	#pragma HLS ARRAY_PARTITION variable=ori_tmp_2 complete dim=0

		g_div MI_tmp_1[aver_BNAD];
		#pragma HLS ARRAY_PARTITION variable=ori_tmp_1 complete dim=0
		g_div MI_tmp_2[aver_BNAD];
		#pragma HLS ARRAY_PARTITION variable=ori_tmp_2 complete dim=0

		stream<attribute_16nt> buf_yunatu[aver_BNAD];
	#pragma HLS STREAM variable=buf_yunatu depth=302  
	#pragma HLS ARRAY_PARTITION variable=buf_yunatu complete dim=0
		stream<g_div> buf_MI[aver_BNAD];
	#pragma HLS STREAM variable=buf_MI depth=302  
	#pragma HLS ARRAY_PARTITION variable=buf_MI complete dim=0
	attribute__t i, j,k;
	mode__2 m;
	for (i = 0; i<H+gr; i++){
		for (j = 0; j<W+gr; j++){
#pragma HLS DEPENDENCE array inter false
#pragma HLS LOOP_FLATTEN off
#pragma HLS PIPELINE
			for(m=0;m<aver_BNAD;m++){
				if(j<W&&i<H){
					g_II[m]>>(line_mem_tmp[m]);
					g_in1[m]>>(ori_tmp_1[m]);
					buf_yunatu[m]<<(ori_tmp_1[m]);
					g_MI[m]>>(MI_tmp_1[m]);
					buf_MI[m]<<(MI_tmp_1[m]);
				}
				if(i<H){
					row_Sum[m]=com_SUM_II(slide_win[m],line_mem_tmp[m]);
					g_windowshift_II(slide_win[m],line_mem_tmp[m],j);
				}
				else{
					row_Sum[m]=(gg_row)0;
				}
				if(j>=gr){
					res_tmp[m]=col_com_SUM_II(line_mem[m],row_Sum[m],j-gr);
					g_col_windowshift_II(line_mem[m],row_Sum[m],i,j-gr);
					if(i>=gr){
						col_result[m]<<(res_tmp[m]);
						buf_yunatu[m]>>(ori_tmp_2[m]);
						g_in_copy[m]<<(ori_tmp_2[m]);
						buf_MI[m]>>(MI_tmp_2[m]);
						g_MI_copy[m]<<MI_tmp_2[m];
						g_MM[m]<<(gg_MM)(MI_tmp_2[m]*MI_tmp_2[m]);
					}
				}
			}
		}
	}
}

void box_divide_II(stream<attribute_16nt> g_in_copy[aver_BNAD],stream<gg_col> img_col_result[aver_BNAD],
		stream<g_div> g_MI_copy[aver_BNAD],stream<gg_MM> g_MM_copy[aver_BNAD],
		stream<attribute_16nt> g_in2[aver_BNAD],
		stream<gg_div> g_C[aver_BNAD],stream<g_div> g_MI_copy1[aver_BNAD]
		,stream<gg_MM> g_MM_copy1[aver_BNAD])
{
	gg_col img_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=img_tmp complete dim=0
	attribute_16nt yunatu[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=yunatu complete dim=0

	g_div MI_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=MI_tmp complete dim=0
	gg_MM MM_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=MM_tmp complete dim=0

	gg_div divide_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=divide_tmp complete dim=0
   g_eps invers_tmp;
   attribute__t i,j;
   mode__2 m;
   invers_tmp = 1.0 / gR2;//æ”?
   for (i=0; i<H; i++) {
	   for(j=0;j<W;j++){
#pragma HLS PIPELINE
		   for(m=0;m<aver_BNAD;m++){
			   img_col_result[m]>>(img_tmp[m]);
			   g_in_copy[m]>>(yunatu[m]);
			   g_MI_copy[m]>>(MI_tmp[m]);
			   g_MM_copy[m]>>(MM_tmp[m]);
			   divide_tmp[m] = img_tmp[m]*invers_tmp; 

			   g_C[m]<<(divide_tmp[m]);
			   g_in2[m]<<(yunatu[m]);
			   g_MI_copy1[m]<<(MI_tmp[m]);
			   g_MM_copy1[m]<<(MM_tmp[m]);
		   }
	 }
   }
}



void guide_stage_1(stream<attribute_16nt> g_in[aver_BNAD]
				,stream<attribute_16nt> g_in1[aver_BNAD],stream<g_2m> g_II[aver_BNAD],stream<g_div> g_MI[aver_BNAD]){
#pragma HLS DATAFLOW

	stream<attribute_16nt> g_in_copy[aver_BNAD];
#pragma HLS STREAM variable=g_in_copy depth=1 dim=1
	stream<g_2m> II_1[aver_BNAD];
#pragma HLS STREAM variable=II_1 depth=1 dim=1

	stream<g_col> img_col_result[aver_BNAD];
#pragma HLS STREAM variable=img_col_result depth=1 dim=1


	boxfilter_sum(g_in,g_in_copy,II_1,img_col_result);
	box_divide(g_in_copy,II_1,img_col_result,g_in1,g_II,g_MI);
}

void guide_stage_2(stream<attribute_16nt> g_in1[aver_BNAD],
		stream<g_2m> g_II[aver_BNAD],stream<g_div> g_MI[aver_BNAD],
		stream<attribute_16nt> g_in2[aver_BNAD],stream<gg_div> g_C[aver_BNAD],
		stream<g_div> g_MI_1[aver_BNAD],stream<gg_MM> g_MM[aver_BNAD]){
#pragma HLS DATAFLOW

	stream<attribute_16nt> g_in_copy[aver_BNAD];
#pragma HLS STREAM variable=g_in_copy depth=1 dim=1

	stream<g_div> g_MI_copy[aver_BNAD];
#pragma HLS STREAM variable=g_MI_copy depth=1 dim=1
	stream<gg_MM> g_MM_copy[aver_BNAD];
#pragma HLS STREAM variable=g_MM_copy depth=1 dim=1

	stream<gg_col> img_col_result[aver_BNAD];
#pragma HLS STREAM variable=img_col_result depth=1 dim=1


	boxfilter_sum_II(g_in1,g_II,g_MI,g_in_copy,img_col_result,g_MI_copy,g_MM_copy);
	box_divide_II(g_in_copy,img_col_result,g_MI_copy,g_MM_copy,g_in2,g_C,g_MI_1,g_MM);
}

void guide_stage_3(stream<attribute_16nt> g_in2[aver_BNAD],stream<gg_div> g_C[aver_BNAD],
		stream<g_div> g_MI_1[aver_BNAD],stream<gg_MM> g_MM[aver_BNAD],
		stream<attribute_16nt> g_in3[aver_BNAD],stream<gg_V> g_V[aver_BNAD],
		stream<gg_V> g_VE[aver_BNAD],stream<g_div> g_MI_2[aver_BNAD]){

	attribute_16nt in2_tmp;
	gg_div C_tmp;
	g_div MI_tmp;
	gg_MM MM_tmp;
	gg_V V_tmp;
	gg_V VE_tmp;
		   attribute__t i,j;
		   mode__2 m;
		   for (i = 0; i<H; i++){
				 for(j=0;j<W;j++){
	#pragma HLS PIPELINE II=1
					 for(m=0;m<aver_BNAD;m++){
						 g_in2[m]>>(in2_tmp);
						 g_C[m]>>(C_tmp);
						 g_MI_1[m]>>(MI_tmp);
						 g_MM[m]>>(MM_tmp);
						 g_in3[m]<<in2_tmp;
						 V_tmp=(gg_V)C_tmp-MM_tmp;
						 VE_tmp=C_tmp-MM_tmp+(gg_V)gEPS;
						 g_V[m]<<V_tmp;
						 g_VE[m]<<VE_tmp;
						 g_MI_2[m]<<MI_tmp;
					 }
				 }
		   }

}

void guide_stage_4(stream<attribute_16nt> g_in3[aver_BNAD],stream<gg_V> g_V[aver_BNAD],
		stream<gg_V> g_VE[aver_BNAD],stream<g_div> g_MI_2[aver_BNAD],
		stream<attribute_16nt> g_in4[aver_BNAD],stream<gg_a> g_a[aver_BNAD],
		stream<g_div> g_MI_3[aver_BNAD]){

	attribute_16nt in3_tmp;
	gg_V V_tmp;
	gg_V VE_tmp;
	g_div MI2_tmp;
	gg_a a_tmp;
		   attribute__t i,j;
		   mode__2 m;
		   for (i = 0; i<H; i++){
				 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
					 for(m=0;m<aver_BNAD;m++){
						 g_in3[m]>>(in3_tmp);
						 g_V[m]>>(V_tmp);
						 g_VE[m]>>(VE_tmp);
						 g_MI_2[m]>>(MI2_tmp);
						 a_tmp=(gg_a)V_tmp/VE_tmp;
						 g_in4[m]<<in3_tmp;
						 g_a[m]<<a_tmp;
						 g_MI_3[m]<<MI2_tmp;
					 }
				 }
		   }

}

void guide_stage_5(stream<attribute_16nt> g_in4[aver_BNAD],stream<gg_a> g_a[aver_BNAD],
		stream<g_div> g_MI_3[aver_BNAD],stream<gg_b> g_b[aver_BNAD],stream<gg_aI> g_aI[aver_BNAD]){

	attribute_16nt in4_tmp;
	gg_a a_tmp;
	gg_b b_tmp;
	g_div MI3_tmp;
	gg_aI aI_tmp;
		   attribute__t i,j;
		   mode__2 m;
		   for (i = 0; i<H; i++){
				 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
					 for(m=0;m<aver_BNAD;m++){
						 g_in4[m]>>(in4_tmp);
						 g_a[m]>>(a_tmp);
						 g_MI_3[m]>>(MI3_tmp);
						 b_tmp=(1-a_tmp)*MI3_tmp;
						 aI_tmp=in4_tmp*a_tmp;
						 g_b[m]<<b_tmp;
						 g_aI[m]<<aI_tmp;
					 }
				 }
		   }

}

void guide_stage_6(stream<gg_b> g_b[aver_BNAD],stream<gg_aI> g_aI[aver_BNAD],stream<gg_O> g_O[aver_BNAD]){

	gg_b b_tmp;
	gg_aI aI_tmp;
	gg_O O_tmp;
		   attribute__t i,j;
		   mode__2 m;
		   for (i = 0; i<H; i++){
				 for(j=0;j<W;j++){
#pragma HLS PIPELINE II=1
					 for(m=0;m<aver_BNAD;m++){
						 g_b[m]>>(b_tmp);
						 g_aI[m]>>(aI_tmp);
						 O_tmp=b_tmp+aI_tmp;
						 g_O[m]<<O_tmp;
					 }
				 }
		   }

}
void guide_filter(stream<attribute_16nt> g_in[aver_BNAD]
				,stream<gg_O> g_O[aver_BNAD]){
#pragma HLS DATAFLOW

	stream<attribute_16nt> g_in1[aver_BNAD];
	#pragma HLS STREAM variable=g_in1 depth=1
	stream<g_2m> g_II[aver_BNAD];
	#pragma HLS STREAM variable=g_II depth=1
	stream<g_div> g_MI[aver_BNAD];
	#pragma HLS STREAM variable=g_MI depth=1

	stream<attribute_16nt> g_in2[aver_BNAD];
	#pragma HLS STREAM variable=g_in2 depth=1
	stream<gg_div> g_C[aver_BNAD];
	#pragma HLS STREAM variable=g_C depth=1
	stream<g_div> g_MI_1[aver_BNAD];
	#pragma HLS STREAM variable=g_MI_1 depth=1
	stream<gg_MM> g_MM[aver_BNAD];
	#pragma HLS STREAM variable=g_MM depth=1

	stream<attribute_16nt> g_in3[aver_BNAD];
	#pragma HLS STREAM variable=g_in3 depth=1
	stream<gg_V> g_V[aver_BNAD];
	#pragma HLS STREAM variable=g_V depth=1
	stream<gg_V> g_VE[aver_BNAD];
	#pragma HLS STREAM variable=g_V depth=1
	stream<g_div> g_MI_2[aver_BNAD];
	#pragma HLS STREAM variable=g_MI_2 depth=1

	stream<attribute_16nt> g_in4[aver_BNAD];
	#pragma HLS STREAM variable=g_in4 depth=1
	stream<gg_a> g_a[aver_BNAD];
	#pragma HLS STREAM variable=g_a depth=1
	stream<g_div> g_MI_3[aver_BNAD];
	#pragma HLS STREAM variable=g_MI_3 depth=1

	stream<gg_b> g_b[aver_BNAD];
	#pragma HLS STREAM variable=g_b depth=1
	stream<gg_aI> g_aI[aver_BNAD];
	#pragma HLS STREAM variable=g_aI depth=1

	guide_stage_1(g_in,g_in1,g_II,g_MI);
	guide_stage_2(g_in1,g_II,g_MI,g_in2,g_C,g_MI_1,g_MM);
	guide_stage_3(g_in2,g_C,g_MI_1,g_MM,g_in3,g_V,g_VE,g_MI_2);
	guide_stage_4(g_in3,g_V,g_VE,g_MI_2,g_in4,g_a,g_MI_3);
	guide_stage_5(g_in4,g_a,g_MI_3,g_b,g_aI);
	guide_stage_6(g_b,g_aI,g_O);

}
void fusion_pro(stream<gg_O> g_O[aver_BNAD],stream<gg_O> &result){

	gg_O O_tmp[aver_BNAD];
#pragma HLS ARRAY_PARTITION variable=O_tmp complete dim=0
	attribute__t i,j;
	mode__2 m;
	gg_O result_tmp;
	gg_O result_tmp1;
		   for (i = 0; i<H; i++){
				 for(j=0;j<W;j++){
#pragma HLS PIPELINE
					 for(m=0;m<aver_BNAD;m++){
						 g_O[m]>>(O_tmp[m]);
				 }
					 result_tmp=O_tmp[0]+O_tmp[1];
					 result_tmp1=result_tmp*(gg_O)0.5;
					 result<<result_tmp1;
				 }
		   }

}
void store(stream<gg_O>&result,gg_O *out_image)
{
	attribute__t i,j;
	int startAddr=0;
	bool flag;
//	data_16nt tempat;
	gg_O tem;
	for(i=0;i<H;i++){
		for(j=0;j<W;j++){
#pragma HLS PIPELINE
			result>>tem;
			out_image[i*W+j]=tem;
		}
	}
}
void parallel_AGD( data_g * image_data, gg_O *out_image)
{
#pragma HLS INTERFACE m_axi depth=1050000 port=image_data offset=slave bundle=gmem0  
#pragma HLS INTERFACE m_axi depth=75000 port=out_image offset=slave bundle=gmem1
#pragma HLS INTERFACE s_axilite port=image_data bundle=control
#pragma HLS INTERFACE s_axilite port=out_image bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

#pragma HLS DATAFLOW
	stream<data_g> ImgdataBW_in[aver_BNAD];
	#pragma HLS STREAM variable=ImgdataBW_in depth=1
	stream<fusion_sum> g_result[aver_BNAD];
	#pragma HLS STREAM variable=g_result depth=1

	stream<attribute_16nt> lineImgMemA[aver_BNAD];
	#pragma HLS STREAM variable=lineImgMemA depth=1
	stream<attribute_16nt> lineImgMemB[aver_BNAD];
	#pragma HLS STREAM variable=lineImgMemB depth=1
	stream<attribute_16nt> lineImgMemC[aver_BNAD];
	#pragma HLS STREAM variable=lineImgMemC depth=1 dim=1
	stream<attribute_16nt> lineImgMemD[aver_BNAD];
	#pragma HLS STREAM variable=lineImgMemD depth=1 dim=1
	stream<attribute_16nt> g_lineImgMemA[aver_BNAD];
	#pragma HLS STREAM variable=g_lineImgMemA depth=1

	stream<gg_O> g_O[aver_BNAD];
	#pragma HLS STREAM variable=g_O depth=1

	stream<gg_O> pro_result("output_stream");
	#pragma HLS STREAM variable=g_result depth=1

	Data_in(image_data, ImgdataBW_in);
	Data_sum(ImgdataBW_in, g_result);
	Data_AVER(g_result, lineImgMemA,lineImgMemB,lineImgMemC,lineImgMemD);
	morph_recon(lineImgMemA,lineImgMemB,lineImgMemC,lineImgMemD,g_lineImgMemA);
	guide_filter(g_lineImgMemA,g_O);
	fusion_pro(g_O,pro_result);
	store(pro_result,out_image);
}


