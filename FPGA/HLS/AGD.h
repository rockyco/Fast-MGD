#include <assert.h>
#include <stdint.h>
#include <iostream>
#include <hls_stream.h>
#include <ap_int.h>
#include <time.h>
#include <ap_fixed.h>
using namespace hls;

#define W 300
#define H 250
#define SIZE W*H
#define ori_BAND 224
#define PINNUM 16
#define BAND 224
#define aver_BNAD 2
#define aver 112
#define iteration 20



#define attributer 4
#define attributeMi 1
#define attributeMa 0
#define attributeR  2*attributer+1

#define gEPS  0.001
#define gr 1
#define gR 2*gr+1
#define gR2 9.0

typedef ap_uint<256> data_g;
typedef ap_ufixed<16,1,AP_RND,AP_SAT>  attribute_16nt;

typedef ap_uint<9> attribute__t;
typedef ap_uint<6> mode__1;
typedef ap_uint<3> mode__2;
typedef ap_uint<17>  fusion_dix;


typedef ap_ufixed<26,8,AP_TRN,AP_WRAP> fusion_sum;
typedef ap_ufixed<20,1,AP_TRN,AP_WRAP> fusion_var;
typedef ap_uint<16> data_16nt;



typedef ap_fixed<20,3,AP_RND,AP_WRAP> g_row;
typedef ap_ufixed<20,4,AP_RND,AP_WRAP> g_col;
typedef ap_ufixed<18,1,AP_RND,AP_WRAP> g_div;
typedef ap_ufixed<16,1,AP_RND,AP_WRAP> g_eps;


typedef ap_ufixed<29,1,AP_RND,AP_WRAP> g_2m;
typedef ap_fixed<25,3,AP_RND,AP_WRAP> gg_row;
typedef ap_ufixed<20,1,AP_RND,AP_WRAP> gg_col;
typedef ap_ufixed<23,1,AP_RND,AP_WRAP> gg_div;

typedef ap_ufixed<23,1,AP_RND,AP_WRAP> gg_MM;
typedef ap_ufixed<23,1,AP_RND,AP_WRAP> gg_V;
typedef ap_ufixed<24,1,AP_RND,AP_WRAP> gg_a;
typedef ap_ufixed<24,1,AP_RND,AP_WRAP> gg_b;
typedef ap_ufixed<24,1,AP_RND,AP_WRAP> gg_aI;
typedef ap_ufixed<32,1,AP_RND,AP_WRAP> gg_O;  //set 32 for synthesis





void parallel_AGD( data_g * image_data, gg_O *out_image);
