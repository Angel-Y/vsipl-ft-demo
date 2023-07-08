//一维复数out-of-place FFT: vsip_ccfftop_f()
//测试程序

#include "vsip.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

//views of blocks
vsip_cvview_f *vector_x;
vsip_cvview_f *vector_y;
vsip_cvview_f *vector_z;

//print first of 10 number in vector
void VU_cvprint_f(vsip_cvview_f* a){
   vsip_length i, n;
   vsip_cscalar_f x;
   n = vsip_cvgetlength_f(a);
   for(i=0; i<n; i++){
        x = vsip_cvget_f(a,i);
        printf("%5.3f+%5.3fi\t",x.r,x.i);
   }
   printf("\n");
   return;
}

//out-of-place fft / ifft plans
vsip_fft_f      *fft_op_plan;
vsip_fft_f      *ifft_op_plan;

int repeat = 10;

//time
struct timeval tv_b_cfar,tv_e_cfar,tv1,tv2;
double time_cfar_used;

//initialize data_x & data_y
void initialize_test_data(vsip_cvview_f *a)
{
	srand(time(NULL));
	vsip_length n = vsip_cvgetlength_f(a);
    vsip_cscalar_f x;
    for(int i=0;i<n;i++)
    {
        x.r=((float)(rand() % 100)) / 100.0f;
        x.i=((float)(rand() % 100)) / 100.0f; 
        // x.r=i+1;
        // x.i=0;
        vsip_cvput_f(a,i,x);
    }
}


void initialize(vsip_length N)
{
    vsip_init(NULL);

    vector_x = vsip_cvcreate_f(N,0);
    vector_y = vsip_cvcreate_f(N,0);
    vector_z = vsip_cvcreate_f(N,0);

	printf("vector_x attribute:\n");
	printf("vector_x->length = %d\n", vsip_cvgetlength_f(vector_x));
	printf("vector_x->offset = %d\n", vsip_cvgetoffset_f(vector_x));
	printf("vector_x->stride = %d\n", vsip_cvgetstride_f(vector_x));

    // initialize_test_data(vector_y);

    //初始化FFT正变换参数
    fft_op_plan = vsip_ccfftop_create_f(N,1.0,-1,1,0);
    if(fft_op_plan==NULL)
    {
        printf("create fft plan failed\n");
    }
    
    //初始化FFT逆变换参数
    ifft_op_plan= vsip_ccfftop_create_f(N,1.0/N,1,1,0);
    if(ifft_op_plan==NULL)
    {
        printf("create ifft plan failed\n");
    }

}


void finalize(void)
{
    
    vsip_cvalldestroy_f(vector_x);
    vsip_cvalldestroy_f(vector_y);
    vsip_cvalldestroy_f(vector_z);

    int s = vsip_fft_destroy_f(fft_op_plan);
    if(s!=0)
    {
        printf("fft_destroy failed\n");
    }
    s = vsip_fft_destroy_f(ifft_op_plan);
    if(s!=0)
    {
        printf("ifft_destroy ok\n");
    }
    
    vsip_finalize((void *)0);

}


int main(int argc,char *argv[])
{
	//fft dimension
	int N = (vsip_length) atoi(argv[1]);

    printf("\n\n---------------------------------------start---------------------------------------------\n");
    printf("fft dimension = %d, repeat = %d \n\n", N, repeat);

//vsipl初始化，初始化向量数据数据，初始化FFT/IFFT所需参数
        initialize(N);
        //初始化进行fft运算的vector_x数据
        initialize_test_data(vector_x);
        printf("vector_x : \n");VU_cvprint_f(vector_x);printf("\n");
		
		printf("vector_x attribute:\n");
		printf("vector_x->length = %d\n", vsip_cvgetlength_f(vector_x));
		printf("vector_x->offset = %d\n", vsip_cvgetoffset_f(vector_x));
		printf("vector_x->stride = %d\n", vsip_cvgetstride_f(vector_x));

//FFT运算            
		//对vector_x作FFT正变换，结果保存在vector_y中
        // gettimeofday(&tv1,0);
        // for (int i = 0; i < repeat; i++)
        // {
            vsip_ccfftop_f( fft_op_plan,vector_x,vector_y);
        // }
        // gettimeofday(&tv2,0);
        // time_cfar_used=(double)(tv2.tv_sec-tv1.tv_sec)*1000000+(double)(tv2.tv_usec-tv1.tv_usec);
        printf("FFT(vector_x) -> vector_y : \n");VU_cvprint_f(vector_y);printf("\n");
        // printf("FFT time:%.2f us\n",time_cfar_used/repeat);


        //对vector_y作FFT逆变换，结果保存在vector_z中
        printf("*************\n");
        // gettimeofday(&tv1,0);
        // for (int i = 0; i < repeat; i++)
        // {
            vsip_ccfftop_f(ifft_op_plan,vector_y,vector_z);
        // }
        // gettimeofday(&tv2,0);
        // time_cfar_used=(double)(tv2.tv_sec-tv1.tv_sec)*1000000+(double)(tv2.tv_usec-tv1.tv_usec);
        printf("FFT(vector_y) -> vector_z : \n");VU_cvprint_f(vector_z);printf("\n");
        // printf("IFFT time:%.2f us\n",time_cfar_used/repeat);printf("\n");
        finalize();
        printf("\n----------------------------------------end---------------------------------------------\n\n");
        return 0;
}

