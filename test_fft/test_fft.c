//一维复数out-of-place FFT: vsip_ccfftop_f()
//测试程序

#include "vsip.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>


vsip_cvview_f *vector_x;
vsip_cvview_f *vector_y;


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


vsip_fft_f      *fft_op_plan;

int repeat = 1;

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


    //初始化FFT正变换参数
    fft_op_plan = vsip_ccfftop_create_f(N,1.0,-1,1,0);
    if(fft_op_plan==NULL)
    {
        printf("create fft plan failed\n");
    }

}


void finalize(void)
{
    
    vsip_cvalldestroy_f(vector_x);
    vsip_cvalldestroy_f(vector_y);

    int s = vsip_fft_destroy_f(fft_op_plan);
    if(s!=0)
    {
        printf("fft_destroy failed\n");
    }
    
    vsip_finalize((void *)0);

}


int main(int argc,char *argv[])
{
	int N = (vsip_length) atoi(argv[1]);
    printf("\n\n---------------------------------------start---------------------------------------------\n");
    printf("fft dimension = %d, repeat = %d \n\n", N, repeat);
	initialize(N);
	//initialize_test_data(vector_x);
	vsip_ccfftop_f( fft_op_plan,vector_x,vector_y);        
    finalize();
    printf("\n----------------------------------------end---------------------------------------------\n\n");
    return 0;
}


