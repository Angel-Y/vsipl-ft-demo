#include<vsip.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define REPEAT_TIMES 500

struct timeval tv_b_cfar,tv_e_cfar;
double time_cfar_used;

int main(){vsip_init((void*)0);    //初始化vsipl库
{
   void VU_vprint_f(vsip_vview_f*);
   
   printf("repeat :%d\n",REPEAT_TIMES);
   /*********初始化向量，赋予随机数********/
   srand(time(NULL));
   for(int i = 128; i <= 128*1024; i *= 2)
   // for(int i = 1; i <= 10; i++)
   {
      //float* in = (float*)malloc(i*sizeof(float));
      //float* out = (float*)malloc(i*sizeof(float));
      //创建视图A、R
      vsip_vview_f *A = vsip_vcreate_f(i,0),
                   *R = vsip_vcreate_f(i,0);
      for(int j=0;j<i;j++){
         float num_random = ((float) (rand() % 1000+1)) / 10.0f;
         //填充数据块
         vsip_vput_f(A,j,num_random);
         //in[j] = num_random;
      }
      
      time_cfar_used = 0.0f;
      for(int j=0;j<REPEAT_TIMES;j++){
         gettimeofday(&tv_b_cfar,0);
         vsip_vsin_f(A,R);
         gettimeofday(&tv_e_cfar,0);
	      time_cfar_used += (double)(tv_e_cfar.tv_sec-tv_b_cfar.tv_sec)*1000000+(double)(tv_e_cfar.tv_usec-tv_b_cfar.tv_usec);
      }
      printf("%d vsip_vsin_f : %.3f us\n",i,time_cfar_used/REPEAT_TIMES);
      printf("-------------------\n");

      //销毁视图A、R
      vsip_valldestroy_f(A);
      vsip_valldestroy_f(R);
   }
} vsip_finalize((void*)0); return 1;}   //释放vsipl库使用资源

void VU_vprint_f(vsip_vview_f* a){
   vsip_length i, n;
   n = vsip_vgetlength_f(a);
   for(i=0; i<n; i++)
      printf("%5.3f\t",vsip_vget_f(a,i));
   printf("\n");
   return;
}
