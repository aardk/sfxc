#include "sfxc_fft_float.h"
#include "utils.h"

#ifdef USE_IPP
#include <ippcore.h>

sfxc_fft_ipp_float::sfxc_fft_ipp_float():buffer(NULL), buffer_r2c(NULL), 
                                         ippspec(NULL), ippspec_r2c(NULL),
                                         specbuffer(NULL), specbuffer_r2c(NULL){
  size = 0;
  order = 0;
}

sfxc_fft_ipp_float::~sfxc_fft_ipp_float(){
  free_buffers();
}

void
sfxc_fft_ipp_float::free_buffers(){
  if(buffer != NULL){
    ippFree(buffer);
    buffer = NULL;
  }
  if(buffer_r2c != NULL){
    ippFree(buffer_r2c);
    buffer_r2c = NULL;
  }
  if(specbuffer != NULL){
    ippFree(specbuffer);
    specbuffer = NULL;
  }
  if(specbuffer_r2c != NULL){
    ippFree(specbuffer_r2c);
    specbuffer_r2c = NULL;
  }
  if(ippspec != NULL){
    ippFree(ippspec);
    ippspec = NULL;
  }
  if(ippspec_r2c != NULL){
    ippFree(ippspec_r2c);
    ippspec_r2c = NULL;
  }
}

void
sfxc_fft_ipp_float::resize(int size_){
  if(size_ == size)
    return;
  free_buffers();
  size = size_;
  order=-1;
  for(int n = size; n > 0; order++)
    n/=2;
}

void
sfxc_fft_ipp_float::alloc(){
  IppStatus status;
  int specbufsize, initbufsize, extbufsize;
  status = ippsFFTGetSize_C_32fc(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                                 &specbufsize, &initbufsize, &extbufsize);
  buffer = (Ipp8u*) ippMalloc(extbufsize);
  Ipp8u *initbuffer = (Ipp8u*) ippMalloc(initbufsize);
  specbuffer = (Ipp8u*) ippMalloc(specbufsize);
  status = ippsFFTInit_C_32fc(&ippspec, order, IPP_FFT_NODIV_BY_ANY, 
                              ippAlgHintFast, specbuffer, initbuffer);
  if(status != ippStsNoErr)
    sfxc_abort("Unable to allocate IPP fft");
  ippFree(initbuffer);
}

void
sfxc_fft_ipp_float::alloc_r2c(){
  IppStatus status;
  int specbufsize, initbufsize, extbufsize;
  status = ippsFFTGetSize_R_32f(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                                 &specbufsize, &initbufsize, &extbufsize);
  buffer_r2c = (Ipp8u*) ippMalloc(extbufsize);
  Ipp8u *initbuffer = (Ipp8u*) ippMalloc(initbufsize);
  specbuffer_r2c = (Ipp8u*) ippMalloc(specbufsize);
  status = ippsFFTInit_R_32f(&ippspec_r2c, order, IPP_FFT_NODIV_BY_ANY, 
                              ippAlgHintFast, specbuffer_r2c, initbuffer);
  if(status != ippStsNoErr)
    sfxc_abort("Unable to allocate IPP rfft");
  ippFree(initbuffer);
}

void
sfxc_fft_ipp_float::fft(const std::complex<float> *in, std::complex<float> *out){
  if(ippspec == NULL)
    alloc();
  if(in == out)
    ippsFFTFwd_CToC_32fc_I((Ipp32fc *) out, ippspec, buffer);
  else
    ippsFFTFwd_CToC_32fc((Ipp32fc *)in, (Ipp32fc *)out, ippspec, buffer);
}  

void
sfxc_fft_ipp_float::rfft(const float *in, std::complex<float> *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(ippspec_r2c == NULL)
    alloc_r2c();
  ippsFFTFwd_RToCCS_32f((Ipp32f *)in, (Ipp32f *)out, ippspec_r2c, buffer_r2c);
} 

void
sfxc_fft_ipp_float::ifft(const std::complex<float> *in, std::complex<float> *out){
  if(ippspec == NULL)
    alloc();
  if(in == out)
    ippsFFTInv_CToC_32fc_I((Ipp32fc *)out, ippspec, buffer);
  else
    ippsFFTInv_CToC_32fc((Ipp32fc *)in, (Ipp32fc *)out, ippspec, buffer);
}

void
sfxc_fft_ipp_float::irfft(const std::complex<float> *in, float *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(ippspec_r2c == NULL)   
    alloc_r2c();
  ippsFFTInv_CCSToR_32f((Ipp32f *)in, (Ipp32f *)out, ippspec_r2c, buffer_r2c);
}
#else // USE FFTW
sfxc_fft_fftw_float::sfxc_fft_fftw_float(){
  plan_forward_set = false;
  plan_backward_set = false;
  plan_forward_I_set = false;
  plan_backward_I_set = false;
  plan_forward_r2c_set = false;
  plan_backward_r2c_set = false;
}

sfxc_fft_fftw_float::~sfxc_fft_fftw_float(){
  free_buffers();
}

void
sfxc_fft_fftw_float::free_buffers(){
  if(plan_forward_set){
    fftwf_destroy_plan(plan_forward);
    plan_forward_set = false;
  }
  if(plan_backward_set){
    fftwf_destroy_plan(plan_backward);
    plan_backward_set = false;
  }
  if(plan_forward_I_set){
    fftwf_destroy_plan(plan_forward_I);
    plan_forward_I_set = false;
  }
  if(plan_backward_I_set){
    fftwf_destroy_plan(plan_backward_I);
    plan_backward_I_set = false;
  }
  if(plan_forward_r2c_set){
    fftwf_destroy_plan(plan_forward_r2c);
    plan_forward_r2c_set = false;
  }
  if(plan_backward_r2c_set){
    fftwf_destroy_plan(plan_backward_r2c);
    plan_backward_r2c_set = false;
  }
}

void
sfxc_fft_fftw_float::resize(int size_){
  free_buffers();
  size = size_;
}

fftwf_plan
sfxc_fft_fftw_float::alloc(int sign, bool inplace){
  fftwf_complex *temp_in = (fftwf_complex *) fftwf_malloc(size * sizeof(fftwf_complex));
  fftwf_complex *temp_out;
  if(inplace)
    temp_out = temp_in;
  else
    temp_out = (fftwf_complex *) fftwf_malloc(size * sizeof(fftwf_complex));
  if((temp_in == NULL) || (temp_out == NULL))
    sfxc_abort("Unable to allocate buffer for fft\n");
//  fftwf_plan plan = fftwf_plan_dft_1d(size, temp_in, temp_out, sign, FFTW_MEASURE);
  fftwf_plan plan = fftwf_plan_dft_1d(size, temp_in, temp_out, sign, FFTW_ESTIMATE);
  fftwf_free(temp_in);
  if(!inplace)
    fftwf_free(temp_out);
  return plan;
}

fftwf_plan
sfxc_fft_fftw_float::
alloc_r2c(int sign){
  float *temp_real = (float *) fftwf_malloc(size * sizeof(float));
  fftwf_complex *temp_complex = (fftwf_complex *)fftwf_malloc(size * sizeof(fftwf_complex));
  if((temp_real == NULL) || (temp_complex == NULL))
    sfxc_abort("Unable to allocate buffer for fft\n");
  fftwf_plan plan;
//  if(sign == FFTW_FORWARD)
//    plan = fftwf_plan_dft_r2c_1d(size, temp_real, temp_complex, FFTW_MEASURE);
//  else
//    plan = fftwf_plan_dft_c2r_1d(size, temp_complex, temp_real, FFTW_MEASURE);
  if(sign == FFTW_FORWARD)
    plan = fftwf_plan_dft_r2c_1d(size, temp_real, temp_complex, FFTW_ESTIMATE);
  else
    plan = fftwf_plan_dft_c2r_1d(size, temp_complex, temp_real, FFTW_ESTIMATE);

  fftwf_free(temp_real);
  fftwf_free(temp_complex);
  return plan;
}

void
sfxc_fft_fftw_float::fft(const std::complex<float> *in, std::complex<float> *out){
  bool inplace = (in == out);
  if((inplace) && (!plan_forward_I_set)){
    plan_forward_I = alloc(FFTW_FORWARD, inplace);
    plan_forward_I_set = true;
  }else if(!plan_forward_set){
    plan_forward = alloc(FFTW_FORWARD, inplace);
    plan_forward_set = true;
  }

  if(inplace)
    fftwf_execute_dft(plan_forward_I, (fftwf_complex *)in, (fftwf_complex *)out);
  else
    fftwf_execute_dft(plan_forward, (fftwf_complex *)in, (fftwf_complex *)out);
}  

void
sfxc_fft_fftw_float::rfft(const float *in, std::complex<float> *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(!plan_forward_r2c_set){
    plan_forward_r2c = alloc_r2c(FFTW_FORWARD);
    plan_forward_r2c_set = true;
  }
  fftwf_execute_dft_r2c(plan_forward_r2c, (float *)in, (fftwf_complex *)out);
} 

void
sfxc_fft_fftw_float::ifft(const std::complex<float> *in, std::complex<float> *out){
  bool inplace = (in == out);
  if((inplace) && (!plan_backward_I_set)){
    plan_backward_I = alloc(FFTW_BACKWARD, inplace);
    plan_backward_I_set = true;
  }else if(!plan_backward_set){
    plan_backward = alloc(FFTW_BACKWARD, inplace);
    plan_backward_set = true;
  }

  if(inplace)
    fftwf_execute_dft(plan_backward_I, (fftwf_complex *)in, (fftwf_complex *)out);
  else
    fftwf_execute_dft(plan_backward, (fftwf_complex *)in, (fftwf_complex *)out);
}

void
sfxc_fft_fftw_float::irfft(const std::complex<float> *in, float *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(!plan_backward_r2c_set){
    plan_backward_r2c = alloc_r2c(FFTW_BACKWARD);
    plan_backward_r2c_set = true;
  }
  fftwf_execute_dft_c2r(plan_backward_r2c, (fftwf_complex *)in, (float *)out);
}
#endif // USE_IPP
