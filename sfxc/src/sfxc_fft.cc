#include "sfxc_fft.h"
#include "utils.h"

#ifdef USE_IPP
#include <ippcore.h>

sfxc_fft_ipp::sfxc_fft_ipp():buffer(NULL), buffer_r2c(NULL), 
                             ippspec(NULL), ippspec_r2c(NULL),
                             specbuffer(NULL), specbuffer_r2c(NULL){
  size = 0;
  order = 0;
}

sfxc_fft_ipp::~sfxc_fft_ipp(){
  free_buffers();
}

void
sfxc_fft_ipp::free_buffers(){
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
}

void
sfxc_fft_ipp::resize(int size_){
  if(size == size_)
    return;
  free_buffers();
  size = size_;
  order=-1;
  for(int n = size; n > 0; order++)
    n/=2;
}

void
sfxc_fft_ipp::alloc(){
  IppStatus status;
  int specbufsize, initbufsize, extbufsize;
  status = ippsFFTGetSize_C_64fc(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                                 &specbufsize, &initbufsize, &extbufsize);
  // Create buffers
  Ipp8u *initbuffer = NULL;
  if (initbufsize > 0)
    initbuffer = (Ipp8u*) ippMalloc(initbufsize);
  if (extbufsize > 0)
    buffer = (Ipp8u*) ippMalloc(extbufsize);
  specbuffer = (Ipp8u*) ippMalloc(specbufsize);
  // Initialize FFT spec
  status = ippsFFTInit_C_64fc(&ippspec, order, IPP_FFT_NODIV_BY_ANY,
                              ippAlgHintFast, specbuffer, initbuffer);
  if(status != ippStsNoErr)
    sfxc_abort("Unable to allocate IPP fft");
  if (initbufsize > 0)
    ippFree(initbuffer);
}

void
sfxc_fft_ipp::alloc_r2c(){
  IppStatus status;
  int specbufsize, initbufsize, extbufsize;
  status = ippsFFTGetSize_R_64f(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                                 &specbufsize, &initbufsize, &extbufsize);
  // Create buffers
  Ipp8u *initbuffer = NULL;
  if (initbufsize > 0)
    initbuffer = (Ipp8u*) ippMalloc(initbufsize);
  if (extbufsize > 0)
    buffer_r2c = (Ipp8u*) ippMalloc(extbufsize);
  specbuffer_r2c = (Ipp8u*) ippMalloc(specbufsize);
  // Initialize FFT spec
  status = ippsFFTInit_R_64f(&ippspec_r2c, order, IPP_FFT_NODIV_BY_ANY,
                              ippAlgHintFast, specbuffer_r2c, initbuffer);
  if (status != ippStsNoErr)
    sfxc_abort("Unable to allocate IPP rfft");
  if (initbufsize > 0)
    ippFree(initbuffer);
}

void
sfxc_fft_ipp::fft(const std::complex<double> *in, std::complex<double> *out){
  if(ippspec == NULL)
    alloc();
  if(in == out)
    ippsFFTFwd_CToC_64fc_I((Ipp64fc *) out, ippspec, buffer);
  else
    ippsFFTFwd_CToC_64fc((Ipp64fc *)in, (Ipp64fc *)out, ippspec, buffer);
}  

void
sfxc_fft_ipp::rfft(const double *in, std::complex<double> *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(ippspec_r2c == NULL)
    alloc_r2c();
  ippsFFTFwd_RToCCS_64f((Ipp64f *)in, (Ipp64f *)out, ippspec_r2c, buffer_r2c);
} 

void
sfxc_fft_ipp::ifft(const std::complex<double> *in, std::complex<double> *out){
  if(ippspec == NULL)
    alloc();
  if(in == out)
    ippsFFTInv_CToC_64fc_I((Ipp64fc *)out, ippspec, buffer);
  else
    ippsFFTInv_CToC_64fc((Ipp64fc *)in, (Ipp64fc *)out, ippspec, buffer);
}

void
sfxc_fft_ipp::irfft(const std::complex<double> *in, double *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(ippspec_r2c == NULL)   
    alloc_r2c();
  ippsFFTInv_CCSToR_64f((Ipp64f *)in, (Ipp64f *)out, ippspec_r2c, buffer_r2c);
}
#else // USE FFTW
sfxc_fft_fftw::sfxc_fft_fftw(){
  plan_forward_set = false;
  plan_backward_set = false;
  plan_forward_I_set = false;
  plan_backward_I_set = false;
  plan_forward_r2c_set = false;
  plan_backward_r2c_set = false;
}

sfxc_fft_fftw::~sfxc_fft_fftw(){
  free_buffers();
}

void
sfxc_fft_fftw::free_buffers(){
  if(plan_forward_set){
    fftw_destroy_plan(plan_forward);
    plan_forward_set = false;
  }
  if(plan_backward_set){
    fftw_destroy_plan(plan_backward);
    plan_backward_set = false;
  }
  if(plan_forward_I_set){
    fftw_destroy_plan(plan_forward_I);
    plan_forward_I_set = false;
  }
  if(plan_backward_I_set){
    fftw_destroy_plan(plan_backward_I);
    plan_backward_I_set = false;
  }
  if(plan_forward_r2c_set){
    fftw_destroy_plan(plan_forward_r2c);
    plan_forward_r2c_set = false;
  }
  if(plan_backward_r2c_set){
    fftw_destroy_plan(plan_backward_r2c);
    plan_backward_r2c_set = false;
  }
}

void
sfxc_fft_fftw::resize(int size_){
  free_buffers();
  size = size_;
}

fftw_plan
sfxc_fft_fftw::alloc(int sign, bool inplace){
  fftw_complex *temp_in = (fftw_complex *) fftw_malloc(size * sizeof(fftw_complex));
  fftw_complex *temp_out;
  if(inplace)
    temp_out = temp_in;
  else
    temp_out = (fftw_complex *) fftw_malloc(size * sizeof(fftw_complex));
  if((temp_in == NULL) || (temp_out == NULL))
    sfxc_abort("Unable to allocate buffer for fft\n");
//  fftw_plan plan = fftw_plan_dft_1d(size, temp_in, temp_out, sign, FFTW_MEASURE);
  fftw_plan plan = fftw_plan_dft_1d(size, temp_in, temp_out, sign, FFTW_ESTIMATE);
  fftw_free(temp_in);
  if(!inplace)
    fftw_free(temp_out);
  return plan;
}

fftw_plan
sfxc_fft_fftw::
alloc_r2c(int sign){
  double *temp_real = (double *) fftw_malloc(size * sizeof(double));
  fftw_complex *temp_complex = (fftw_complex *)fftw_malloc(size * sizeof(fftw_complex));
  if((temp_real == NULL) || (temp_complex == NULL))
    sfxc_abort("Unable to allocate buffer for fft\n");
  fftw_plan plan;
//  if(sign == FFTW_FORWARD)
//    plan = fftw_plan_dft_r2c_1d(size, temp_real, temp_complex, FFTW_MEASURE);
//  else
//    plan = fftw_plan_dft_c2r_1d(size, temp_complex, temp_real, FFTW_MEASURE);
  if(sign == FFTW_FORWARD)
    plan = fftw_plan_dft_r2c_1d(size, temp_real, temp_complex, FFTW_ESTIMATE);
  else
    plan = fftw_plan_dft_c2r_1d(size, temp_complex, temp_real, FFTW_ESTIMATE);

  fftw_free(temp_real);
  fftw_free(temp_complex);
  return plan;
}

void
sfxc_fft_fftw::fft(const std::complex<double> *in, std::complex<double> *out){
  bool inplace = (in == out);
  if((inplace) && (!plan_forward_I_set)){
    plan_forward_I = alloc(FFTW_FORWARD, inplace);
    plan_forward_I_set = true;
  }else if(!plan_forward_set){
    plan_forward = alloc(FFTW_FORWARD, inplace);
    plan_forward_set = true;
  }

  if(inplace)
    fftw_execute_dft(plan_forward_I, (fftw_complex *)in, (fftw_complex *)out);
  else
    fftw_execute_dft(plan_forward, (fftw_complex *)in, (fftw_complex *)out);
}  

void
sfxc_fft_fftw::rfft(const double *in, std::complex<double> *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(!plan_forward_r2c_set){
    plan_forward_r2c = alloc_r2c(FFTW_FORWARD);
    plan_forward_r2c_set = true;
  }
  fftw_execute_dft_r2c(plan_forward_r2c, (double *)in, (fftw_complex *)out);
} 

void
sfxc_fft_fftw::ifft(const std::complex<double> *in, std::complex<double> *out){
  bool inplace = (in == out);
  if((inplace) && (!plan_backward_I_set)){
    plan_backward_I = alloc(FFTW_BACKWARD, inplace);
    plan_backward_I_set = true;
  }else if(!plan_backward_set){
    plan_backward = alloc(FFTW_BACKWARD, inplace);
    plan_backward_set = true;
  }

  if(inplace)
    fftw_execute_dft(plan_backward_I, (fftw_complex *)in, (fftw_complex *)out);
  else
    fftw_execute_dft(plan_backward, (fftw_complex *)in, (fftw_complex *)out);
}

void
sfxc_fft_fftw::irfft(const std::complex<double> *in, double *out){
  SFXC_ASSERT((void *)in != (void *)out);
  if(!plan_backward_r2c_set){
    plan_backward_r2c = alloc_r2c(FFTW_BACKWARD);
    plan_backward_r2c_set = true;
  }
  fftw_execute_dft_c2r(plan_backward_r2c, (fftw_complex *)in, (double *)out);
}

#endif // USE_IPP
