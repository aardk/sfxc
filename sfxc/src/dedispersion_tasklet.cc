#include "dedispersion_tasklet.h"

Dedispersion_tasklet::Dedispersion_tasklet() {
}

Dedispersion_tasklet::~Dedispersion_tasklet(){
}

bool
Dedispersion_tasklet::do_task(){
  bool done_work = false;
  for (size_t i = 0; i < dedispersion_modules.size(); i++) {
    if (dedispersion_modules[i] != Coherent_dedispersion_ptr()) {
      if (dedispersion_modules[i]->has_work()) {
        dedispersion_modules[i]->do_task();
        done_work = true;
      }
    }
  }
  return done_work;
}

void
Dedispersion_tasklet::empty_output_queue(){
  for (size_t i = 0; i < dedispersion_modules.size(); i++) {
    if (dedispersion_modules[i] != Coherent_dedispersion_ptr()) {
      dedispersion_modules[i]->empty_output_queue();
    }
  }
}

void
Dedispersion_tasklet::create_dedispersion_filter(){
  filter.resize(fft_size_dedispersion + 1);
  double dnu = channel_bw / fft_size_dedispersion;
  double f0 = channel_freq + sideband * channel_bw / 2;
  double arg0 = 0;//channel_offset * channel_freq;
  for (int i = 0; i < fft_size_dedispersion + 1; i++) {
    double f = sideband*dnu*(i - fft_size_dedispersion / 2);
    double arg1 = sideband*f*f*DM / (2.41e-10*f0*f0*(f+f0)); // For f, f0 in [MHz] 
    double arg = arg0 + arg1;
    arg = -2 * M_PI * (arg - floor(arg)); 
    filter[i] = std::complex<FLOAT>(cos(arg) / fft_size_dedispersion, 
                                    sin(arg) / fft_size_dedispersion);
    std::cerr.precision(8);
    if(RANK_OF_NODE == -10)
      std::cerr << "filter["<<i<<"]="<<filter[i] 
                << ", f= " << f << ", f0=" << f0
                << ", dnu = "<< dnu<<", DM="<< DM<< "\n";
  }
}

Dedispersion_tasklet::Delay_queue_ptr
Dedispersion_tasklet::get_output_buffer(int stream_nr) {
  return dedispersion_modules[stream_nr]->get_output_buffer();
}

void Dedispersion_tasklet::connect_to(Delay_queue_ptr buffer, int stream_nr) {
  if (dedispersion_modules.size() <= stream_nr)
    dedispersion_modules.resize(stream_nr+1);
  dedispersion_modules[stream_nr] = Coherent_dedispersion_ptr(new Coherent_dedispersion(stream_nr, 
                                                  fft, fft_cor, filter, dedispersion_buffer, zeropad_buffer));
  dedispersion_modules[stream_nr]->connect_to(buffer);
}

void 
Dedispersion_tasklet::set_parameters(const Correlation_parameters &parameters, Pulsar &pulsar){
  total_input_fft = 0;

  DM = pulsar.polyco_params[0].DM;
  channel_freq = parameters.channel_freq / 1000000.; // In MHz
  channel_bw = parameters.bandwidth / 1000000.; // In MHz
  sideband = (parameters.sideband == 'L') ? -1 : 1;
  fft_size_dedispersion = parameters.fft_size_dedispersion;
  fft_size_correlation = parameters.fft_size_correlation;
  start_time = parameters.slice_start;
  stop_time = parameters.slice_start + parameters.slice_time; 

  dedispersion_buffer.resize(fft_size_dedispersion + 1);
  zeropad_buffer.resize(2 * fft_size_correlation);
  // zero padding
  SFXC_ZERO_F(&zeropad_buffer[fft_size_correlation], fft_size_correlation);
  // Initialize the FFT's
  fft.resize(2 * fft_size_dedispersion);
  fft_cor.resize(2 * fft_size_correlation);
  create_dedispersion_filter();

  // Apply parameters to modules
  for (size_t i = 0; i < dedispersion_modules.size(); i++) {
    if (dedispersion_modules[i] != Coherent_dedispersion_ptr()) {
      // We share the dedispersion filter and some buffers between dedispersion_modules
      // to reduce memory usage, which can become enormous at P band frequencies. 
      dedispersion_modules[i]->set_parameters(parameters);
    }
  }
}
