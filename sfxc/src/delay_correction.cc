#include "delay_correction.h"
#include "sfxc_math.h"
#include "config.h"

Delay_correction::Delay_correction(int stream_nr_)
    : output_buffer(Output_buffer_ptr(new Output_buffer())),
      output_memory_pool(32),current_time(-1),
      stream_nr(stream_nr_), stream_idx(-1)
{
}

Delay_correction::~Delay_correction() {
#if PRINT_TIMER
  int N = number_channels();
  int numiterations = total_ffts;
  double time = delay_timer.measured_time()*1000000;
  PROGRESS_MSG("MFlops: " << 5.0*N*log2(N) * numiterations / (1.0*time));
#endif
}

void Delay_correction::do_task() {
  SFXC_ASSERT(has_work());
  Input_buffer_element input = input_buffer->front_and_pop();
  int nbuffer=input->nfft;
  current_fft+=nbuffer;
  // Allocate output buffer
  int output_stride =  fft_cor_size()/2 + 4; // there are fft_size+1 points and each fft should be 16 bytes alligned
  Output_buffer_element cur_output = output_memory_pool.allocate();
  cur_output->stride = output_stride;
  int window_func = correlation_parameters.window;
  int nfft_cor;
  if (window_func == SFXC_WINDOW_PFB) {
    int nb = (nbuffer * fft_size() + tbuf_end - tbuf_start ) / (fft_rot_size() / 2);
    nfft_cor = std::max(0, nb - 2*SFXC_NTAPS + 1);
  } else
    nfft_cor = (nbuffer * fft_size() + tbuf_end - tbuf_start) / (fft_rot_size() / 2);
  // The windowing touches each block twice, so the last block needs to be preserved
  if ((window_func != SFXC_WINDOW_NONE) && (window_func != SFXC_WINDOW_PFB))
    nfft_cor -= 1;
  if (cur_output->data.size() != nfft_cor * output_stride)
    cur_output->data.resize(nfft_cor * output_stride);
#ifndef DUMMY_CORRELATION
  size_t tbuf_size = time_buffer.size();
  for(int buf=0;buf<nbuffer;buf++) {
    double delay = get_delay(current_time + fft_length/2);
    double delay_in_samples = delay*sample_rate();
    int integer_delay = (int)std::floor(delay_in_samples+.5);

    // Output is in frequency_buffer
    fractional_bit_shift(&input->data[buf * fft_size()],
                         integer_delay,
                         delay_in_samples - integer_delay);

    // Input is from frequency_buffer
    fringe_stopping(&time_buffer[tbuf_end%tbuf_size]);
    tbuf_end += fft_size();

    current_time.inc_samples(fft_size());
    total_ffts++;
  }
  SFXC_ASSERT(tbuf_end - tbuf_start <= tbuf_size);

  size_t nsamp_per_window;
  if (window_func == SFXC_WINDOW_NONE) 
    nsamp_per_window = fft_rot_size() / 2;
  else if (window_func == SFXC_WINDOW_PFB) 
    nsamp_per_window = fft_rot_size() * SFXC_NTAPS;
  else
    nsamp_per_window = fft_rot_size();

  for(int i=0; i<nfft_cor; i++){
    // apply window function
    size_t eob = tbuf_size - tbuf_start%tbuf_size; // how many samples to end of buffer
    size_t nsamp = std::min(eob, nsamp_per_window);
    SFXC_MUL_F(&time_buffer[tbuf_start%tbuf_size], &window[0], &temp_buffer[0], nsamp);
    if(nsamp < nsamp_per_window)
      SFXC_MUL_F(&time_buffer[0], &window[nsamp], &temp_buffer[nsamp], nsamp_per_window - nsamp);
    // Flip sideband if needed
    if (correlation_parameters.sideband != correlation_parameters.station_streams[stream_idx].sideband)
      SFXC_MUL_F(&temp_buffer[0], &flip[0], &temp_buffer[0], nsamp_per_window);
    // When SFXC_WINDOW_NONE is set we zeropad
    if (window_func == SFXC_WINDOW_NONE)
      memset(&temp_buffer[nsamp_per_window], 0, nsamp_per_window*sizeof(FLOAT));
    else if (window_func == SFXC_WINDOW_PFB) {
      for (int j=1; j<SFXC_NTAPS; j++) {
        SFXC_ADD_F_I(&temp_buffer[j * fft_rot_size()], &temp_buffer[0], fft_rot_size());
      }
    }

    tbuf_start += fft_rot_size()/2;
    SFXC_ASSERT(tbuf_start <= tbuf_end);
    // Do the final fft from time to frequency
    fft_t2f_cor.rfft(&temp_buffer[0], &temp_fft_buffer[temp_fft_offset]);
    memcpy(&cur_output->data[i * output_stride], &temp_fft_buffer[output_offset], output_stride * sizeof(std::complex<FLOAT>));
  }
#endif // DUMMY_CORRELATION
  if(nfft_cor > 0){
    output_buffer->push(cur_output);
  }
}

void Delay_correction::fractional_bit_shift(FLOAT *input,
    int integer_shift,
    double fractional_delay) {
  // 3) execute the complex to complex FFT, from Time to Frequency domain
  //    input: sls. output sls_freq
  fft_t2f.rfft(&input[0], &frequency_buffer[0]);
  total_ffts++;

  // Element 0 and (fft_size() / 2) are real numbers
  frequency_buffer[0] *= 0.5;
  frequency_buffer[fft_size() / 2] *= 0.5; // Nyquist frequency

  // 4c) zero the unused subband (?)
  SFXC_ZERO_FC(&frequency_buffer[(fft_size() / 2) + 1], (fft_size() / 2) - 1);

  // 5a)calculate the fract bit shift (=phase corrections in freq domain)
  // the following should be double
  const double dfr  = (double)sample_rate() / fft_size(); // delta frequency
  const double tmp1 = -2.0*M_PI*fractional_delay/sample_rate();
  const double tmp2 = M_PI*(integer_shift & (4*oversamp - 1))/(2*oversamp);
  const double constant_term = tmp2 -tmp1 * (bandwidth() / 2);
  const double linear_term = tmp1*dfr;

  // 5b)apply phase correction in frequency range
  const int size = (fft_size() / 2) + 1;

  double phi = constant_term;
  // in the loop we calculate sin(phi) and cos(phi) with phi=contant_term + i*linear_term
  // This we do efficiently using a recurrence relation
  // sin(t+delta)=sin(t)-[a*sin(t)-b*cos(t)] ; cos(t+delta)=cos(t)-[a*cos(t)+b*sin(t)]
  // a=2*sin^2(delta/2) ; b=sin(delta)
  double temp=sin(linear_term/2);
  double a=2*temp*temp,b=sin(linear_term);
  double cos_phi, sin_phi;
#ifdef HAVE_SINCOS
  sincos(phi, &sin_phi, &cos_phi);
#else
  sin_phi = sin(phi);
  cos_phi = cos(phi);
#endif
  for (int i = 0; i < size; i++) {
    // the following should be double
    exp_array[i] = std::complex<FLOAT>(cos_phi,-sin_phi);
    // Compute sin_phi=sin(phi); cos_phi = cos(phi);
    temp=sin_phi-(a*sin_phi-b*cos_phi);
    cos_phi=cos_phi-(a*cos_phi+b*sin_phi);
    sin_phi=temp;
  }
  SFXC_MUL_FC_I(&exp_array[0], &frequency_buffer[0], size);

  // 6a)execute the complex to complex FFT, from Frequency to Time domain
  //    input: sls_freq. output sls
  fft_f2t.ifft(&frequency_buffer[0], &frequency_buffer[0]);

  total_ffts++;
}

void Delay_correction::fringe_stopping(FLOAT output[]) {
  const double mult_factor_phi = -sideband() * 2.0 * M_PI;
  const double center_freq = channel_freq() + sideband() * (bandwidth() / 2) + LO_offset;

  double phi, delta_phi, sin_phi, cos_phi;
  double lo_phase = start_phase + LO_offset*current_time.diff(correlation_parameters.stream_start);
  phi = center_freq * get_delay(current_time) + lo_phase + get_phase(current_time) / (2 * M_PI);
  double floor_phi = std::floor(phi);
  phi = mult_factor_phi*(phi-floor_phi);

  { // compute delta_phi
    SFXC_ASSERT(((int64_t)fft_size() * 1000000) % sample_rate() == 0);
    double phi_end = center_freq * get_delay(current_time + fft_length) + 
                     lo_phase + fft_length.get_time()*LO_offset +
                     get_phase(current_time + fft_length) / (2 * M_PI);
    phi_end = mult_factor_phi*(phi_end-floor_phi);

    delta_phi = (phi_end - phi) / fft_size();
  }

  // We use a constant amplitude factor over the fft
  double amplitude = get_amplitude(current_time + fft_length/2);
  // We perform a recursion for the (co)sines similar to what is done in the fractional bitshift
  double temp=sin(delta_phi/2);
  double a=2*temp*temp,b=sin(delta_phi);
#ifdef HAVE_SINCOS
  sincos(phi, &sin_phi, &cos_phi);
  sin_phi *= amplitude;
  cos_phi *= amplitude;
#else
  sin_phi = amplitude * sin(phi);
  cos_phi = amplitude * cos(phi);
#endif

  for (size_t i = 0; i < fft_size(); i++) {
    // Compute sin_phi=sin(phi); cos_phi = cos(phi);
    // 7)subtract dopplers and put real part in Bufs for the current segment
    output[i] =
      frequency_buffer[i].real()*cos_phi + frequency_buffer[i].imag()*sin_phi;

    // Compute sin_phi=sin(phi); cos_phi = cos(phi);
    temp=sin_phi-(a*sin_phi-b*cos_phi);
    cos_phi=cos_phi-(a*cos_phi+b*sin_phi);
    sin_phi=temp;
  }
}

void
Delay_correction::set_parameters(const Correlation_parameters &parameters, Delay_table_akima &delays) {
  stream_idx = 0;
  while ((stream_idx < parameters.station_streams.size()) &&
         (parameters.station_streams[stream_idx].station_stream != stream_nr))
    stream_idx++;
  if (stream_idx == parameters.station_streams.size()) {
    // Data stream is not participating in current time slice
    return;
  }

  delay_table = delays;
  bits_per_sample = parameters.station_streams[stream_idx].bits_per_sample;
  correlation_parameters = parameters;
  oversamp = sample_rate() / (2 * bandwidth());

  current_time = parameters.stream_start;
  current_time.set_sample_rate(sample_rate());
  SFXC_ASSERT(((int64_t)fft_size() * 1000000) % sample_rate() == 0);
  fft_length = Time((double)fft_size() / (sample_rate() / 1000000));

  // We use a ringbuffer to implement the window overlap and PFB.
  // That buffer needs to be large enough to prevent samples that we
  // still need to be overwritten by the buffers we receive from the
  // bit2float module.  The following calculation over-estimates the
  // size a little bit, but that is ok.
  size_t nfft_min = std::max(fft_rot_size() / fft_size(), (size_t)1);
  size_t nfft_max = nfft_min * SFXC_NTAPS +
    (std::max(CORRELATOR_BUFFER_SIZE / fft_size(), nfft_min) *
     sample_rate()) / correlation_parameters.sample_rate;
  time_buffer.resize(nfft_max * fft_size());

  exp_array.resize(fft_size());
  frequency_buffer.resize(fft_size());
  if (parameters.window == SFXC_WINDOW_PFB)
    temp_buffer.resize(fft_rot_size() * SFXC_NTAPS);
  else
    temp_buffer.resize(fft_rot_size());

  if (fft_cor_size() > fft_rot_size())
    temp_fft_buffer.resize(fft_cor_size()/2 + 4);
  else
    temp_fft_buffer.resize(fft_rot_size()/2 + 4);
  memset(&temp_fft_buffer[0], 0, temp_fft_buffer.size() * sizeof(temp_fft_buffer[0]));

  fft_t2f.resize(fft_size());
  fft_f2t.resize(fft_size());
  fft_t2f_cor.resize(fft_rot_size());
  create_window();
  create_flip();

  // Calculate the offset into temp_fft_buffer where we can find the
  // spectral points that we want to correlate.
  int64_t delta, freq = channel_freq();
  if (parameters.sideband != parameters.station_streams[stream_idx].sideband)
    freq += sideband() * bandwidth();
  if (parameters.sideband == 'L')
    delta = freq - parameters.channel_freq;
  else
    delta = parameters.channel_freq - freq;
  if (delta > 0){
    output_offset = delta * fft_cor_size() / parameters.sample_rate;
    temp_fft_offset = 0;
  } else {
    output_offset = 0;
    temp_fft_offset = -delta * fft_cor_size() / parameters.sample_rate;
  }

  SFXC_ASSERT(parameters.fft_size_correlation >= parameters.fft_size_delaycor);
  n_ffts_per_integration =
    (int64_t) parameters.station_streams[stream_idx].sample_rate * parameters.slice_size /
     ((int64_t) parameters.sample_rate * parameters.fft_size_delaycor);

  LO_offset = parameters.station_streams[stream_idx].LO_offset;
  double dt = current_time.diff(parameters.experiment_start);
  // Compute start phase of LO offset with maximum numerical precision
  start_phase = LO_offset*(dt-trunc(dt)) + (LO_offset-trunc(LO_offset))*trunc(dt);
  start_phase = start_phase - floor(start_phase);

  extra_delay = parameters.station_streams[stream_idx].extra_delay;

  current_fft = 0;
  tbuf_start = 0;
  tbuf_end = 0;
}

void Delay_correction::connect_to(Input_buffer_ptr new_input_buffer) {
  SFXC_ASSERT(input_buffer == Input_buffer_ptr());
  input_buffer = new_input_buffer;
}

double Delay_correction::get_delay(Time time) {
  return delay_table.delay(time) + extra_delay;
}

double Delay_correction::get_phase(Time time) {
  return delay_table.phase(time);
}

double Delay_correction::get_amplitude(Time time) {
  return delay_table.amplitude(time);
}

bool Delay_correction::has_work() {
  if (input_buffer->empty())
    return false;
  if (output_memory_pool.empty())
    return false;
  if (n_ffts_per_integration == current_fft)
    return false;
  SFXC_ASSERT((current_fft<=n_ffts_per_integration)&&(current_fft>=0))
  return true;
}

Delay_correction::Output_buffer_ptr
Delay_correction::get_output_buffer() {
  SFXC_ASSERT(output_buffer != Output_buffer_ptr());
  return output_buffer;
}

int Delay_correction::sideband() {
  return (correlation_parameters.station_streams[stream_idx].sideband == 'L' ? -1 : 1);
}

void 
Delay_correction::create_window(){
  const int n = (correlation_parameters.window == SFXC_WINDOW_PFB) ? 
                SFXC_NTAPS * fft_rot_size() : fft_rot_size();
  window.resize(n);
  switch(correlation_parameters.window){
  case SFXC_WINDOW_NONE:
    //  Identical to the case without windowing
    for (int i=0; i<n/2; i++)
      window[i] = 1;
    for (int i = n/2; i < n; i++)
      window[i] = 0;
    break;
  case SFXC_WINDOW_RECT:
    // rectangular window (including zero padding)
    for (int i=0; i<n/4; i++)
      window[i] = 0;
    for (int i = n/4; i < 3*n/4; i++)
      window[i] = 1;
    for (int i = 3*n/4 ; i < n ; i++)
      window[i] = 0;
    break;
  case SFXC_WINDOW_COS:
    // Cosine window
    for (int i=0; i<n; i++){
      window[i] = sin(M_PI * i /(n-1));
    }
    break;
  case SFXC_WINDOW_PFB:
    // Hann windowed SINC
    for (int i=0; i<n; i++){
      double x = M_PI * (i-n/2) / fft_rot_size();
      FLOAT f = (i == n/2) ? 1 : sin(x) / x;
      window[i] = 0.5 * (1 - cos(2*M_PI*i/(n-1))) * f;
    }
    break;
  case SFXC_WINDOW_HAMMING:
    // Hamming window
    for (int i=0; i<n; i++){
      window[i] = 0.54 - 0.46 * cos(2*M_PI*i/(n-1));
    }
    break;
  case SFXC_WINDOW_HANN:
    // Hann window
    for (int i=0; i<n; i++){
      window[i] = 0.5 * (1 - cos(2*M_PI*i/(n-1)));
    }
    break;
  default:
    sfxc_abort("Invalid windowing function");
  }
}

// It is possible to flip the sidebandedness of a subband by flipping
// the sign of every other sample in the time domain.  This function
// constructs a vector to do this.
void
Delay_correction::create_flip() {
  const int n = (correlation_parameters.window == SFXC_WINDOW_PFB)?
                fft_rot_size() * SFXC_NTAPS : fft_rot_size();
  flip.resize(n);
  for (int i = 0; i < n; i++)
    flip[i] = ((i % 2) == 0) ? 1 : -1;
}

void Delay_correction::get_state(std::ostream &out) {
  out << "\t\t{\n"
      << "\t\t\"stream_nr\": " << stream_nr << ",\n"
      << "\t\t\"memory_pool_free\": " <<  output_memory_pool.number_free_element() << ",\n"
      << "\t\t\"current_time\": \"" << current_time.date_string(6) << "\",\n"
      << "\t\t\"n_input_buffer\": " << input_buffer->size() << ",\n"
      << "\t\t\"LO_offset\": " << LO_offset << ",\n"
      << "\t\t\"current_fft\": "<< current_fft << ",\n"
      << "\t\t\"n_ffts_per_integration\": " << n_ffts_per_integration << "\n"
      << "\t\t}";
}
