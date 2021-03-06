#include <math.h>
#include "bit2float_worker.h"

const FLOAT sample_value_ms[] = {
                                  -7, -2, 2, 7
                                };

const FLOAT sample_value_m[]  = {
                                  -5, 5
                                };

Bit2float_worker::Bit2float_worker(int stream_nr_, bit_statistics_ptr statistics_)
  : output_buffer_(new Output_queue()),
    fft_size(-1),
    bits_per_sample(-1),
    sample_rate(-1),
    tsys_freq(80),
    memory_pool_(32),
    stream_nr(stream_nr_),
    n_ffts_per_integration(0), current_fft(0), state(IDLE), statistics(statistics_)
    /**/
{
  SFXC_ASSERT(!memory_pool_.empty());
  // Lookup tables used in the bit2float conversion
  for (int i=0; i<256; i++) {
    lookup_table[i][0] = sample_value_ms[i & 3];
    lookup_table[i][1] = sample_value_ms[(i>>2) & 3];
    lookup_table[i][2] = sample_value_ms[(i>>4) & 3];
    lookup_table[i][3] = sample_value_ms[(i>>6) & 3];
    for (int j=0; j<8 ; j++)
      lookup_table_1bit[i][j] = sample_value_m[(i>>j) & 1];
  }
}

int64_t
Bit2float_worker::do_task() {
  std::vector<unsigned char> &inp_data = input_buffer_->data;
  int64_t samples_written = 0;

  size_t inp_size = inp_data.size();

  Output_pool_data &out_frame = out_element.data();
  size_t output_buffer_size = out_frame.data.size();

  while (out_index != output_buffer_size) {
    uint64_t read = input_buffer_->read;
    uint64_t write = input_buffer_->write;

    SFXC_ASSERT(sample_in_byte < 8);
    SFXC_ASSERT(read <= write);
    switch(state) {
    case IDLE:
    {
      if ((write - read) < 3){
        return samples_written;
      }

      uint8_t header = inp_data[read++ % inp_size];
      switch(header) {
      case HEADER_DATA:
        bytes_left = inp_data[read++ % inp_size];
        bytes_left |= (inp_data[read++ % inp_size] << 8);
        SFXC_ASSERT(bytes_left > 0);
        state = SEND_DATA;
        break;
      case HEADER_DELAY:
      {
        int8_t new_delay = inp_data[read++ % inp_size];
        if (cur_delay < 0) {
          // the initial delay at the beginning of the integration
          sample_in_byte = new_delay;
        } else {
          int diff = new_delay - cur_delay;
          if ((diff == -1) || (diff > 1)) {
            // insert a sample into the stream
            memset(&out_frame.data[out_index], 0, sizeof(FLOAT));
            out_index++;
	  } else {
            // remove a sample from the stream
            SFXC_ASSERT(sample_in_byte == 0);
            sample_in_byte = 1;
          }
        }
        cur_delay = new_delay;
        break;
      }
      case HEADER_INVALID:{
        invalid_left = inp_data[read++ % inp_size];
        invalid_left |= (inp_data[read++ % inp_size] << 8);
        statistics->inc_invalid(invalid_left);
        int start = current_fft * fft_size;
        invalid.push_back((Invalid){start + out_index, invalid_left});
        state = SEND_INVALID;
        break;
      }
      default:
        SFXC_ASSERT_MSG(false, "Read invalid header from input buffer");
      }
      break;
    }
    case SEND_INVALID:
    {
      size_t output_buffer_free = output_buffer_size - out_index;
      size_t invalid_to_write = std::min(output_buffer_free, invalid_left);
      // check if there was a delay change and we should remove samples from the stream
      if (sample_in_byte > 0) {
        // NB after receiving the first delay we can have sample_in_byte > 1
        size_t samples_to_drop = std::min(sample_in_byte, invalid_to_write);
        sample_in_byte -= samples_to_drop;
        invalid_to_write -= samples_to_drop;
        invalid_left -= samples_to_drop;
      }

      if (invalid_to_write > 0) {
        memset(&out_frame.data[out_index], 0, invalid_to_write * sizeof(FLOAT));
        invalid_left -= invalid_to_write;
        out_index += invalid_to_write;
        samples_written += invalid_to_write;
      }
      if (invalid_left == 0)
        state = IDLE;
      break;
    }
    case SEND_DATA:
    {
      if (read == write) {
        input_buffer_->read = read;
        return samples_written;
      }

      SFXC_ASSERT((write - read) > 0);

      int samples_per_byte = (8 / bits_per_sample);
      size_t output_buffer_free = output_buffer_size - out_index;
      size_t bytes_in_input_buffer = std::min((size_t)(write - read), bytes_left);
      size_t samp_to_write = std::min(bytes_in_input_buffer * samples_per_byte - sample_in_byte,
                                      output_buffer_free);
      SFXC_ASSERT(bytes_left > 0);
      SFXC_ASSERT(bytes_in_input_buffer > 0);
      SFXC_ASSERT(output_buffer_free > 0);

      bytes_left -= (samp_to_write + sample_in_byte) / samples_per_byte;
      sample_in_byte = bit2float(&out_frame.data[out_index], sample_in_byte, samp_to_write, &read);
      out_index += samp_to_write;
      if (bytes_left == 0) {
        SFXC_ASSERT(sample_in_byte == 0);
        state = IDLE;
      }
      samples_written += samp_to_write;
      break;
    }
    case PURGE_STREAM:
      size_t bytes_to_advance = std::min((size_t)(write - read), bytes_left);
      read += bytes_to_advance;
      bytes_left -= bytes_to_advance;
      if (bytes_left == 0) {
        if ((write - read) > 0) {
          switch(inp_data[read % inp_size]) {
          case HEADER_ENDSTREAM:
            read += 1;
            state = IDLE;
            break;
          case HEADER_DELAY:
            if ((write - read) >= 2)
              read += 2;
            break;
          case HEADER_INVALID:
            if ((write - read) >= 3)
              read += 3;
            break;
          case HEADER_DATA:
            if ((write - read) >= 3) {
              read += 1;
              bytes_left = inp_data[read++ % inp_size];
              bytes_left |= (inp_data[read++ % inp_size] << 8);
              SFXC_ASSERT(bytes_left > 0);
            }
            break;
          }
        }
      }
      break;
    }
    SFXC_ASSERT(read <= write);
    input_buffer_->read = read;
  }
  if (out_index == output_buffer_size) {
    output_buffer_->push(out_element);
    current_fft += out_element.data().nfft;
    if (current_fft == n_ffts_per_integration)
      state = PURGE_STREAM;
    else
      allocate_element();
  }

  return samples_written;
}


bool
Bit2float_worker::has_work() {
  if( queue_.size() > 0) {
    set_parameters();
  }

  if (sample_rate <= 0)
    return false;

  if (current_fft == n_ffts_per_integration)
    return false;

  if((input_buffer_->bytes_read() == 0)&&(state!=SEND_INVALID))
    return false;

//  if (memory_pool_.empty())
  if (memory_pool_.number_free_element()<2)
    return false;

  return true;
}

void
Bit2float_worker::
connect_to(Input_buffer_ptr buffer) {
  input_buffer_ = buffer;
}

Bit2float_worker::Output_queue_ptr
Bit2float_worker::
get_output_buffer() {
  return output_buffer_;
}

void
Bit2float_worker::
set_new_parameters(const Correlation_parameters &parameters, Delay_table_akima &delay_table) {
  int stream_idx = 0;
  while ((stream_idx < parameters.station_streams.size()) &&
         (parameters.station_streams[stream_idx].station_stream != stream_nr))
    stream_idx++;
  if (stream_idx == parameters.station_streams.size()) {
    // Data stream is not participating in current time slice
    return;
  }
  Bit2float_parameters new_parameters;
  new_parameters.stream_start = parameters.stream_start;
  new_parameters.bits_per_sample = parameters.station_streams[stream_idx].bits_per_sample;
  new_parameters.sample_rate = parameters.station_streams[stream_idx].sample_rate;
  new_parameters.base_sample_rate = parameters.sample_rate;
  new_parameters.tsys_freq = parameters.station_streams[stream_idx].tsys_freq;
  new_parameters.fft_size_delaycor = parameters.fft_size_delaycor;
  new_parameters.fft_size_correlation = parameters.fft_size_correlation;
  int nfft_min = std::max(parameters.fft_size_correlation / parameters.fft_size_delaycor, 1);
  new_parameters.n_ffts_per_buffer =
    (std::max(CORRELATOR_BUFFER_SIZE / parameters.fft_size_delaycor, nfft_min) *
     parameters.station_streams[stream_idx].sample_rate) / parameters.sample_rate;
  new_parameters.delay_in_samples = 
    delay_table.delay(parameters.stream_start) * new_parameters.sample_rate;

  new_parameters.n_ffts_per_integration =
     (int64_t) new_parameters.sample_rate * parameters.slice_size / 
     ((int64_t) new_parameters.base_sample_rate * parameters.fft_size_delaycor);

  queue_.push(new_parameters);
}

void
Bit2float_worker::
set_parameters() {
  Bit2float_parameters &new_parameters = queue_.front();
  bits_per_sample = new_parameters.bits_per_sample;
  sample_rate = new_parameters.sample_rate;
  base_sample_rate = new_parameters.base_sample_rate;
  tsys_freq = new_parameters.tsys_freq;

  fft_size = new_parameters.fft_size_delaycor;
  SFXC_ASSERT(((int64_t)fft_size * 1000000) % sample_rate == 0);
  nfft_max = new_parameters.n_ffts_per_buffer;
  n_ffts_per_integration = new_parameters.n_ffts_per_integration;
  statistics->reset_statistics(bits_per_sample, sample_rate, base_sample_rate);

  current_fft = 0;
  invalid_left = 0;
  sample_in_byte = 0;
  cur_delay = -1;
  allocate_element();
  invalid.resize(0);
  if (bits_per_sample == 2) {
    tsys_count = std::floor(new_parameters.delay_in_samples + 0.5);
    // Calculate number of samples within the full cycle.
    tsys_count %= (sample_rate / tsys_freq);
    while (tsys_count < 0)
      tsys_count += sample_rate / tsys_freq;
    tsys_on = (tsys_count < (sample_rate / (2 * tsys_freq)));
    // Calculate number of samples within this falf-cycle.
    tsys_count %= (sample_rate / (2 * tsys_freq));
    // Calculate number of samples left for this half-cycle.
    tsys_count = ((sample_rate / (2 * tsys_freq)) - tsys_count);
    // Convert to bytes.  
    tsys_count /= 4;
  }
  queue_.pop();
}

// Empty the input queue, called from the destructor of Input_node
void Bit2float_worker::empty_input_queue() {
  input_buffer_->read = 0;
  input_buffer_->write = 0;
}

Bit2float_worker_sptr
Bit2float_worker::new_sptr(int stream_nr_, bit_statistics_ptr statistics_){
  return Bit2float_worker_sptr(new Bit2float_worker(stream_nr_, statistics_));
}

int 
Bit2float_worker::bit2float(FLOAT *output, int start, int nsamples, uint64_t *readp) {
  std::vector<unsigned char> &input_data = input_buffer_->data;
  int dsize = input_data.size();
  uint64_t read = *readp;
  // avoid the overhead of the shared pointer by dereferencing it
  bit_statistics *stats = statistics.get();

  int iout = 0;
  if (bits_per_sample == 2) {
    // Write the first byte
    int samp_to_write = std::min(nsamples, 4 - start);
    memcpy((char*)&output[iout],
           &lookup_table[(int)input_data[read % dsize]][start],
           samp_to_write * sizeof(FLOAT));
    stats->inc_counter(input_data[read % dsize], tsys_on);
    if (--tsys_count == 0) {
      tsys_count = (sample_rate / (2 * tsys_freq)) / 4;
      tsys_on = !tsys_on;
    }
    nsamples -= samp_to_write;
    iout += samp_to_write;
    if (samp_to_write + start == 4)
      read++;
    else
      return start + samp_to_write;

    // Write the main bulk of the data
    int nbytes = nsamples / 4;
    while (nbytes > 0) {
      int index = read % dsize;
      int towrite = std::min(nbytes, dsize-index);
      for (int end = index+towrite; index < end; index++) {
        memcpy((char*)&output[iout],
               &lookup_table[(int)input_data[index]][0],
               4 * sizeof(FLOAT));
        stats->inc_counter(input_data[index], tsys_on);
	if (--tsys_count == 0) {
	  tsys_count = (sample_rate / (2 * tsys_freq)) / 4;
	  tsys_on = !tsys_on;
	}
        iout += 4;
      }
      nbytes -= towrite;
      nsamples -= towrite * 4;
      read += towrite;
    }
    // Write the final byte
    if (nsamples > 0) {
      memcpy((char*)&output[iout],
             &lookup_table[(int)input_data[read % dsize]][0],
             nsamples * sizeof(FLOAT));
      *readp = read;
      return nsamples;
    }
  } else { // 1 bit samples
    SFXC_ASSERT(bits_per_sample == 1);
    // Write the first byte
    int samp_to_write = std::min(nsamples, 8 - start);
    memcpy((char*)&output[iout],
           &lookup_table_1bit[(int)input_data[read % dsize]][start],
           samp_to_write * sizeof(FLOAT));
    stats->inc_counter(input_data[read % dsize], true);
    nsamples -= samp_to_write;
    iout += samp_to_write;
    if (samp_to_write+start == 8)
      read++;
    else
      return start+samp_to_write;

    // Write the main bulk of the data
    int nbytes = nsamples / 8;
    while (nbytes>0) {
      int index = read % dsize;
      int towrite = std::min(nbytes, dsize-index);
      for (int end = index+towrite; index < end; index++) {
        memcpy((char*)&output[iout],
               &lookup_table_1bit[(int)input_data[index]][0],
               8 * sizeof(FLOAT));
        stats->inc_counter(input_data[index], true);
        iout += 8;
      }
      nbytes -= towrite;
      nsamples -= towrite * 8;
      read += towrite;
    }
    // Write the final byte
    if (nsamples>0) {
      memcpy((char*)&output[iout],
             &lookup_table_1bit[(int)input_data[read % dsize]][0],
             nsamples * sizeof(FLOAT));
      *readp = read;
      return nsamples;
    }
  }

  *readp = read;
  return 0;
}

void
Bit2float_worker::allocate_element(){
  int nfft = std::min(nfft_max, n_ffts_per_integration - current_fft);
  int nsamples = nfft * fft_size;
  out_element = memory_pool_.allocate();
  
  if(out_element.data().data.size() != nsamples)
    out_element.data().data.resize(nsamples);
  out_element.data().nfft = nfft;
  out_index=0;
}

std::vector<Bit2float_worker::Invalid> *
Bit2float_worker::get_invalid(){
  return &invalid;
}

void Bit2float_worker::get_state(std::ostream &out) {
  out << "\t\t{\n"
      << "\t\t\"memory_pool_free\": " << memory_pool_.number_free_element() << ",\n"
      << "\t\t\"current_fft\": " << current_fft << ",\n"
      << "\t\t\"n_ffts_per_integration\": " << n_ffts_per_integration << ",\n"
      << "\t\t\"state\": ";
  switch(state) {
  case IDLE:
    out << "\"IDLE\"\n";
    break;
  case SEND_INVALID:
    out << "\"SEND_INVALID\"\n";
    break;
  case SEND_DATA:
    out << "\"SEND_DATA\"\n";
    break;
  case PURGE_STREAM:
    out << "\"PURGE_STREAM\"\n";
    break;
  }
  out << "\t\t}";
}
