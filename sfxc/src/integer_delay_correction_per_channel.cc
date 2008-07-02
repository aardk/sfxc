#include "integer_delay_correction_per_channel.h"

Integer_delay_correction_per_channel::
Integer_delay_correction_per_channel()
    : output_buffer_(new Output_buffer()),
    nr_output_bytes(-1),
    bits_per_sample(-1),
    sample_rate(-1),
    _current_time(-1),
    integration_time(-1),
    current_delay(1,1),
    memory_pool_(10)
    /**/
{
  SFXC_ASSERT(!memory_pool_.empty());
}

void
Integer_delay_correction_per_channel::do_task() {
  SFXC_ASSERT(has_work());
  SFXC_ASSERT(current_delay.first <= 0);

  // Acquire the input and output buffers
  Input_buffer_element input_element = input_buffer_->front();
  Output_buffer_element output_element;
  output_element.delay = (char)current_delay.second;

  // Compute the offset in bytes
  // If byte_offset < -nr_output_bytes
  //   The requested output lies before the input data, send invalid data
  // If byte_offset < 0
  //   The requested output lies partially before the input data,
  //   send invalid data for the first part and then the data from the input
  // If byte_offset > 0 && byte_offset + nr_output_bytes < input_element.size()
  //   the data requested is in the input file and lies within one data block,
  //   this is the normal case
  // If byte_offset > 0 && byte_offset + nr_output_bytes >= input_element.size()
  //   the data requested is in the input file and lies within two data block,
  //   - send the data from the first block,
  //   - send a message to release the first block,
  //   - send the data from the second block

  int byte_offset =
    (_current_time-input_element.start_time) *
    sample_rate*bits_per_sample/8/1000000 +
    current_delay.first;

  if (byte_offset >= 0) {
    // we have valid data

    if (byte_offset >=
        (int)input_element.channel_data.data().data.size()) {
      // This can happen when we go to a next integration slice
      // And the integer delay changes at the same time
      // Release the current block
      input_buffer_->pop();
      return;
    }

    // Default case with normal data
    SFXC_ASSERT (byte_offset >= 0);
    int input_data_size = input_element.channel_data.data().data.size();
    if ((byte_offset + nr_output_bytes) < input_data_size) {
      // Normal case where all data lies within one input block

      { // Set invalid samples
        output_element.invalid_samples_begin = 0;
        output_element.nr_invalid_samples = 0;
        // Check whether we have invalid samples
        if ((byte_offset*8/bits_per_sample+output_element.delay) <
            (input_element.invalid_samples_begin + input_element.nr_invalid_samples)) {
          // nr of invalid samples from the beginning of the fft
          int nr_invalid_samples =
            /* stop  */ (input_element.invalid_samples_begin +
                         input_element.nr_invalid_samples) -
            /* start */ (byte_offset * 8 / bits_per_sample + output_element.delay);
          output_element.nr_invalid_samples =
            std::min(nr_output_bytes*8/bits_per_sample, nr_invalid_samples);
        }
        SFXC_ASSERT(output_element.nr_invalid_samples >= 0);
      }

      // Send data
      output_element.channel_data = input_element.channel_data;
      output_element.first_byte = byte_offset;
      output_element.nr_bytes = nr_output_bytes+1;
      output_buffer_->push(output_element);

    } else {
      // Case where the data lies within two input blocks

      int bytes_in_current_block = input_data_size-byte_offset;
      // Send first block of data
      output_element.channel_data = input_element.channel_data;
      output_element.first_byte = byte_offset;
      output_element.nr_bytes = bytes_in_current_block;

      // Get the second block of data
      input_buffer_->pop();
      SFXC_ASSERT(!input_buffer_->empty());
      input_element = input_buffer_->front();

      { // Set invalid samples
        output_element.invalid_samples_begin = 0;
        output_element.nr_invalid_samples = 0;

        // Simplification because the only invalid samples are missing data
        // or mark5-headers
        SFXC_ASSERT(input_element.invalid_samples_begin == 0);

        // start of the invalid data is the begin of the second data block
        output_element.invalid_samples_begin =
          bytes_in_current_block*8/bits_per_sample - output_element.delay;

        int samples_in_second_block =
          (nr_output_bytes-bytes_in_current_block)*8/bits_per_sample + output_element.delay;

        output_element.nr_invalid_samples =
          std::min(samples_in_second_block, input_element.nr_invalid_samples);
      }
      output_buffer_->push(output_element);

      // Send second block of data
      // Don't send the delay again:
      output_element.delay = -1;
      output_element.channel_data = input_element.channel_data;
      output_element.first_byte = 0;
      output_element.nr_bytes = nr_output_bytes+1-bytes_in_current_block;
      output_buffer_->push(output_element);
    }
  } else if (byte_offset < -(nr_output_bytes+1)) {
    // Completely random data

    // Do not do the bit offset
    output_element.delay = char(0);
    output_element.invalid_samples_begin = 0;
    output_element.nr_invalid_samples = nr_output_bytes*8/bits_per_sample;

    // Send random data
    output_element.channel_data = allocate_random_element();
    output_element.first_byte = 0;
    output_element.nr_bytes = nr_output_bytes+1;
    output_buffer_->push(output_element);
  } else {
    SFXC_ASSERT(byte_offset < 0);

    // Partially random data
    output_element.invalid_samples_begin = 0;
    // Either all data is random, or part of the second block contains valid data
    output_element.nr_invalid_samples =
      std::min(nr_output_bytes*8/bits_per_sample,
               // random data in the first block
               -byte_offset*8/bits_per_sample - output_element.delay +
               // random data in the second block
               input_element.nr_invalid_samples);

    // Send random data
    output_element.channel_data = allocate_random_element();
    output_element.first_byte = 0;
    output_element.nr_bytes = -byte_offset;
    output_buffer_->push(output_element);

    // Send real data
    // Don't send the delay again:
    output_element.delay = -1;
    output_element.channel_data = input_element.channel_data;
    output_element.first_byte = 0;
    output_element.nr_bytes = nr_output_bytes+1+byte_offset;
    output_buffer_->push(output_element);

  }

  // Increase the time to the beginning of the next fft
  _current_time += delta_time;
  current_delay = get_delay(_current_time);


  // Check whether the next fft crosses an integration border
  // Then we set the time to the beginning of the next integration slice
  if ((_current_time/integration_time) !=
      ((_current_time+delta_time-1)/integration_time)) {
    // Continue with the next integration slice
    _current_time =
      ((_current_time+delta_time-1)/integration_time)*integration_time;
    current_delay = get_delay(_current_time);
  }
}


void
Integer_delay_correction_per_channel::fetch_next_time_interval() {
  SFXC_ASSERT(delay_table.initialised());

  // We retreive the current interval
  current_interval_ = intervals_.front_and_pop();
  SFXC_ASSERT( !current_interval_.empty() );

  _current_time = current_interval_.start_time_;
  current_delay = get_delay(_current_time);

  //DEBUG_MSG(__PRETTY_FUNCTION__<<":: set interval["<<current_interval_.start_time_ <<":"<< current_interval_.stop_time_<<"]")
}

void
Integer_delay_correction_per_channel::add_time_interval(
	uint64_t start, uint64_t stop) {
  SFXC_ASSERT( start < stop );
  intervals_.push( Time_interval(start, stop) );
}

bool
Integer_delay_correction_per_channel::has_work() {
  SFXC_ASSERT(output_buffer_ != Output_buffer_ptr());

  if ( _current_time >= current_interval_.stop_time_ ) {
    if ( intervals_.empty() ) return false;
    fetch_next_time_interval();
  }

  if (sample_rate <= 0)
    return false;
  if (current_delay.first > 0)
    return false;
  if (input_buffer_ == Input_buffer_ptr())
    return false;
  if (input_buffer_->empty())
    return false;
  if (memory_pool_.empty())
    return false;

  // Check whether we cross a block boundary:
  Input_buffer_element &input_element = input_buffer_->front();
  int byte_offset =
    (_current_time-input_element.start_time)*sample_rate*bits_per_sample/8/1000000 +
    current_delay.first;
  if (size_t(byte_offset + nr_output_bytes +1) >=
      input_element.channel_data.data().data.size()) {
    if (input_buffer_->size() < 2)
      return false;
  }

  return true;
}

void
Integer_delay_correction_per_channel::
connect_to(Input_buffer_ptr buffer) {
  input_buffer_ = buffer;
}

Integer_delay_correction_per_channel::Output_buffer_ptr
Integer_delay_correction_per_channel::
get_output_buffer() {
  return output_buffer_;
}

double
Integer_delay_correction_per_channel::delay(int64_t time) {
  return delay_table.delay(time+delta_time/2);
}


Integer_delay_correction_per_channel::Delay_type
Integer_delay_correction_per_channel::get_delay(int64_t time) {
  SFXC_ASSERT(delay_table.initialised());
  SFXC_ASSERT(delta_time > 0);
  SFXC_ASSERT(delta_time%2 == 0);
  double delay_ = delay(time);
  int delay_in_samples = (int)std::floor(delay_*sample_rate+.5);

  // All because modulo doesn't work for negative values
  // delay_in_bytes = std::floor(delay_in_samples/(8./bits_per_sample))
  SFXC_ASSERT(delay_in_samples < 0);
  int delay_in_bytes = -((-delay_in_samples)/(8/bits_per_sample))-1;
  int delay_in_remaining_samples =
    delay_in_samples-delay_in_bytes*(8/bits_per_sample);
  if (delay_in_remaining_samples*bits_per_sample == 8) {
    delay_in_bytes++;
    delay_in_remaining_samples = 0;
  }

  SFXC_ASSERT((delay_in_bytes <= 0) &&
         (delay_in_remaining_samples < 8));
  SFXC_ASSERT((delay_in_bytes*(8/bits_per_sample) +
          delay_in_remaining_samples) ==
         delay_in_samples);

  return Delay_type(delay_in_bytes, delay_in_remaining_samples);
}

void
Integer_delay_correction_per_channel::
set_delay_table(Delay_table_akima &table) {
  delay_table = table;
}

void
Integer_delay_correction_per_channel::
set_parameters(const Input_node_parameters &parameters,
               int node_nr) {
  bits_per_sample = parameters.bits_per_sample();

  SFXC_ASSERT((parameters.number_channels*bits_per_sample)%8 == 0);
  // The offset is not counted
  nr_output_bytes = parameters.number_channels*bits_per_sample/8;
  sample_rate = parameters.sample_rate();

  SFXC_ASSERT(((nr_output_bytes*(8/bits_per_sample))*1000000) % sample_rate== 0);
  delta_time = (nr_output_bytes*(8/bits_per_sample))*1000000/sample_rate;
  integration_time = parameters.integr_time*1000;

  nr_bytes_per_integration_slice =
    Control_parameters::nr_bytes_per_integration_slice_input_node_to_correlator_node
    (parameters.integr_time,
     sample_rate,
     bits_per_sample,
     parameters.number_channels);


}

int
Integer_delay_correction_per_channel::
bytes_of_output() {
  return nr_bytes_per_integration_slice;
}


Integer_delay_correction_per_channel::Input_data_block
Integer_delay_correction_per_channel::
allocate_random_element() {
  Input_data_block result = memory_pool_.allocate();
  if (result.data().data.size() != (size_t)nr_output_bytes+1) {
    result.data().data.resize(nr_output_bytes+1);
#ifdef SFXC_INVALIDATE_SAMPLES
#ifdef SFXC_CHECK_INVALID_SAMPLES

    for (int i=0; i<nr_output_bytes+1; i++) {
      result.data().data[i] = INVALID_PATTERN;
    }
#endif
#else
    // Randomize data
    for (int i=0; i<nr_output_bytes+1; i++) {
      result.data().data[i] = (char)park_miller_random();
    }
#endif

  }

  return result;
}

// Empty the input queue, called from the destructor of Input_node
void Integer_delay_correction_per_channel::empty_input_queue() {
  while (!input_buffer_->empty()) {
    input_buffer_->pop();
  }
}

