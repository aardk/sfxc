/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *            Aard Keimpema <keimpema@jive.nl>, 2008
 *
 * $Id:$
 *
 */

#include <sched.h>
#include "input_node_tasklet.h"
#include "utils.h"
#include "mark5a_reader.h"
#include "mark5b_reader.h"
#include "vlba_reader.h"
#include "vdif_reader.h"
#include "monitor.h"

typedef Input_node_types::Data_memory_pool  Data_memory_pool;
typedef shared_ptr<Data_memory_pool>        Data_memory_pool_ptr;
typedef shared_ptr<Data_reader>             Data_reader_ptr;

Input_node_tasklet *
get_input_node_tasklet_mark5a(Data_reader_ptr reader, Data_memory_pool_ptr memory_pool_,
                              Time ref_date) {
  shared_ptr<Mark5a_reader> mark5a_reader_ptr =
    shared_ptr<Mark5a_reader>( new Mark5a_reader(reader, ref_date) );

  return new Input_node_tasklet(mark5a_reader_ptr, memory_pool_);
}

Input_node_tasklet *
get_input_node_tasklet_vlba(Data_reader_ptr reader, Data_memory_pool_ptr memory_pool_,
                            Time ref_date) {
  // Maximal buffer size
  Input_data_format_reader::Data_frame data;
  data.buffer = memory_pool_->allocate();

  shared_ptr<VLBA_reader> vlba_reader_ptr =
    shared_ptr<VLBA_reader>( new VLBA_reader(reader, ref_date));

  return new Input_node_tasklet(vlba_reader_ptr, memory_pool_);
}

Input_node_tasklet *
get_input_node_tasklet_mark5b(Data_reader_ptr reader, Data_memory_pool_ptr memory_pool_,
                              Time ref_date) {
  typedef shared_ptr<Input_data_format_reader> Input_reader_ptr;
  Input_data_format_reader::Data_frame   data;
  data.buffer = memory_pool_->allocate();

  Input_reader_ptr   mark5b_reader_ptr(new Mark5b_reader(reader, data, ref_date));

  return new Input_node_tasklet(mark5b_reader_ptr, memory_pool_);
}

Input_node_tasklet *
get_input_node_tasklet_vdif(Data_reader_ptr reader, Data_memory_pool_ptr memory_pool_,
                            Time ref_date) {
  typedef shared_ptr<VDIF_reader> Input_reader_ptr;
  Input_data_format_reader::Data_frame   data;
  data.buffer = memory_pool_->allocate();

  Input_reader_ptr  vdif_reader_ptr(new VDIF_reader(reader, data, ref_date));
  return new Input_node_tasklet(vdif_reader_ptr, memory_pool_);
}


Input_node_tasklet *
get_input_node_tasklet(shared_ptr<Data_reader> reader,
                       TRANSPORT_TYPE type, Time ref_date) {
  SFXC_ASSERT(type != UNINITIALISED);
  shared_ptr<Data_memory_pool> memory_pool_(new Data_memory_pool(2 * 32 * 64));

  if (type == MARK5A) {
    return get_input_node_tasklet_mark5a(reader, memory_pool_, ref_date);
  }
  if (type == VLBA) {
    return get_input_node_tasklet_vlba(reader, memory_pool_, ref_date);
  }
  if (type == MARK5B) {
    return get_input_node_tasklet_mark5b(reader, memory_pool_, ref_date);
  }
  if (type == VDIF) {
    return get_input_node_tasklet_vdif(reader, memory_pool_, ref_date);
  }
  return NULL;
}



Input_node_tasklet::
Input_node_tasklet(Input_reader_ptr_ reader_ptr, Data_memory_pool_ptr memory_pool_)
    : reader_(reader_ptr, memory_pool_), delay_pool(10) 
{
  last_duration_ = 0;
  initialized = false;
}


void
Input_node_tasklet::
add_time_interval(Time &start_time, Time &stop_time, Time &leave_time) {
  akima_delays = delay_table.create_akima_spline(start_time, stop_time - start_time);
  /// Create a list of integer delay changes for the data writers
  Delay_memory_pool_element delay_list = delay_pool.allocate();
  delay_list.data().resize(0);
  delay_list.data().push_back(get_delay(start_time));

  start_time.set_sample_rate(sample_rate);
  stop_time.set_sample_rate(sample_rate);
  int64_t nsamples = stop_time.diff_samples(start_time);
  get_delays(start_time, nsamples, delay_list.data());
  data_writer_.add_delay(delay_list);
  data_writer_.add_time_interval(start_time, stop_time);

  // A new interval is added to the mark5 reader-tasklet 
  // We adjust the start and stop times to take into account the integer delay
  Time tbh = reader_.get_data_reader()->time_between_headers();
  Time delay_start = Time(akima_delays.delay(start_time)*1e6) - overlap_time;
  Time delay_stop = Time(akima_delays.delay(stop_time)*1e6) + overlap_time * 3;
  delay_start += min_extra_delay;
  delay_stop += max_extra_delay;
  int32_t start_frames = (int32_t) std::floor(delay_start/tbh);
  int32_t stop_frames = (int32_t) std::ceil(delay_stop/tbh);
  Time start_time_reader = start_time + tbh * start_frames;
  Time stop_time_reader = stop_time + tbh * stop_frames;
  Time leave_time_reader = leave_time + tbh * stop_frames;
  reader_.add_time_interval(start_time_reader, stop_time_reader, leave_time_reader);
}

void Input_node_tasklet::initialise(int num_tracks)
{
  if (reader_.get_data_reader()->get_transport_type() == VDIF && num_tracks == 0)
    channel_extractor_= Channel_extractor_tasklet_ptr( 
	new Channel_extractor_tasklet_VDIF(reader_.get_data_reader()));
  else
    channel_extractor_= Channel_extractor_tasklet_ptr( 
        new Channel_extractor_tasklet(reader_.get_data_reader()));

  channel_extractor_->connect_to(reader_.get_output_buffer());
  initialized = true;
}

Input_node_tasklet::~Input_node_tasklet() {
  if (channel_extractor_)
    channel_extractor_->empty_input_queue();
  data_writer_.empty_input_queue();

	PROGRESS_MSG( "Total duration:" << rttimer_processing_.measured_time() << " sec" );
	PROGRESS_MSG( "      reading:" << toMB(reader_.get_num_processed_bytes())/rttimer_processing_.measured_time() << " MB/s" );
	PROGRESS_MSG( "  channelizer:" << toMB(channel_extractor_->get_num_processed_bytes())/rttimer_processing_.measured_time() << " MB/s duration:" << channel_extractor_->get_sec() );
	PROGRESS_MSG( "      writing:" << toMB(data_writer_.get_num_processed_bytes())/data_writer_.get_sec() << " MB/s duration:" << data_writer_.get_sec() );
}


void
Input_node_tasklet::wait_termination() {
  /// Block until all the thread into the pool terminates.
  wait( pool_ );
}

void
Input_node_tasklet::start_tasklets() {
	rttimer_processing_.start();
  pool_.register_thread( channel_extractor_->start() );
  pool_.register_thread( reader_.start() );
}

void
Input_node_tasklet::stop_tasklets() {
  reader_.stop();
  if (channel_extractor_)
    channel_extractor_->stop();
  data_writer_.stop_threads();
  rttimer_processing_.stop();
}

void Input_node_tasklet::set_delay_table(Delay_table &table) {
  delay_table.add_scans(table);
}

void
Input_node_tasklet::
set_parameters(const Input_node_parameters &input_node_param,
               int station_number) {
  reader_.set_parameters(input_node_param);
  if(!initialized)
    initialise(input_node_param.n_tracks);

  channel_extractor_->set_parameters(input_node_param);

  size_t number_frequency_channels = input_node_param.channels.size();

  sample_rate=input_node_param.sample_rate();
  bits_per_sample=input_node_param.bits_per_sample();

  for (size_t i=0; i < number_frequency_channels; i++)
    data_writer_.add_channel();

  for (size_t i=0; i < number_frequency_channels; i++) {
    data_writer_.connect_to(i, channel_extractor_->get_output_buffer(i) );
    data_writer_.set_parameters(i, input_node_param, station_number);
  }

  overlap_time = input_node_param.overlap_time;
  int min_samples = 0, max_samples = 0;
  for (size_t i=0; i < number_frequency_channels; i++) {
    min_samples = 
      MIN(input_node_param.channels[i].extra_delay_in_samples, min_samples);
    max_samples =
      MAX(input_node_param.channels[i].extra_delay_in_samples, max_samples);
  }
  min_extra_delay = Time((min_samples * 1e6) / sample_rate);
  max_extra_delay = Time((max_samples * 1e6) / sample_rate);
}


Time
Input_node_tasklet::
get_current_time() {
  // Current time in [ms], if the delay correction hasn't progressed as far as
  // the reader we return the current time position of the delay correction
  Time reader_time = reader_.get_current_time();

  Time writer_time = data_writer_.get_current_time();
  if (writer_time < reader_time)
    return writer_time;
  else
    return reader_time;
}

void
Input_node_tasklet::add_data_writer(size_t i, Data_writer_sptr data_writer,
				    Time slice_start, Time slice_stop,
				    int64_t slice_samples) {
  /// Add a new timeslice to stream to the given data_writer into the
  /// data_writer queue.
  data_writer_.add_timeslice_to_stream(i, data_writer, slice_start,
				       slice_stop, slice_samples);
}

void
Input_node_tasklet::get_delays(Time start_time, int64_t nsamples, std::vector<Delay> &delay_list)
{
  Time stop_time = start_time;
  stop_time.inc_samples(nsamples);

  SFXC_ASSERT(stop_time>start_time);
  Delay dstart = get_delay(start_time);
  Delay dstop = get_delay(stop_time);
  bool delay_different=((dstop.bytes != dstart.bytes)||
                        (dstop.remaining_samples != dstart.remaining_samples));

  if(delay_different){
    if(nsamples == 1){
      delay_list.push_back(dstop);
    }else{
      get_delays(start_time, nsamples / 2, delay_list);
      start_time.inc_samples(nsamples / 2);
      get_delays(start_time, nsamples - (nsamples / 2), delay_list);
    }
  }else{
    bool rate_different = (akima_delays.rate(start_time) * akima_delays.rate(stop_time)) < 0;
    if((rate_different) && (nsamples >= 3)){
      get_delays(start_time, nsamples / 2, delay_list);
      start_time.inc_samples(nsamples / 2);
      get_delays(start_time, nsamples - (nsamples / 2), delay_list);
    }
  }
}

Delay
Input_node_tasklet::get_delay(Time time) {
  SFXC_ASSERT(delay_table.initialised());
  double delay_ = akima_delays.delay(time);
  int32_t delay_in_samples = (int32_t) std::floor(delay_*sample_rate+.5);

  int32_t delay_in_bytes = (int) floor((delay_in_samples-1)*1./(8/bits_per_sample));
  int32_t delay_in_remaining_samples =
                    delay_in_samples-delay_in_bytes*(8/bits_per_sample);
  if (delay_in_remaining_samples*bits_per_sample == 8) {
    delay_in_bytes++;
    delay_in_remaining_samples = 0;
  }

  SFXC_ASSERT((delay_in_bytes*(8/bits_per_sample) +
          delay_in_remaining_samples) ==
         delay_in_samples);
  return (Delay){time, delay_in_bytes, delay_in_remaining_samples};
}

void Input_node_tasklet::get_state(std::ostream &out) {
  out << "\t\"input_node_tasklet\": {\n" 
      << "\t\t\"initialized\": " << std::boolalpha << initialized <<",\n"
      << "\t\t\"delay_pool_free\": " << delay_pool.number_free_element() << "\n"
      << "\t},\n";
  reader_.get_state(out);
  if (initialized)
    channel_extractor_->get_state(out);
  data_writer_.get_state(out);
}
