#ifndef INPUT_NODE_TASKLET_H
#define INPUT_NODE_TASKLET_H

#include <queue>

#if __cplusplus >= 201103L
#include <memory>
using std::shared_ptr;
#else
#include <tr1/memory>
using std::tr1::shared_ptr;
#endif

#include "tasklet/tasklet.h"
#include "thread.h"

#include "data_reader.h"
#include "data_writer.h"
#include "delay_table_akima.h"
#include "control_parameters.h"

#include "input_data_format_reader_tasklet.h"
#include "channel_extractor_tasklet.h"
#include "channel_extractor_tasklet_vdif.h"

#include "input_node_data_writer_tasklet.h"
#include "correlator_time.h"

#include "rttimer.h"

// for RUNTIME_STATISTIC
#include "monitor.h"

class Input_node_tasklet {
public:
  typedef Input_data_format_reader              Input_reader_;
  typedef shared_ptr<Input_reader_>             Input_reader_ptr_;
  typedef Input_data_format_reader_tasklet      Input_reader_tasklet_;
  typedef Input_node_types::Data_memory_pool    Data_memory_pool;
  typedef shared_ptr<Data_memory_pool>          Data_memory_pool_ptr;
  typedef shared_ptr<Channel_extractor_tasklet> Channel_extractor_tasklet_ptr;

  // The mark5a-reader and the first data block
  Input_node_tasklet(Input_reader_ptr_ input_reader_ptr, Data_memory_pool_ptr memory_pool_);

  ~Input_node_tasklet();

  void initialise(int num_tracks);
  void start_tasklets();
  void stop_tasklets();
  void wait_termination();

  /// Write state of input node tasklet for debug purposes
  void get_state(std::ostream &out);

  /// Sets a new time interval for which it should output data
  /// (typically the duration of a scan, or part thereof).
  /// \param start_time in milliseconds
  /// \param stop_time in milliseconds
  /// It is not possible to go back in time, as the data might be
  /// streamed in to the input node and not be buffered anymore.
  /// Time is in milliseconds
  void add_time_interval(Time &start_time, Time &stop_time, Time &leave_time);

  // Inherited from Input_node_tasklet
  void set_delay_table(Delay_table &delay);
  void set_parameters(const Input_node_parameters &input_node_param,
                      int station_number);


  /// Returns the current time in microseconds
  Time get_current_time();

  /// Sets the output writer for channel i
  void add_data_writer(size_t i, Data_writer_sptr data_writer,
		       Time slice_start, Time slice_stop,
		       int64_t slice_samples);

  /// Compute a list of delays, note we only store the times(+delay) 
  /// where the integer delay changes
  void get_delays(Time start_time, int64_t nsamples, std::vector<Delay> &delay_list);
  /// Calculates the delay at time
  Delay get_delay(Time time);
private:
  /// All the thread created in this class are stored in the thread pool
  ThreadPool pool_;

  /// Memory pool to store the delays for each time interval
  Delay_memory_pool delay_pool;

  /// we need one thread for reading
  Input_reader_tasklet_            reader_;

  /// We need one thread for the allocation
  Channel_extractor_tasklet_ptr    channel_extractor_;

  /// We need one thread for the writing
  Input_node_data_writer_tasklet   data_writer_;

  Timer rttimer_processing_;
  double last_duration_;

  Delay_table delay_table;
  Delay_table_akima akima_delays;

  bool initialized;
  uint64_t sample_rate;
  int bits_per_sample;
  int64_t size_slice; // Number of samples for one integration slice
  Time overlap_time;  // Size of buffer additional data needed for dedispersion filter
  Time min_extra_delay, max_extra_delay;
};


/** Returns an input_node_tasklet for the data reader.
 * It determines the number of tracks from the data
 **/
Input_node_tasklet *
get_input_node_tasklet(shared_ptr<Data_reader> reader,
                       TRANSPORT_TYPE type, Time ref_date);

#endif // INPUT_NODE_TASKLET_H
