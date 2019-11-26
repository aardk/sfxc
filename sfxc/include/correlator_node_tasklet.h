/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id: correlator_node.h 2206 2019-08-27 12:04:06Z kettenis $
 *
 */

#ifndef CORRELATOR_NODE_TASKLET_H
#define CORRELATOR_NODE_TASKLET_H
#include <queue>
#include <set>
#include "multiple_data_readers_controller.h"
#include "single_data_writer_controller.h"
#include "control_parameters.h"
#include "correlator_node_data_reader_tasklet.h"
#include "correlator_node_bit2float_tasklet.h"
#include "log_writer_mpi.h"
#include "correlation_core.h"
#include "correlation_core_pulsar.h"
#include "delay_correction.h"
#include <tasklet/tasklet_manager.h>
#include "timer.h"
#include "thread.h"

#include "monitor.h"
#include "eventor_poll.h"

/**
 *  The correlation_node_tasklet implements the main loop of the correlation.
 **/
class Correlator_node_tasklet : public Thread {
friend class Correlator_node_tasklet_controller;
public:
  typedef shared_ptr<Correlator_node_data_reader_tasklet> Bit_sample_reader_ptr;
  typedef shared_ptr<Data_reader>          Data_reader_ptr;
  typedef shared_ptr<Data_writer>          Data_writer_ptr;
  typedef shared_ptr<Delay_correction>     Delay_correction_ptr;

  bool has_requested;

  /// The states of the correlator_node.
  enum Status {
    // Initialise the Correlate node
    STOPPED=0,
    // The node is correlating
    CORRELATING
  };

  Correlator_node_tasklet(int nr_corr_node, bool psr_binning, bool phased_array);
  ~Correlator_node_tasklet();

  /// The the main_loop of the correlator node.
  void do_execute();

  /// Terminate the main_loop.
  void terminate();

  /// Callback function for adding a data_reader:
  void hook_added_data_reader(size_t reader, Data_reader_ptr data_reader);

  /// Callback function for adding a data_writer:
  void hook_added_data_writer(size_t writer, Data_writer_ptr data_writer);

  void add_delay_table(Delay_table &table, int sn1, int sn2);
  void add_uvw_table(Uvw_model &table, int sn1);

  void output_node_set_timeslice(int slice_nr, int stream_nr, int band, int accum,
				 int bytes, int nbins);

  void add_new_slice(const Correlation_parameters &parameters);
  void add_source_list(const std::map<std::string, int> &sources);

  void set_parameters(const Correlation_parameters &parameters);

  int get_correlate_node_number();


class Reader_thread : public Thread {
    std::vector< Bit_sample_reader_ptr >        bit_sample_readers_;

    struct job {
      std::vector<int> bits_per_sample;
      std::vector<int> stream_list;
      int station_streams_size;
    };

    Threadsafe_queue<struct job> queue_;
    bool readers_active_;
    /// Time spend in waiting for new slice
    Timer timer_waiting_;

    /// Time spend in waiting for data to arrive
    Timer timer_breading_;

    /// Time spend in really reading the data
    Timer timer_reading_;

  public:
    std::vector< Bit_sample_reader_ptr >& bit_sample_readers() {
      return bit_sample_readers_;
    }

    ~Reader_thread() {
      double total_duration = timer_waiting_.measured_time()+timer_reading_.measured_time();
      double ratio1 = ((100.0*timer_waiting_.measured_time())/total_duration);
      double ratio2 = ((100.0*timer_reading_.measured_time())/total_duration);
      PROGRESS_MSG( "reading ratio:(waiting: "<< ratio1 <<"%, reading:"<< ratio2 <<"%)" );
    }

  class Listener : public FdEventListener {
      Bit_sample_reader_ptr& reader_;

    public:
      Listener(Bit_sample_reader_ptr& ptr ) : reader_(ptr) {};

      void on_event(short event) {
        if ( reader_->has_work() )
          reader_->do_task();
      };

      void on_error(short event) {};
    };

    Eventor_poll eventsrc_;

    void do_execute() {
      for (unsigned int i=0;i<bit_sample_readers_.size();i++) {
        eventsrc_.add_listener( POLLIN,
                                bit_sample_readers_[i]->get_fd(),
                                new Listener( bit_sample_readers_[i] ) );
      }

      //eventsrc_.randomize();

      isrunning_ = true;
      DEBUG_MSG("reading thread started ! n_readers = " << bit_sample_readers_.size());
      try {
        while ( isrunning_ ) {
					if ( readers_active_ == false ) {
//            timer_waiting_.resume();
            fetch_new_time_slice();
//            timer_waiting_.stop();
          } else {
            timer_reading_.resume();
            /// Wait something happens.
            eventsrc_.wait_until_any_event();
            readers_active_ = false;
            for ( unsigned int i= 0;(i<bit_sample_readers_.size())&&(!readers_active_);i++) {
             readers_active_=bit_sample_readers_[i]->active();
            }
            timer_reading_.stop();
          }
        }
      } catch (QueueClosedException& exp) {
        DEBUG_MSG(" : The queue is closed !");
      }
    }

    void stop() {
      isrunning_ = false;
      queue_.close();
    }

    void fetch_new_time_slice() {
      // Immediately start prefetching the data:
      struct job jb = queue_.front();
      DEBUG_MSG("New input fetched:" << jb.station_streams_size);

      readers_active_=false;
      for (size_t i=0; i<bit_sample_readers_.size(); i++) {
        SFXC_ASSERT(bit_sample_readers_[i] !=
                    Bit_sample_reader_ptr());
        if (jb.stream_list[i] >= 0) {
          bit_sample_readers_[i]->set_parameters();
          readers_active_=true;
        }
      }
      queue_.pop();
    }

    /// This function add a new timeslice to read...
    void add_time_slice_to_read(const Correlation_parameters& parameters) {
      struct job jb;
      // First create a list of input streams
      jb.stream_list.resize(bit_sample_readers_.size());
      jb.station_streams_size = parameters.station_streams.size();
      for(int i = 0; i < bit_sample_readers_.size(); i++)
        jb.stream_list[i] = -1;
      for(int i = 0; i < jb.station_streams_size; i++)
        jb.stream_list[parameters.station_streams[i].station_stream] = i;

      //DEBUG_MSG("Add A Time slice:" << jb.station_streams_size );
      queue_.push(jb);
    }
  };

private:
  Reader_thread reader_thread_;
  Correlator_node_bit2float_tasklet bit2float_thread_;
  /// We need one thread for the integer delay correction
  ThreadPool threadpool_;
  void start_threads();
  void stop_threads();

  /// Main "usefull" function in which the real correlation computation is
  /// done.
  void correlate();

  bool pulsar_binning; // Set to true if pulsar binning is enabled
  bool phased_array; // Set to true if in phased array mode

  /// State variables:
  Status status;

  /// Number of the correlator node
  int nr_corr_node;
  bool isinitialized_;

  std::vector< Delay_correction_ptr >         delay_modules;
  Correlation_core                            *correlation_core, *correlation_core_normal;
  Correlation_core_pulsar                     *correlation_core_pulsar;

  Threadsafe_queue<Correlation_parameters>    integration_slices_queue;
  std::vector<int>                            delay_index;
  std::vector<Delay_table>                    delay_tables;
  std::vector<Uvw_model>                      uvw_tables;

  Timer bit_sample_reader_timer_, bits_to_float_timer_, delay_timer_, correlation_timer_;

#ifdef RUNTIME_STATISTIC
  QOS_MonitorSpeed reader_state_;
  QOS_MonitorSpeed delaycorrection_state_;
  QOS_MonitorSpeed correlation_state_;
  QOS_MonitorSpeed dotask_state_;
#endif //RUNTIME_STATISTIC
};

#endif // CORRELATOR_NODE_H
