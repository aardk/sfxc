/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#ifndef OUTPUT_NODE_H
#define OUTPUT_NODE_H

#include "node.h"
#include "multiple_data_readers_controller.h"
#include "multiple_data_writers_controller.h"
#include "output_header.h"

#include <memory_pool.h>

#include <fstream>
#include <map>
#include <queue>

class Output_node;

/// Controller for output node specific commands
class Output_node_controller : public Controller {
public:
  Output_node_controller(Output_node &node);

  Process_event_status process_event(MPI_Status &status);

private:
  Output_node &node;

  std::ofstream phasecal_file;
  std::ofstream tsys_file;
};

class Stream_param {
 public:
  Stream_param(int stream_, int band_, bool accum_) {
    stream = stream_;
    band = band_;
    accum = accum_;
  };

  int stream;
  int band;
  bool accum;
};

/**
 * The output node will receive a message from the controller node where
 * to store the data and it allows connections from the correlate node to
 * be opened. The node sorts the received data from the correlate nodes
 * and stores it for further processing. The output node has to make the
 * received data available to the user and it should be archived in a
 * proper way.
 *
 * \ingroup Node
 **/
class Output_node : public Node {
public:
  // Input types
  typedef Multiple_data_readers_controller::value_type input_value_type;
  typedef std::map<int, Stream_param>            Input_stream_order_map;
  typedef Input_stream_order_map::value_type     Input_stream_order_map_value;

  /**
   * Manages the input from one correlator node. The input stream is
   * used to store the data from one correlator node. The data is read
   * from the reader and the slice_size queue contains the size of the
   * subsequent time slices. This has to be a queue (not a single
   * slice) as the output node might have several slices from
   * correlator that it has to output in the case that the output node
   * becomes too slow to store the data. In practice this is not the case
   **/
  class Input_stream {
  public:
    Input_stream(boost::shared_ptr<Data_reader> reader);

    /** Fills the buffer with as much data as possible and returns the number of
     * bytes written.
     **/
    int read_bytes(std::vector<char> &buffer);
    /** returns whether we reached the end of the current time slice
     **/
    bool end_of_slice();

    /** sets the length of a new time slice
     **/
    void set_length_time_slice(int64_t nBytes, int nbins);

    /** Goto the next data slice **/
    void goto_next_slice(int &new_slice_size, int &new_nbins);
    struct Slice{
      int64_t nBytes;
      int32_t nBins;
    };
  private:
    // Data_reader from which the input data can be read
    boost::shared_ptr<Data_reader> reader;
    // list with sizes of the time slices
    std::queue<Slice> slice_size;
    // read offset within slice
    size_t offset;
  };

  Output_node(int rank, Log_writer *writer, int buffer_size = 10);
  Output_node(int rank, int buffer_size = 10);
  void initialise();

  ~Output_node();

  void start();
  void terminate();

  enum STATUS {
    STOPPED=0,
    START_NEW_SLICE,
    READ_INPUT,
    ACCUMULATE_INPUT,
    WRITE_OUTPUT,
    END_SLICE,
    END_NODE
  };

  /**
   * This function is called once at the beginning of the correlation
   * to write the global header. Make sure it is called before any
   * correlator node sends data, otherwise the output correlation file
   * might not start with the global header.
   **/
  void write_global_header(const Output_header_global &global_header);
  /**
   * Notifies the output node that there is a block of data arriving
   * from a correlator node.
   **/
  void set_order_of_input_stream(int stream, int order, int band, int accum,
				 size_t size, int nbins);

  /**
   * This function sets the total number of time slices so that the
   * output node knows when it's done writing data. The number of
   * slices is sent from the manager node.
   **/
  void set_number_of_time_slices(int n_time_slices);

  // Callback functions:
  void hook_added_data_reader(size_t reader);
  void hook_added_data_writer(size_t writer);

private:

  /**
   * Function that writes a bit of data.
   * Returns whether it wrote something
   **/
  bool write_output(int nBytes);

  /// The number of output files we are writing to
  int n_data_writers;

  /// Data buffer into which the data from the correlator node is read
  std::vector<char>                   input_buffer;
  size_t			      total_bytes_read;

  std::vector<std::vector<char> >     accum_buffer;
  std::vector<int>		      integration;

  // Controllers:
  Output_node_controller              output_node_ctrl;
  Multiple_data_readers_controller    data_readers_ctrl;
  Multiple_data_writers_controller    data_writer_ctrl;

  STATUS                              status;
  // Ordered map of the input streams, based on the order (sequence
  // number of the data blocks)
  Input_stream_order_map           input_streams_order;
  // One input stream for every correlate node
  std::vector<Input_stream *>         input_streams;

  int32_t curr_slice, number_of_time_slices, curr_stream, curr_slice_size;
  int32_t curr_band;
  bool finalize_integration;
  int32_t total_bytes_written, number_of_bins;
  int32_t current_output_file;
  // In one read operation we might only partially receive the output file nr, 
  // this tracks the number of received bytes.
  int output_file_index;

  uint32_t number_channels;
};

#endif // OUTPUT_NODE_H
