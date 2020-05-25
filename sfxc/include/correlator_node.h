/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#ifndef CORRELATOR_NODE_H
#define CORRELATOR_NODE_H
#include "node.h"
#include "multiple_data_readers_controller.h"
#include "single_data_writer_controller.h"
#include "control_parameters.h"
#include "correlator_node_tasklet.h"
#include "uvw_model.h"
#include "delay_table_akima.h"

// Declare the correlator controller:
class Correlator_node;

/**
 * Correlator_node_controller processes specific signals for the Correlator node.
 **/
class Correlator_node_controller : public Controller {
public:
  Correlator_node_controller(Correlator_node &node);
  ~Correlator_node_controller();

  Process_event_status process_event(MPI_Status &status);

private:
  Correlator_node &node;
};

/**
 * A correlate node will initialize the correlation process and connect
 * to the output node. It can receive messages from a data node asking to
 * open an input connection and from the controller node to process a
 * time slice. After the slice is processed the node will send a message
 * to the controller node saying it is available for a next job.
 *
 * \ingroup Node
 **/
class Correlator_node : public Node {
friend class Correlator_node_controller;
public:
  /// The states of the correlator_node.
  enum Status {
    CORRELATING = 0,
    END_NODE
  };

  Correlator_node(int rank, int nr_corr_node, bool psr_binning, bool phased_array);
  ~Correlator_node();

  /// The the main_loop of the correlator node.
  void start();

  /// Terminate the main_loop.
  void terminate();

  /// Callback function for adding a data_reader:
  void hook_added_data_reader(size_t reader);

  /// Callback function for adding a data_writer:
  void hook_added_data_writer(size_t writer);

  void get_state(std::ostream &out);

  void add_delay_table(Delay_table &table, int sn1, int sn2);
  void add_uvw_table(Uvw_model &table, int sn1);

  void output_node_set_timeslice(int slice_nr, int stream_nr, int band, int accum,
				 int bytes, int nbins);

  void receive_parameters(const Correlation_parameters &parameters);
  void add_source_list(const std::map<std::string, int> &sources);
  void set_parameters();
  

private:
  void start_threads();
  void stop_threads();

  /// Main loop just processes mpi messages 
  void main_loop();

  /// Everything related to the actual correlation is implemented in the tasklet
  Correlator_node_tasklet          tasklet;
  Correlator_node_controller       correlator_node_ctrl;

  /// The correlator node is connected to each of the input nodes.
  Multiple_data_readers_controller data_readers_ctrl;

  Single_data_writer_controller    data_writer_ctrl;

  // Contains all timing/binning parameters relating to any pulsar in the current experiment
  Pulsar_parameters pulsar_parameters; 
  Mask_parameters mask_parameters;
  
  bool pulsar_binning; // Set to true if pulsar binning is enabled

  /// State variables:
  Status status;
};

#endif // CORRELATOR_NODE_H
