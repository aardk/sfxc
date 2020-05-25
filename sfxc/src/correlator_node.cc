/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *            Aard Keimpema <Keimpema@JIVE.nl>, 2008
 *
 * $Id$
 *
 */

#include "correlator_node.h"
#include "utils.h"
#ifdef USE_IPP
#include <ippcore.h>
#endif
Correlator_node::Correlator_node(int rank, int nr_corr_node, bool pulsar_binning_, bool phased_array_)
    : Node(rank),
    correlator_node_ctrl(*this),
    data_readers_ctrl(*this),
    data_writer_ctrl(*this),
    status(CORRELATING),
    pulsar_binning(pulsar_binning_),
    pulsar_parameters(get_log_writer()),
    tasklet(nr_corr_node, pulsar_binning_, phased_array_) {
  #ifdef USE_IPP
  ippSetNumThreads(1);
  #endif
  get_log_writer()(1) << "Correlator_node(" << nr_corr_node << ")" << std::endl;
  add_controller(&correlator_node_ctrl);
  add_controller(&data_readers_ctrl);
  add_controller(&data_writer_ctrl);


  int32_t msg;
  MPI_Send(&msg, 1, MPI_INT32,
           RANK_MANAGER_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD);

  start_threads();
}

Correlator_node::~Correlator_node() {
}

void Correlator_node::start_threads() {
  tasklet.start();
}

void Correlator_node::stop_threads() {
  tasklet.terminate();
  wait(tasklet);
}

void Correlator_node::start() {
  /// We enter the main loop of the coorelator node.
  //DEBUG_MSG("START MAIN LOOp !");
  main_loop();
}

void Correlator_node::terminate() {
  DEBUG_MSG("Correlator node received terminate signal.");
  status = END_NODE;
}

void Correlator_node::main_loop() {
  while ( status != END_NODE ) {
    if (check_and_process_waiting_message() == NO_MESSAGE)
      usleep(100);
  }
  stop_threads();
}

void Correlator_node::add_delay_table(Delay_table &table, int sn1, int sn2) {
  tasklet.add_delay_table(table, sn1, sn2);
}

void Correlator_node::add_uvw_table(Uvw_model &table, int sn) {
  tasklet.add_uvw_table(table, sn);
}

void Correlator_node::hook_added_data_reader(size_t stream_nr) {
  tasklet.hook_added_data_reader(stream_nr, data_readers_ctrl.get_data_reader(stream_nr));
}


void Correlator_node::hook_added_data_writer(size_t i) {
  SFXC_ASSERT(i == 0);
  tasklet.hook_added_data_writer(i, data_writer_ctrl.get_data_writer(0));
}

void
Correlator_node::add_source_list(const std::map<std::string, int> &sources) {
  tasklet.add_source_list(sources);
}

void
Correlator_node::receive_parameters(const Correlation_parameters &parameters) {
  tasklet.add_new_slice(parameters);
}

void Correlator_node::get_state(std::ostream &out) {
  out << "{\n"
      << "\t\"rank\": " << RANK_OF_NODE << ",\n"
      << "\t\"host\": \"" << HOSTNAME_OF_NODE << "\",\n"
      << "\t\"id\": \"" << ID_OF_NODE << "\",\n"
      << "\t\"now\": \"" << Time::now() << "\",\n";
  tasklet.get_state(out);
  out << "}";
}
