/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#include "log_node.h"
#include "correlator_time.h"
#include "utils.h"
#include "types.h"
#include "log_writer.h"
#include "log_writer_cout.h"
#include "log_writer_file.h"

#include <iostream>

Log_node::Log_node(int rank, int nNodes)
    : Node(rank), log_node_ctrl(*this, nNodes), running(true) {
  get_log_writer()(1) << "Log_node()" << std::endl;
  add_controller(&log_node_ctrl);

  int32_t msg;
  MPI_Send(&msg, 1, MPI_INT32,
           RANK_MANAGER_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD);
}

Log_node::Log_node(int rank, int nNodes, Log_writer *writer)
    : Node(rank, writer), log_node_ctrl(*this, nNodes) {
  get_log_writer()(1) << "Log_node()" << std::endl;
  add_controller(&log_node_ctrl);

  int32_t msg;
  MPI_Send(&msg, 1, MPI_INT32,
           RANK_MANAGER_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD);
 }

Log_node::~Log_node() {
  check_and_process_waiting_message();
}

void Log_node::start() {
  while ((!log_node_ctrl.ready()) && running) {
    check_and_process_message();
  }
  //DEBUG_MSG("Log node is terminate.");
}

void Log_node::terminate() {
  running=false;
  DEBUG_MSG("Log node terminate.");
}

void Log_node::hook_added_data_reader(size_t reader) {}
void Log_node::hook_added_data_writer(size_t writer) {}

void Log_node::get_state(std::ostream &out) {
  out << "{\n"
      << "\t\"rank\": " << RANK_OF_NODE << ",\n"
      << "\t\"host\": \"" << HOSTNAME_OF_NODE << "\",\n"
      << "\t\"id\": \"" << ID_OF_NODE << "\",\n"
      << "\t\"now\": \"" << Time::now() << "\"\n"
      << "}";
}
