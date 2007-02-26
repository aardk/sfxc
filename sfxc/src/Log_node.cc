/* Author(s): Nico Kruithof, 2007
 * 
 * $Id$
 */

#include <Log_node.h>

#include <types.h>
#include <Log_writer.h>
#include <Log_writer_cout.h>
#include <Log_writer_file.h>

#include <iostream>
#include <assert.h>

Log_node::Log_node(int rank, int nNodes) 
  : Node(rank), log_node_ctrl(*this, nNodes)
{
  add_controller(&log_node_ctrl);

  get_log_writer().message(0, "Log_node: ready");

  INT32 msg;
  MPI_Send(&msg, 1, MPI_INT32, 
           RANK_MANAGER_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD);
}

void Log_node::start() {
  while (!log_node_ctrl.ready()) {
    check_and_process_message();
  }
}

void Log_node::hook_added_data_reader(int reader) {
}
void Log_node::hook_added_data_writer(int writer) {
}
