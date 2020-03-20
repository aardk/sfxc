/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#include "correlator_node.h"
#include "data_reader_file.h"
#include "data_reader_tcp.h"

#include "data_writer_file.h"
#include "data_writer_tcp.h"
#include "delay_table_akima.h"
#include "uvw_model.h"

#include "utils.h"

#include "mpi_transfer.h"

#include <signal.h>
#include <iostream>


Correlator_node_controller::Correlator_node_controller(Correlator_node &node)
    : Controller(node), node(node) {}

Correlator_node_controller::~Correlator_node_controller() {
}

Controller::Process_event_status
Correlator_node_controller::process_event(MPI_Status &status) {
  switch (status.MPI_TAG) {
  case MPI_TAG_DELAY_TABLE: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      Delay_table table;
      int sn[2];
      MPI_Transfer::receive_bcast(status, table, sn);
      if(sn[1] >= 0)
        node.add_delay_table(table, sn[0], sn[1]);
      else
        node.add_delay_table(table, sn[0], sn[0]);
      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_BP_TABLE: 
      // Passthrough
  case MPI_TAG_CL_TABLE: {
      MPI_Status status2;
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      int size;
      MPI_Get_elements(&status, MPI_CHAR, &size);
      char table[size];
      MPI_Recv(&table[0], size, MPI_CHAR, status.MPI_SOURCE,
               status.MPI_TAG, MPI_COMM_WORLD, &status2);
      if (status.MPI_TAG == MPI_TAG_BP_TABLE){
        node.correlation_core_normal->add_bp_table(table);
        if(node.correlator_node_type == CORRELATOR_NODE_PULSAR_BINNING)
          node.correlation_core_pulsar->add_bp_table(table);
      }else{
        node.correlation_core_normal->add_cl_table(table);
        if(node.correlator_node_type == CORRELATOR_NODE_PULSAR_BINNING)
          node.correlation_core_pulsar->add_cl_table(table);
      }
      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_UVW_TABLE: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      Uvw_model table;
      int sn;
      MPI_Transfer::receive_bcast(status, table, sn);
      node.add_uvw_table(table, sn);
      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_CORR_PARAMETERS: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      Correlation_parameters parameters;
      MPI_Transfer::receive(status, parameters);
      parameters.mask_parameters = &node.mask_parameters;
      node.receive_parameters(parameters);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_PULSAR_PARAMETERS: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      MPI_Transfer::receive(status, node.pulsar_parameters);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_MASK_PARAMETERS: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      Mask_parameters parameters;
      MPI_Transfer::receive_bcast(status, node.mask_parameters);
      return PROCESS_EVENT_STATUS_SUCCEEDED;
  }
  case MPI_TAG_SOURCE_LIST:{
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      std::map<std::string, int> sources;
      MPI_Transfer::receive(status, sources);
      node.correlation_core_normal->add_source_list(sources);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  }
  return PROCESS_EVENT_STATUS_UNKNOWN;
}
