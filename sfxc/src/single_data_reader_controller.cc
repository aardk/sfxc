/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#include "single_data_reader_controller.h"
#include "sfxc_mpi.h"
#include "data_reader_factory.h"
#include "data_reader_tcp.h"

Single_data_reader_controller::
Single_data_reader_controller(Node &node)
  : Controller(node) {}

Single_data_reader_controller::Process_event_status
Single_data_reader_controller::process_event(MPI_Status &status) {
  MPI_Status status2;
  switch (status.MPI_TAG) {
  case MPI_TAG_ADD_DATA_READER: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;

      int size;
      MPI_Get_elements(&status, MPI_CHAR, &size);
      SFXC_ASSERT(size > sizeof(int32_t));
      char msg[size];
      char *p = msg;
      MPI_Recv(&msg, size, MPI_CHAR, status.MPI_SOURCE,
               status.MPI_TAG, MPI_COMM_WORLD, &status2);
      SFXC_ASSERT(status.MPI_SOURCE == status2.MPI_SOURCE);
      SFXC_ASSERT(status.MPI_TAG == status2.MPI_TAG);

      int32_t stream_nr;
      memcpy(&stream_nr, p, sizeof(int32_t));
      size -= sizeof(int32_t);
      p += sizeof(int32_t);

      // Make sure the array is null-terminated
      if (size > 0)
	p[size - 1] = 0;

      std::vector<std::string> sources;
      while (size > 0) {
	sources.push_back(p);
	size -= strlen(p) + 1;
	p += strlen(p) + 1;
      }
      SFXC_ASSERT(sources.size() > 0);

      shared_ptr<Data_reader> reader(Data_reader_factory::get_reader(sources));

      set_data_reader(stream_nr, reader);

      MPI_Send(&stream_nr, 1, MPI_INT32,
               status.MPI_SOURCE, MPI_TAG_CONNECTION_ESTABLISHED,
               MPI_COMM_WORLD);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }

  case MPI_TAG_ADD_DATA_READER_TCP2: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;

      int size;
      MPI_Get_elements(&status, MPI_INT64, &size);
      SFXC_ASSERT(size >= 3); // stream_nr, [ip-addr]+, port
      uint64_t ip_addr[size];
      MPI_Recv(&ip_addr, size, MPI_INT64, status.MPI_SOURCE,
               status.MPI_TAG, MPI_COMM_WORLD, &status2);

      SFXC_ASSERT(status.MPI_SOURCE == status2.MPI_SOURCE);
      SFXC_ASSERT(status.MPI_TAG == status2.MPI_TAG);

      int32_t stream_nr = ip_addr[0];
      uint64_t port = ip_addr[size-1];

      shared_ptr<Data_reader>
      reader(new Data_reader_tcp(ip_addr+1, size-2, port));
      set_data_reader(stream_nr, reader);

      MPI_Send(&stream_nr, 1, MPI_INT32,
               status.MPI_SOURCE, MPI_TAG_CONNECTION_ESTABLISHED,
               MPI_COMM_WORLD);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  }
  return PROCESS_EVENT_STATUS_UNKNOWN;
}

shared_ptr<Data_reader>
Single_data_reader_controller::get_data_reader(int reader) {
  SFXC_ASSERT(reader == 0);
  SFXC_ASSERT(reader_ != Data_reader_ptr() );

  return reader_;
}

void
Single_data_reader_controller::
set_data_reader(int stream_nr, Data_reader_ptr reader) {
  SFXC_ASSERT(stream_nr == 0);
  SFXC_ASSERT(reader_ == Data_reader_ptr());

  reader_ = reader;
  node.hook_added_data_reader(0);
}
