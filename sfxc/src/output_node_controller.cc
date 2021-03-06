/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

// Output_node_controller is defined in Output_node.h:
#include "output_node.h"

#include "utils.h"
#include "data_writer_file.h"
#include "data_reader_tcp.h"
#include "tcp_connection.h"
#include "correlator_time.h"

#include <iostream>


Output_node_controller::Output_node_controller(Output_node &node)
  : Controller(node), node(node)
{
}

Controller::Process_event_status
Output_node_controller::process_event(MPI_Status &status) {
  MPI_Status status2;
  switch (status.MPI_TAG) {
  case MPI_TAG_OUTPUT_NODE_GLOBAL_HEADER: {
      int size;
      MPI_Get_elements(&status, MPI_CHAR, &size);
      SFXC_ASSERT(size >= sizeof(Output_header_global));
      Output_header_global *global_header;
      global_header = (Output_header_global *)malloc(size);
      SFXC_ASSERT(global_header != NULL);
      MPI_Recv((char *)global_header, size, MPI_CHAR, status.MPI_SOURCE,
               status.MPI_TAG, MPI_COMM_WORLD, &status2);
      SFXC_ASSERT(global_header->header_size == size);
      node.write_global_header(*global_header);

      struct Output_header_phasecal phasecal_header;
      phasecal_header.header_size = sizeof(phasecal_header);
      memcpy(&phasecal_header.experiment, global_header->experiment,
	     sizeof(phasecal_header.experiment));
      phasecal_header.output_format_version =
	global_header->output_format_version;
      phasecal_header.correlator_version = global_header->correlator_version;
      memcpy(&phasecal_header.correlator_branch,
	     global_header->correlator_branch,
	     sizeof(phasecal_header.correlator_branch));
      phasecal_header.job_nr = global_header->job_nr;
      phasecal_header.subjob_nr = global_header->subjob_nr;
      phasecal_file.write((char *)&phasecal_header, sizeof(phasecal_header));

      if (tsys_file.is_open()) {
	struct Output_header_tsys tsys_header;
	tsys_header.header_size = sizeof(tsys_header);
	memcpy(&tsys_header.experiment, global_header->experiment,
	       sizeof(tsys_header.experiment));
	tsys_header.output_format_version =
	  global_header->output_format_version;
	tsys_header.correlator_version = global_header->correlator_version;
	memcpy(&tsys_header.correlator_branch,
	       global_header->correlator_branch,
	       sizeof(tsys_header.correlator_branch));
	tsys_header.job_nr = global_header->job_nr;
	tsys_header.subjob_nr = global_header->subjob_nr;
	tsys_file.write((char *)&tsys_header, sizeof(tsys_header));
      }

      free(global_header);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_OUTPUT_STREAM_SLICE_SET_ORDER: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      int32_t param[6]; // stream, order, band, accum, size (in bytes), n_bins
      MPI_Recv(&param, 6, MPI_INT32, status.MPI_SOURCE,
               status.MPI_TAG, MPI_COMM_WORLD, &status2);

      SFXC_ASSERT(status.MPI_SOURCE == status2.MPI_SOURCE);
      SFXC_ASSERT(status.MPI_TAG == status2.MPI_TAG);

      // Create an output buffer:
      node.set_order_of_input_stream(param[0], param[1], param[2], param[3],
				     param[4], param[5]);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_OUTPUT_NODE_CORRELATION_READY: {
      get_log_writer()(3) << print_MPI_TAG(status.MPI_TAG) << std::endl;
      int32_t nr_of_time_slices;
      MPI_Recv(&nr_of_time_slices, 1, MPI_INT32, status.MPI_SOURCE,
               status.MPI_TAG, MPI_COMM_WORLD, &status2);
      SFXC_ASSERT(nr_of_time_slices >= 0);
      node.set_number_of_time_slices(nr_of_time_slices);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_OUTPUT_NODE_SET_PHASECAL_FILE: {
      int len;
      MPI_Get_elements(&status, MPI_CHAR, &len);
      SFXC_ASSERT(len > 0);

      char filename[len];
      MPI_Recv(&filename, len, MPI_CHAR, status.MPI_SOURCE,
	       status.MPI_TAG, MPI_COMM_WORLD, &status2);
      SFXC_ASSERT(filename[len - 1] == 0);
      SFXC_ASSERT(strncmp(filename, "file://", 7) == 0);
      phasecal_file.open(filename + 7, std::ios::out | std::ios::trunc | std::ios::binary);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_OUTPUT_NODE_WRITE_PHASECAL: {
      int len;
      MPI_Get_elements(&status, MPI_CHAR, &len);

      char msg[len];
      MPI_Recv(&msg, len, MPI_CHAR, status.MPI_SOURCE,
	       status.MPI_TAG, MPI_COMM_WORLD, &status2);

      int pos = 0;
      uint8_t station_number, frequency_number, sideband, polarisation;
      MPI_Unpack(msg, len, &pos, &station_number, 1, MPI_UINT8, MPI_COMM_WORLD);
      MPI_Unpack(msg, len, &pos, &frequency_number, 1, MPI_UINT8, MPI_COMM_WORLD);
      MPI_Unpack(msg, len, &pos, &sideband, 1, MPI_UINT8, MPI_COMM_WORLD);
      MPI_Unpack(msg, len, &pos, &polarisation, 1, MPI_UINT8, MPI_COMM_WORLD);

      Time start_time;
      uint64_t start_time_ticks;
      MPI_Unpack(msg, len, &pos, &start_time_ticks, 1, MPI_INT64, MPI_COMM_WORLD);
      start_time.set_clock_ticks(start_time_ticks);
      uint32_t mjd = start_time.get_mjd();
      uint32_t secs = start_time.get_time();

      Time integration_time;
      uint64_t integration_time_ticks;
      MPI_Unpack(msg, len, &pos, &integration_time_ticks, 1, MPI_INT64, MPI_COMM_WORLD);
      integration_time.set_clock_ticks(integration_time_ticks);
      uint32_t integration_time_secs = integration_time.get_time();

      uint32_t num_samples;
      MPI_Unpack(msg, len, &pos, &num_samples, 1, MPI_INT32, MPI_COMM_WORLD);
      std::vector<int32_t> samples;
      samples.resize(num_samples);
      MPI_Unpack(msg, len, &pos, &samples[0], num_samples, MPI_INT32, MPI_COMM_WORLD);

      SFXC_ASSERT(phasecal_file.is_open());
      phasecal_file.write((char *)&station_number, sizeof(station_number));
      phasecal_file.write((char *)&frequency_number, sizeof(frequency_number));
      phasecal_file.write((char *)&sideband, sizeof(sideband));
      phasecal_file.write((char *)&polarisation, sizeof(polarisation));
      phasecal_file.write((char *)&mjd, sizeof(mjd));
      phasecal_file.write((char *)&secs, sizeof(secs));
      phasecal_file.write((char *)&integration_time_secs, sizeof(integration_time_secs));
      phasecal_file.write((char *)&num_samples, sizeof(num_samples));
      phasecal_file.write((char *)&samples[0], num_samples * sizeof(samples[0]));

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_OUTPUT_NODE_SET_TSYS_FILE: {
      int len;
      MPI_Get_elements(&status, MPI_CHAR, &len);
      SFXC_ASSERT(len > 0);

      char filename[len];
      MPI_Recv(&filename, len, MPI_CHAR, status.MPI_SOURCE,
	       status.MPI_TAG, MPI_COMM_WORLD, &status2);
      SFXC_ASSERT(filename[len - 1] == 0);
      SFXC_ASSERT(strncmp(filename, "file://", 7) == 0);
      tsys_file.open(filename + 7, std::ios::out | std::ios::trunc | std::ios::binary);

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  case MPI_TAG_OUTPUT_NODE_WRITE_TSYS: {
      int len;
      MPI_Get_elements(&status, MPI_CHAR, &len);

      char msg[len];
      MPI_Recv(&msg, len, MPI_CHAR, status.MPI_SOURCE,
	       status.MPI_TAG, MPI_COMM_WORLD, &status2);

      int pos = 0;
      uint8_t station_number, frequency_number, sideband, polarisation;
      MPI_Unpack(msg, len, &pos, &station_number, 1, MPI_UINT8, MPI_COMM_WORLD);
      MPI_Unpack(msg, len, &pos, &frequency_number, 1, MPI_UINT8, MPI_COMM_WORLD);
      MPI_Unpack(msg, len, &pos, &sideband, 1, MPI_UINT8, MPI_COMM_WORLD);
      MPI_Unpack(msg, len, &pos, &polarisation, 1, MPI_UINT8, MPI_COMM_WORLD);

      Time start_time;
      uint64_t start_time_ticks;
      MPI_Unpack(msg, len, &pos, &start_time_ticks, 1, MPI_INT64, MPI_COMM_WORLD);
      start_time.set_clock_ticks(start_time_ticks);
      uint32_t mjd = start_time.get_mjd();
      uint32_t secs = start_time.get_time();

      uint64_t tsys[4];
      MPI_Unpack(msg, len, &pos, &tsys[0], 4, MPI_INT64, MPI_COMM_WORLD);

      if (tsys_file.is_open()) {
	tsys_file.write((char *)&station_number, sizeof(station_number));
	tsys_file.write((char *)&frequency_number, sizeof(frequency_number));
	tsys_file.write((char *)&sideband, sizeof(sideband));
	tsys_file.write((char *)&polarisation, sizeof(polarisation));
	tsys_file.write((char *)&mjd, sizeof(mjd));
	tsys_file.write((char *)&secs, sizeof(secs));
	tsys_file.write((char *)&tsys[0], sizeof(tsys));
      }

      return PROCESS_EVENT_STATUS_SUCCEEDED;
    }
  }
  return PROCESS_EVENT_STATUS_UNKNOWN;
}
