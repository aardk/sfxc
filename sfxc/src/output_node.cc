/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */
#include "mpi_transfer.h"
#include "output_header.h"
#include "output_node.h"
#include "utils.h"

#include <iostream>

Output_node::Output_node(int rank, int size)
    : Node(rank),
    output_node_ctrl(*this),
    data_readers_ctrl(*this),
    data_writer_ctrl(*this),
    status(STOPPED), n_data_writers(0), output_file_index(0),
      curr_slice(0), number_of_time_slices(-1), curr_band(-1), curr_stream(-1), current_output_file(-1) {
  initialise();
}

Output_node::Output_node(int rank, Log_writer *writer, int size)
    : Node(rank, writer),
    output_node_ctrl(*this),
    data_readers_ctrl(*this),
    data_writer_ctrl(*this),
    status(STOPPED), n_data_writers(0), output_file_index(0),
      curr_slice(0), number_of_time_slices(-1), curr_band(-1), curr_stream(-1), current_output_file(-1) {
  initialise();
}

void Output_node::initialise() {
  add_controller(&data_readers_ctrl);
  add_controller(&output_node_ctrl);
  add_controller(&data_writer_ctrl);

  /// initialize and retreive the listening addresses/port
  std::vector<uint64_t> addrs;
  data_readers_ctrl.get_listening_ip( addrs );

  /// send this to the manager NODE
  MPI_Transfer::send_ip_address(addrs, RANK_MANAGER_NODE);


  int32_t msg=0;
  MPI_Send(&msg, 1, MPI_INT32,
           RANK_MANAGER_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD);
}

Output_node::~Output_node() {
  if (!get_assertion_raised()) {
    DEBUG_MSG("Output node finished");

    // empty the input buffers to the output
    SFXC_ASSERT(status == END_NODE);
    SFXC_ASSERT(input_streams_order.empty());
  }
}

void Output_node::terminate() {
  DEBUG_MSG("Output node terminate");
  status = END_NODE;
}

void Output_node::start() {
  while (status != END_NODE) {
    switch (status) {

    case STOPPED: {
        SFXC_ASSERT(curr_stream == -1);
        // blocking:
        check_and_process_message();

        if (curr_slice == number_of_time_slices) {
          status = END_NODE;
        } else if (!input_streams_order.empty()) {
          if (input_streams_order.begin()->first == curr_slice) {
            status = START_NEW_SLICE;
          }
        }
        break;
      }
    case START_NEW_SLICE: {
        SFXC_ASSERT(curr_stream == -1);

        SFXC_ASSERT(!input_streams_order.empty());
        SFXC_ASSERT(input_streams_order.begin()->first == curr_slice);
        curr_stream = input_streams_order.begin()->second.stream;
        curr_band = input_streams_order.begin()->second.band;
        finalize_integration = !input_streams_order.begin()->second.accum;
        SFXC_ASSERT(curr_stream >= 0);
        input_streams_order.erase(input_streams_order.begin());
        input_streams[curr_stream]->goto_next_slice(curr_slice_size, number_of_bins);
        total_bytes_written = 0;
        if (input_buffer.size() < (number_of_bins * curr_slice_size))
          input_buffer.resize(number_of_bins * curr_slice_size);
	if (curr_band >= accum_buffer.size()) {
	  accum_buffer.resize(curr_band + 1);
	  integration.resize(curr_band + 1, -1);
	}
	if (accum_buffer[curr_band].size() < (number_of_bins * curr_slice_size))
	  accum_buffer[curr_band].resize(number_of_bins * curr_slice_size);
	total_bytes_read = 0;
        status = READ_INPUT;
        break;
      }
    case READ_INPUT: {
        process_all_waiting_messages();

        int bytes_read = input_streams[curr_stream]->read_bytes(input_buffer);
        if (bytes_read <= 0) {
          usleep(100);
	} else {
	  total_bytes_read += bytes_read;
	}

        // Check whether we arrived at the end of the slice
        if (input_streams[curr_stream]->end_of_slice()) {
	  SFXC_ASSERT(total_bytes_read == number_of_bins * curr_slice_size);
          status = ACCUMULATE_INPUT;
        }
	break;
      }
    case ACCUMULATE_INPUT: {
        Output_header_timeslice *timeslice =
	  (Output_header_timeslice *)&input_buffer[4];
	size_t stats_offset = 4 + sizeof(Output_header_timeslice) +
	  timeslice->number_uvw_coordinates * sizeof(Output_uvw_coordinates);
	size_t data_offset = stats_offset + 
	  timeslice->number_statistics * sizeof(Output_header_bitstatistics);

	if (integration[curr_band] != timeslice->integration_slice) {
	  integration[curr_band] = timeslice->integration_slice;

	  // Initialize metadata
	  accum_buffer[curr_band] = input_buffer;

	  // Initialize visibilities if have more than one integration
	  // slice per integeration
	  if (!finalize_integration) {
	    for (int i = 0; i < timeslice->number_baselines; i++) {
	      Output_header_baseline *baseline =
		(Output_header_baseline *)&input_buffer[data_offset];
	      data_offset += sizeof(Output_header_baseline);
	      float *input = (float *)&input_buffer[data_offset];
	      float *accum = (float *)&accum_buffer[curr_band][data_offset];

	      // Complex is simply a pair of doubles.
	      for (int j = 0; j < 2 * number_channels; j++)
		accum[j] = input[j] * baseline->weight;;

	      data_offset += 2 * number_channels * sizeof(float);
	    }
	  }
	} else {
	  // Accumulate statistics
	  for (int i = 0; i < timeslice->number_statistics; i++) {
	    Output_header_bitstatistics *input =
	      (Output_header_bitstatistics *)&input_buffer[stats_offset];
	    Output_header_bitstatistics *accum =
	      (Output_header_bitstatistics *)&accum_buffer[curr_band][stats_offset];

	    for (int j = 0; j < sizeof(input->levels) / sizeof(input->levels[0]); j++)
	      accum->levels[j] += input->levels[j];
	    accum->n_invalid += input->n_invalid;
	    stats_offset += sizeof(Output_header_bitstatistics);
	  }

	  // Accumulate visibilities
	  for (int i = 0; i < timeslice->number_baselines; i++) {
	    Output_header_baseline *input_baseline =
	      (Output_header_baseline *)&input_buffer[data_offset];
	    Output_header_baseline *accum_baseline =
	      (Output_header_baseline *)&accum_buffer[curr_band][data_offset];
	    data_offset += sizeof(Output_header_baseline);
	    float *input = (float *)&input_buffer[data_offset];
	    float *accum = (float *)&accum_buffer[curr_band][data_offset];

	    // Accumulate visibility weights
	    accum_baseline->weight += input_baseline->weight;

	    // Complex is simply a pair of doubles.
	    for (int j = 0; j < 2 * number_channels; j++) {
	      accum[j] += input[j] * input_baseline->weight;
	      if (finalize_integration && accum_baseline->weight != 0)
		accum[j] /= accum_baseline->weight;
	    }

	    data_offset += 2 * number_channels * sizeof(float);
	  }
	}

	if (finalize_integration)
	  status = WRITE_OUTPUT;
	else
	  status = END_SLICE;
	break;
      }
    case WRITE_OUTPUT: {
        write_output(number_of_bins * curr_slice_size);
        total_bytes_written += number_of_bins * curr_slice_size;
        status = END_SLICE;
        break;
      }
    case END_SLICE: {
        curr_stream = -1;
        curr_slice ++;
        if (curr_slice == number_of_time_slices) {
          status = END_NODE;
        } else if (input_streams_order.empty()) {
          status = STOPPED;
        } else if (input_streams_order.begin()->first != curr_slice) {
          status = STOPPED;
        } else {
          status = START_NEW_SLICE;
        }
        break;
      }
    case END_NODE: {
        // For completeness sake
        break;
      }
    }
  }

  DEBUG_MSG("Shutting down !");
  data_readers_ctrl.stop();
  ///DEBUG_MSG("WANT TO SHUT DOWN THE WRITER !");

  // End the node;
  int32_t msg=0;
  MPI_Send(&msg, 1, MPI_INT32,
           RANK_MANAGER_NODE, MPI_TAG_OUTPUT_NODE_FINISHED, MPI_COMM_WORLD);
}

void
Output_node::
write_global_header(const Output_header_global &global_header) {
  int nbytes = global_header.header_size;
  for(int i=0;i<n_data_writers;i++)
    data_writer_ctrl.get_data_writer(i)->put_bytes(nbytes, (char *)&global_header);

  number_channels = (global_header.number_channels + 1);
}

void
Output_node::
set_order_of_input_stream(int stream, int order, int band, int accum, size_t size,
			  int nbins) {
  SFXC_ASSERT(stream >= 0);

  SFXC_ASSERT(stream < (int)input_streams.size());
  // Check that the weight does not exist yet:
  SFXC_ASSERT(input_streams_order.find(order) == input_streams_order.end());
  // Check that the ordering is right (not before the current element):
  if (!input_streams_order.empty()) {
    SFXC_ASSERT(order >= curr_slice);
  }

  // Add the stream to the queue:
  input_streams_order.insert(Input_stream_order_map_value(order, Stream_param(stream, band, accum)));
  input_streams[stream]->set_length_time_slice(size, nbins);

  SFXC_ASSERT(status != END_NODE);
}

bool Output_node::write_output(int nBytes) {
  if (nBytes <= 0)
    return false;

  SFXC_ASSERT(curr_stream >= 0);
  SFXC_ASSERT(input_streams[curr_stream] != NULL);

  int nbytes_per_file = curr_slice_size;
  int bytes_written=0;
  int index_in_file = total_bytes_written%nbytes_per_file;
  while(bytes_written<nBytes){
    // Get index of next output file
    if(output_file_index < 4){
      char *current_output_file_ptr = (char *)&current_output_file;
      int to_read = std::min(4 - output_file_index, nBytes - bytes_written);
      memcpy(&current_output_file_ptr[output_file_index], &accum_buffer[curr_band][bytes_written], to_read);
      output_file_index += to_read;
      bytes_written += to_read;
      index_in_file += to_read;
    }
    // Write the data
    int to_write = std::min(nbytes_per_file-index_in_file, nBytes-bytes_written);
//    std::cout << "current_output_file = " << current_output_file <<"\n";
    data_writer_ctrl.get_data_writer(current_output_file)->put_bytes(to_write, &accum_buffer[curr_band][bytes_written]);
    bytes_written += to_write;
    index_in_file += to_write;
    if(index_in_file >= nbytes_per_file){
      output_file_index = 0;
      index_in_file=0;
    }
  }
  SFXC_ASSERT(bytes_written==nBytes);
  return true;
}

void Output_node::hook_added_data_reader(size_t reader) {
  // Create an output buffer:
  data_readers_ctrl.enable_buffering(reader);

  // Create the data_stream:
  if (input_streams.size() <= reader) {
    input_streams.resize(reader+1, NULL);
  }

  input_streams[reader] =
    new Input_stream(data_readers_ctrl.get_data_reader(reader));
}

void Output_node::hook_added_data_writer(size_t writer) {
  n_data_writers++;
}

void
Output_node::set_number_of_time_slices(int n_time_slices) {
  number_of_time_slices = n_time_slices;
}


/*
 *  Input_stream
 */

Output_node::Input_stream::Input_stream(shared_ptr<Data_reader> reader)
    : reader(reader) {
  reader->set_size_dataslice(0);
}

int
Output_node::Input_stream::read_bytes(std::vector<char> &buffer) {
  SFXC_ASSERT(reader != shared_ptr<Data_reader>());
  size_t nBytes = std::min(buffer.size(),
                           (size_t)reader->get_size_dataslice());
  nBytes = reader->get_bytes(nBytes, &buffer[offset]);
  if (nBytes > 0)
    offset += nBytes;
  return nBytes;
}


bool
Output_node::Input_stream::end_of_slice() {
  return reader->end_of_dataslice();
}

void
Output_node::Input_stream::set_length_time_slice(int64_t nBytes, int nBins) {
  Slice new_slice={nBytes,nBins};
  slice_size.push(new_slice);
}

void
Output_node::Input_stream::goto_next_slice(int &new_slice_size, int &new_nbins) {
  SFXC_ASSERT(!slice_size.empty());
  SFXC_ASSERT(slice_size.front().nBytes > 0);
  SFXC_ASSERT(slice_size.front().nBins > 0);
  SFXC_ASSERT(reader->end_of_dataslice());
  new_slice_size = slice_size.front().nBytes;
  new_nbins = slice_size.front().nBins;
  reader->set_size_dataslice(new_nbins * new_slice_size);
  offset = 0;

  SFXC_ASSERT(reader->get_size_dataslice() > 0);
  slice_size.pop();
}

void Output_node::get_state(std::ostream &out) {
   out << "{\n"
      << "\t\"rank\": " << RANK_OF_NODE << ",\n"
      << "\t\"host\": \"" << HOSTNAME_OF_NODE << "\",\n"
      << "\t\"id\": \"" << ID_OF_NODE << "\",\n"
      << "\t\"now\": \"" << Time::now() << "\",\n"
      << "\t\"state\": ";
  switch (status) {
    case STOPPED:
     out << "\"STOPPED\",\n";
     break;
    case READ_INPUT:
     out << "\"READ_INPUT\",\n";
     break;
    case ACCUMULATE_INPUT:
     out << "\"ACCUMULATE_INPUT\",\n";
     break;
    case WRITE_OUTPUT:
     out << "\"WRITE_OUTPUT\",\n";
     break;
    case END_SLICE:
     out << "\"END_SLICE\",\n";
     break;
    case END_NODE:
     out << "\"END_NODE\",\n";
     break;
  }
  out << "\t\"curr_stream\": " << curr_stream << ",\n"
      << "\t\"curr_band\": " << curr_band << ",\n"
      << "\t\"curr_slice\": " << curr_slice << ",\n"
      << "\t\"slice_size\": " << curr_slice_size << ",\n"
      << "\t\"number_of_bins\": " << number_of_bins << ",\n"
      << "\t\"nr_time_slices\": " << number_of_time_slices << ",\n";
  data_readers_ctrl.get_state(out);
  out << "}";
}
