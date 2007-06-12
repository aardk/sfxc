/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 * 
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 * 
 * $Id$
 *
 *  Tests reading a file from disk and then writing it back using a Data_node
 */

#include <types.h>
#include <Input_node.h>
#include <Output_node.h>
#include <Log_node.h>
#include <Log_writer_cout.h>

#include <fstream>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "constPrms.h"
#include "runPrms.h"
#include "genPrms.h"
#include "staPrms.h"
#include "genFunctions.h"
#include "InData.h"
#include "delayTable.h"
#include "MPI_Transfer.h"


#include <Node.h>
#include <Data_reader2buffer.h>
#include <TCP_Connection.h>
#include <Buffer2data_writer.h>
#include <Data_writer.h>
#include <Data_writer_file.h>
#include <Data_reader_tcp.h>
#include <utils.h>

int input_node = 2;
int output_node = 3;
char output_file_to_disk[] = "output_file_to_disk.dat";


RunP RunPrms;
GenP GenPrms;
StaP StaPrms[NstationsMax];

// MPI
int numtasks, rank;


void wait_for_setting_up_channel() {
  MPI_Status status;
  INT64 channel;
  MPI_Recv(&channel, 1, MPI_INT64, MPI_ANY_SOURCE,
           MPI_TAG_INPUT_CONNECTION_ESTABLISHED, MPI_COMM_WORLD, &status);
}

void test_writing_data_to_file(char *control_file,
                               char *output_file) {
  MPI_Status status;
  
  if (rank == RANK_MANAGER_NODE) {
    int msg=0;
    // Log node:
    MPI_Send(&msg, 1, MPI_INT32, 
             RANK_LOG_NODE, MPI_TAG_SET_LOG_NODE, MPI_COMM_WORLD);
    MPI_Send(&msg, 1, MPI_INT32, 
             RANK_LOG_NODE, MPI_TAG_LOG_NODE_SET_OUTPUT_COUT, MPI_COMM_WORLD);
    MPI_Recv(&msg, 1, MPI_INT32, 
             RANK_LOG_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD, &status);

    // Initialise control parameters
    Log_writer_mpi log_writer(rank);
    initialise_control(control_file, log_writer, RunPrms, GenPrms, StaPrms);


    // Input_node node:
    MPI_Send(&msg, 1, MPI_INT32, 
             input_node, MPI_TAG_SET_INPUT_NODE, MPI_COMM_WORLD);

    MPI_Transfer mpi_transfer;
    mpi_transfer.send_general_parameters(input_node);

    MPI_Recv(&msg, 1, MPI_INT32, 
             input_node, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD, &status);
        
    // Connect the input and output:
    // strlen+1 so that \0 gets transmitted as well
    MPI_Send(StaPrms[0].get_mk4file(), strlen(StaPrms[0].get_mk4file())+1, 
             MPI_CHAR, 
             input_node, MPI_TAG_SET_DATA_READER_FILE, MPI_COMM_WORLD);
    wait_for_setting_up_channel();
    
    { // send output file to input node
      int size = strlen(output_file)+1+sizeof(int), position=0;
      char buffer[size];
      MPI_Pack(&output_node, 1, MPI_INT32, 
               buffer, size, &position, MPI_COMM_WORLD);
      MPI_Pack(output_file, strlen(output_file)+1, MPI_CHAR, 
               buffer, size, &position, MPI_COMM_WORLD);
    
      MPI_Send(buffer, size, MPI_PACKED, 
               input_node, MPI_TAG_ADD_DATA_WRITER_FILE, MPI_COMM_WORLD);
      wait_for_setting_up_channel();
    }

    // set priorities:
    INT64 priority_in[] = {output_node,0,0};
    MPI_Send(&priority_in, 2, MPI_INT64, 
             input_node, MPI_TAG_INPUT_STREAM_SET_PRIORITY, MPI_COMM_WORLD);

    // Wait for data nodes to finish
    MPI_Status status;
    int i;
    MPI_Recv(&i, 1, MPI_INT32, MPI_ANY_SOURCE,
             MPI_TAG_DATASTREAM_EMPTY, MPI_COMM_WORLD, &status);
    assert(status.MPI_TAG == MPI_TAG_DATASTREAM_EMPTY);
    assert(status.MPI_SOURCE == input_node);
   
    MPI_Send(&rank, 1, MPI_INT, 
             RANK_LOG_NODE, MPI_TAG_LOG_MESSAGES_ENDED, MPI_COMM_WORLD);

  } else {
    // Don't use the output node:
    if (rank == output_node) {
      MPI_Send(&rank, 1, MPI_INT, 
               RANK_LOG_NODE, MPI_TAG_LOG_MESSAGES_ENDED, MPI_COMM_WORLD);
      return;
    }
    
    INT32 msg;
    MPI_Recv(&msg, 1, MPI_INT32, 
             RANK_MANAGER_NODE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    switch (status.MPI_TAG) {
    case MPI_TAG_SET_LOG_NODE: {
      assert (RANK_LOG_NODE == rank);
      Log_node log_node(rank,numtasks);
      log_node.start();
      break;
    }
    case MPI_TAG_SET_INPUT_NODE: {
      assert(rank == input_node);
      Input_node node(rank);
      node.start();
      break;
    }
    default: {
      assert(false);
      break;
    }
    }
  }
}

void test_writing_data_to_output_node_using_TCP(char *control_file,
                                                char *output_file) {
  MPI_Status status;
  
  if (rank == RANK_MANAGER_NODE) {
    int msg=0;
    // Log node:
    MPI_Send(&msg, 1, MPI_INT32, 
             RANK_LOG_NODE, MPI_TAG_SET_LOG_NODE, MPI_COMM_WORLD);
    MPI_Send(&msg, 1, MPI_INT32, 
             RANK_LOG_NODE, MPI_TAG_LOG_NODE_SET_OUTPUT_COUT, MPI_COMM_WORLD);
    MPI_Recv(&msg, 1, MPI_INT32, 
             RANK_LOG_NODE, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD, &status);

    // Initialise control parameters
    Log_writer_mpi log_writer(rank);
    initialise_control(control_file, log_writer, RunPrms, GenPrms, StaPrms);

    // Output_node node:
    MPI_Send(&msg, 1, MPI_INT32, 
             output_node, MPI_TAG_SET_OUTPUT_NODE, MPI_COMM_WORLD);
    MPI_Recv(&msg, 1, MPI_INT32, 
             output_node, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD, &status);
    MPI_Send(output_file, strlen(output_file)+1, MPI_CHAR, 
             output_node, MPI_TAG_SET_DATA_WRITER_FILE, MPI_COMM_WORLD);
    wait_for_setting_up_channel();

    // Input_node node:
    MPI_Send(&msg, 1, MPI_INT32, 
             input_node, MPI_TAG_SET_INPUT_NODE, MPI_COMM_WORLD);
    MPI_Transfer mpi_transfer;
    mpi_transfer.send_general_parameters(input_node);
    MPI_Recv(&msg, 1, MPI_INT32, 
             input_node, MPI_TAG_NODE_INITIALISED, MPI_COMM_WORLD, &status);
    MPI_Send(StaPrms[0].get_mk4file(), strlen(StaPrms[0].get_mk4file())+1, 
             MPI_CHAR, 
             input_node, MPI_TAG_SET_DATA_READER_FILE, MPI_COMM_WORLD);
    wait_for_setting_up_channel();

    // Connect input and output
    // first two numbers are arbitrary
    int writer_stream_nr = 5;
    int reader_stream_nr = 7;
    INT32 ranks[] = { writer_stream_nr, reader_stream_nr, output_node};
    MPI_Send(ranks, 3, MPI_INT32, 
             input_node,
             MPI_TAG_ADD_OUTPUT_CONNECTION_MULTIPLE_INPUT_TCP, MPI_COMM_WORLD);
    wait_for_setting_up_channel();
    
    // set priorities:
    {
      INT64 priority_in[] = {writer_stream_nr,0,0};
      MPI_Send(&priority_in, 3, MPI_INT64, 
               input_node, MPI_TAG_INPUT_STREAM_SET_PRIORITY, MPI_COMM_WORLD);
    }
    {
      INT64 priority[] = {reader_stream_nr, 0};
      MPI_Send(&priority, 2, MPI_INT64, 
               output_node, MPI_TAG_OUTPUT_STREAM_SET_PRIORITY, MPI_COMM_WORLD);
    }

    // Wait for data nodes to finish
    MPI_Status status;
    int i;
    MPI_Recv(&i, 1, MPI_INT32, MPI_ANY_SOURCE,
             MPI_TAG_DATASTREAM_EMPTY, MPI_COMM_WORLD, &status);
    assert(status.MPI_TAG == MPI_TAG_DATASTREAM_EMPTY);
    assert(status.MPI_SOURCE == input_node);
   
    int slicenr = 0;
    MPI_Send(&slicenr, 1, MPI_INT32, output_node, 
             MPI_TAG_OUTPUT_NODE_CORRELATION_READY, MPI_COMM_WORLD);
    MPI_Send(&slicenr, 1, MPI_INT32, output_node, 
             MPI_TAG_CORRELATION_READY, MPI_COMM_WORLD);

    MPI_Send(&rank, 1, MPI_INT, 
             RANK_LOG_NODE, MPI_TAG_LOG_MESSAGES_ENDED, MPI_COMM_WORLD);

  } else {
    INT32 msg;
    MPI_Recv(&msg, 1, MPI_INT32, 
             RANK_MANAGER_NODE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    switch (status.MPI_TAG) {
    case MPI_TAG_SET_LOG_NODE: {
      assert (RANK_LOG_NODE == rank);
      Log_node log_node(rank,numtasks);
      log_node.start();
      break;
    }
    case MPI_TAG_SET_INPUT_NODE: {
      assert(rank == input_node);
      Input_node node(rank);
      node.start();
      break;
    }
    case MPI_TAG_SET_OUTPUT_NODE: {
      assert(rank == output_node);
      Output_node node(rank);
      node.start();
      break;
    }
    default: {
      assert(false);
      break;
    }
    }
  }
}
int main(int argc, char *argv[]) {
  assert(input_node+output_node+RANK_MANAGER_NODE+RANK_LOG_NODE == 6);

  char *control_file, *output_directory;



  //initialisation
  int stat = MPI_Init(&argc,&argv);
  if (stat != MPI_SUCCESS) {
    std::cout << "Error starting MPI program. Terminating.\n";
    MPI_Abort(MPI_COMM_WORLD, stat);
  }

  // get the number of tasks set at commandline (= number of processors)
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  // get the ID (rank) of the task, fist rank=0, second rank=1 etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  std::cout << "#" << rank << " pid = " << getpid() << std::endl;

  // 0: controller
  // 1: log node
  // 2: input node
  // 3: output node
  assert(numtasks == 4);
  
  ///////////////////////////
  //  The real work
  ///////////////////////////
  assert(argc == 3);
  control_file = argv[1];
  output_directory = argv[2];

  // First test
  char output_file1[strlen(output_directory)+strlen(output_file_to_disk)+2];
  char output_file2[strlen(output_directory)+strlen(output_file_to_disk)+2];
  sprintf(output_file1, "%s/%s.1", output_directory, output_file_to_disk);
  sprintf(output_file2, "%s/%s.2", output_directory, output_file_to_disk);

  test_writing_data_to_file(control_file, output_file1);
  test_writing_data_to_file(control_file, output_file2);
  
  
  if (rank == input_node) {
    sync();
    std::stringstream cmd; 
    cmd << "cmp " << output_file1 << " " << output_file2;
    if (system(cmd.str().c_str())) {
      std::cout << "cmp file output failed" << std::endl;
      return 1;
    }
  }

  /////////////////////
  MPI_Barrier( MPI_COMM_WORLD );
  /////////////////////

  // Second test
  test_writing_data_to_output_node_using_TCP(control_file, output_file1);
  
  {
    sync();
    sleep(1);
    std::stringstream cmd; cmd << "cmp " << output_file1 << " " << output_file2;
    if (system(cmd.str().c_str())) {
      std::cout << "cmp tcp output failed" << std::endl;
      return 1;
    }
  }

  //close the mpi stuff
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();

  return 0;
}
