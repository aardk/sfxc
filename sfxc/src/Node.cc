/*
CVS keywords
$Author$
$Date$
$Name$
$Revision$
$Source$

Author     : NGH Kruithof
StartDate  : 20061101
Last change: 20061124
*/

#include <Node.h>
#include <iostream>
#include <assert.h>

Node::Node(int rank) : rank(rank), log_writer(0) {
  log_writer.set_mpilevel(10);
  log_writer.set_rank(rank);
}

void 
Node::add_controller(Controller *controller) {
  controllers.push_back(controller);
}

void Node::start() {
  while (true) {
    MESSAGE_RESULT result = check_and_process_messages();
    
    // NGHK: Make an enum for the result types:
    switch (result) {
      case MESSAGE_PROCESSED: {
        break;
      }
      case TERMINATE_NODE: {
        return;
      }
      case NO_MESSAGE: {
        assert(false);
        return;
      }
      case ERROR_IN_PROCESSING: {
        log_writer.error("Error, failed to process message");
        break;
      }
      case MESSAGE_UNKNOWN: {
        break;
      }
    }
  }
}

Node::MESSAGE_RESULT Node::check_and_process_waiting_messages() {
    MPI_Status status;
    int result;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &result, &status);
    if (result) {
      return process_event(status);
    }
    return NO_MESSAGE;
}

Node::MESSAGE_RESULT Node::check_and_process_messages() {
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    return process_event(status);
}

Node::MESSAGE_RESULT Node::process_event(MPI_Status &status) {
  if (status.MPI_TAG == MPI_TAG_CORRELATION_READY) {
    MPI_Status status2; int msg;
    MPI_Recv(&msg, 1, MPI_INT, status.MPI_SOURCE,
             status.MPI_TAG, MPI_COMM_WORLD, &status2);
    
    log_writer.MPI(0, "MPI_TAG_CORRELATION_READY <===");
    return TERMINATE_NODE;
  }
  for (Controller_iterator it = controllers.begin();
       it != controllers.end();
       it++) {
    Controller::Process_event_status result = (*it)->process_event(status);
    switch (result) {
    case Controller::PROCESS_EVENT_STATUS_SUCCEEDED: // Processing succeeded
      {
        return MESSAGE_PROCESSED;
      }
    case Controller::PROCESS_EVENT_STATUS_UNKNOWN: // Unknown command, try next controller
      {
        continue;
        break;
      }
    case Controller::PROCESS_EVENT_STATUS_FAILED: // Processing failed
      {
        log_writer(0) << "Error in processing tag:" << status.MPI_TAG << std::endl;
        return ERROR_IN_PROCESSING;
      }
    }
  }
  {
    char msg[80];
    snprintf(msg, 80, "Unknown event %d", status.MPI_TAG);
    log_writer.error(msg);
  }
  
  // Remove event:  
  int size;
  MPI_Get_elements(&status, MPI_CHAR, &size);
  assert(size >= 0);
  char msg[size];
  MPI_Status status2;
  MPI_Recv(&msg, size, MPI_CHAR, status.MPI_SOURCE,
           status.MPI_TAG, MPI_COMM_WORLD, &status2);
  
  return MESSAGE_UNKNOWN;
}

