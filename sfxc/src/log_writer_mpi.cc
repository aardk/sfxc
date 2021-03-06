/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#include "log_writer_mpi.h"
#include "sfxc_mpi.h"
#include "utils.h"

#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <cstring>

class Log_writer_mpi_buffer : public Log_writer_buffer {
public:
  Log_writer_mpi_buffer(int rank,
                        int messagelevel,
                        int buffer_size=160);
  ~Log_writer_mpi_buffer();

private:
  void put_buffer(void);
  void put_char(int);

protected:
  int overflow(int);
  int sync();

  int rank;
};

Log_writer_mpi::Log_writer_mpi(int rank,
                               int messagelevel)
    : Log_writer(new Log_writer_mpi_buffer(rank, messagelevel)) {}

// Buffer
Log_writer_mpi_buffer::Log_writer_mpi_buffer(int rank,
    int message_level,
    int buffer_size)
    : Log_writer_buffer(message_level, buffer_size),
    rank(rank) {}

Log_writer_mpi_buffer::~Log_writer_mpi_buffer() {
  sync();
}

int Log_writer_mpi_buffer::overflow(int c) {
  put_buffer();

  if (c != EOF) {
    if (pbase() == epptr()) {
      put_char(c);
    } else {
      sputc(c);
    }
  }

  return 0;
}

int Log_writer_mpi_buffer::sync() {
  put_buffer();
  return 0;
}

void Log_writer_mpi_buffer::put_char(int chr) {
  if (current_level <= max_level) {
    std::cout << char(chr);
  }
}
void Log_writer_mpi_buffer::put_buffer() {
  if (pbase() != pptr()) {
    int     len = (pptr() - pbase());
    SFXC_ASSERT(len > 0);
    char    *buffer = new char[len + 19 + 1];

    if (current_level <= max_level) {
      struct timeval time_struct;
      gettimeofday(&time_struct,NULL);
      // Round to the nearest millisecond
      suseconds_t ms = (time_struct.tv_usec + 500) / 1000;
      time_struct.tv_sec += ms / 1000;
      time_struct.tv_usec = (ms % 1000) * 1000;
      struct tm *tm_struct = localtime(&time_struct.tv_sec);
      snprintf(buffer, 20, "%02dh%02dm%02d.%03lds, %02d, ",
               tm_struct->tm_hour, tm_struct->tm_min, tm_struct->tm_sec,
               (time_struct.tv_usec+500)/1000, rank);
      SFXC_ASSERT(strlen(buffer) == 19);
      strncpy(buffer + 19, pbase(), len);
      buffer[len + 19] = '\0';

      // If MT_MPI is defined then acquire  the mutex.
      // otherwise do nothing.
      IF_MT_MPI_ENABLED( RAIIMutex mutex(g_mpi_thebig_mutex) );

      MPI_Send(buffer, len+20+1, MPI_CHAR, RANK_LOG_NODE, MPI_TAG_LOG_MESSAGE, MPI_COMM_WORLD);
    }

    setp(pbase(), epptr());
    delete [] buffer;
  }
}
