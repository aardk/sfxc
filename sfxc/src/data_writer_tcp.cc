/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 *
 */

#include "data_writer_tcp.h"
#include "tcp_connection.h"
#include "utils.h"

#include <iostream>
// defines send:
#include <sys/types.h>
#include <sys/socket.h>
#include <errno.h>

#include <sys/poll.h>

Data_writer_tcp::Data_writer_tcp()
    : socket(-1) {}

void Data_writer_tcp::open_connection(TCP_Connection &tcp_connection) {
  socket = tcp_connection.open_connection();
  SFXC_ASSERT_MSG(socket > 0,
                  "Couldn't connect nodes, please kill other sfxc: 'killall -9 sfxc'");
}

Data_writer_tcp::~Data_writer_tcp() {
  close(socket);
}

size_t
Data_writer_tcp::do_put_bytes(size_t nBytes, const char *buff) {
  if (socket <= 0) return 0;
  SFXC_ASSERT(nBytes > 0);
  size_t bytes_written = 0;
  while (bytes_written != nBytes) {
    ssize_t result = write(socket, buff+bytes_written, nBytes-bytes_written);

    if (result <= 0) {
      return bytes_written;
    }
    bytes_written += result;
  }
  SFXC_ASSERT(bytes_written == nBytes);
  return bytes_written;
}

bool Data_writer_tcp::can_write() {
//   struct pollfd {
//     int fd;           /* file descriptor */
//     short events;     /* requested events */
//     short revents;    /* returned events */
//   };
  pollfd fds[1];
  fds[0].fd = socket;
  fds[0].events = POLLOUT;

  int ret = poll(fds, 1, /*timeout in miliseconds*/ 0);
  if (ret > 0) {
    return ((fds[0].revents & POLLOUT) != 0);
  }
  return false;
}

