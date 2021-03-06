/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id: tcp_connection.h 412 2007-12-05 12:13:20Z kruithof $
 *
 */

#ifndef TCP_CONNECTION_H_
#define TCP_CONNECTION_H_

#include <vector>
#include <string>

#include <stdio.h>
#include <stdint.h>
#include <string.h>

class TCP_Connection {
public:
  TCP_Connection(bool verbose = false);
  virtual ~TCP_Connection();

  /// Open a port on the server side
  bool open_port(unsigned short int port,
                 int connections);
  /// Open a connection on the server size
  unsigned int open_connection();

  /// Client side connect
  int do_connect(uint64_t, unsigned short int port);
  int do_connect(const char *hostname, unsigned short int port);

  void get_ip_addresses(std::vector<uint64_t> &addr);
  void get_ip_addresses(std::vector<std::string> &addr);

  bool is_localhost(uint64_t ip_addr) {
    return (ip_addr == 0x100007f);
  }
  bool is_localhost(const char *hostname) {
    return (strcmp(hostname, "127.0.0.1") == 0);
  }
  std::string ip_addr(uint64_t ip) {
    char addr[16];
    sprintf(addr, "%d.%d.%d.%d",
            (int)ip&255, (int)(ip>>8)&255, (int)(ip>>16)&255, (int)ip>>24);
    return std::string(addr);
  }
  std::string get_address() {
    return connection_addr;
  }

  int get_port();
private:
  bool verbose;
  int connection_socket;
  int port_nr;
  char connection_addr[17];
};

#endif /*TCP_CONNECTION_H_*/
