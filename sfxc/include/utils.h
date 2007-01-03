#ifndef UTILS_H
#define UTILS_H

#include <types.h>
#include <list>
#include <string>

/// Interface_pair is a pair of two strings 
/// containing the interface name and its ip address;
typedef std::pair<std::string, std::string> Interface_pair;

/** 
 * get_ip_address retrieves the ip-addresses of the current host.
 * @returns a list of Interface_pair.
 **/ 
void get_ip_address(std::list<Interface_pair> &addresses,
                    bool IPv4_only = true);
          
/** 
 * get the time in micro seconds from an array (year, day, hour, minute, second).
 * @returns a list of Interface_pair.
 **/ 
INT64 get_us_time(int time[]);


/** Initialises the global control files, this should be removed at some point.
  **/
int initialise_control(char *filename);
                    
/** Constructs an array containing all control data for a correlate node (Gen and Run).
 *  @returns the size of the array
  **/
void send_control_data(int rank);

/** Stores all control data for a correlate node (Gen and Run).
  **/
void receive_control_data(MPI_Status &status);
                    
                    
void send_station_control_data(int rank, int station);
void receive_station_control_data(MPI_Status &status);
                    
#endif // UTILS_H
