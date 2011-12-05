/* Copyright (c) 2007 Joint Institute for VLBI in Europe (Netherlands)
 * All rights reserved.
 *
 * Author(s): Nico Kruithof <Kruithof@JIVE.nl>, 2007
 *
 * $Id$
 */
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <ifaddrs.h>
#include <iostream>

#include "data_reader_file.h"
#include "exception_common.h"
#include "utils.h"
#ifdef USE_MPI
#ifdef FLOAT
#undef FLOAT
#endif
#include "sfxc_mpi.h"
int RANK_OF_NODE = -1; // Rank of the current node
#else
int RANK_OF_NODE = getpid(); // Rank of the current node
#endif


void abort_sfxc_assertion(const char *file, int line, const char* message) 
{
  std::cout << "#" << RANK_OF_NODE << " "
  << file << ", l" << line
  << ", Assertion failed: " << message << std::endl;

#ifdef USE_MPI
  int numtasks;
  // get the number of tasks set at commandline (= number of processors)
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  for (int i=0; i<numtasks; i++) {
    if (i!=RANK_OF_NODE) {
      int32_t msg=1; // 1 means error
      MPI_Send(&msg, 1, MPI_INT32, i,
               MPI_TAG_END_NODE, MPI_COMM_WORLD);
    }
  }

  // Close this node
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
#endif

  exit(-1);
}

void sfxc_abort(const char *msg){
  // Note : default argument msg=""
  if(strlen(msg)>0)
    std::cout << "Node #" << RANK_OF_NODE << " fatal error : " << msg << "\n";
  else
    std::cout << "Node #" << RANK_OF_NODE << " caused termination of all processes\n";


#ifdef USE_MPI
  // Kill all other nodes 
  int numtasks;
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  for (int i=0; i<numtasks; i++) {
    if (i!=RANK_OF_NODE) {
      int32_t msg=1; // 1 means error
      MPI_Send(&msg, 1, MPI_INT32, i,
               MPI_TAG_END_NODE, MPI_COMM_WORLD);
    }
  }

  // Close this node
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
#endif

  exit(-1);
}

bool isPower2(int val){
  return ((val > 0) && ((val & (~val + 1)) == val));
}

int64_t get_us_time(int time[]) {
  int64_t result = 0;
  // time[0] is year
  result = time[2]; // hours
  result = time[3] +   60* result; // minutes
  result = time[4] +   60* result; // seconds
  result = 0       + 1000* result; // milisecs
  result = 0       + 1000* result; // microsecs

  return result;
}

void get_ip_address(std::list<Interface_pair> &addresses,
                    bool IPv4_only) {
  struct ifaddrs *ifa = NULL;

  if (getifaddrs (&ifa) < 0) {
    perror ("getifaddrs");
    return;
  }

  for (; ifa; ifa = ifa->ifa_next) {
    char ip[ 200 ];
    socklen_t salen;

    if (ifa->ifa_addr->sa_family == AF_INET) {
      // IP v4
      salen = sizeof (struct sockaddr_in);
    } else if (!IPv4_only && (ifa->ifa_addr->sa_family == AF_INET6)) {
      // IP v6
      salen = sizeof (struct sockaddr_in6);
    } else {
      continue;
    }

    if (getnameinfo (ifa->ifa_addr, salen,
                     ip, sizeof (ip), NULL, 0, NI_NUMERICHOST) < 0) {
      perror ("getnameinfo");
      continue;
    }
    addresses.push_back(Interface_pair(std::string(ifa->ifa_name),
                                       std::string(ip)));
  }

  freeifaddrs (ifa);
}

//*****************************************************************************
//  irbit2: random seeding
//  See Numerical Recipes
//  primitive polynomial mod 2 of order n produces 2^n - 1 random bits
//*****************************************************************************
unsigned long iseed = 42;
void set_seed_1_bit(unsigned long seed_) {
  SFXC_ASSERT(seed_ != 0);
  iseed = seed_;
}

int irbit2() {
#define IB1 1
  //  #define IB2 2
#define IB4 8
  //  #define IB5 16
#define IB6 32
  //  #define IB18 131072
#define IB30 536870912
#define MASK (IB1+IB4+IB6)

  if (iseed & IB30) {
    iseed=((iseed ^ MASK) << 1) | IB1;
    return 1;
  } else {
    iseed <<= 1;
    return 0;
  }
#undef MASK
#undef IB30
  //  #undef IB18
#undef IB6
  //  #undef IB5
#undef IB4
  //  #undef IB2
#undef IB1
}

long park_miller_seed = 42;
void park_miller_set_seed(unsigned long seed_) {
  SFXC_ASSERT(seed_ != 0);
  park_miller_seed = seed_;
}

// Generated 31 random bits at a time
long unsigned int park_miller_random() {
  long unsigned int hi, lo;

  lo = 16807 * (park_miller_seed & 0xFFFF);
  hi = 16807 * (park_miller_seed >> 16);

  lo += (hi & 0x7FFF) << 16;
  lo += hi >> 15;

  if (lo > 0x7FFFFFFF)
    lo -= 0x7FFFFFFF;

  return (park_miller_seed = (long)lo);
}

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
bool directory_exist(const std::string& path) {
  struct stat buf;
  if ( stat(path.c_str(), &buf) == 0) {
    return S_ISDIR( buf.st_mode );
  }

  if ( errno == ENOENT ) {
    return false;
  }
  MTHROW(std::string("Unable to test if the directory exist: ")+strerror(errno));
}

void create_directory(const std::string& path) {
  if ( mkdir(path.c_str(), S_IRWXU) < 0 ) {
    MTHROW(std::string("Unable to create a directory :")+path+" because of:"+strerror(errno));
  }
}

std::string itoa (int32_t n) {
  char s[17];
  std::string u;

  if (n < 0) { //turns n positive
    n = (-1 * n);
    u = "-"; //adds '-' on result string
  }

  int i=0; //s counter

  do {
    s[i++]= n%10 + '0'; //conversion of each digit of n to char
    n -= n%10; //update n value
  }

  while ((n /= 10) > 0);

  for (int j = i-1; j >= 0; j--) {
    u += s[j]; //building our string number
  }

  return u;
}

int mjd(int day, int month, int year)
// Calculate the modified julian day, formula taken from the all knowing wikipedia
{
  int a = (14-month)/12;
  int y = year + 4800 - a;
  int m = month + 12*a - 3;
  int jdn = day + ((153*m+2)/5) + 365*y + (y/4) - (y/100) + (y/400) - 32045;
  return jdn - 2400000.5;
}
