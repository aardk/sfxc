#include <utils.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <ifaddrs.h>

#include <assert.h>

INT64 get_us_time(int time[]) {
  INT64 result = 0;
  // time[0] is year
  result = time[1];                // days
  result = time[2] +   24* result; // hours
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


#include <genFunctions.h>
#include <iostream>
using namespace std;
//constants
#include "constPrms.h"

//class and function definitions
#include "runPrms.h"
#include "genPrms.h"
#include "staPrms.h"
#include "genFunctions.h"
#include "InData.h"
#include "ProcessData.h"

#include <Data_reader_file.h>

#define SEED 10

//global variables
//declaration and default settings run parameters
extern RunP RunPrms;
//declaration and default settings general parameters
extern GenP GenPrms;
//station parameters class, declaration and default settings
extern StaP StaPrms[NstationsMax];
// used for randomising numbers for Headers in Mk4 file
extern UINT32 seed;
//declarations for offsets
extern INT64 sliceStartByte[NstationsMax][NprocessesMax];
extern INT64 sliceStartTime [NprocessesMax];
extern INT64 sliceStopTime  [NprocessesMax];
extern INT64 sliceTime;


/** Initialises the global control files, this should be removed at some point.
 **/
int
initialise_control(char *filename) {
  int    i, Nstations;
  
//  int status, numtasks, rank;

  // seed the random number generator (global variable!)
  seed = (UINT32) time((time_t *)NULL);
  seed = SEED;
  std::cout << "seed: " << seed << std::endl;

  //parse control file for run parameters
  if (RunPrms.parse_ctrlFile(filename) != 0) {
    cerr << "ERROR: Control file "<< filename <<", program aborted.\n";
    return -1;
  }

  //show version information and control file info
  if (RunPrms.get_messagelvl() > 0)
    cout << "\nSource " << __FILE__ << " compiled at: "
         << __DATE__ << " " <<__TIME__ << "\n\n"
         << "Control file name "  <<  filename << "\n\n";
  
  //check control parameters, optionally show them
  if (RunPrms.check_params() != 0) {
    cerr << "ERROR: Run control parameter, program aborted.\n";
    return -1;
  }
  
  if (RunPrms.get_interactive() && RunPrms.get_messagelvl()> 0)
    askContinue();

  //parse control file for general parameters
  if (GenPrms.parse_ctrlFile(filename) != 0) {
    cerr << "ERROR: Control file "<< filename <<", program aborted.\n";
    return -1;
  }

  //check general control parameters, optionally show them
  if (GenPrms.check_params() != 0) {
    cerr << "ERROR: General control parameter, program aborted.\n";
    return -1;
  }
  
  if (RunPrms.get_interactive() && RunPrms.get_messagelvl()> 0)
    askContinue();

  //get the number of stations
  Nstations = GenPrms.get_nstations();
  
  //parse the control file for all station parameters
  for (i=0; i<Nstations; i++)
    if (StaPrms[i].parse_ctrlFile(filename,i) != 0 ) {
      cerr << "ERROR: Control file "<< filename <<", program aborted.\n";
      return -1;
    }
    
  //check station control parameters, optionally show them
  for (i=0; i<Nstations; i++){
    if (StaPrms[i].check_params() != 0 ) {
      cerr << "ERROR: Station control parameter, program aborted.\n";
      return -1;
    }
    if (RunPrms.get_interactive() && RunPrms.get_messagelvl()> 0)
      askContinue();
  }
  
  return 0;  
}

void send_control_data(int rank) {
}

void receive_control_data(MPI_Status &status) {
}