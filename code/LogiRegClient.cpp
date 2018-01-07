#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"
#include "gwasiter.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <iostream>
#include <sstream>

using namespace NTL;
using namespace std;

int main(int argc, char** argv) {
  if (argc < 3) {
    cout << "Usage: LogiRegClient party_id param_file" << endl;
    return 1;
  }

  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 2) {
    cout << "Error: party_id should be 0, 1, or 2" << endl;
    return 1;
  }

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));

  /* Initialize MPC environment */
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }

  bool success = logireg_protocol(mpc, pid);

  // This is here just to keep P0 online until the end for data transfer
  // In practice, P0 would send data in advance before each phase and go offline
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  if (success) {
    cout << "Protocol successfully completed" << endl;
    return 0;
  } else {
    cout << "Protocol abnormally terminated" << endl;
    return 1;
  }
}
