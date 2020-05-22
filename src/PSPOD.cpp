/**
 * @file PSPOD.cpp
 * @author Weiqi Shen (weiqishen1994@ufl.edu)
 * @brief Code entrance
 * @version 0.1
 * @date 2019-01-28
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <time.h>
#include "global.h"
#include "funcs.h"

using namespace std;

int main(int argc, char *argv[])
{
  //Parsing command line options
  if (argc < 2)
  {
    cout << "For help use -h,--help" << endl;
    return 0;
  }
  if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) //help
  {
    cout << "PSPOD <input_file>" << endl;
    cout << "For more help, go to Github repository." << endl;
    return 0;
  }

  //read input
  run_input.setup(argv[1]);

  clock_t t = clock();
  //start main program
  switch (run_input.task)
  {
  case SNAPSHOT_POD:
    CalculateSpatialPOD();
    break;
  case SPECTRAL_POD:
    CalculateSpectralPOD();
    break;
  case AZIMUTHAL_SPOD:
    CalculateAzimuthalSpectralPOD();
    break;
  default:
    Fatal_Error("Feature not implemented");
  }

  t = clock() - t;
  cout << "Execution time= " << (double)t / ((double)CLOCKS_PER_SEC) << " seconds" << endl;
  return 0;
}
