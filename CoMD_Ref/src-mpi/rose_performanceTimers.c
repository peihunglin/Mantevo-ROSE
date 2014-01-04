/// \file
/// Performance timer functions.
///
/// Use the timer functionality to collect timing and number of calls
/// information for chosen computations (such as force calls) and
/// communication (such as sends, receives, reductions).  Timing results
/// are reported at the end of the run showing overall timings and
/// statistics of timings across ranks.
///
/// A new timer can be added as follows:
/// -# add new handle to the TimerHandle in performanceTimers.h
/// -# provide a corresponding name in timerName
///
/// Note that the order of the handles and names must be the
/// same. This order also determines the order in which the timers are
/// printed. Names can contain leading spaces to show a hierarchical
/// ordering.  Timers with zero calls are omitted from the report.
///
/// Raw timer data is obtained from the getTime() and getTick()
/// functions.  The supplied portable versions of these functions can be
/// replaced with platform specific versions for improved accuracy or
/// lower latency.
/// \see TimerHandle
/// \see getTime
/// \see getTick
///
#include "performanceTimers.h"
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "performanceTimers.h"
#include "mytype.h"
#include "parallel.h"
#include "yamlOutput.h"
static uint64_t getTime();
static double getTick();
static void timerStats();
/// You must add timer name in same order as enum in .h file.
/// Leading spaces can be specified to show a hierarchy of timers.
char *timerName[numberOfTimers] = {("total"), ("loop"), ("timestep"), ("  position"), ("  velocity"), ("  redistribute"), ("    atomHalo"), ("  force"), ("    eamHalo"), ("commHalo"), ("commReduce")};
/// Timer data collected.  Also facilitates computing averages and
/// statistics.
typedef struct TimersSt {
//!< call start time
uint64_t start;
//!< current total time
uint64_t total;
//!< current call count
uint64_t count;
//!< lap time
uint64_t elapsed;
//!< rank with min value
int minRank;
//!< rank with max value
int maxRank;
//!< min over ranks
double minValue;
//!< max over ranks
double maxValue;
//!< average over ranks
double average;
//!< stdev across ranks
double stdev;}Timers;
/// Global timing data collected.
typedef struct TimerGlobalSt {
//!< average time (us) per atom per rank 
double atomRate;
//!< average time (us) per atom
double atomAllRate;
//!< average atoms per time (us)
double atomsPerUSec;}TimerGlobal;
static Timers perfTimer[numberOfTimers];
static TimerGlobal perfGlobal;

void profileStart(const enum TimerHandle handle)
{
  perfTimer[handle] . start = getTime();
}

void profileStop(const enum TimerHandle handle)
{
  perfTimer[handle] . count += 1;
  uint64_t delta = getTime() - perfTimer[handle] . start;
  perfTimer[handle] . total += delta;
  perfTimer[handle] . elapsed += delta;
}
/// \details
/// Return elapsed time (in seconds) since last call with this handle
/// and clear for next lap.

double getElapsedTime(const enum TimerHandle handle)
{
  double etime = getTick() * ((double )perfTimer[handle] . elapsed);
  perfTimer[handle] . elapsed = 0;
  return etime;
}
/// \details
/// The report contains two blocks.  The upper block is performance
/// information for the printRank.  The lower block is statistical
/// information over all ranks.

void printPerformanceResults(int nGlobalAtoms,int printRate)
{
// Collect timer statistics overall and across ranks
  timerStats();
  if (!printRank()) {
    return ;
  }
// only print timers with non-zero values.
  double tick = getTick();
  double loopTime = perfTimer[loopTimer] . total * tick;
  fprintf(stdout,"\n\nTimings for Rank %d\n",getMyRank());
  fprintf(stdout,"        Timer        # Calls    Avg/Call (s)   Total (s)    %% Loop\n");
  fprintf(stdout,"___________________________________________________________________\n");
  for (int ii = 0; ii < numberOfTimers; ++ii) {
    double totalTime = perfTimer[ii] . total * tick;
    if (perfTimer[ii] . count > 0) {
      fprintf(stdout,"%-16s%12lu     %8.4f      %8.4f    %8.2f\n",timerName[ii],perfTimer[ii] . count,totalTime / ((double )perfTimer[ii] . count),totalTime,totalTime / loopTime * 100.0);
    }
  }
  fprintf(stdout,"\nTiming Statistics Across %d Ranks:\n",getNRanks());
  fprintf(stdout,"        Timer        Rank: Min(s)       Rank: Max(s)      Avg(s)    Stdev(s)\n");
  fprintf(stdout,"_____________________________________________________________________________\n");
  for (int ii = 0; ii < numberOfTimers; ++ii) {
    if (perfTimer[ii] . count > 0) {
      fprintf(stdout,"%-16s%6d:%10.4f  %6d:%10.4f  %10.4f  %10.4f\n",timerName[ii],perfTimer[ii] . minRank,perfTimer[ii] . minValue * tick,perfTimer[ii] . maxRank,perfTimer[ii] . maxValue * tick,perfTimer[ii] . average * tick,perfTimer[ii] . stdev * tick);
    }
  }
  double atomsPerTask = nGlobalAtoms / ((real_t )(getNRanks()));
  perfGlobal . atomRate = perfTimer[timestepTimer] . average * tick * 1e6 / (atomsPerTask * perfTimer[timestepTimer] . count * printRate);
  perfGlobal . atomAllRate = perfTimer[timestepTimer] . average * tick * 1e6 / (nGlobalAtoms * perfTimer[timestepTimer] . count * printRate);
  perfGlobal . atomsPerUSec = 1.0 / perfGlobal . atomAllRate;
  fprintf(stdout,"\n---------------------------------------------------\n");
  fprintf(stdout," Average atom update rate:     %6.2f us/atom/task\n",perfGlobal . atomRate);
  fprintf(stdout,"---------------------------------------------------\n\n");
  fprintf(stdout,"\n---------------------------------------------------\n");
  fprintf(stdout," Average all atom update rate: %6.2f us/atom\n",perfGlobal . atomAllRate);
  fprintf(stdout,"---------------------------------------------------\n\n");
  fprintf(stdout,"\n---------------------------------------------------\n");
  fprintf(stdout," Average atom rate:            %6.2f atoms/us\n",perfGlobal . atomsPerUSec);
  fprintf(stdout,"---------------------------------------------------\n\n");
}

void printPerformanceResultsYaml(FILE *file)
{
  if (!printRank()) {
    return ;
  }
  double tick = getTick();
  double loopTime = perfTimer[loopTimer] . total * tick;
  fprintf(file,"\nPerformance Results:\n");
  fprintf(file,"  TotalRanks: %d\n",getNRanks());
  fprintf(file,"  ReportingTimeUnits: seconds\n");
  fprintf(file,"Performance Results For Rank %d:\n",getMyRank());
  for (int ii = 0; ii < numberOfTimers; ii++) {
    if (perfTimer[ii] . count > 0) {
      double totalTime = perfTimer[ii] . total * tick;
      fprintf(file,"  Timer: %s\n",timerName[ii]);
      fprintf(file,"    CallCount:  %lu\n",perfTimer[ii] . count);
      fprintf(file,"    AvgPerCall: %8.4f\n",totalTime / ((double )perfTimer[ii] . count));
      fprintf(file,"    Total:      %8.4f\n",totalTime);
      fprintf(file,"    PercentLoop: %8.2f\n",totalTime / loopTime * 100);
    }
  }
  fprintf(file,"Performance Results Across Ranks:\n");
  for (int ii = 0; ii < numberOfTimers; ii++) {
    if (perfTimer[ii] . count > 0) {
      fprintf(file,"  Timer: %s\n",timerName[ii]);
      fprintf(file,"    MinRank: %d\n",perfTimer[ii] . minRank);
      fprintf(file,"    MinTime: %8.4f\n",perfTimer[ii] . minValue * tick);
      fprintf(file,"    MaxRank: %d\n",perfTimer[ii] . maxRank);
      fprintf(file,"    MaxTime: %8.4f\n",perfTimer[ii] . maxValue * tick);
      fprintf(file,"    AvgTime: %8.4f\n",perfTimer[ii] . average * tick);
      fprintf(file,"    StdevTime: %8.4f\n",perfTimer[ii] . stdev * tick);
    }
  }
  fprintf(file,"Performance Global Update Rates:\n");
  fprintf(file,"  AtomUpdateRate:\n");
  fprintf(file,"    AverageRate: %6.2f\n",perfGlobal . atomRate);
  fprintf(file,"    Units: us/atom/task\n");
  fprintf(file,"  AllAtomUpdateRate:\n");
  fprintf(file,"    AverageRate: %6.2f\n",perfGlobal . atomAllRate);
  fprintf(file,"    Units: us/atom\n");
  fprintf(file,"  AtomRate:\n");
  fprintf(file,"    AverageRate: %6.2f\n",perfGlobal . atomsPerUSec);
  fprintf(file,"    Units: atoms/us\n");
  fprintf(file,"\n");
}
/// Returns current time as a 64-bit integer.  This portable version
/// returns the number of microseconds since mindight, Jamuary 1, 1970.
/// Hence, timing data will have a resolution of 1 microsecond.
/// Platforms with access to calls with lower latency or higher
/// resolution (such as a cycle counter) may wish to replace this
/// implementation and change the conversion factor in getTick as
/// appropriate.
/// \see getTick for the conversion factor between the integer time
/// units of this function and seconds.

static uint64_t getTime()
{
  struct timeval ptime;
  uint64_t t = 0;
  gettimeofday(&ptime,((struct timezone *)((void *)0)));
  t = ((uint64_t )1000000) * ((uint64_t )ptime . tv_sec) + ((uint64_t )ptime . tv_usec);
  return t;
}
/// Returns the factor for converting the integer time reported by
/// getTime into seconds.  The portable getTime returns values in units
/// of microseconds so the conversion is simply 1e-6.
/// \see getTime

static double getTick()
{
  double seconds_per_cycle = 1.0e-6;
  return seconds_per_cycle;
}
/// Collect timer statistics across ranks.

static void timerStats()
{
  double sendBuf[numberOfTimers];
  double recvBuf[numberOfTimers];
// Determine average of each timer across ranks
  for (int ii = 0; ii < numberOfTimers; ii++) 
    sendBuf[ii] = ((double )perfTimer[ii] . total);
  addDoubleParallel(sendBuf,recvBuf,numberOfTimers);
  for (int ii = 0; ii < numberOfTimers; ii++) 
    perfTimer[ii] . average = recvBuf[ii] / ((double )(getNRanks()));
// Determine min and max across ranks and which rank
  RankReduceData reduceSendBuf[numberOfTimers];
  RankReduceData reduceRecvBuf[numberOfTimers];
  for (int ii = 0; ii < numberOfTimers; ii++) {
    reduceSendBuf[ii] . val = ((double )perfTimer[ii] . total);
    reduceSendBuf[ii] . rank = getMyRank();
  }
  minRankDoubleParallel(reduceSendBuf,reduceRecvBuf,numberOfTimers);
  for (int ii = 0; ii < numberOfTimers; ii++) {
    perfTimer[ii] . minValue = reduceRecvBuf[ii] . val;
    perfTimer[ii] . minRank = reduceRecvBuf[ii] . rank;
  }
  maxRankDoubleParallel(reduceSendBuf,reduceRecvBuf,numberOfTimers);
  for (int ii = 0; ii < numberOfTimers; ii++) {
    perfTimer[ii] . maxValue = reduceRecvBuf[ii] . val;
    perfTimer[ii] . maxRank = reduceRecvBuf[ii] . rank;
  }
// Determine standard deviation
  for (int ii = 0; ii < numberOfTimers; ii++) {
    double temp = ((double )perfTimer[ii] . total) - perfTimer[ii] . average;
    sendBuf[ii] = temp * temp;
  }
  addDoubleParallel(sendBuf,recvBuf,numberOfTimers);
  for (int ii = 0; ii < numberOfTimers; ii++) {
    perfTimer[ii] . stdev = sqrt(recvBuf[ii] / ((double )(getNRanks())));
  }
}