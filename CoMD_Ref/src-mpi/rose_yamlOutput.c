/// \file
/// Write simulation information in YAML format.
///
/// Information regarding platform, run parameters, performance, etc.,
/// are written to a file whose name is generated from the CoMDVariant
/// and the time of the run.  This provides a simple mechanism to track
/// and compare performance etc.
///
/// There are much more sophisticated libraries and routines available
/// to handle YAML, but this simple implemenation handles everything we
/// really need.
#include "yamlOutput.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "CoMD_info.h"
#include "mytype.h"
#include "parallel.h"
FILE *yamlFile = ((void *)0);
static const char *CoMDVersion = "1.1";
static const char *CoMDVariant = "CoMD-mpi";

static void getTimeString(char *timestring)
{
  time_t rawtime;
  struct tm *timeinfo;
  time(&rawtime);
  timeinfo = localtime((&rawtime));
  sprintf(timestring,"%4d-%02i-%02d, %02d:%02d:%02d",timeinfo -> tm_year + 1900,timeinfo -> tm_mon + 1,timeinfo -> tm_mday,timeinfo -> tm_hour,timeinfo -> tm_min,timeinfo -> tm_sec);
}

void yamlBegin()
{
  if (!printRank()) {
    return ;
  }
  char filename[64];
  time_t rawtime;
  time(&rawtime);
  struct tm *ptm = localtime((&rawtime));
  char sdate[25];
//use tm_mon+1 because tm_mon is 0 .. 11 instead of 1 .. 12
  sprintf(sdate,"%04d:%02d:%02d-%02d:%02d:%02d",ptm -> tm_year + 1900,ptm -> tm_mon + 1,ptm -> tm_mday,ptm -> tm_hour,ptm -> tm_min,ptm -> tm_sec);
  sprintf(filename,"%s.%s.yaml",CoMDVariant,sdate);
  yamlFile = fopen(filename,"w");
}

void yamlAppInfo(FILE *file)
{
  if (!printRank()) {
    return ;
  }
  printSeparator(file);
  fprintf(file,"Mini-Application Name    : %s\n",CoMDVariant);
  fprintf(file,"Mini-Application Version : %s\n",CoMDVersion);
  fprintf(file,"Platform:\n");
  fprintf(file,"  hostname: %s\n","hudson-rose-33.llnl.gov");
  fprintf(file,"  kernel name: %s\n","'Linux'");
  fprintf(file,"  kernel release: %s\n","'2.6.32-358.23.2.el6.x86_64'");
  fprintf(file,"  processor: %s\n","'x86_64'");
  fprintf(file,"Build:\n");
  fprintf(file,"  CC: %s\n","'/home/lin32/Development/projects/rose/master/install/bin/identityTranslator'");
  fprintf(file,"  compiler version: %s\n","'ROSE (pre-release beta version: 0.9.6a)'");
  fprintf(file,"  CFLAGS: %s\n","'-I/home/lin32/Development/opt/OpenMPI/include -pthread -L/home/lin32/Development/opt/OpenMPI/lib  -std=c99 -fopenmp -rose:OpenMP -DDOUBLE -DDO_MPI -g -O5   '");
  fprintf(file,"  LDFLAGS: %s\n","' -lmpi -lm '");
  fprintf(file,"  using MPI: %s\n",(builtWithMpi()?"true" : "false"));
  fprintf(file,"  Threading: none\n");
  fprintf(file,"  Double Precision: %s\n","true");
  char timestring[32];
  getTimeString(timestring);
  fprintf(file,"Run Date/Time: %s\n",timestring);
  fprintf(file,"\n");
  fflush(file);
}

void yamlEnd()
{
  if (!printRank()) {
    return ;
  }
  fclose(yamlFile);
}

void printSeparator(FILE *file)
{
//fprintf(file,"=========================================================================\n");
  fprintf(file,"\n");
}
