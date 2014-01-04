/// \file
/// Wrappers for MPI functions.  This should be the only compilation 
/// unit in the code that directly calls MPI functions.  To build a pure
/// serial version of the code with no MPI, do not define DO_MPI.  If
/// DO_MPI is not defined then all MPI functionality is replaced with
/// equivalent single task behavior.
#include "parallel.h"
#ifdef DO_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>
static int myRank = 0;
static int nRanks = 1;
#ifdef DO_MPI
#ifdef SINGLE
#define REAL_MPI_TYPE MPI_FLOAT
#else
#define REAL_MPI_TYPE MPI_DOUBLE
#endif
#endif

int getNRanks()
{
  return nRanks;
}

int getMyRank()
{
  return myRank;
}
/// \details
/// For now this is just a check for rank 0 but in principle it could be
/// more complex.  It is also possible to suppress practically all
/// output by causing this function to return 0 for all ranks.

int printRank()
{
  if (myRank == 0) {
    return 1;
  }
  return 0;
}

void timestampBarrier(const char *msg)
{
  barrierParallel();
  if (!printRank()) {
    return ;
  }
  time_t t = time(((void *)0));
  char *timeString = ctime((&t));
// clobber newline
  timeString[24] = '\0';
  fprintf(stdout,"%s: %s\n",timeString,msg);
  fflush(stdout);
}

void initParallel(int *argc,char ***argv)
{
#ifdef DO_MPI
  MPI_Init(argc,argv);
  MPI_Comm_rank(((MPI_Comm )((void *)(&ompi_mpi_comm_world))),&myRank);
  MPI_Comm_size(((MPI_Comm )((void *)(&ompi_mpi_comm_world))),&nRanks);
#endif
}

void destroyParallel()
{
#ifdef DO_MPI
  MPI_Finalize();
#endif
}

void barrierParallel()
{
#ifdef DO_MPI
  MPI_Barrier(((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#endif
}
/// \param [in]  sendBuf Data to send.
/// \param [in]  sendLen Number of bytes to send.
/// \param [in]  dest    Rank in MPI_COMM_WORLD where data will be sent.
/// \param [out] recvBuf Received data.
/// \param [in]  recvLen Maximum number of bytes to receive.
/// \param [in]  source  Rank in MPI_COMM_WORLD from which to receive.
/// \return Number of bytes received.

int sendReceiveParallel(void *sendBuf,int sendLen,int dest,void *recvBuf,int recvLen,int source)
{
#ifdef DO_MPI
  int bytesReceived;
  MPI_Status status;
  MPI_Sendrecv(sendBuf,sendLen,((MPI_Datatype )((void *)(&ompi_mpi_byte))),dest,0,recvBuf,recvLen,((MPI_Datatype )((void *)(&ompi_mpi_byte))),source,0,((MPI_Comm )((void *)(&ompi_mpi_comm_world))),&status);
  MPI_Get_count(&status,((MPI_Datatype )((void *)(&ompi_mpi_byte))),&bytesReceived);
  return bytesReceived;
#else
#endif
}

void addIntParallel(int *sendBuf,int *recvBuf,int count)
{
#ifdef DO_MPI
  MPI_Allreduce(sendBuf,recvBuf,count,((MPI_Datatype )((void *)(&ompi_mpi_int))),((MPI_Op )((void *)(&ompi_mpi_op_sum))),((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#else
#endif
}

void addRealParallel(real_t *sendBuf,real_t *recvBuf,int count)
{
#ifdef DO_MPI
  MPI_Allreduce(sendBuf,recvBuf,count,((MPI_Datatype )((void *)(&ompi_mpi_double))),((MPI_Op )((void *)(&ompi_mpi_op_sum))),((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#else
#endif
}

void addDoubleParallel(double *sendBuf,double *recvBuf,int count)
{
#ifdef DO_MPI
  MPI_Allreduce(sendBuf,recvBuf,count,((MPI_Datatype )((void *)(&ompi_mpi_double))),((MPI_Op )((void *)(&ompi_mpi_op_sum))),((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#else
#endif
}

void maxIntParallel(int *sendBuf,int *recvBuf,int count)
{
#ifdef DO_MPI
  MPI_Allreduce(sendBuf,recvBuf,count,((MPI_Datatype )((void *)(&ompi_mpi_int))),((MPI_Op )((void *)(&ompi_mpi_op_max))),((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#else
#endif
}

void minRankDoubleParallel(RankReduceData *sendBuf,RankReduceData *recvBuf,int count)
{
#ifdef DO_MPI
  MPI_Allreduce(sendBuf,recvBuf,count,((MPI_Datatype )((void *)(&ompi_mpi_double_int))),((MPI_Op )((void *)(&ompi_mpi_op_minloc))),((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#else
#endif
}

void maxRankDoubleParallel(RankReduceData *sendBuf,RankReduceData *recvBuf,int count)
{
#ifdef DO_MPI
  MPI_Allreduce(sendBuf,recvBuf,count,((MPI_Datatype )((void *)(&ompi_mpi_double_int))),((MPI_Op )((void *)(&ompi_mpi_op_maxloc))),((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#else
#endif
}
/// \param [in] count Length of buf in bytes.

void bcastParallel(void *buf,int count,int root)
{
#ifdef DO_MPI
  MPI_Bcast(buf,count,((MPI_Datatype )((void *)(&ompi_mpi_byte))),root,((MPI_Comm )((void *)(&ompi_mpi_comm_world))));
#endif
}

int builtWithMpi()
{
#ifdef DO_MPI
  return 1;
#else
#endif
}
