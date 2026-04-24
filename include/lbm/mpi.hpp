#pragma once

#if defined(TOP_ENABLE_MPI)

#include <mpi.h>

#else

#include <chrono>
#include <cstring>

typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Request;
typedef struct MPI_Status {
  int dummy;
} MPI_Status;

#define MPI_COMM_WORLD 0
#define MPI_DOUBLE 0
#define MPI_STATUS_IGNORE nullptr

inline int MPI_Init(int*, char***) {
  return 0;
}

inline int MPI_Finalize() {
  return 0;
}

inline int MPI_Comm_rank(MPI_Comm, int* rank) {
  *rank = 0;
  return 0;
}

inline int MPI_Comm_size(MPI_Comm, int* size) {
  *size = 1;
  return 0;
}

inline int MPI_Barrier(MPI_Comm) {
  return 0;
}

inline double MPI_Wtime() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

inline int MPI_Send(const void*, int, MPI_Datatype, int, int, MPI_Comm) {
  return 0;
}

inline int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status*) {
  return 0;
}

inline int MPI_Sendrecv(
  const void* sendbuf,
  int sendcount,
  MPI_Datatype,
  int,
  int,
  void* recvbuf,
  int recvcount,
  MPI_Datatype,
  int,
  int,
  MPI_Comm,
  MPI_Status*
) {
  if (sendbuf != nullptr && recvbuf != nullptr && sendcount == recvcount) {
    std::memcpy(recvbuf, sendbuf, static_cast<size_t>(sendcount) * sizeof(double));
  }
  return 0;
}

inline int MPI_Sendrecv_replace(void* buf, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status*) {
  (void)buf;
  return 0;
}

#endif
