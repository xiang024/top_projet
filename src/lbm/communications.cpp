#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <lbm/communications.hpp>
#include <lbm/physics.hpp>

void save_frame(FILE* fp, const Mesh* mesh) {
  lbm_file_entry_t buffer[WRITE_BUFFER_ENTRIES];
  size_t cnt = 0;

  for (size_t i = 1; i < mesh->width - 1; i++) {
    for (size_t j = 1; j < mesh->height - 1; j++) {
      const lbm_mesh_cell_t cell = Mesh_get_cell(mesh, i, j);
      const double density      = get_cell_density(cell);
      Vector v;
      get_cell_velocity(v, cell, density);
      const double norm = std::sqrt(get_vect_norm_2(v, v));

      buffer[cnt].rho = (float)density;
      buffer[cnt].v   = (float)norm;
      cnt++;
      assert(cnt <= WRITE_BUFFER_ENTRIES);

      if (cnt == WRITE_BUFFER_ENTRIES) {
        fwrite(buffer, sizeof(lbm_file_entry_t), cnt, fp);
        cnt = 0;
      }
    }
  }

  if (cnt != 0) {
    fwrite(buffer, sizeof(lbm_file_entry_t), cnt, fp);
  }
}

static int helper_get_rank_id(int nb_x, int nb_y, int rank_x, int rank_y) {
  if (rank_x < 0 || rank_x >= nb_x || rank_y < 0 || rank_y >= nb_y) {
    return -1;
  }
  return rank_x + rank_y * nb_x;
}

static void choose_domain_split(int comm_size, uint32_t width, uint32_t height, int* nb_x, int* nb_y) {
  // Prefer an x-only split when possible: it preserves the file order and avoids top/bottom halo traffic.
  if (comm_size > 0 && width % (uint32_t)comm_size == 0) {
    *nb_x = comm_size;
    *nb_y = 1;
    return;
  }

  int best_x      = 1;
  int best_y      = comm_size;
  double best_key = 1.0e300;

  for (int x = 1; x <= comm_size; ++x) {
    if (comm_size % x != 0) {
      continue;
    }
    const int y = comm_size / x;
    if (width % (uint32_t)x != 0 || height % (uint32_t)y != 0) {
      continue;
    }

    const double local_w = (double)width / (double)x;
    const double local_h = (double)height / (double)y;
    const double key     = 2.0 * (local_w + local_h);
    if (key < best_key) {
      best_key = key;
      best_x   = x;
      best_y   = y;
    }
  }

  *nb_x = best_x;
  *nb_y = best_y;
}

void lbm_comm_print(const lbm_comm_t* mesh_comm) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  static bool first_call = true;
  if (first_call && rank == RANK_MASTER) {
    first_call = false;
    fprintf(
      stderr,
      "%4s| %8s %8s %8s %8s | %12s %12s %12s %12s | %6s %6s | %6s %6s\n",
      "RANK",
      "TOP",
      "BOTTOM",
      "LEFT",
      "RIGHT",
      "TOP LEFT",
      "TOP RIGHT",
      "BOTTOM LEFT",
      "BOTTOM RIGHT",
      "POS X",
      "POS Y",
      "DIM X",
      "DIM Y"
    );
  }

  // One barrier only for readable startup printing, not in the halo hot path.
  MPI_Barrier(MPI_COMM_WORLD);
  fprintf(
    stderr,
    "%4d| %7d  %7d  %7d  %7d  | %11d  %11d  %11d  %11d  | %5d  %5d  | %5d  %5d \n",
    rank,
    mesh_comm->top_id,
    mesh_comm->bottom_id,
    mesh_comm->left_id,
    mesh_comm->right_id,
    mesh_comm->corner_id[CORNER_TOP_LEFT],
    mesh_comm->corner_id[CORNER_TOP_RIGHT],
    mesh_comm->corner_id[CORNER_BOTTOM_LEFT],
    mesh_comm->corner_id[CORNER_BOTTOM_RIGHT],
    mesh_comm->x,
    mesh_comm->y,
    mesh_comm->width,
    mesh_comm->height
  );
}

void lbm_comm_init(lbm_comm_t* mesh_comm, int rank, int comm_size, uint32_t width, uint32_t height) {
  int nb_x = 1;
  int nb_y = 1;
  choose_domain_split(comm_size, width, height, &nb_x, &nb_y);

  assert(nb_x * nb_y == comm_size);
  if (width % (uint32_t)nb_x != 0 || height % (uint32_t)nb_y != 0) {
    fatal("Can't get a valid 2D cut for current problem size and number of processes.");
  }

  const int rank_x = rank % nb_x;
  const int rank_y = rank / nb_x;

  mesh_comm->nb_x = nb_x;
  mesh_comm->nb_y = nb_y;

  mesh_comm->width  = width / nb_x + 2;
  mesh_comm->height = height / nb_y + 2;

  mesh_comm->x = rank_x * width / nb_x;
  mesh_comm->y = rank_y * height / nb_y;

  mesh_comm->left_id   = helper_get_rank_id(nb_x, nb_y, rank_x - 1, rank_y);
  mesh_comm->right_id  = helper_get_rank_id(nb_x, nb_y, rank_x + 1, rank_y);
  mesh_comm->top_id    = helper_get_rank_id(nb_x, nb_y, rank_x, rank_y - 1);
  mesh_comm->bottom_id = helper_get_rank_id(nb_x, nb_y, rank_x, rank_y + 1);

  mesh_comm->corner_id[CORNER_TOP_LEFT]     = helper_get_rank_id(nb_x, nb_y, rank_x - 1, rank_y - 1);
  mesh_comm->corner_id[CORNER_TOP_RIGHT]    = helper_get_rank_id(nb_x, nb_y, rank_x + 1, rank_y - 1);
  mesh_comm->corner_id[CORNER_BOTTOM_LEFT]  = helper_get_rank_id(nb_x, nb_y, rank_x - 1, rank_y + 1);
  mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] = helper_get_rank_id(nb_x, nb_y, rank_x + 1, rank_y + 1);

  MPI_Type_vector(mesh_comm->width - 2, DIRECTIONS, mesh_comm->height * DIRECTIONS, MPI_DOUBLE, &mesh_comm->row_type);
  MPI_Type_commit(&mesh_comm->row_type);
  mesh_comm->buffer = NULL;

  lbm_comm_print(mesh_comm);
}

void lbm_comm_release(lbm_comm_t* mesh_comm) {
  mesh_comm->x         = 0;
  mesh_comm->y         = 0;
  mesh_comm->width     = 0;
  mesh_comm->height    = 0;
  mesh_comm->right_id  = -1;
  mesh_comm->left_id   = -1;
  mesh_comm->top_id    = -1;
  mesh_comm->bottom_id = -1;
  if (mesh_comm->row_type != MPI_DATATYPE_NULL) {
    MPI_Type_free(&mesh_comm->row_type);
    mesh_comm->row_type = MPI_DATATYPE_NULL;
  }
  if (mesh_comm->buffer != NULL) {
    free(mesh_comm->buffer);
    mesh_comm->buffer = NULL;
  }
}

static inline int mpi_rank_or_null(int rank) {
  return rank < 0 ? MPI_PROC_NULL : rank;
}

static bool has_any_neighbor(const lbm_comm_t* mesh) {
  return mesh->left_id >= 0 || mesh->right_id >= 0 || mesh->top_id >= 0 || mesh->bottom_id >= 0
         || mesh->corner_id[CORNER_TOP_LEFT] >= 0 || mesh->corner_id[CORNER_TOP_RIGHT] >= 0
         || mesh->corner_id[CORNER_BOTTOM_LEFT] >= 0 || mesh->corner_id[CORNER_BOTTOM_RIGHT] >= 0;
}

static void exchange_column(Mesh* mesh, int send_to, int recv_from, uint32_t send_x, uint32_t recv_x, int tag) {
  MPI_Sendrecv(
    Mesh_get_cell(mesh, send_x, 1),
    (mesh->height - 2) * DIRECTIONS,
    MPI_DOUBLE,
    mpi_rank_or_null(send_to),
    tag,
    Mesh_get_cell(mesh, recv_x, 1),
    (mesh->height - 2) * DIRECTIONS,
    MPI_DOUBLE,
    mpi_rank_or_null(recv_from),
    tag,
    MPI_COMM_WORLD,
    MPI_STATUS_IGNORE
  );
}

static void exchange_row(
  const lbm_comm_t* mesh_comm,
  Mesh* mesh,
  int send_to,
  int recv_from,
  uint32_t send_y,
  uint32_t recv_y,
  int tag
) {
  MPI_Sendrecv(
    Mesh_get_cell(mesh, 1, send_y),
    1,
    mesh_comm->row_type,
    mpi_rank_or_null(send_to),
    tag,
    Mesh_get_cell(mesh, 1, recv_y),
    1,
    mesh_comm->row_type,
    mpi_rank_or_null(recv_from),
    tag,
    MPI_COMM_WORLD,
    MPI_STATUS_IGNORE
  );
}

static void exchange_corner(
  Mesh* mesh,
  int send_to,
  int recv_from,
  uint32_t send_x,
  uint32_t send_y,
  uint32_t recv_x,
  uint32_t recv_y,
  int tag
) {
  MPI_Sendrecv(
    Mesh_get_cell(mesh, send_x, send_y),
    DIRECTIONS,
    MPI_DOUBLE,
    mpi_rank_or_null(send_to),
    tag,
    Mesh_get_cell(mesh, recv_x, recv_y),
    DIRECTIONS,
    MPI_DOUBLE,
    mpi_rank_or_null(recv_from),
    tag,
    MPI_COMM_WORLD,
    MPI_STATUS_IGNORE
  );
}

void lbm_comm_halo_exchange(lbm_comm_t* mesh, Mesh* mesh_to_process) {
  if (!has_any_neighbor(mesh)) {
    return;
  }

  exchange_column(mesh_to_process, mesh->right_id, mesh->left_id, mesh->width - 2, 0, 10);
  exchange_column(mesh_to_process, mesh->left_id, mesh->right_id, 1, mesh->width - 1, 11);

  exchange_row(mesh, mesh_to_process, mesh->bottom_id, mesh->top_id, mesh->height - 2, 0, 20);
  exchange_row(mesh, mesh_to_process, mesh->top_id, mesh->bottom_id, 1, mesh->height - 1, 21);

  exchange_corner(
    mesh_to_process,
    mesh->corner_id[CORNER_TOP_LEFT],
    mesh->corner_id[CORNER_BOTTOM_RIGHT],
    1,
    1,
    mesh->width - 1,
    mesh->height - 1,
    30
  );
  exchange_corner(
    mesh_to_process,
    mesh->corner_id[CORNER_TOP_RIGHT],
    mesh->corner_id[CORNER_BOTTOM_LEFT],
    mesh->width - 2,
    1,
    0,
    mesh->height - 1,
    31
  );
  exchange_corner(
    mesh_to_process,
    mesh->corner_id[CORNER_BOTTOM_LEFT],
    mesh->corner_id[CORNER_TOP_RIGHT],
    1,
    mesh->height - 2,
    mesh->width - 1,
    0,
    32
  );
  exchange_corner(
    mesh_to_process,
    mesh->corner_id[CORNER_BOTTOM_RIGHT],
    mesh->corner_id[CORNER_TOP_LEFT],
    mesh->width - 2,
    mesh->height - 2,
    0,
    0,
    33
  );
}

void save_frame_all_domain(FILE* fp, Mesh* source_mesh, Mesh* temp) {
  int comm_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (comm_size <= 1) {
    save_frame(fp, source_mesh);
    return;
  }

  if (rank == RANK_MASTER) {
    save_frame(fp, source_mesh);
    for (int i = 1; i < comm_size; i++) {
      MPI_Status status;
      MPI_Recv(
        temp->cells,
        source_mesh->width * source_mesh->height * DIRECTIONS,
        MPI_DOUBLE,
        i,
        0,
        MPI_COMM_WORLD,
        &status
      );
      save_frame(fp, temp);
    }
  } else {
    MPI_Send(
      source_mesh->cells,
      source_mesh->width * source_mesh->height * DIRECTIONS,
      MPI_DOUBLE,
      RANK_MASTER,
      0,
      MPI_COMM_WORLD
    );
  }
}
