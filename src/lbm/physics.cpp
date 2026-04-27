#include <lbm/physics.hpp>

#include <cassert>
#include <cstdlib>

#include <omp.h>

#include <lbm/communications.hpp>
#include <lbm/config.hpp>
#include <lbm/structures.hpp>

#if DIRECTIONS == 9 && DIMENSIONS == 2
/// Definition of the 9 base vectors used to discretize the directions on each mesh.
const Vector direction_matrix[DIRECTIONS] = {
  // clang-format off
  {+0.0, +0.0},
  {+1.0, +0.0}, {+0.0, +1.0}, {-1.0, +0.0}, {+0.0, -1.0},
  {+1.0, +1.0}, {-1.0, +1.0}, {-1.0, -1.0}, {+1.0, -1.0},
  // clang-format on
};
#else
#error Need to define adapted direction matrix.
#endif

#if DIRECTIONS == 9
/// Weigths used to compensate the differences in lenght of the 9 directional vectors.
const double equil_weight[DIRECTIONS] = {
  // clang-format off
  4.0 / 9.0,
  1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
  1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
  // clang-format on
};

/// Opposite directions for bounce back implementation.
const int opposite_of[DIRECTIONS] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
#else
#error Need to define adapted equilibrium distribution function
#endif

double get_vect_norm_2(Vector const a, Vector const b) {
  double res = 0.0;
  for (size_t k = 0; k < DIMENSIONS; k++) {
    res += a[k] * b[k];
  }
  return res;
}

double get_cell_density(const lbm_mesh_cell_t cell) {
  assert(cell != NULL);
  double res = 0.0;
  for (size_t k = 0; k < DIRECTIONS; k++) {
    res += cell[k];
  }
  return res;
}

void get_cell_velocity(Vector v, const lbm_mesh_cell_t cell, double cell_density) {
  assert(v != NULL);
  assert(cell != NULL);

  for (size_t d = 0; d < DIMENSIONS; d++) {
    v[d] = 0.0;
    for (size_t k = 0; k < DIRECTIONS; k++) {
      v[d] += cell[k] * direction_matrix[k][d];
    }
    v[d] /= cell_density;
  }
}

double compute_equilibrium_profile(Vector velocity, double density, int direction) {
  const double v2 = get_vect_norm_2(velocity, velocity);
  const double p  = get_vect_norm_2(direction_matrix[direction], velocity);
  const double p2 = p * p;

  double f_eq = 1.0 + (3.0 * p) + ((9.0 / 2.0) * p2) - ((3.0 / 2.0) * v2);
  f_eq *= equil_weight[direction] * density;

  return f_eq;
}

#ifndef LBM_FAST_COLLISION
#define LBM_FAST_COLLISION 0
#endif

#if LBM_FAST_COLLISION
static inline double relax_population(double f, double feq) {
  return f - RELAX_PARAMETER * (f - feq);
}
#endif

void compute_cell_collision(lbm_mesh_cell_t cell_out, const lbm_mesh_cell_t cell_in) {
#if LBM_FAST_COLLISION
  // D2Q9 fast path. Enable only after validating the chosen reference checksum.
  const double f0 = cell_in[0];
  const double f1 = cell_in[1];
  const double f2 = cell_in[2];
  const double f3 = cell_in[3];
  const double f4 = cell_in[4];
  const double f5 = cell_in[5];
  const double f6 = cell_in[6];
  const double f7 = cell_in[7];
  const double f8 = cell_in[8];

  const double rho     = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
  const double inv_rho = 1.0 / rho;
  const double ux      = (f1 - f3 + f5 - f6 - f7 + f8) * inv_rho;
  const double uy      = (f2 - f4 + f5 + f6 - f7 - f8) * inv_rho;
  const double u2      = ux * ux + uy * uy;
  const double common  = 1.0 - 1.5 * u2;

  const double w0 = 4.0 / 9.0;
  const double ws = 1.0 / 9.0;
  const double wd = 1.0 / 36.0;

  const double c1 = ux;
  const double c2 = uy;
  const double c3 = -ux;
  const double c4 = -uy;
  const double c5 = ux + uy;
  const double c6 = -ux + uy;
  const double c7 = -ux - uy;
  const double c8 = ux - uy;

  cell_out[0] = relax_population(f0, w0 * rho * common);
  cell_out[1] = relax_population(f1, ws * rho * (common + 3.0 * c1 + 4.5 * c1 * c1));
  cell_out[2] = relax_population(f2, ws * rho * (common + 3.0 * c2 + 4.5 * c2 * c2));
  cell_out[3] = relax_population(f3, ws * rho * (common + 3.0 * c3 + 4.5 * c3 * c3));
  cell_out[4] = relax_population(f4, ws * rho * (common + 3.0 * c4 + 4.5 * c4 * c4));
  cell_out[5] = relax_population(f5, wd * rho * (common + 3.0 * c5 + 4.5 * c5 * c5));
  cell_out[6] = relax_population(f6, wd * rho * (common + 3.0 * c6 + 4.5 * c6 * c6));
  cell_out[7] = relax_population(f7, wd * rho * (common + 3.0 * c7 + 4.5 * c7 * c7));
  cell_out[8] = relax_population(f8, wd * rho * (common + 3.0 * c8 + 4.5 * c8 * c8));
#else
  // Reference-safe version: keeps the original operation order.
  const double density = get_cell_density(cell_in);
  Vector v;
  get_cell_velocity(v, cell_in, density);

  for (size_t k = 0; k < DIRECTIONS; k++) {
    const double f_eq = compute_equilibrium_profile(v, density, k);
    cell_out[k]       = cell_in[k] - RELAX_PARAMETER * (cell_in[k] - f_eq);
  }
#endif
}

void compute_bounce_back(lbm_mesh_cell_t cell) {
  const double f1 = cell[1];
  const double f2 = cell[2];
  const double f5 = cell[5];
  const double f6 = cell[6];

  cell[1] = cell[3];
  cell[2] = cell[4];
  cell[3] = f1;
  cell[4] = f2;
  cell[5] = cell[7];
  cell[6] = cell[8];
  cell[7] = f5;
  cell[8] = f6;
}

double helper_compute_poiseuille(const size_t i, const size_t size) {
  const double y = (double)(i - 1);
  const double L = (double)(size - 1);
  return 4.0 * INFLOW_MAX_VELOCITY / (L * L) * (L * y - y * y);
}

void compute_inflow_zou_he_poiseuille_distr(const Mesh* mesh, lbm_mesh_cell_t cell, size_t id_y) {
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif
  (void)mesh;

  // Use the global height: id_y is already a global Y coordinate.
  const double v   = helper_compute_poiseuille(id_y, MESH_HEIGHT);
  const double rho = (cell[0] + cell[2] + cell[4] + 2.0 * (cell[3] + cell[6] + cell[7])) / (1.0 - v);
  const double rv  = rho * v;

  cell[1] = cell[3] + (2.0 / 3.0) * rv;
  cell[5] = cell[7] - 0.5 * (cell[2] - cell[4]) + (1.0 / 6.0) * rv;
  cell[8] = cell[6] + 0.5 * (cell[2] - cell[4]) + (1.0 / 6.0) * rv;
}

void compute_outflow_zou_he_const_density(lbm_mesh_cell_t cell) {
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif

  const double rho = 1.0;
  const double v   = -1.0 + (cell[0] + cell[2] + cell[4] + 2.0 * (cell[1] + cell[5] + cell[8])) / rho;

  cell[3] = cell[1] - (2.0 / 3.0) * rho * v;
  cell[7] = cell[5] + 0.5 * (cell[2] - cell[4]) - (1.0 / 6.0) * rho * v;
  cell[6] = cell[8] + 0.5 * (cell[4] - cell[2]) - (1.0 / 6.0) * rho * v;
}

static void build_special_cell_cache(lbm_mesh_type_t* mesh_type) {
  size_t count = 0;

  // Keep the exact same domain as the previous special_cells() loop.
  for (size_t i = 1; i < mesh_type->width - 1; i++) {
    for (size_t j = 1; j < mesh_type->height - 1; j++) {
      if (*(lbm_cell_type_t_get_cell(mesh_type, i, j)) != CELL_FUILD) {
        count++;
      }
    }
  }

  free(mesh_type->special_x);
  free(mesh_type->special_y);
  free(mesh_type->special_type);

  mesh_type->special_x     = NULL;
  mesh_type->special_y     = NULL;
  mesh_type->special_type  = NULL;
  mesh_type->special_count = count;
  mesh_type->special_ready = 1;

  if (count == 0) {
    return;
  }

  mesh_type->special_x = static_cast<uint32_t*>(malloc(count * sizeof(uint32_t)));
  mesh_type->special_y = static_cast<uint32_t*>(malloc(count * sizeof(uint32_t)));
  mesh_type->special_type = static_cast<lbm_cell_type_t*>(malloc(count * sizeof(lbm_cell_type_t)));

  if (mesh_type->special_x == NULL || mesh_type->special_y == NULL || mesh_type->special_type == NULL) {
    perror("malloc");
    abort();
  }

  size_t p = 0;
  for (size_t i = 1; i < mesh_type->width - 1; i++) {
    for (size_t j = 1; j < mesh_type->height - 1; j++) {
      const lbm_cell_type_t type = *(lbm_cell_type_t_get_cell(mesh_type, i, j));
      if (type != CELL_FUILD) {
        mesh_type->special_x[p] = (uint32_t)i;
        mesh_type->special_y[p] = (uint32_t)j;
        mesh_type->special_type[p] = type;
        p++;
      }
    }
  }
}

void special_cells(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  if (!mesh_type->special_ready) {
    build_special_cell_cache(mesh_type);
  }

  // The cached list is usually much smaller than the full mesh.
  // Keep it serial to avoid creating an OpenMP team at every time step.
  for (size_t p = 0; p < mesh_type->special_count; p++) {
    const uint32_t i = mesh_type->special_x[p];
    const uint32_t j = mesh_type->special_y[p];

    switch (mesh_type->special_type[p]) {
    case CELL_FUILD:
      break;
    case CELL_BOUNCE_BACK:
      compute_bounce_back(Mesh_get_cell(mesh, i, j));
      break;
    case CELL_LEFT_IN:
      compute_inflow_zou_he_poiseuille_distr(mesh, Mesh_get_cell(mesh, i, j), j + mesh_comm->y);
      break;
    case CELL_RIGHT_OUT:
      compute_outflow_zou_he_const_density(Mesh_get_cell(mesh, i, j));
      break;
    }
  }
}

void collision(Mesh* mesh_out, const Mesh* mesh_in) {
#if DIRECTIONS != 9 || DIMENSIONS != 2
#error collision() fast loop is specialized for D2Q9.
#endif
  assert(mesh_in->width == mesh_out->width);
  assert(mesh_in->height == mesh_out->height);

  const size_t width  = mesh_in->width;
  const size_t height = mesh_in->height;
  const size_t stride_x = height * DIRECTIONS;

  double* __restrict__ const out_cells = mesh_out->cells;
  const double* __restrict__ const in_cells = mesh_in->cells;

  const double omega = RELAX_PARAMETER;
  const double w0 = 4.0 / 9.0;
  const double ws = 1.0 / 9.0;
  const double wd = 1.0 / 36.0;

#pragma omp parallel for collapse(2) schedule(static)
  for (size_t i = 1; i < width - 1; i++) {
    for (size_t j = 1; j < height - 1; j++) {
      const size_t id = i * stride_x + j * DIRECTIONS;

      const double f0 = in_cells[id + 0];
      const double f1 = in_cells[id + 1];
      const double f2 = in_cells[id + 2];
      const double f3 = in_cells[id + 3];
      const double f4 = in_cells[id + 4];
      const double f5 = in_cells[id + 5];
      const double f6 = in_cells[id + 6];
      const double f7 = in_cells[id + 7];
      const double f8 = in_cells[id + 8];

      const double rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
      const double inv_rho = 1.0 / rho;
      const double ux = (f1 - f3 + f5 - f6 - f7 + f8) * inv_rho;
      const double uy = (f2 - f4 + f5 + f6 - f7 - f8) * inv_rho;
      const double u2 = ux * ux + uy * uy;
      const double common = 1.0 - 1.5 * u2;

      const double c1 = ux;
      const double c2 = uy;
      const double c3 = -ux;
      const double c4 = -uy;
      const double c5 = ux + uy;
      const double c6 = -ux + uy;
      const double c7 = -ux - uy;
      const double c8 = ux - uy;

      const double feq0 = w0 * rho * common;
      const double feq1 = ws * rho * (common + 3.0 * c1 + 4.5 * c1 * c1);
      const double feq2 = ws * rho * (common + 3.0 * c2 + 4.5 * c2 * c2);
      const double feq3 = ws * rho * (common + 3.0 * c3 + 4.5 * c3 * c3);
      const double feq4 = ws * rho * (common + 3.0 * c4 + 4.5 * c4 * c4);
      const double feq5 = wd * rho * (common + 3.0 * c5 + 4.5 * c5 * c5);
      const double feq6 = wd * rho * (common + 3.0 * c6 + 4.5 * c6 * c6);
      const double feq7 = wd * rho * (common + 3.0 * c7 + 4.5 * c7 * c7);
      const double feq8 = wd * rho * (common + 3.0 * c8 + 4.5 * c8 * c8);

      out_cells[id + 0] = f0 - omega * (f0 - feq0);
      out_cells[id + 1] = f1 - omega * (f1 - feq1);
      out_cells[id + 2] = f2 - omega * (f2 - feq2);
      out_cells[id + 3] = f3 - omega * (f3 - feq3);
      out_cells[id + 4] = f4 - omega * (f4 - feq4);
      out_cells[id + 5] = f5 - omega * (f5 - feq5);
      out_cells[id + 6] = f6 - omega * (f6 - feq6);
      out_cells[id + 7] = f7 - omega * (f7 - feq7);
      out_cells[id + 8] = f8 - omega * (f8 - feq8);
    }
  }
}

void propagation(Mesh* mesh_out, const Mesh* mesh_in) {
#if DIRECTIONS != 9 || DIMENSIONS != 2
#error propagation() is specialized for D2Q9.
#endif

  const size_t width = mesh_out->width;
  const size_t height = mesh_in->height;
  const size_t stride_x = height * DIRECTIONS;
  double* __restrict__ const out_cells = mesh_out->cells;
  const double* __restrict__ const in_cells = mesh_in->cells;

#pragma omp parallel for collapse(2) schedule(static)
  for (size_t i = 1; i < width - 1; i++) {
    for (size_t j = 1; j < height - 1; j++) {
      const size_t id = i * stride_x + j * DIRECTIONS;

      out_cells[id + 0] = in_cells[id + 0];
      out_cells[id + 1] = in_cells[id - stride_x + 1];
      out_cells[id + 2] = in_cells[id - DIRECTIONS + 2];
      out_cells[id + 3] = in_cells[id + stride_x + 3];
      out_cells[id + 4] = in_cells[id + DIRECTIONS + 4];
      out_cells[id + 5] = in_cells[id - stride_x - DIRECTIONS + 5];
      out_cells[id + 6] = in_cells[id + stride_x - DIRECTIONS + 6];
      out_cells[id + 7] = in_cells[id + stride_x + DIRECTIONS + 7];
      out_cells[id + 8] = in_cells[id - stride_x + DIRECTIONS + 8];
    }
  }
}
