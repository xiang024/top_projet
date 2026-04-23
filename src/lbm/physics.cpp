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

/// Opposite directions for bounce back implementation
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
  const double eps = 1e-12;

  // Loop on all dimensions
  for (size_t d = 0; d < DIMENSIONS; d++) {
    v[d] = 0.0;

    // Sum all directions
    for (size_t k = 0; k < DIRECTIONS; k++) {
      v[d] += cell[k] * direction_matrix[k][d];
    }

    // Normalize
    v[d] = (cell_density > eps) ? (v[d] / cell_density) : 0.0;
  }
}

double compute_equilibrium_profile(Vector velocity, double density, int direction) {
  const double v2 = get_vect_norm_2(velocity, velocity);

  // Compute `e_i * v_i / c`
  const double p  = get_vect_norm_2(direction_matrix[direction], velocity);
  const double p2 = p * p;

  // Terms without density and direction weight
  double f_eq = 1.0 + (3.0 * p) + ((9.0 / 2.0) * p2) - ((3.0 / 2.0) * v2);

  // Multiply everything by the density and direction weight
  f_eq *= equil_weight[direction] * density;

  return f_eq;
}

void compute_cell_collision(lbm_mesh_cell_t cell_out, const lbm_mesh_cell_t cell_in) {
  const double eps = 1e-12;

  // Compute macroscopic values directly in the hot path to reduce helper-call overhead.
  const double density = cell_in[0] + cell_in[1] + cell_in[2] + cell_in[3] + cell_in[4]
                       + cell_in[5] + cell_in[6] + cell_in[7] + cell_in[8];

  Vector v = {0.0, 0.0};
  if (density > eps) {
    v[0] = (cell_in[1] - cell_in[3] + cell_in[5] - cell_in[6] - cell_in[7] + cell_in[8]) / density;
    v[1] = (cell_in[2] - cell_in[4] + cell_in[5] + cell_in[6] - cell_in[7] - cell_in[8]) / density;
  }

  // Loop on microscopic directions
  #pragma omp simd
  for (size_t k = 0; k < DIRECTIONS; k++) {
    // Compute f at equilibrium
    double f_eq = compute_equilibrium_profile(v, density, k);
    // Compute f_out
    cell_out[k] = cell_in[k] - RELAX_PARAMETER * (cell_in[k] - f_eq);
  }
}

void compute_bounce_back(lbm_mesh_cell_t cell) {
  double tmp[DIRECTIONS];
  for (size_t k = 0; k < DIRECTIONS; k++) {
    tmp[k] = cell[opposite_of[k]];
  }
  for (size_t k = 0; k < DIRECTIONS; k++) {
    cell[k] = tmp[k];
  }
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

  // Set macroscopic fluid info
  // Poiseuille distribution on X and null on Y
  // We just want the norm, so `v = v_x`
  const double v = helper_compute_poiseuille(id_y, MESH_HEIGHT);


  // Compute rho from U and inner flow on surface
  const double rho = (cell[0] + cell[2] + cell[4] + 2 * (cell[3] + cell[6] + cell[7])) / (1.0 - v);

  // Now compute unknown microscopic values
  cell[1] = cell[3]; // + (2.0/3.0) * density * v_y <--- no velocity on Y so v_y = 0
  cell[5] = cell[7] - (1.0 / 2.0) * (cell[2] - cell[4])
            + (1.0 / 6.0) * (rho * v); // + (1.0/2.0) * rho * v_y    <--- no velocity on Y so v_y = 0
  cell[8] = cell[6] + (1.0 / 2.0) * (cell[2] - cell[4])
            + (1.0 / 6.0) * (rho * v); //- (1.0/2.0) * rho * v_y    <--- no velocity on Y so v_y = 0

  // No need to copy already known one as the value will be "loss" in the wall at propagatation time
}

void compute_outflow_zou_he_const_density(lbm_mesh_cell_t cell) {
#if DIRECTIONS != 9
#error Implemented only for 9 directions
#endif

  double const rho = 1.0;
  // Compute macroscopic velocity depending on inner flow going onto the wall
  const double v = -1.0 + (1.0 / rho) * (cell[0] + cell[2] + cell[4] + 2 * (cell[1] + cell[5] + cell[8]));

  // Now can compute unknown microscopic values
  cell[3] = cell[1] - (2.0 / 3.0) * rho * v;
  cell[7] = cell[5]
            + (1.0 / 2.0) * (cell[2] - cell[4])
            // - (1.0/2.0) * (rho * v_y)    <--- no velocity on Y so v_y = 0
            - (1.0 / 6.0) * (rho * v);
  cell[6] = cell[8]
            + (1.0 / 2.0) * (cell[4] - cell[2])
            // + (1.0/2.0) * (rho * v_y)    <--- no velocity on Y so v_y = 0
            - (1.0 / 6.0) * (rho * v);
}

void special_cells(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  // Handle left inflow boundary in a dedicated linear pass.
  for (size_t j = 1; j < mesh->height - 1; j++) {
    if (*(lbm_cell_type_t_get_cell(mesh_type, 1, j)) == CELL_LEFT_IN) {
      compute_inflow_zou_he_poiseuille_distr(mesh, Mesh_get_cell(mesh, 1, j), j + mesh_comm->y);
    }
  }

  // Handle right outflow boundary in a dedicated linear pass.
  for (size_t j = 1; j < mesh->height - 1; j++) {
    if (*(lbm_cell_type_t_get_cell(mesh_type, mesh->width - 2, j)) == CELL_RIGHT_OUT) {
      compute_outflow_zou_he_const_density(Mesh_get_cell(mesh, mesh->width - 2, j));
    }
  }

  // Loop on all inner cells for bounce-back cells (obstacle and top/bottom walls).
  #pragma omp parallel for schedule(static)
  for (size_t i = 1; i < mesh->width - 1; i++) {
    for (size_t j = 1; j < mesh->height - 1; j++) {
      if (*(lbm_cell_type_t_get_cell(mesh_type, i, j)) == CELL_BOUNCE_BACK) {
        compute_bounce_back(Mesh_get_cell(mesh, i, j));
      }
    }
  }
}

void collision(Mesh* mesh_out, const Mesh* mesh_in) {
  assert(mesh_in->width == mesh_out->width);
  assert(mesh_in->height == mesh_out->height);

  // Loop on all inner cells
  #pragma omp parallel for schedule(static)
  for (size_t i = 1; i < mesh_in->width - 1; i++) {
    for (size_t j = 1; j < mesh_in->height - 1; j++) {
      lbm_mesh_cell_t cell_out       = Mesh_get_cell(mesh_out, i, j);
      const lbm_mesh_cell_t cell_in  = Mesh_get_cell(mesh_in, i, j);
      compute_cell_collision(cell_out, cell_in);
    }
  }
}

void propagation(Mesh* mesh_out, const Mesh* mesh_in) {
  const ssize_t width  = static_cast<ssize_t>(mesh_out->width);
  const ssize_t height = static_cast<ssize_t>(mesh_out->height);

  // Process each direction on its valid source rectangle to avoid bounds checks
  // in the inner loops.
  #pragma omp parallel
{
  for (size_t k = 0; k < DIRECTIONS; k++) {
    const ssize_t dx = static_cast<ssize_t>(direction_matrix[k][0]);
    const ssize_t dy = static_cast<ssize_t>(direction_matrix[k][1]);

    const ssize_t i_begin = (dx < 0) ? -dx : 0;
    const ssize_t i_end   = (dx > 0) ? (width - dx) : width;
    const ssize_t j_begin = (dy < 0) ? -dy : 0;
    const ssize_t j_end   = (dy > 0) ? (height - dy) : height;
    #pragma omp for schedule(static) nowait
    for (ssize_t i = i_begin; i < i_end; i++) {
      for (ssize_t j = j_begin; j < j_end; j++) {
        Mesh_get_cell(mesh_out, i + dx, j + dy)[k] = Mesh_get_cell(mesh_in, i, j)[k];
      }
    }
  }
}
}
