#include <lbm/initialization.hpp>

#include <cassert>
#include <cstddef>

#include <mpi.h>

#include <lbm/physics.hpp>

void init_cond_velocity_0_density_1(Mesh* mesh) {
  assert(mesh != NULL);

  for (size_t i = 0; i < mesh->width; i++) {
    for (size_t j = 0; j < mesh->height; j++) {
      for (size_t k = 0; k < DIRECTIONS; k++) {
        Mesh_get_cell(mesh, i, j)[k] = equil_weight[k];
      }
    }
  }
}

void setup_init_state_circle_obstacle(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  for (size_t j = mesh_comm->y; j < mesh->height + mesh_comm->y; j++) {
    for (size_t i = mesh_comm->x; i < mesh->width + mesh_comm->x; i++) {
      if (((i - OBSTACLE_X) * (i - OBSTACLE_X)) + ((j - OBSTACLE_Y) * (j - OBSTACLE_Y)) <=
          OBSTACLE_R * OBSTACLE_R) {
        *(lbm_cell_type_t_get_cell(mesh_type, i - mesh_comm->x, j - mesh_comm->y)) = CELL_BOUNCE_BACK;
      }
    }
  }
}

void setup_init_state_global_poiseuille_profile(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  Vector v         = {0.0, 0.0};
  const double rho = 1.0;

  for (size_t i = 0; i < mesh->width; i++) {
    const size_t global_i = i + mesh_comm->x;

    for (size_t j = 0; j < mesh->height; j++) {
      const size_t global_j = j + mesh_comm->y;
      v[0]                  = helper_compute_poiseuille(global_j, MESH_HEIGHT);

      lbm_mesh_cell_t cell = Mesh_get_cell(mesh, i, j);
      *(lbm_cell_type_t_get_cell(mesh_type, i, j)) = CELL_FUILD;

      for (size_t k = 0; k < DIRECTIONS; k++) {
        cell[k] = compute_equilibrium_profile(v, rho, k);

        // Use global x coordinate: otherwise MPI ranks after rank 0 initialize their
        // local left ghost/interior columns as if they were the physical inflow.
        if (global_i > 1) {
          cell[k] = equil_weight[k];
        }
      }
    }
  }
}

void setup_init_state_border(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  Vector v         = {0.0, 0.0};
  const double rho = 1.0;

  if (mesh_comm->left_id == -1) {
    for (size_t j = 1; j < mesh->height - 1; j++) {
      *(lbm_cell_type_t_get_cell(mesh_type, 0, j)) = CELL_LEFT_IN;
    }
  }

  if (mesh_comm->right_id == -1) {
    for (size_t j = 1; j < mesh->height - 1; j++) {
      *(lbm_cell_type_t_get_cell(mesh_type, mesh->width - 1, j)) = CELL_RIGHT_OUT;
    }
  }

  if (mesh_comm->top_id == -1) {
    for (size_t i = 0; i < mesh->width; i++) {
      for (size_t k = 0; k < DIRECTIONS; k++) {
        Mesh_get_cell(mesh, i, 0)[k] = compute_equilibrium_profile(v, rho, k);
        *(lbm_cell_type_t_get_cell(mesh_type, i, 0)) = CELL_BOUNCE_BACK;
      }
    }
  }

  if (mesh_comm->bottom_id == -1) {
    for (size_t i = 0; i < mesh->width; i++) {
      for (size_t k = 0; k < DIRECTIONS; k++) {
        Mesh_get_cell(mesh, i, mesh->height - 1)[k] = compute_equilibrium_profile(v, rho, k);
        *(lbm_cell_type_t_get_cell(mesh_type, i, mesh->height - 1)) = CELL_BOUNCE_BACK;
      }
    }
  }
}

void setup_init_state(Mesh* mesh, lbm_mesh_type_t* mesh_type, const lbm_comm_t* mesh_comm) {
  setup_init_state_global_poiseuille_profile(mesh, mesh_type, mesh_comm);
  setup_init_state_border(mesh, mesh_type, mesh_comm);
  setup_init_state_circle_obstacle(mesh, mesh_type, mesh_comm);
}