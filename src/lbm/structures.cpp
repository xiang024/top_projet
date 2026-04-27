#include <lbm/structures.hpp>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

void Mesh_init(Mesh* mesh, uint32_t width, uint32_t height) {
  mesh->width  = width;
  mesh->height = height;

  const size_t count = (size_t)width * (size_t)height * (size_t)DIRECTIONS;
  mesh->cells = static_cast<double*>(aligned_alloc(64, ((count * sizeof(double) + 63) / 64) * 64));
  if (mesh->cells == NULL) {
    perror("aligned_alloc");
    abort();
  }
}

void Mesh_release(Mesh* mesh) {
  mesh->width  = 0;
  mesh->height = 0;
  free(mesh->cells);
  mesh->cells = NULL;
}

void lbm_mesh_type_t_init(lbm_mesh_type_t* meshtype, uint32_t width, uint32_t height) {
  meshtype->width  = width;
  meshtype->height = height;

  meshtype->types = static_cast<lbm_cell_type_t*>(malloc((width + 2) * (size_t)height * sizeof(lbm_cell_type_t)));
  if (meshtype->types == NULL) {
    perror("malloc");
    abort();
  }

  meshtype->special_x = NULL;
  meshtype->special_y = NULL;
  meshtype->special_type = NULL;
  meshtype->special_count = 0;
  meshtype->special_ready = 0;
}

void lbm_mesh_type_t_release(lbm_mesh_type_t* mesh) {
  mesh->width  = 0;
  mesh->height = 0;
  free(mesh->types);
  free(mesh->special_x);
  free(mesh->special_y);
  free(mesh->special_type);
  mesh->types = NULL;
  mesh->special_x = NULL;
  mesh->special_y = NULL;
  mesh->special_type = NULL;
  mesh->special_count = 0;
  mesh->special_ready = 0;
}

void fatal(const char* message) {
  fprintf(stderr, "FATAL ERROR : %s\n", message);
  abort();
}
