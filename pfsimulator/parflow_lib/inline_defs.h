#ifndef _INLINE_DEFS_H
#define _INLINE_DEFS_H

#include "inline_saturation.h"
#include "inline_phase_source.h"
#include "inline_phase_density.h"
#include "inline_phase_rel_perm.h"
#include "inline_bc_pressure.h"
#include "inline_richards_internal.h"
#include "inline_overland_flow.h"
#include "inline_overland_flow_diffusive.h"
#include "inline_overland_flow_kinematic.h"

#define InlinePFVConstInit(c, z)                \
  {                                             \
    Grid       *grid = VectorGrid(z);           \
    Subgrid    *subgrid;                        \
    Subvector  *z_sub;                          \
    double     *zp;                             \
    int ix, iy, iz;                             \
    int nx, ny, nz;                             \
    int nx_z, ny_z, nz_z;                       \
    int sg, i, j, k, i_z;                       \
    ForSubgridI(sg, GridSubgrids(grid))         \
    {                                           \
      subgrid = GridSubgrid(grid, sg);          \
      z_sub = VectorSubvector(z, sg);           \
      ix = SubgridIX(subgrid);                  \
      iy = SubgridIY(subgrid);                  \
      iz = SubgridIZ(subgrid);                  \
      nx = SubgridNX(subgrid);                  \
      ny = SubgridNY(subgrid);                  \
      nz = SubgridNZ(subgrid);                  \
      nx_z = SubvectorNX(z_sub);                \
      ny_z = SubvectorNY(z_sub);                \
      nz_z = SubvectorNZ(z_sub);                \
      zp = SubvectorElt(z_sub, ix, iy, iz);     \
      i_z = 0;                                  \
      BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,  \
                i_z, nx_z, ny_z, nz_z, 1, 1, 1,   \
      {                                           \
        zp[i_z] = c;                              \
      });                                         \
    }                                             \
  }

#define InlineInitVector(vector, value)         \
  {                                            \
    Grid *grid = VectorGrid(vector);            \
    Subvector *v_sub;                           \
    double *vp;                                 \
    Subgrid *subgrid;                           \
    int ix;                                     \
    int iy;                                     \
    int iz;                                     \
    int nx;                                     \
    int ny;                                     \
    int nz;                                     \
    int nx_v;                                   \
    int ny_v;                                   \
    int nz_v;                                   \
    int i_s;                                    \
    int iv;                                     \
    ForSubgridI(i_s, GridSubgrids(grid))        \
    {                                           \
      subgrid = GridSubgrid(grid, i_s);         \
      v_sub = VectorSubvector(vector, i_s);     \
      ix = SubgridIX(subgrid);                  \
      iy = SubgridIY(subgrid);                  \
      iz = SubgridIZ(subgrid);                  \
                                                \
      nx = SubgridNX(subgrid);                  \
      ny = SubgridNY(subgrid);                  \
      nz = SubgridNZ(subgrid);                  \
                                                \
      nx_v = SubvectorNX(v_sub);                \
      ny_v = SubvectorNY(v_sub);                \
      nz_v = SubvectorNZ(v_sub);                \
                                                \
      vp = SubvectorElt(v_sub, ix, iy, iz);     \
                                                \
      iv = 0;                                   \
      BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,  \
                iv, nx_v, ny_v, nz_v, 1, 1, 1,    \
      {                                           \
        vp[iv] = value;                           \
      });                                         \
    }                                             \
  }

#define InlineInitMatrix(matrix, value)         \
  {                                             \
    Grid *grid = MatrixGrid(matrix);            \
    Submatrix *A_sub;                           \
    double *Ap;                                 \
    int im;                                     \
    SubgridArray *subgrids;                     \
    Subgrid *subgrid;                           \
    Stencil *stencil;                            \
    int is;                                     \
    int s;                                      \
    int ix;                                     \
    int iy;                                     \
    int iz;                                     \
    int nx;                                     \
    int ny;                                     \
    int nz;                                     \
    int nx_m;                                   \
    int ny_m;                                   \
    int nz_m;                                   \
    subgrids = GridSubgrids(grid);                  \
    ForSubgridI(is, subgrids)                       \
    {                                               \
      subgrid = SubgridArraySubgrid(subgrids, is);  \
      A_sub = MatrixSubmatrix(matrix, is);          \
      ix = SubgridIX(subgrid);                      \
      iy = SubgridIY(subgrid);                      \
      iz = SubgridIZ(subgrid);                      \
                                                    \
      nx = SubgridNX(subgrid);                      \
      ny = SubgridNY(subgrid);                      \
      nz = SubgridNZ(subgrid);                      \
                                                    \
      nx_m = SubmatrixNX(A_sub);                    \
      ny_m = SubmatrixNY(A_sub);                    \
      nz_m = SubmatrixNZ(A_sub);                    \
                                                    \
      stencil = MatrixStencil(matrix);              \
      for (s = 0; s < StencilSize(stencil); s++)    \
      {                                             \
        Ap = SubmatrixElt(A_sub, s, ix, iy, iz);    \
                                                    \
        im = 0;                                     \
        BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,  \
                  im, nx_m, ny_m, nz_m, 1, 1, 1,    \
        {                                           \
          Ap[im] = value;                           \
        });                                         \
      }                                             \
    }                                               \
  }

#define InlineNewVectorType(grid, nc, num_ghost, vtype, var)            \
  {                                                                     \
    Vector *new_vector;                                                 \
    new_vector = ctalloc(Vector, 1);                                    \
    (new_vector->subvectors) = ctalloc(Subvector*, GridNumSubgrids(grid)); \
    int data_size = 0;                                                  \
    int i;                                                              \
    int n;                                                              \
    VectorDataSpace(new_vector) = NewSubgridArray();                    \
    ForSubgridI(i, GridSubgrids(grid))                                  \
    {                                                                   \
      Subvector *new_sub = ctalloc(Subvector, 1);                       \
      Subgrid *subgrid = GridSubgrid(grid, i);                          \
      SubvectorDataSpace(new_sub)                                       \
        = NewSubgrid(SubgridIX(subgrid) - num_ghost,                    \
                     SubgridIY(subgrid) - num_ghost,                    \
                     SubgridIZ(subgrid) - num_ghost,                    \
                     SubgridNX(subgrid) + 2 * num_ghost,                \
                     SubgridNY(subgrid) + 2 * num_ghost,                \
                     SubgridNZ(subgrid) + 2 * num_ghost,                \
                     SubgridRX(subgrid),                                \
                     SubgridRY(subgrid),                                \
                     SubgridRZ(subgrid),                                \
                     SubgridProcess(subgrid));                          \
      AppendSubgrid(SubvectorDataSpace(new_sub), VectorDataSpace(new_vector)); \
      n = SubvectorNX(new_sub) * SubvectorNY(new_sub) * SubvectorNZ(new_sub); \
      data_size = n;                                                    \
      VectorSubvector(new_vector, i) = new_sub;                         \
    }                                                                   \
    (new_vector->data_size) = data_size;                                \
    VectorGrid(new_vector) = grid;                                      \
    VectorSize(new_vector) = GridSize(grid);                            \
                                                                        \
    /* Allocate vector data */                                          \
    Grid *vgrid = VectorGrid(new_vector);                               \
    for (i = 0; i < 10; i++)                                            \
      FreeCommPkg(VectorCommPkg(new_vector, i));                        \
                                                                        \
    ForSubgridI(i, GridSubgrids(vgrid))                                 \
    {                                                                   \
      Subvector *subvector = VectorSubvector(new_vector, i);            \
      data_size = SubvectorNX(subvector) * SubvectorNY(subvector) * SubvectorNZ(subvector); \
      SubvectorDataSize(subvector) = data_size;                         \
      double *data = amps_CTAlloc(double, data_size);                   \
      VectorSubvector(new_vector, i)->allocated = TRUE;                 \
      SubvectorData(VectorSubvector(new_vector, i)) = data;             \
    }                                                                   \
    for (i = 0; i < 10; i++)                                            \
      VectorCommPkg(new_vector, i) = NewVectorCommPkg(new_vector, GridComputePkg(vgrid, i)); \
                                                                        \
    new_vector->type = vector_non_samrai;                               \
    var = new_vector;                                                   \
  }


#endif // _INLINE_DEFS_H
