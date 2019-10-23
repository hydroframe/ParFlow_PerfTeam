#ifndef _INLINE_DEFS_H
#define _INLINE_DEFS_H

#include "inline_saturation.h"

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


#define PHASE_DENSITY_MODULE(module, pd_phase, pressure, density, pressure_der, density_der, ftype) \
  {                                                                     \
    GetModulePublicXtra(PhaseDensity, module, public_xtra);             \
    int phase = pd_phase;                                               \
    Vector *phase_pressure = pressure;                                  \
    Vector *density_v = density;                                        \
    double *pressure_d = pressure_der;                                  \
    double *density_d = density_der;                                    \
    int fcn = ftype;                                                    \
    int sg;                                                             \
    switch ((public_xtra->type[phase]))                                 \
    {                                                                   \
      case 0:                                                           \
      {                                                                 \
        double constant;                                                \
        GetDummyType(PhaseDensity, 0, (public_xtra->data[phase]), dummy); \
        constant = (dummy->constant);                                   \
                                                                        \
        if (density_v != NULL)                                          \
        {                                                               \
          Grid *pd_grid = VectorGrid(density_v);                        \
          ForSubgridI(sg, GridSubgrids(pd_grid))                        \
          {                                                             \
            subgrid = GridSubgrid(pd_grid, sg);                         \
                                                                        \
            d_sub = VectorSubvector(density_v, sg);                     \
                                                                        \
            ix = SubgridIX(subgrid) - 1;                                \
            iy = SubgridIY(subgrid) - 1;                                \
            iz = SubgridIZ(subgrid) - 1;                                \
                                                                        \
            nx = SubgridNX(subgrid) + 2;                                \
            ny = SubgridNY(subgrid) + 2;                                \
            nz = SubgridNZ(subgrid) + 2;                                \
                                                                        \
            nx_d = SubvectorNX(d_sub);                                  \
            ny_d = SubvectorNY(d_sub);                                  \
            nz_d = SubvectorNZ(d_sub);                                  \
                                                                        \
            dp = SubvectorElt(d_sub, ix, iy, iz);                       \
                                                                        \
            id = 0;                                                     \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,                \
                        id, nx_d, ny_d, nz_d, 1, 1, 1,                  \
              {                                                         \
                dp[id] = constant;                                      \
              });                                                       \
            }                                                           \
            else   /* fcn = CALCDER */                                  \
            {                                                           \
              BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz,                \
                        id, nx_d, ny_d, nz_d, 1, 1, 1,                  \
              {                                                         \
                dp[id] = 0.0;                                           \
              });                                                       \
            }   /* End if fcn */                                        \
          }    /* End subgrid loop */                                   \
        }      /* End if density_v is not NULL */                       \
        else                                                            \
        {                                                               \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            (*density_d) = constant;                                    \
          }                                                             \
          else  /* fcn = CALCDER */                                     \
          {                                                             \
            (*density_d) = 0.0;                                         \
          }                                                             \
        }      /* End else */                                           \
        break;                                                          \
      }                                                                 \
      case 1:                                                           \
      {                                                                 \
        double ref;                                                     \
        double comp;                                                    \
        GetDummyType(PhaseDensity, 1, (public_xtra->data[phase]), dummy); \
        ref = (dummy->reference_density);                               \
        comp = (dummy->compressibility_constant);                       \
                                                                        \
        if (density_v != NULL)                                          \
        {                                                               \
          Grid *pd_grid = VectorGrid(density_v);                        \
          ForSubgridI(sg, GridSubgrids(pd_grid))                        \
          {                                                             \
            subgrid = GridSubgrid(pd_grid, sg);                         \
                                                                        \
            p_sub = VectorSubvector(phase_pressure, sg);                \
            d_sub = VectorSubvector(density_v, sg);                     \
                                                                        \
            ix = SubgridIX(subgrid) - 1;                                \
            iy = SubgridIY(subgrid) - 1;                                \
            iz = SubgridIZ(subgrid) - 1;                                \
                                                                        \
            nx = SubgridNX(subgrid) + 2;                                \
            ny = SubgridNY(subgrid) + 2;                                \
            nz = SubgridNZ(subgrid) + 2;                                \
                                                                        \
            nx_p = SubvectorNX(p_sub);                                  \
            ny_p = SubvectorNY(p_sub);                                  \
            nz_p = SubvectorNZ(p_sub);                                  \
                                                                        \
            nx_d = SubvectorNX(d_sub);                                  \
            ny_d = SubvectorNY(d_sub);                                  \
            nz_d = SubvectorNZ(d_sub);                                  \
                                                                        \
            pp = SubvectorElt(p_sub, ix, iy, iz);                       \
            dp = SubvectorElt(d_sub, ix, iy, iz);                       \
                                                                        \
            ip = 0;                                                     \
            id = 0;                                                     \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,                \
                        ip, nx_p, ny_p, nz_p, 1, 1, 1,                  \
                        id, nx_d, ny_d, nz_d, 1, 1, 1,                  \
              {                                                         \
                dp[id] = ref * exp(pp[ip] * comp);                      \
              });                                                       \
            }                                                           \
            else          /* fcn = CALCDER */                           \
            {                                                           \
              BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,                \
                        ip, nx_p, ny_p, nz_p, 1, 1, 1,                  \
                        id, nx_d, ny_d, nz_d, 1, 1, 1,                  \
              {                                                         \
                dp[id] = comp * ref * exp(pp[ip] * comp);               \
              });                                                       \
            }                                                           \
          }                                                             \
        }                                                               \
        else                                                            \
        {                                                               \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            (*density_d) = ref * exp((*pressure_d) * comp);             \
          }                                                             \
          else                                                          \
          {                                                             \
            (*density_d) = comp * ref * exp((*pressure_d) * comp);      \
          }                                                             \
        }                                                               \
                                                                        \
        break;                                                          \
      }        /* End case 1 */                                         \
    }                                                                   \
  }

// proble, problem_data, vector f, matrix a, double time, vector pressure, int fcn
//#define RICHARDS_BC_INTERNAL_MODULE(vector, matrix, time, pressure, ftype)

#endif // _INLINE_DEFS_H
