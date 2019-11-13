#ifndef _INLINE_PHASE_DENSITY_H
#define _INLINE_PHASE_DENSITY_H

#define PHASE_DENSITY_MODULE(module, pd_phase, pressure, density, pressure_der, density_der, ftype) \
  {                                                                     \
    GetModulePublicXtra(PhaseDensity, module, public_xtra);             \
    int phase = pd_phase;                                               \
    Vector *phase_pressure = pressure;                                  \
    Vector *density_v = density;                                        \
    double *pressure_d = pressure_der;                                  \
    double *density_d = density_der;                                    \
    int fcn = ftype;                                                    \
    /* This is annoying but necessary to inline this body into other macro loops */ \
    int sg;                                                             \
    int id;                                                             \
    int nx_d;                                                           \
    int ny_d;                                                           \
    int nz_d;                                                           \
    int nx_p;                                                           \
    int ny_p;                                                           \
    int nz_p;                                                           \
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

#endif // _INLINE_PHASE_DENSITY_H
