#ifndef _INLINE_RICHARDS_INTERNAL_H
#define _INLINE_RICHARDS_INTERNAL_H

#define RICHARDS_INTERNAL_MODULE(module, problem, problem_data,\
                                 vector, matrix, time, pressure, ftype)\
  {                                                                     \
    GetModulePublicXtra(RichardsBCInternal, module, public_xtra);       \
    Vector *f = vector;                                                 \
    Matrix *A = matrix;                                                 \
    int fcn = ftype;                                                    \
    WellData         *well_data = ProblemDataWellData(problem_data);    \
    WellDataPhysical *well_data_physical;                               \
    WellDataValue    *well_data_value;                                  \
                                                                        \
    TimeCycleData    *time_cycle_data;                                  \
                                                                        \
    int num_conditions = (public_xtra->num_conditions);                 \
    int num_wells, total_num;                                           \
                                                                        \
    Grid             *grid = VectorGrid(pressure);                      \
                                                                        \
    SubgridArray     *internal_bc_subgrids = NULL;                      \
                                                                        \
    Subgrid          *subgrid, *subgrid_ind, *new_subgrid;              \
                                                                        \
    Subvector        *p_sub;                                            \
                                                                        \
                                                                        \
    double           *pp;                                               \
    double           *internal_bc_conditions = NULL;                    \
                                                                        \
    double dx, dy, dz;                                                  \
    double value;                                                       \
                                                                        \
    int ix, iy, iz;                                                     \
    int nx, ny, nz;                                                     \
    int rx, ry, rz;                                                     \
    int process;                                                        \
                                                                        \
    int i, j, k;                                                        \
    int grid_index, well, index;                                        \
    int cycle_number, interval_number;                                  \
    int ip, im;                                                         \
                                                                        \
                                                                        \
    if (num_conditions > 0)                                             \
    {                                                                   \
      internal_bc_subgrids = NewSubgridArray();                         \
      internal_bc_conditions = ctalloc(double, num_conditions);         \
                                                                        \
      for (i = 0; i < num_conditions; i++)                              \
      {                                                                 \
        switch ((public_xtra->type[i]))                                 \
        {                                                               \
          case 0:                                                       \
          {                                                             \
            GetDummyType(RichardsBCInternal, 0, (public_xtra->data[i]), dummy0); \
                                                                        \
            ix = IndexSpaceX((dummy0->xlocation), 0);                   \
            iy = IndexSpaceY((dummy0->ylocation), 0);                   \
            iz = IndexSpaceZ((dummy0->zlocation), 0);                   \
                                                                        \
            nx = 1;                                                     \
            ny = 1;                                                     \
            nz = 1;                                                     \
                                                                        \
            rx = 0;                                                     \
            ry = 0;                                                     \
            rz = 0;                                                     \
                                                                        \
            process = amps_Rank(amps_CommWorld);                        \
                                                                        \
            new_subgrid = NewSubgrid(ix, iy, iz,                        \
                                     nx, ny, nz,                        \
                                     rx, ry, rz,                        \
                                     process);                          \
                                                                        \
            AppendSubgrid(new_subgrid, internal_bc_subgrids);           \
                                                                        \
            internal_bc_conditions[i] = (dummy0->value);                \
                                                                        \
            break;                                                      \
          }                                                             \
        }                                                               \
      }                                                                 \
    }                                                                   \
                                                                        \
                                                                        \
    num_wells = WellDataNumPressWells(well_data);                       \
    total_num = num_conditions + num_wells;                             \
                                                                        \
    if ((num_conditions > 0) || (num_wells > 0))                        \
    {                                                                   \
      /* Set explicit pressure assignments*/                            \
                                                                        \
      for (grid_index = 0; grid_index < GridNumSubgrids(grid); grid_index++) \
      {                                                                 \
        subgrid = GridSubgrid(grid, grid_index);                        \
                                                                        \
        p_sub = VectorSubvector(pressure, grid_index);                  \
        pp = SubvectorData(p_sub);                                      \
                                                                        \
                                                                        \
        for (index = 0; index < total_num; index++)                     \
        {                                                               \
          if (index < num_conditions)                                   \
          {                                                             \
            subgrid_ind = SubgridArraySubgrid(internal_bc_subgrids, index); \
            value = internal_bc_conditions[index];                      \
          }                                                             \
          else                                                          \
          {                                                             \
            well = index - num_conditions;                              \
            time_cycle_data = WellDataTimeCycleData(well_data);         \
            well_data_physical = WellDataPressWellPhysical(well_data, well); \
            cycle_number = WellDataPhysicalCycleNumber(well_data_physical); \
            interval_number =                                           \
                              TimeCycleDataComputeIntervalNumber(problem, time, \
                                                                 time_cycle_data, \
                                                                 cycle_number); \
            well_data_value =                                           \
                              WellDataPressWellIntervalValue(well_data, well, \
                                                             interval_number); \
            subgrid_ind = WellDataPhysicalSubgrid(well_data_physical);  \
            value = WellDataValuePhaseValue(well_data_value, 0);        \
          }                                                             \
                                                                        \
          ix = SubgridIX(subgrid_ind);                                  \
          iy = SubgridIY(subgrid_ind);                                  \
          iz = SubgridIZ(subgrid_ind);                                  \
                                                                        \
          nx = SubgridNX(subgrid_ind);                                  \
          ny = SubgridNY(subgrid_ind);                                  \
          nz = SubgridNZ(subgrid_ind);                                  \
                                                                        \
          dx = SubgridDX(subgrid_ind);                                  \
          dy = SubgridDY(subgrid_ind);                                  \
          dz = SubgridDZ(subgrid_ind);                                  \
                                                                        \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            Subvector *f_sub = VectorSubvector(f, grid_index);          \
            double *fp = SubvectorData(f_sub);                          \
                                                                        \
            BoxLoopI0(i, j, k,                                          \
                      ix, iy, iz, nx, ny, nz,                           \
            {                                                           \
              /* Need to check if i,j,k is part of this subgrid or not */ \
              if (((i >= SubgridIX(subgrid)) &&                         \
                   (i < SubgridIX(subgrid) + SubgridNX(subgrid))) &&    \
                  ((j >= SubgridIY(subgrid)) &&                         \
                   (j < SubgridIY(subgrid) + SubgridNY(subgrid))) &&    \
                  ((k >= SubgridIZ(subgrid)) &&                         \
                   (k < SubgridIZ(subgrid) + SubgridNZ(subgrid))))      \
              {                                                         \
                ip = SubvectorEltIndex(f_sub, i, j, k);                 \
                fp[ip] = pp[ip] - value;                                \
              }                                                         \
            });                                                         \
          }                                                             \
          else if (fcn == CALCDER)                                      \
          {                                                             \
            Submatrix        *A_sub = MatrixSubmatrix(A, grid_index);   \
                                                                        \
            double *cp = SubmatrixStencilData(A_sub, 0);                \
            double *wp = SubmatrixStencilData(A_sub, 1);                \
            double *ep = SubmatrixStencilData(A_sub, 2);                \
            double *sp = SubmatrixStencilData(A_sub, 3);                \
            double *np = SubmatrixStencilData(A_sub, 4);                \
            double *lp = SubmatrixStencilData(A_sub, 5);                \
            double *up = SubmatrixStencilData(A_sub, 6);                \
                                                                        \
            BoxLoopI0(i, j, k,                                          \
                      ix, iy, iz, nx, ny, nz,                           \
            {                                                           \
              /* Need to check if i,j,k is part of this subgrid or not */ \
              if (((i >= SubgridIX(subgrid)) &&                         \
                   (i < SubgridIX(subgrid) + SubgridNX(subgrid))) &&    \
                  ((j >= SubgridIY(subgrid)) &&                         \
                   (j < SubgridIY(subgrid) + SubgridNY(subgrid))) &&    \
                  ((k >= SubgridIZ(subgrid)) &&                         \
                   (k < SubgridIZ(subgrid) + SubgridNZ(subgrid))))      \
              {                                                         \
                im = SubmatrixEltIndex(A_sub, i, j, k);                 \
                cp[im] = 1.0;                                           \
                wp[im] = 0.0;                                           \
                ep[im] = 0.0;                                           \
                sp[im] = 0.0;                                           \
                np[im] = 0.0;                                           \
                lp[im] = 0.0;                                           \
                up[im] = 0.0;                                           \
              }                                                         \
            });                                                         \
          }                                                             \
        }           /* End loop over conditions */                      \
      }             /* End loop over processor subgrids */              \
                                                                        \
                                                                        \
      if (num_conditions > 0)                                           \
      {                                                                 \
        FreeSubgridArray(internal_bc_subgrids);                         \
      }                                                                 \
      tfree(internal_bc_conditions);                                    \
    }               /* End if have well or internal pressure conditions */ \
  }

#endif // _INLINE_RICHARDS_INTERNAL_H
