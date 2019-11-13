#ifndef _INLINE_PHASE_SOURCE_H
#define _INLINE_PHASE_SOURCE_H

// Assumes Problem* problem and ProblemData* problem_data are in scope already
#define PHASE_SOURCE_MODULE(_module, _source, _phase, _time)            \
  {                                                                     \
    GetModulePublicXtra(PhaseSource, _module, public_xtra);             \
    Vector *phase_source = _source;                                     \
    int phase = _phase;                                                 \
    double time = _time;                                                \
                                                                        \
    WellData         *well_data = ProblemDataWellData(problem_data);    \
    WellDataPhysical *well_data_physical;                               \
    WellDataValue    *well_data_value;                                  \
                                                                        \
    TimeCycleData    *time_cycle_data;                                  \
                                                                        \
    Vector           *perm_x = ProblemDataPermeabilityX(problem_data);  \
    Vector           *perm_y = ProblemDataPermeabilityY(problem_data);  \
    Vector           *perm_z = ProblemDataPermeabilityZ(problem_data);  \
                                                                        \
    Grid             *grid = VectorGrid(phase_source);                  \
                                                                        \
    SubgridArray     *subgrids = GridSubgrids(grid);                    \
                                                                        \
    Subgrid          *subgrid, *well_subgrid, *tmp_subgrid;             \
    Subvector        *px_sub, *py_sub, *pz_sub, *ps_sub;                \
                                                                        \
    double           *data, *px, *py, *pz;                              \
                                                                        \
    int ix, iy, iz;                                                     \
    int nx, ny, nz;                                                     \
    int r;                                                              \
                                                                        \
    int nx_p, ny_p, nz_p;                                               \
    int nx_ps, ny_ps, nz_ps;                                            \
                                                                        \
    int is, i, j, k, ip, ips;                                           \
                                                                        \
    /* Locals associated with wells */                                  \
    int well;                                                           \
    int cycle_number, interval_number;                                  \
    double volume, flux, well_value;                                    \
    double weight = -FLT_MAX;                                           \
    double area_x, area_y, area_z, area_sum;                            \
    double avg_x, avg_y, avg_z;                                         \
    double dx, dy, dz;                                                  \
                                                                        \
    /*----------------------------------------------------------------------- \
     * Put in any user defined sources for this phase                   \
     *-----------------------------------------------------------------------*/ \
                                                                        \
    InitVectorAll(phase_source, 0.0);                                   \
                                                                        \
    switch ((public_xtra->type[phase]))                                 \
    {                                                                   \
      case 0:                                                           \
      {                                                                 \
        int num_regions;                                                \
        int     *region_indices;                                        \
        double  *values;                                                \
                                                                        \
        GrGeomSolid  *gr_solid;                                         \
        double value;                                                   \
        int ir;                                                         \
                                                                        \
        GetDummyType(PhaseSource, 0, (public_xtra->data[phase]), dummy0); \
                                                                        \
        num_regions = (dummy0->num_regions);                            \
        region_indices = (dummy0->region_indices);                      \
        values = (dummy0->values);                                      \
                                                                        \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
          value = values[ir];                                           \
                                                                        \
          ForSubgridI(is, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, is);                \
            ps_sub = VectorSubvector(phase_source, is);                 \
                                                                        \
            ix = SubgridIX(subgrid);                                    \
            iy = SubgridIY(subgrid);                                    \
            iz = SubgridIZ(subgrid);                                    \
                                                                        \
            nx = SubgridNX(subgrid);                                    \
            ny = SubgridNY(subgrid);                                    \
            nz = SubgridNZ(subgrid);                                    \
                                                                        \
            /* RDF: assume resolution is the same in all 3 directions */ \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            data = SubvectorData(ps_sub);                               \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
                                                                        \
              data[ips] = value;                                        \
            });                                                         \
          }                                                             \
        }                                                               \
                                                                        \
        break;                                                          \
      }  /* End case 0 */                                               \
                                                                        \
      case 1:                                                           \
      {                                                                 \
        GrGeomSolid  *gr_domain;                                        \
        double x, y, z;                                                 \
        int function_type;                                              \
                                                                        \
        GetDummyType(PhaseSource, 1, (public_xtra->data[phase]), dummy1);\
                                                                        \
        function_type = (dummy1->function_type);                        \
                                                                        \
        gr_domain = ProblemDataGrDomain(problem_data);                  \
                                                                        \
        ForSubgridI(is, subgrids)                                       \
        {                                                               \
          subgrid = SubgridArraySubgrid(subgrids, is);                  \
          ps_sub = VectorSubvector(phase_source, is);                   \
                                                                        \
          ix = SubgridIX(subgrid);                                      \
          iy = SubgridIY(subgrid);                                      \
          iz = SubgridIZ(subgrid);                                      \
                                                                        \
          nx = SubgridNX(subgrid);                                      \
          ny = SubgridNY(subgrid);                                      \
          nz = SubgridNZ(subgrid);                                      \
                                                                        \
          /* RDF: assume resolution is the same in all 3 directions */  \
          r = SubgridRX(subgrid);                                       \
                                                                        \
          data = SubvectorData(ps_sub);                                 \
                                                                        \
          switch (function_type)                                        \
          {                                                             \
            case 1: /* p= x */                                          \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                x = RealSpaceX(i, SubgridRX(subgrid));                  \
                /* nonlinear case -div(p grad p) = f */                 \
                data[ips] = -1.0;                                       \
              });                                                       \
              break;                                                    \
            } /* End case 1 */                                          \
                                                                        \
            case 2: /* p= x+y+z */                                      \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                /* nonlinear case -div(p grad p) = f */                 \
                data[ips] = -3.0;                                       \
              });                                                       \
              break;                                                    \
            } /* End case 2 */                                          \
                                                                        \
            case 3: /* p= x^3y^2 + sinxy + 1 */                         \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                x = RealSpaceX(i, SubgridRX(subgrid));                  \
                y = RealSpaceY(j, SubgridRY(subgrid));                  \
                /* nonlinear case -div(p grad p) = f */                 \
                data[ips] = -pow((3 * x * x * y * y + y * cos(x * y)), 2) - pow((2 * x * x * x * y + x * cos(x * y)), 2) - (x * x * x * y * y + sin(x * y) + 1) * (6 * x * y * y + 2 * x * x * x - (x * x + y * y) * sin(x * y)); \
              });                                                       \
              break;                                                    \
            } /* End case 3 */                                          \
                                                                        \
            case 4: /* f for p = x^3y^4 + x^2 + sinxy cosy + 1 */       \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                x = RealSpaceX(i, SubgridRX(subgrid));                  \
                y = RealSpaceY(j, SubgridRY(subgrid));                  \
                z = RealSpaceZ(k, SubgridRZ(subgrid));                  \
                                                                        \
                data[ips] = -pow(3 * x * x * pow(y, 4) + 2 * x + y * cos(x * y) * cos(y), 2) - pow(4 * x * x * x * y * y * y + x * cos(x * y) * cos(y) - sin(x * y) * sin(y), 2) - (x * x * x * pow(y, 4) + x * x + sin(x * y) * cos(y) + 1) * (6 * x * pow(y, 4) + 2 - (x * x + y * y + 1) * sin(x * y) * cos(y) + 12 * x * x * x * y * y - 2 * x * cos(x * y) * sin(y)); \
              });                                                       \
              break;                                                    \
            } /* End case 4 */                                          \
                                                                        \
            case 5: /* f = xyz-y^2z^2t^2-x^2z^2t^2-x^2y^2t^2 (p=xyzt+1)*/ \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                x = RealSpaceX(i, SubgridRX(subgrid));                  \
                y = RealSpaceY(j, SubgridRY(subgrid));                  \
                z = RealSpaceZ(k, SubgridRZ(subgrid));                  \
                                                                        \
                data[ips] = x * y * z - time * time * (y * y * z * z + x * x * z * z + x * x * y * y); \
              });                                                       \
              break;                                                    \
            } /* End case 5 */                                          \
                                                                        \
            case 6: /* f = xyz-y^2z^2t^2-2x^2z^2t^2-3x^2y^2t^2 (p=xyzt+1, \
                     *                                          K=(1; 2; 3) )*/ \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                x = RealSpaceX(i, SubgridRX(subgrid));                  \
                y = RealSpaceY(j, SubgridRY(subgrid));                  \
                z = RealSpaceZ(k, SubgridRZ(subgrid));                  \
                                                                        \
                data[ips] = x * y * z                                   \
                            - time * time * (y * y * z * z + x * x * z * z * 2.0 + x * x * y * y * 3.0); \
              });                                                       \
              break;                                                    \
            } /* End case 6 */                                          \
          }   /* End switch statement on function_types */              \
        }     /* End subgrid loop */                                    \
                                                                        \
        break;                                                          \
      }  /* End case 1 for input types */                               \
    }    /* End switch statement for input types */                     \
                                                                        \
    /*----------------------------------------------------------------------- \
     * Put in any flux wells from the well package                      \
     *-----------------------------------------------------------------------*/ \
                                                                        \
    if (WellDataNumFluxWells(well_data) > 0)                            \
    {                                                                   \
      time_cycle_data = WellDataTimeCycleData(well_data);               \
                                                                        \
      for (well = 0; well < WellDataNumFluxWells(well_data); well++)    \
      {                                                                 \
        well_data_physical = WellDataFluxWellPhysical(well_data, well); \
        cycle_number = WellDataPhysicalCycleNumber(well_data_physical); \
                                                                        \
        interval_number = TimeCycleDataComputeIntervalNumber(problem, time, time_cycle_data, cycle_number); \
                                                                        \
        well_data_value = WellDataFluxWellIntervalValue(well_data, well, interval_number); \
                                                                        \
        well_subgrid = WellDataPhysicalSubgrid(well_data_physical);     \
                                                                        \
        well_value = 0.0;                                               \
        if (WellDataPhysicalAction(well_data_physical) == INJECTION_WELL) \
        {                                                               \
          well_value = WellDataValuePhaseValue(well_data_value, phase); \
        }                                                               \
        else if (WellDataPhysicalAction(well_data_physical)             \
                 == EXTRACTION_WELL)                                    \
                                                                        \
        {                                                               \
          well_value = -WellDataValuePhaseValue(well_data_value, phase); \
        }                                                               \
                                                                        \
        volume = WellDataPhysicalSize(well_data_physical);              \
        flux = well_value / volume;                                     \
                                                                        \
        avg_x = WellDataPhysicalAveragePermeabilityX(well_data_physical); \
        avg_y = WellDataPhysicalAveragePermeabilityY(well_data_physical); \
        avg_z = WellDataPhysicalAveragePermeabilityZ(well_data_physical); \
                                                                        \
        ForSubgridI(is, subgrids)                                       \
        {                                                               \
          subgrid = SubgridArraySubgrid(subgrids, is);                  \
                                                                        \
          px_sub = VectorSubvector(perm_x, is);                         \
          py_sub = VectorSubvector(perm_y, is);                         \
          pz_sub = VectorSubvector(perm_z, is);                         \
                                                                        \
          ps_sub = VectorSubvector(phase_source, is);                   \
                                                                        \
          nx_p = SubvectorNX(ps_sub);                                   \
          ny_p = SubvectorNY(ps_sub);                                   \
          nz_p = SubvectorNZ(ps_sub);                                   \
                                                                        \
          nx_ps = SubvectorNX(ps_sub);                                  \
          ny_ps = SubvectorNY(ps_sub);                                  \
          nz_ps = SubvectorNZ(ps_sub);                                  \
                                                                        \
          /*  Get the intersection of the well with the subgrid  */     \
          if ((tmp_subgrid = IntersectSubgrids(subgrid, well_subgrid))) \
          {                                                             \
            /*  If an intersection;  loop over it, and insert value  */ \
            ix = SubgridIX(tmp_subgrid);                                \
            iy = SubgridIY(tmp_subgrid);                                \
            iz = SubgridIZ(tmp_subgrid);                                \
                                                                        \
            nx = SubgridNX(tmp_subgrid);                                \
            ny = SubgridNY(tmp_subgrid);                                \
            nz = SubgridNZ(tmp_subgrid);                                \
                                                                        \
            dx = SubgridDX(tmp_subgrid);                                \
            dy = SubgridDY(tmp_subgrid);                                \
            dz = SubgridDZ(tmp_subgrid);                                \
                                                                        \
            area_x = dy * dz;                                           \
            area_y = dx * dz;                                           \
            area_z = dx * dy;                                           \
            area_sum = area_x + area_y + area_z;                        \
                                                                        \
            px = SubvectorElt(px_sub, ix, iy, iz);                      \
            py = SubvectorElt(py_sub, ix, iy, iz);                      \
            pz = SubvectorElt(pz_sub, ix, iy, iz);                      \
                                                                        \
            data = SubvectorElt(ps_sub, ix, iy, iz);                    \
                                                                        \
            ip = 0;                                                     \
            ips = 0;                                                    \
            BoxLoopI2(i, j, k, ix, iy, iz, nx, ny, nz,                  \
                      ip, nx_p, ny_p, nz_p, 1, 1, 1,                    \
                      ips, nx_ps, ny_ps, nz_ps, 1, 1, 1,                \
            {                                                           \
              if (WellDataPhysicalMethod(well_data_physical)            \
                  == FLUX_STANDARD)                                     \
              {                                                         \
                weight = 1.0;                                           \
              }                                                         \
              else if (WellDataPhysicalMethod(well_data_physical)       \
                       == FLUX_WEIGHTED)                                \
              {                                                         \
                weight = (px[ip] / avg_x) * (area_x / area_sum)         \
                         + (py[ip] / avg_y) * (area_y / area_sum)       \
                         + (pz[ip] / avg_z) * (area_z / area_sum);      \
              }                                                         \
              else if (WellDataPhysicalMethod(well_data_physical)       \
                       == FLUX_PATTERNED)                               \
              {                                                         \
                weight = 0.0;                                           \
              }                                                         \
              data[ips] += weight * flux;                               \
            });                                                         \
                                                                        \
            /* done with this temporay subgrid */                       \
            FreeSubgrid(tmp_subgrid);                                   \
          }                                                             \
        }                                                               \
      }                                                                 \
    }  /* End well data */                                              \
}

#endif // _INLINE_PHASE_SOURCE_H
