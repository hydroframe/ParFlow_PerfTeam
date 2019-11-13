#ifndef _INLINE_BC_PRESSURE_H
#define _INLINE_BC_PRESSURE_H

#define BC_PRESSURE_MODULE(_module, bc_struct, grid, gr_domain, time)   \
  {                                                                     \
    GetModulePublicXtra(BCPressure, _module, public_xtra);              \
    GetModuleInstanceXtra(BCPressure, _module, instance_xtra);          \
                                                                        \
    PFModule       *phase_density = (instance_xtra->phase_density);     \
                                                                        \
    BCPressureData *bc_pressure_data = ProblemDataBCPressureData(problem_data); \
                                                                        \
    TimeCycleData  *time_cycle_data;                                    \
                                                                        \
    int num_phases = (public_xtra->num_phases);                         \
                                                                        \
    Problem        *problem = (instance_xtra->problem);                 \
                                                                        \
    SubgridArray   *subgrids = GridSubgrids(grid);                      \
                                                                        \
    Subgrid        *subgrid;                                            \
                                                                        \
    Vector      *z_mult = ProblemDataZmult(problem_data);               \
    Vector      *rsz = ProblemDataRealSpaceZ(problem_data);             \
    Subvector   *z_mult_sub;                                            \
    Subvector   *rsz_sub;                                               \
    double      *z_mult_dat;                                            \
    double      *rsz_dat;                                               \
                                                                        \
    double       ***values;                                             \
                                                                        \
    double         *patch_values;                                       \
    int patch_values_size;                                              \
                                                                        \
    int            *fdir;                                               \
                                                                        \
    int num_patches;                                                    \
    int ipatch, is, i, j, k, ival, ips, phase;                          \
    int cycle_number, interval_number;                                  \
                                                                        \
    num_patches = BCPressureDataNumPatches(bc_pressure_data);           \
                                                                        \
    if (num_patches > 0)                                                \
    {                                                                   \
      time_cycle_data = BCPressureDataTimeCycleData(bc_pressure_data);  \
                                                                        \
                                                                        \
      bc_struct = NewBCStruct(subgrids, gr_domain,                      \
                              num_patches,                              \
                              BCPressureDataPatchIndexes(bc_pressure_data), \
                              BCPressureDataBCTypes(bc_pressure_data),  \
                              NULL);                                    \
                                                                        \
                                                                        \
      values = ctalloc(double **, num_patches);                         \
      BCStructValues(bc_struct) = values;                               \
                                                                        \
      for (ipatch = 0; ipatch < num_patches; ipatch++)                  \
      {                                                                 \
        values[ipatch] = ctalloc(double *, SubgridArraySize(subgrids)); \
                                                                        \
        cycle_number = BCPressureDataCycleNumber(bc_pressure_data, ipatch); \
        interval_number = TimeCycleDataComputeIntervalNumber(           \
                                                             problem, time, time_cycle_data, cycle_number); \
                                                                        \
        switch (BCPressureDataType(bc_pressure_data, ipatch))           \
        {                                                               \
          case DirEquilRefPatch:                                        \
          {                                                             \
                                                                        \
            GeomSolid       *ref_solid;                                 \
            double z, dz2, dtmp;                                        \
            double offset, interface_press, interface_den;              \
            double ref_den, ref_press, nonlin_resid;                    \
            double density_der, density, fcn_val;                       \
            double height;                                              \
            double gravity = -ProblemGravity(problem);                  \
                                                                        \
            int ref_patch;                                              \
            int max_its = 10;                                           \
            int iterations;                                             \
            int ix, iy, iz, nx, ny, nz, r, iel;                         \
                                                                        \
            double         **elevations;                                \
                                                                        \
            GetBCPressureTypeStruct(DirEquilRefPatch, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
            if (instance_xtra->elevations == NULL)                      \
            {                                                           \
              instance_xtra->elevations = ctalloc(double **, num_patches); \
              instance_xtra->problem_data = problem_data;               \
              instance_xtra->grid = grid;                               \
            }                                                           \
                                                                        \
            if (instance_xtra->elevations[ipatch] == NULL)              \
            {                                                           \
              ref_solid = ProblemDataSolid(problem_data,                \
                                           DirEquilRefPatchRefSolid(interval_data)); \
              ref_patch = DirEquilRefPatchRefPatch(interval_data);      \
                                                                        \
              /* Calculate elevations at (x,y) points on reference patch. */ \
              instance_xtra->elevations[ipatch] = CalcElevations(ref_solid, ref_patch, subgrids, problem_data); \
            }                                                           \
                                                                        \
            elevations = instance_xtra->elevations[ipatch];             \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              z_mult_sub = VectorSubvector(z_mult, is);                 \
              rsz_sub = VectorSubvector(rsz, is);                       \
              z_mult_dat = SubvectorData(z_mult_sub);                   \
              rsz_dat = SubvectorData(rsz_sub);                         \
                                                                        \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              ix = SubgridIX(subgrid);                                  \
              iy = SubgridIY(subgrid);                                  \
              iz = SubgridIZ(subgrid);                                  \
                                                                        \
              nx = SubgridNX(subgrid);                                  \
              ny = SubgridNY(subgrid);                                  \
              nz = SubgridNZ(subgrid);                                  \
                                                                        \
              /* RDF: assume resolution is the same in all 3 directions */ \
              r = SubgridRX(subgrid);                                   \
                                                                        \
              dz2 = SubgridDZ(subgrid) * 0.5;                           \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct,         \
                                ipatch, is,                             \
              {                                                         \
                ref_press = DirEquilRefPatchValue(interval_data);       \
                /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, &ref_press, &ref_den, CALCFCN) */ \
                PFModuleInvokeType(PhaseDensityInvoke, phase_density,   \
                                   (0, NULL, NULL, &ref_press, &ref_den, \
                                    CALCFCN));                          \
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);           \
                z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];     \
                iel = (i - ix) + (j - iy) * nx;                         \
                fcn_val = 0.0;                                          \
                nonlin_resid = 1.0;                                     \
                iterations = -1;                                        \
                                                                        \
                                                                        \
                while ((nonlin_resid > 1.0E-6) && (iterations < max_its)) \
                {                                                       \
                  if (iterations > -1)                                  \
                  {                                                     \
                    /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, */ \
                    /*                      &patch_values[ival], &density_der, */ \
                    /*                      CALCDER); */                \
                    PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                       (0, NULL, NULL, &patch_values[ival], \
                                        &density_der, CALCDER));        \
                    dtmp = 1.0 - 0.5 * density_der * gravity            \
                           * (z - elevations[is][iel]);                 \
                    patch_values[ival] = patch_values[ival] - fcn_val / dtmp; \
                  }                                                     \
                  else                                                  \
                  {                                                     \
                    patch_values[ival] = ref_press;                     \
                  }                                                     \
                                                                        \
                  /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, */ \
                  /*                      &patch_values[ival], &density, */ \
                  /*                      CALCFCN); */                  \
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                     (0, NULL, NULL, &patch_values[ival], \
                                      &density, CALCFCN));              \
                                                                        \
                  fcn_val = patch_values[ival] - ref_press              \
                            - 0.5 * (density + ref_den) * gravity       \
                            * (z - elevations[is][iel]);                \
                  nonlin_resid = fabs(fcn_val);                         \
                                                                        \
                  iterations++;                                         \
                }            /* End of while loop */                    \
                                                                        \
                                                                        \
                                                                        \
                for (phase = 1; phase < num_phases; phase++)            \
                {                                                       \
                  interface_press = DirEquilRefPatchValueAtInterface(   \
                                                                     interval_data, phase); \
                  /* PHASE_DENSITY_MODULE(phase_density, (phase-1), NULL, NULL, */ \
                  /*                      &interface_press, &interface_den, */ \
                  /*                      CALCFCN); */                  \
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                     (phase - 1, NULL, NULL, &interface_press, \
                                      &interface_den, CALCFCN));        \
                                                                        \
                  offset = (interface_press - ref_press)                \
                           / (0.5 * (interface_den + ref_den) * gravity); \
                  ref_press = interface_press;                          \
                                                                        \
                  /* PHASE_DENSITY_MODULE(phase_density, phase, NULL, NULL, */ \
                  /*                      &ref_press, &ref_den, */      \
                  /*                      CALCFCN); */                  \
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                     (phase, NULL, NULL, &ref_press, &ref_den, \
                                      CALCFCN));                        \
                                                                        \
                                                                        \
                  if (patch_values[ival] < interface_press)             \
                  {                                                     \
                    height = elevations[is][iel];                       \
                    nonlin_resid = 1.0;                                 \
                    iterations = -1;                                    \
                    while ((nonlin_resid > 1.0E-6) && (iterations < max_its)) \
                    {                                                   \
                      if (iterations > -1)                              \
                      {                                                 \
                        /* PHASE_DENSITY_MODULE(phase_density, phase, NULL, NULL, */ \
                        /*                      &patch_values[ival], &density_der, */ \
                        /*                      CALCDER); */            \
                        PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                           (phase, NULL, NULL, &patch_values[ival], \
                                            &density_der, CALCDER));    \
                                                                        \
                        dtmp = 1.0 - 0.5 * density_der * gravity        \
                               * (z - height);                          \
                        patch_values[ival] = patch_values[ival]         \
                                             - fcn_val / dtmp;          \
                      }                                                 \
                      else                                              \
                      {                                                 \
                        height = height + offset;                       \
                        patch_values[ival] = ref_press;                 \
                      }                                                 \
                                                                        \
                      /* PHASE_DENSITY_MODULE(phase_density, phase, NULL, NULL, */ \
                      /*                      &patch_values[ival], &density, */ \
                      /*                      CALCFCN); */              \
                      PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                         (phase, NULL, NULL,            \
                                          &patch_values[ival], &density, \
                                          CALCFCN));                    \
                                                                        \
                      fcn_val = patch_values[ival] - ref_press          \
                                - 0.5 * (density + ref_den)             \
                                * gravity * (z - height);               \
                      nonlin_resid = fabs(fcn_val);                     \
                                                                        \
                      iterations++;                                     \
                    }              /* End of while loop */              \
                  }                /* End if above interface */         \
                }                  /* End phase loop */                 \
              });                  /* End BCStructPatchLoop body */     \
            }                      /* End subgrid loop */               \
                                                                        \
                                                                        \
            break;                                                      \
          } /* End DirEquilRefPatch */                                  \
                                                                        \
          case DirEquilPLinear:                                         \
          {                                                             \
                                                                        \
            int num_points;                                             \
            int ip;                                                     \
                                                                        \
            double x, y, z, dx2, dy2, dz2;                              \
            double unitx, unity, line_min, line_length, xy, slope;      \
                                                                        \
            double dtmp, offset, interface_press, interface_den;        \
            double ref_den, ref_press, nonlin_resid;                    \
            double density_der, density, fcn_val;                       \
            double height;                                              \
            double gravity = -ProblemGravity(problem);                  \
                                                                        \
            int max_its = 10;                                           \
            int iterations;                                             \
                                                                        \
            GetBCPressureTypeStruct(DirEquilPLinear, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              z_mult_sub = VectorSubvector(z_mult, is);                 \
              rsz_sub = VectorSubvector(rsz, is);                       \
              z_mult_dat = SubvectorData(z_mult_sub);                   \
              rsz_dat = SubvectorData(rsz_sub);                         \
                                                                        \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              dx2 = SubgridDX(subgrid) / 2.0;                           \
              dy2 = SubgridDY(subgrid) / 2.0;                           \
              dz2 = SubgridDZ(subgrid) / 2.0;                           \
                                                                        \
              /* compute unit direction vector for piecewise linear line */ \
              unitx = DirEquilPLinearXUpper(interval_data)              \
                      - DirEquilPLinearXLower(interval_data);           \
              unity = DirEquilPLinearYUpper(interval_data)              \
                      - DirEquilPLinearYLower(interval_data);           \
              line_length = sqrt(unitx * unitx + unity * unity);        \
              unitx /= line_length;                                     \
              unity /= line_length;                                     \
              line_min = DirEquilPLinearXLower(interval_data) * unitx   \
                         + DirEquilPLinearYLower(interval_data) * unity; \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2;  \
                y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2;  \
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);           \
                z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];     \
                                                                        \
                /* project center of BC face onto piecewise line */     \
                xy = (x * unitx + y * unity - line_min) / line_length;  \
                                                                        \
                /* find two neighboring points */                       \
                ip = 1;                                                 \
                num_points = DirEquilPLinearNumPoints(interval_data);   \
                for (; ip < (num_points - 1); ip++)                     \
                {                                                       \
                  if (xy < DirEquilPLinearPoint(interval_data, ip))     \
                    break;                                              \
                }                                                       \
                                                                        \
                /* compute the slope */                                 \
                slope = ((DirEquilPLinearValue(interval_data, ip)       \
                          - DirEquilPLinearValue(interval_data, (ip - 1))) \
                         / (DirEquilPLinearPoint(interval_data, ip)     \
                            - DirEquilPLinearPoint(interval_data, (ip - 1)))); \
                                                                        \
                ref_press = DirEquilPLinearValue(interval_data, ip - 1) \
                            + slope * (xy - DirEquilPLinearPoint(interval_data, ip - 1)); \
                                                                        \
                /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, */ \
                /*                      &ref_press, &ref_den, */        \
                /*                      CALCFCN); */                    \
                PFModuleInvokeType(PhaseDensityInvoke, phase_density,   \
                                   (0, NULL, NULL, &ref_press, &ref_den, \
                                    CALCFCN));                          \
                                                                        \
                fcn_val = 0.0;                                          \
                nonlin_resid = 1.0;                                     \
                iterations = -1;                                        \
                                                                        \
                                                                        \
                while ((nonlin_resid > 1.0E-6) && (iterations < max_its)) \
                {                                                       \
                  if (iterations > -1)                                  \
                  {                                                     \
                    /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, */ \
                    /*                      &patch_values[ival], &density_der, */ \
                    /*                      CALCDER); */                \
                    PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                       (0, NULL, NULL, &patch_values[ival], \
                                        &density_der, CALCDER));        \
                                                                        \
                    dtmp = 1.0 - 0.5 * density_der * gravity * z;       \
                    patch_values[ival] = patch_values[ival] - fcn_val / dtmp; \
                  }                                                     \
                  else                                                  \
                  {                                                     \
                    patch_values[ival] = ref_press;                     \
                  }                                                     \
                                                                        \
                  /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, */ \
                  /*                      &patch_values[ival], &density, */ \
                  /*                      CALCFCN); */                  \
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                     (0, NULL, NULL, &patch_values[ival], \
                                      &density, CALCFCN));              \
                                                                        \
                  fcn_val = patch_values[ival] - ref_press              \
                            - 0.5 * (density + ref_den) * gravity * z;  \
                  nonlin_resid = fabs(fcn_val);                         \
                                                                        \
                  iterations++;                                         \
                }            /* End of while loop */                    \
                                                                        \
                                                                        \
                for (phase = 1; phase < num_phases; phase++)            \
                {                                                       \
                  interface_press = DirEquilPLinearValueAtInterface(    \
                                                                    interval_data, phase); \
                                                                        \
                                                                        \
                  /* PHASE_DENSITY_MODULE(phase_density, (phase-1), NULL, NULL, */ \
                  /*                      &interface_press, &interface_den, */ \
                  /*                      CALCFCN); */                  \
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                     (phase - 1, NULL, NULL, &interface_press, \
                                      &interface_den, CALCFCN));        \
                                                                        \
                  offset = (interface_press - ref_press)                \
                           / (0.5 * (interface_den + ref_den) * gravity); \
                  ref_press = interface_press;                          \
                                                                        \
                  /* PHASE_DENSITY_MODULE(phase_density, phase, NULL, NULL, */ \
                  /*                      &ref_press, &ref_den, */      \
                  /*                      CALCFCN); */                  \
                  PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                     (phase, NULL, NULL, &ref_press, &ref_den, \
                                      CALCFCN));                        \
                                                                        \
                                                                        \
                  if (patch_values[ival] < interface_press)             \
                  {                                                     \
                    height = 0.0;                                       \
                    nonlin_resid = 1.0;                                 \
                    iterations = -1;                                    \
                    while ((nonlin_resid > 1.0E-6) && (iterations < max_its)) \
                    {                                                   \
                      if (iterations > -1)                              \
                      {                                                 \
                                                                        \
                        /* PHASE_DENSITY_MODULE(phase_density, phase, NULL, NULL, */ \
                        /*                      &patch_values[ival], &density_der, */ \
                        /*                      CALCDER); */            \
                        PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                           (phase, NULL, NULL, &patch_values[ival], \
                                            &density_der, CALCDER));    \
                                                                        \
                        dtmp = 1.0 - 0.5 * density_der * gravity * (z - height); \
                        patch_values[ival] = patch_values[ival]         \
                                             - fcn_val / dtmp;          \
                      }                                                 \
                      else                                              \
                      {                                                 \
                        height = height + offset;                       \
                        patch_values[ival] = ref_press;                 \
                      }                                                 \
                                                                        \
                      /* PHASE_DENSITY_MODULE(phase_density, phase, NULL, NULL, */ \
                      /*                      &patch_values[ival], &density, */ \
                      /*                      CALCFCN); */              \
                      PFModuleInvokeType(PhaseDensityInvoke, phase_density, \
                                         (phase, NULL, NULL,            \
                                          &patch_values[ival], &density, \
                                          CALCFCN));                    \
                                                                        \
                      fcn_val = patch_values[ival] - ref_press          \
                                - 0.5 * (density + ref_den) * gravity   \
                                * (z - height);                         \
                      nonlin_resid = fabs(fcn_val);                     \
                                                                        \
                      iterations++;                                     \
                    }              /* End of while loop */              \
                  }                /* End if above interface */         \
                }                  /* End phase loop */                 \
              });                  /* End BCStructPatchLoop body */     \
            }                                                           \
            break;                                                      \
          } /* End DirEquilPLinear */                                   \
                                                                        \
          case FluxConst:                                               \
          {                                                             \
            /* Constant flux rate value on patch */                     \
            double flux;                                                \
                                                                        \
            GetBCPressureTypeStruct(FluxConst, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
            flux = FluxConstValue(interval_data);                       \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values[ival] = flux;                              \
              });                                                       \
            }       /* End subgrid loop */                              \
            break;                                                      \
          } /* End FluxConst */                                         \
                                                                        \
          case FluxVolumetric:                                          \
          {                                                             \
            /* Constant volumetric flux value on patch */               \
            double dx, dy, dz;                                          \
            double area, volumetric_flux;                               \
                                                                        \
            GetBCPressureTypeStruct(FluxVolumetric, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              z_mult_sub = VectorSubvector(z_mult, is);                 \
              z_mult_dat = SubvectorData(z_mult_sub);                   \
                                                                        \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              dx = SubgridDX(subgrid);                                  \
              dy = SubgridDY(subgrid);                                  \
              dz = SubgridDZ(subgrid);                                  \
                                                                        \
              area = 0.0;                                               \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);           \
                /* primary direction x */                               \
                if (fdir[0])                                            \
                {                                                       \
                  area += dy * dz * z_mult_dat[ips];                    \
                }                                                       \
                /* primary direction y */                               \
                else if (fdir[1])                                       \
                {                                                       \
                  area += dx * dz * z_mult_dat[ips];                    \
                }                                                       \
                /* primary direction z */                               \
                else if (fdir[2])                                       \
                {                                                       \
                  area += dx * dy;                                      \
                }                                                       \
              });                                                       \
                                                                        \
              if (area > 0.0)                                           \
              {                                                         \
                volumetric_flux = FluxVolumetricValue(interval_data)    \
                                  / area;                               \
                BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                {                                                       \
                  patch_values[ival] = volumetric_flux;                 \
                });                                                     \
              }                                                         \
            }           /* End subgrid loop */                          \
            break;                                                      \
          } /* End FluxVolumetric */                                    \
                                                                        \
          case PressureFile:                                            \
          {                                                             \
            Vector          *tmp_vector;                                \
            Subvector       *subvector;                                 \
            char            *filename;                                  \
            double          *tmpp;                                      \
            int itmp;                                                   \
            double z, dz2;                                              \
            double density, dtmp;                                       \
                                                                        \
            double gravity = ProblemGravity(problem);                   \
                                                                        \
            /* Calculate density using dtmp as dummy argument. */       \
            dtmp = 0.0;                                                 \
                                                                        \
            /* PHASE_DENSITY_MODULE(phase_density, 0, NULL, NULL, */    \
            /*                      &dtmp, &density, */                 \
            /*                      CALCFCN); */                        \
            PFModuleInvokeType(PhaseDensityInvoke, phase_density,       \
                               (0, NULL, NULL, &dtmp, &density, CALCFCN)); \
                                                                        \
            GetBCPressureTypeStruct(PressureFile, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              z_mult_sub = VectorSubvector(z_mult, is);                 \
              rsz_sub = VectorSubvector(rsz, is);                       \
              z_mult_dat = SubvectorData(z_mult_sub);                   \
              rsz_dat = SubvectorData(rsz_sub);                         \
                                                                        \
              dz2 = SubgridDZ(subgrid) / 2.0;                           \
                                                                        \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              tmp_vector = NewVectorType(grid, 1, 0, vector_cell_centered); \
                                                                        \
              filename = PressureFileName(interval_data);               \
              ReadPFBinary(filename, tmp_vector);                       \
                                                                        \
              subvector = VectorSubvector(tmp_vector, is);              \
                                                                        \
              tmpp = SubvectorData(subvector);                          \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                ips = SubvectorEltIndex(z_mult_sub, i, j, k);           \
                z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips];     \
                itmp = SubvectorEltIndex(subvector, i, j, k);           \
                                                                        \
                patch_values[ival] = tmpp[itmp];     /*- density*gravity*z;*/ \
                /*last part taken out, very likely to be a bug)*/       \
              });                                                       \
                                                                        \
              FreeVector(tmp_vector);                                   \
            }             /* End subgrid loop */                        \
            break;                                                      \
          } /* End PressureFile */                                      \
                                                                        \
          case FluxFile:                                                \
          {                                                             \
            /* Read input fluxes from file (temporary) */               \
            Vector          *tmp_vector;                                \
            Subvector       *subvector;                                 \
            char            *filename;                                  \
            double          *tmpp;                                      \
            int itmp;                                                   \
                                                                        \
            GetBCPressureTypeStruct(FluxFile, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              tmp_vector = NewVectorType(grid, 1, 0, vector_cell_centered); \
                                                                        \
              filename = FluxFileName(interval_data);                   \
              ReadPFBinary(filename, tmp_vector);                       \
                                                                        \
              subvector = VectorSubvector(tmp_vector, is);              \
                                                                        \
              tmpp = SubvectorData(subvector);                          \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                itmp = SubvectorEltIndex(subvector, i, j, k);           \
                                                                        \
                patch_values[ival] = tmpp[itmp];                        \
              });                                                       \
                                                                        \
              FreeVector(tmp_vector);                                   \
            }         /* End subgrid loop */                            \
            break;                                                      \
          } /* End FluxFile */                                          \
                                                                        \
          case ExactSolution:                                           \
          {                                                             \
            /* Calculate pressure based on pre-defined functions */     \
            double x, y, z, dx2, dy2, dz2;                              \
            int fcn_type;                                               \
                                                                        \
            GetBCPressureTypeStruct(ExactSolution, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              z_mult_sub = VectorSubvector(z_mult, is);                 \
              rsz_sub = VectorSubvector(rsz, is);                       \
              z_mult_dat = SubvectorData(z_mult_sub);                   \
              rsz_dat = SubvectorData(rsz_sub);                         \
                                                                        \
              /* compute patch_values_size */                           \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              dx2 = SubgridDX(subgrid) / 2.0;                           \
              dy2 = SubgridDY(subgrid) / 2.0;                           \
              dz2 = SubgridDZ(subgrid) / 2.0;                           \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              fcn_type = ExactSolutionFunctionType(interval_data);      \
                                                                        \
              switch (fcn_type)                                         \
              {                                                         \
                case 1:  /* p = x */                                    \
                {                                                       \
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                  {                                                     \
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2; \
                                                                        \
                    patch_values[ival] = x;                             \
                  });                                                   \
                                                                        \
                  break;                                                \
                }     /* End case 1 */                                  \
                                                                        \
                case 2:  /* p = x + y + z */                            \
                {                                                       \
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                  {                                                     \
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2; \
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2; \
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);       \
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips]; \
                    patch_values[ival] = x + y + z;                     \
                  });                                                   \
                                                                        \
                  break;                                                \
                }     /* End case 2 */                                  \
                                                                        \
                case 3:  /* p = x^3y^2 + sinxy + 1*/                    \
                {                                                       \
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                  {                                                     \
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2; \
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2; \
                                                                        \
                    patch_values[ival] = x * x * x * y * y + sin(x * y) + 1; \
                  });                                                   \
                  break;                                                \
                }     /* End case 3 */                                  \
                                                                        \
                case 4:  /* p = x^3 y^4 + x^2 + sinxy cosy + 1 */       \
                {                                                       \
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                  {                                                     \
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2; \
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2; \
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);       \
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips]; \
                    patch_values[ival] = pow(x, 3) * pow(y, 4) + x * x + sin(x * y) * cos(y) + 1; \
                  });                                                   \
                  break;                                                \
                }     /* End case 4 */                                  \
                                                                        \
                case 5:  /* p = xyzt +1 */                              \
                {                                                       \
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                  {                                                     \
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2; \
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2; \
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);       \
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips]; \
                    patch_values[ival] = x * y * z * time + 1;          \
                  });                                                   \
                  break;                                                \
                }     /* End case 5 */                                  \
                                                                        \
                case 6:  /* p = xyzt +1 */                              \
                {                                                       \
                  BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
                  {                                                     \
                    x = RealSpaceX(i, SubgridRX(subgrid)) + fdir[0] * dx2; \
                    y = RealSpaceY(j, SubgridRY(subgrid)) + fdir[1] * dy2; \
                    ips = SubvectorEltIndex(z_mult_sub, i, j, k);       \
                    z = rsz_dat[ips] + fdir[2] * dz2 * z_mult_dat[ips]; \
                    patch_values[ival] = x * y * z * time + 1;          \
                  });                                                   \
                  break;                                                \
                }     /* End case 5 */                                  \
              }       /* End switch */                                  \
            }         /* End subgrid loop */                            \
            break;                                                      \
          } /* End ExactSolution */                                     \
                                                                        \
          case OverlandFlow:                                            \
          {                                                             \
            /* Constant "rainfall" rate value on patch */               \
            double flux;                                                \
                                                                        \
            GetBCPressureTypeStruct(FluxConst, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
                                                                        \
            flux = OverlandFlowValue(interval_data);                    \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values[ival] = flux;                              \
              });                                                       \
            }       /* End subgrid loop */                              \
            break;                                                      \
          } /* End OverlandFlow */                                      \
                                                                        \
          case OverlandFlowPFB:                                         \
          {                                                             \
            /* Read input fluxes from file (overland) */                \
            Vector          *tmp_vector;                                \
            Subvector       *subvector;                                 \
            /* double          *data; */                                \
            char            *filename;                                  \
            double          *tmpp;                                      \
            int itmp;                                                   \
            double dtmp;                                                \
                                                                        \
            GetBCPressureTypeStruct(OverlandFlowPFB, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
                                                                        \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              /* compute patch_values_size (this isn't really needed yet) */ \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              tmp_vector = NewVectorType(grid, 1, 0, vector_cell_centered); \
              /*data = ctalloc(double, SizeOfVector(tmp_vector));*/     \
              /* SetTempVectorData(tmp_vector, data); */                \
                                                                        \
              printf("reading overland file \n");                       \
              filename = OverlandFlowPFBFileName(interval_data);        \
              ReadPFBinary(filename, tmp_vector);                       \
                                                                        \
              subvector = VectorSubvector(tmp_vector, is);              \
                                                                        \
              tmpp = SubvectorData(subvector);                          \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                itmp = SubvectorEltIndex(subvector, i, j, k);           \
                                                                        \
                patch_values[ival] = tmpp[itmp];                        \
              });                                                       \
                                                                        \
              /* tfree(VectorData(tmp_vector)); */                      \
              FreeVector(tmp_vector);                                   \
            }              /* End subgrid loop */                       \
            break;                                                      \
          } /* End OverlandFlowPFB */                                   \
                                                                        \
          case SeepageFace:                                             \
          {                                                             \
            GetBCPressureTypeStruct(SeepageFace, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
            double flux;                                                \
                                                                        \
            flux = SeepageFaceValue(interval_data);                     \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values[ival] = flux;                              \
              });                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
          } /* End SeepageFace */                                       \
                                                                        \
          case OverlandKinematic:                                       \
          {                                                             \
            GetBCPressureTypeStruct(OverlandKinematic, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
            double flux;                                                \
                                                                        \
            flux = OverlandKinematicValue(interval_data);               \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values[ival] = flux;                              \
              });                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
          } /* End OverlandKinematic */                                 \
                                                                        \
          case OverlandDiffusive:                                       \
          {                                                             \
            GetBCPressureTypeStruct(OverlandDiffusive, interval_data, bc_pressure_data, \
                                    ipatch, interval_number);           \
            double flux;                                                \
                                                                        \
            flux = OverlandDiffusiveValue(interval_data);               \
            ForSubgridI(is, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, is);              \
                                                                        \
              patch_values_size = 0;                                    \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values_size++;                                    \
              });                                                       \
                                                                        \
              patch_values = ctalloc(double, patch_values_size);        \
              values[ipatch][is] = patch_values;                        \
                                                                        \
              BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is, \
              {                                                         \
                patch_values[ival] = flux;                              \
              });                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
          } /* End OverlandDiffusive */                                 \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }

#endif // _INLINE_BC_PRESSURE_H
