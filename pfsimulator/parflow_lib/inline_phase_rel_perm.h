#ifndef _INLINE_PHASE_REL_PERM_H
#define _INLINE_PHASE_REL_PERM_H

#include "problem_phase_rel_perm.h"

#define PHASERELPERM_MODULE(module, phase_rel_perm, phase_pressure, phase_density, gravity, problem_data, fcn) \
  {                                                                     \
    GetModulePublicXtra(PhaseRelPerm, module, public_xtra);             \
    Grid          *grid = VectorGrid(phase_rel_perm);                   \
    GrGeomSolid   *gr_solid;                                            \
    Subvector     *pr_sub;                                              \
    Subvector     *pp_sub;                                              \
    Subvector     *pd_sub;                                              \
    Subvector     *n_values_sub;                                        \
    Subvector     *alpha_values_sub;                                    \
                                                                        \
    double        *prdat, *ppdat, *pddat;                               \
    double        *n_values_dat, *alpha_values_dat;                     \
                                                                        \
    SubgridArray  *subgrids = GridSubgrids(grid);                       \
                                                                        \
    Subgrid       *subgrid;                                             \
                                                                        \
    int sg;                                                             \
                                                                        \
    int ix, iy, iz;                                                     \
    int nx, ny, nz;                                                     \
                                                                        \
    int i, j, k, r, ipr, ipp, ipd;                                      \
                                                                        \
    int n_index, alpha_index;                                           \
                                                                        \
    int num_regions, *region_indices;                                   \
    int ir, *fdir;                                                      \
                                                                        \
    BeginTiming(public_xtra->time_index);                               \
                                                                        \
    /* Initialize relative permeabilities to 0.0 */                     \
    InitVectorAll(phase_rel_perm, 0.0);                                 \
                                                                        \
    switch ((public_xtra->type))                                        \
    {                                                                   \
      case 0:   /* Constant relative permeability within regions */     \
      {                                                                 \
        double  *values;                                                \
        int ir;                                                         \
                                                                        \
        GetDummyType(PhaseRelPerm, 0, (public_xtra->data), dummy0);     \
                                                                        \
        num_regions = (dummy0->num_regions);                            \
        region_indices = (dummy0->region_indices);                      \
        values = (dummy0->values);                                      \
                                                                        \
        /* Compute rel perms for Dirichlet boundary conditions */       \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
                                                                        \
            ix = SubgridIX(subgrid);                                    \
            iy = SubgridIY(subgrid);                                    \
            iz = SubgridIZ(subgrid);                                    \
                                                                        \
            nx = SubgridNX(subgrid);                                    \
            ny = SubgridNY(subgrid);                                    \
            nz = SubgridNZ(subgrid);                                    \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                prdat[ipr] = values[ir];                                \
              });                                                       \
            }                                                           \
            else    /* fcn = CALCDER */                                 \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                prdat[ipr] = 0.0;                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End loop over regions */                           \
                                                                        \
        /* Compute rel perms inside regions */                          \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
                                                                        \
            ix = SubgridIX(subgrid) - 1;                                \
            iy = SubgridIY(subgrid) - 1;                                \
            iz = SubgridIZ(subgrid) - 1;                                \
                                                                        \
            nx = SubgridNX(subgrid) + 2;                                \
            ny = SubgridNY(subgrid) + 2;                                \
            nz = SubgridNZ(subgrid) + 2;                                \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                prdat[ipr] = values[ir];                                \
              });                                                       \
            }                                                           \
            else    /* fcn = CALCDER */                                 \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                prdat[ipr] = 0.0;                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End loop over regions */                           \
        break;                                                          \
      }        /* End case 0 */                                         \
                                                                        \
      case 1: /* Van Genuchten relative permeability */                 \
      {                                                                 \
        int data_from_file;                                             \
        double  *alphas, *ns, head;                                     \
        double alpha, n, m, opahn, ahnm1, coeff;                        \
                                                                        \
        Vector  *n_values, *alpha_values;                               \
                                                                        \
        GetDummyType(PhaseRelPerm, 1, (public_xtra->data), dummy1);     \
                                                                        \
        num_regions = (dummy1->num_regions);                            \
        region_indices = (dummy1->region_indices);                      \
        alphas = (dummy1->alphas);                                      \
        ns = (dummy1->ns);                                              \
        data_from_file = (dummy1->data_from_file);                      \
                                                                        \
        /* Compute rel perms for Dirichlet boundary conditions */       \
        if (data_from_file == 0)  /* alphas and ns given by region */   \
        {                                                               \
          for (ir = 0; ir < num_regions; ir++)                          \
          {                                                             \
            gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
            ForSubgridI(sg, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, sg);              \
                                                                        \
              pr_sub = VectorSubvector(phase_rel_perm, sg);             \
              pp_sub = VectorSubvector(phase_pressure, sg);             \
              pd_sub = VectorSubvector(phase_density, sg);              \
                                                                        \
              ix = SubgridIX(subgrid);                                  \
              iy = SubgridIY(subgrid);                                  \
              iz = SubgridIZ(subgrid);                                  \
                                                                        \
              nx = SubgridNX(subgrid);                                  \
              ny = SubgridNY(subgrid);                                  \
              nz = SubgridNZ(subgrid);                                  \
                                                                        \
              r = SubgridRX(subgrid);                                   \
                                                                        \
              prdat = SubvectorData(pr_sub);                            \
              ppdat = SubvectorData(pp_sub);                            \
              pddat = SubvectorData(pd_sub);                            \
                                                                        \
                                                                        \
              if (fcn == CALCFCN)                                       \
              {                                                         \
                if (dummy1->lookup_tables[ir])                          \
                {                                                       \
                  switch (dummy1->lookup_tables[ir]->interpolation_method) \
                  {                                                     \
                                                                        \
                    case 0:                                             \
                    {                                                   \
                      GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz, \
                                     nx, ny, nz,                        \
                      {                                                 \
                        /* Table Lookup */                              \
                                                                        \
                        ipr = SubvectorEltIndex(pr_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipp = SubvectorEltIndex(pp_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipd = SubvectorEltIndex(pd_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 1.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                                                                        \
                          prdat[ipr] = VanGLookupSpline(head,           \
                                                        dummy1->lookup_tables[ir], \
                                                        CALCFCN);       \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                                                                        \
                                                                        \
                    case 1:                                             \
                    {                                                   \
                      int pt = 0;                                       \
                      PhaseRelPermVanGTable *lookup_table = dummy1->lookup_tables[ir]; \
                      double interval = lookup_table->interval;         \
                      double min_pressure_head = lookup_table->min_pressure_head; \
                      int num_sample_points = lookup_table->num_sample_points; \
                      int max = num_sample_points + 1;                  \
                                                                        \
                      GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz, \
                                     nx, ny, nz,                        \
                      {                                                 \
                        /* Table Lookup */                              \
                                                                        \
                        ipr = SubvectorEltIndex(pr_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipp = SubvectorEltIndex(pp_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipd = SubvectorEltIndex(pd_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 1.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                                                                        \
                          if (head < fabs(min_pressure_head))           \
                          {                                             \
                            pt = (int)floor(head / interval);           \
                            assert(pt < max);                           \
                                                                        \
                            prdat[ipr] = lookup_table->a[pt] + lookup_table->slope[pt] * \
                                         (head - lookup_table->x[pt]);  \
                          }                                             \
                          else                                          \
                          {                                             \
                            prdat[ipr] = 0.0;                           \
                          }                                             \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                  }                                                     \
                }                                                       \
                else                                                    \
                {                                                       \
                  /* Compute VanG curve */                              \
                                                                        \
                  GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz, \
                                 nx, ny, nz,                            \
                  {                                                     \
                    ipr = SubvectorEltIndex(pr_sub,                     \
                                            i + fdir[0], j + fdir[1], k + fdir[2]); \
                    ipp = SubvectorEltIndex(pp_sub,                     \
                                            i + fdir[0], j + fdir[1], k + fdir[2]); \
                    ipd = SubvectorEltIndex(pd_sub,                     \
                                            i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                    if (ppdat[ipp] >= 0.0)                              \
                      prdat[ipr] = 1.0;                                 \
                    else                                                \
                    {                                                   \
                      alpha = alphas[ir];                               \
                      n = ns[ir];                                       \
                      m = 1.0e0 - (1.0e0 / n);                          \
                                                                        \
                      head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                      opahn = 1.0 + pow(alpha * head, n);               \
                      ahnm1 = pow(alpha * head, n - 1);                 \
                      prdat[ipr] = pow(1.0 - ahnm1 / (pow(opahn, m)), 2) \
                                   / pow(opahn, (m / 2));               \
                    }                                                   \
                  });                                                   \
                }                                                       \
              }                                                         \
              else  /* fcn = CALCDER */                                 \
              {                                                         \
                if (dummy1->lookup_tables[ir])                          \
                {                                                       \
                  switch (dummy1->lookup_tables[ir]->interpolation_method) \
                  {                                                     \
                                                                        \
                    case 0:                                             \
                    {                                                   \
                      GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz, \
                                     nx, ny, nz,                        \
                      {                                                 \
                        /* Table Lookup */                              \
                                                                        \
                        ipr = SubvectorEltIndex(pr_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipp = SubvectorEltIndex(pp_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipd = SubvectorEltIndex(pd_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 0.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                          prdat[ipr] = VanGLookupSpline(head,           \
                                                        dummy1->lookup_tables[ir], \
                                                        CALCDER);       \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                                                                        \
                                                                        \
                    case 1:                                             \
                    {                                                   \
                      int pt = 0;                                       \
                      PhaseRelPermVanGTable *lookup_table = dummy1->lookup_tables[ir]; \
                      double interval = lookup_table->interval;         \
                      double min_pressure_head = lookup_table->min_pressure_head; \
                      int num_sample_points = lookup_table->num_sample_points; \
                      int max = num_sample_points + 1;                  \
                                                                        \
                      GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz, \
                                     nx, ny, nz,                        \
                      {                                                 \
                        /* Table Lookup */                              \
                                                                        \
                        ipr = SubvectorEltIndex(pr_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipp = SubvectorEltIndex(pp_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                        ipd = SubvectorEltIndex(pd_sub,                 \
                                                i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 0.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                          if (ppdat[ipp] >= 0.0)                        \
                            prdat[ipr] = 1.0;                           \
                          else                                          \
                          {                                             \
                            head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                                                                        \
                            if (head < fabs(min_pressure_head))         \
                            {                                           \
                              pt = (int)floor(head / interval);         \
                              assert(pt < max);                         \
                                                                        \
                              prdat[ipr] = lookup_table->a_der[pt] + lookup_table->slope_der[pt] * \
                                           (head - lookup_table->x[pt]); \
                            }                                           \
                            else                                        \
                            {                                           \
                              prdat[ipr] = 0.0;                         \
                            }                                           \
                          }                                             \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                  }                                                     \
                }                                                       \
                else                                                    \
                {                                                       \
                  /* Compute VanG curve */                              \
                                                                        \
                  GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz, \
                                 nx, ny, nz,                            \
                  {                                                     \
                    ipr = SubvectorEltIndex(pr_sub,                     \
                                            i + fdir[0], j + fdir[1], k + fdir[2]); \
                    ipp = SubvectorEltIndex(pp_sub,                     \
                                            i + fdir[0], j + fdir[1], k + fdir[2]); \
                    ipd = SubvectorEltIndex(pd_sub,                     \
                                            i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                    if (ppdat[ipp] >= 0.0)                              \
                      prdat[ipr] = 0.0;                                 \
                    else                                                \
                    {                                                   \
                      alpha = alphas[ir];                               \
                      n = ns[ir];                                       \
                      m = 1.0e0 - (1.0e0 / n);                          \
                                                                        \
                      head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                      opahn = 1.0 + pow(alpha * head, n);               \
                      ahnm1 = pow(alpha * head, n - 1);                 \
                      coeff = 1.0 - ahnm1 * pow(opahn, -m);             \
                                                                        \
                      prdat[ipr] = 2.0 * (coeff / (pow(opahn, (m / 2)))) \
                                   * ((n - 1) * pow(alpha * head, n - 2) * alpha \
                                      * pow(opahn, -m)                  \
                                      - ahnm1 * m * pow(opahn, -(m + 1)) * n * alpha * ahnm1) \
                                   + pow(coeff, 2) * (m / 2) * pow(opahn, (-(m + 2) / 2)) \
                                   * n * alpha * ahnm1;                 \
                    }                                                   \
                  });                                                   \
                }                                                       \
              }   /* End else clause */                                 \
            }     /* End subgrid loop */                                \
          }       /* End loop over regions */                           \
        }         /* End if data not from file */                       \
                                                                        \
        else if (data_from_file == 1)  /* ns and alphas from pfb file */ \
        {                                                               \
          gr_solid = ProblemDataGrDomain(problem_data);                 \
          n_values = dummy1->n_values;                                  \
          alpha_values = dummy1->alpha_values;                          \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
                                                                        \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
            pp_sub = VectorSubvector(phase_pressure, sg);               \
            pd_sub = VectorSubvector(phase_density, sg);                \
                                                                        \
            n_values_sub = VectorSubvector(n_values, sg);               \
            alpha_values_sub = VectorSubvector(alpha_values, sg);       \
                                                                        \
            ix = SubgridIX(subgrid);                                    \
            iy = SubgridIY(subgrid);                                    \
            iz = SubgridIZ(subgrid);                                    \
                                                                        \
            nx = SubgridNX(subgrid);                                    \
            ny = SubgridNY(subgrid);                                    \
            nz = SubgridNZ(subgrid);                                    \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
            pddat = SubvectorData(pd_sub);                              \
                                                                        \
            n_values_dat = SubvectorData(n_values_sub);                 \
            alpha_values_dat = SubvectorData(alpha_values_sub);         \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipp = SubvectorEltIndex(pp_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipd = SubvectorEltIndex(pd_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                n_index = SubvectorEltIndex(n_values_sub, i, j, k);     \
                alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k); \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 1.0;                                     \
                else                                                    \
                {                                                       \
                  alpha = alpha_values_dat[alpha_index];                \
                  n = n_values_dat[n_index];                            \
                  m = 1.0e0 - (1.0e0 / n);                              \
                                                                        \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  opahn = 1.0 + pow(alpha * head, n);                   \
                  ahnm1 = pow(alpha * head, n - 1);                     \
                  prdat[ipr] = pow(1.0 - ahnm1 / (pow(opahn, m)), 2)    \
                               / pow(opahn, (m / 2));                   \
                }                                                       \
              });                                                       \
            }                                                           \
            else    /* fcn = CALCDER */                                 \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipp = SubvectorEltIndex(pp_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipd = SubvectorEltIndex(pd_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                n_index = SubvectorEltIndex(n_values_sub, i, j, k);     \
                alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k); \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  alpha = alpha_values_dat[alpha_index];                \
                  n = n_values_dat[n_index];                            \
                  m = 1.0e0 - (1.0e0 / n);                              \
                                                                        \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  opahn = 1.0 + pow(alpha * head, n);                   \
                  ahnm1 = pow(alpha * head, n - 1);                     \
                  coeff = 1.0 - ahnm1 * pow(opahn, -m);                 \
                                                                        \
                  prdat[ipr] = 2.0 * (coeff / (pow(opahn, (m / 2))))    \
                               * ((n - 1) * pow(alpha * head, n - 2) * alpha \
                                  * pow(opahn, -m)                      \
                                  - ahnm1 * m * pow(opahn, -(m + 1)) * n * alpha * ahnm1) \
                               + pow(coeff, 2) * (m / 2) * pow(opahn, (-(m + 2) / 2)) \
                               * n * alpha * ahnm1;                     \
                }                                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End if data_from_file */                           \
                                                                        \
        /* Compute rel. perms. on interior */                           \
        if (data_from_file == 0)  /* alphas and ns given by region */   \
        {                                                               \
          for (ir = 0; ir < num_regions; ir++)                          \
          {                                                             \
            gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
            ForSubgridI(sg, subgrids)                                   \
            {                                                           \
              subgrid = SubgridArraySubgrid(subgrids, sg);              \
                                                                        \
              pr_sub = VectorSubvector(phase_rel_perm, sg);             \
              pp_sub = VectorSubvector(phase_pressure, sg);             \
              pd_sub = VectorSubvector(phase_density, sg);              \
                                                                        \
              ix = SubgridIX(subgrid) - 1;                              \
              iy = SubgridIY(subgrid) - 1;                              \
              iz = SubgridIZ(subgrid) - 1;                              \
                                                                        \
              nx = SubgridNX(subgrid) + 2;                              \
              ny = SubgridNY(subgrid) + 2;                              \
              nz = SubgridNZ(subgrid) + 2;                              \
                                                                        \
              r = SubgridRX(subgrid);                                   \
                                                                        \
              prdat = SubvectorData(pr_sub);                            \
              ppdat = SubvectorData(pp_sub);                            \
              pddat = SubvectorData(pd_sub);                            \
                                                                        \
              if (fcn == CALCFCN)                                       \
              {                                                         \
                if (dummy1->lookup_tables[ir])                          \
                {                                                       \
                  switch (dummy1->lookup_tables[ir]->interpolation_method) \
                  {                                                     \
                                                                        \
                    case 0:                                             \
                    {                                                   \
                      GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
                      {                                                 \
                        /* Table Lookup */                              \
                        ipr = SubvectorEltIndex(pr_sub, i, j, k);       \
                        ipp = SubvectorEltIndex(pp_sub, i, j, k);       \
                        ipd = SubvectorEltIndex(pd_sub, i, j, k);       \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 1.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                          prdat[ipr] = VanGLookupSpline(head,           \
                                                        dummy1->lookup_tables[ir], \
                                                        CALCFCN);       \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                                                                        \
                    case 1:                                             \
                    {                                                   \
                      int pt = 0;                                       \
                      PhaseRelPermVanGTable *lookup_table = dummy1->lookup_tables[ir]; \
                      double interval = lookup_table->interval;         \
                      double min_pressure_head = lookup_table->min_pressure_head; \
                      int num_sample_points = lookup_table->num_sample_points; \
                      int max = num_sample_points + 1;                  \
                                                                        \
                      GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
                      {                                                 \
                        /* Table Lookup */                              \
                        ipr = SubvectorEltIndex(pr_sub, i, j, k);       \
                        ipp = SubvectorEltIndex(pp_sub, i, j, k);       \
                        ipd = SubvectorEltIndex(pd_sub, i, j, k);       \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 1.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                          if (ppdat[ipp] >= 0.0)                        \
                            prdat[ipr] = 1.0;                           \
                          else                                          \
                          {                                             \
                            head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                                                                        \
                            if (head < fabs(min_pressure_head))         \
                            {                                           \
                              pt = (int)floor(head / interval);         \
                              assert(pt < max);                         \
                                                                        \
                              prdat[ipr] = lookup_table->a[pt] + lookup_table->slope[pt] * \
                                           (head - lookup_table->x[pt]); \
                            }                                           \
                            else                                        \
                            {                                           \
                              prdat[ipr] = 0.0;                         \
                            }                                           \
                          }                                             \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                  }                                                     \
                }                                                       \
                else                                                    \
                {                                                       \
                  /* Compute VanG curve */                              \
                                                                        \
                  GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
                  {                                                     \
                    ipr = SubvectorEltIndex(pr_sub, i, j, k);           \
                    ipp = SubvectorEltIndex(pp_sub, i, j, k);           \
                    ipd = SubvectorEltIndex(pd_sub, i, j, k);           \
                                                                        \
                    if (ppdat[ipp] >= 0.0)                              \
                      prdat[ipr] = 1.0;                                 \
                    else                                                \
                    {                                                   \
                      alpha = alphas[ir];                               \
                      n = ns[ir];                                       \
                      m = 1.0e0 - (1.0e0 / n);                          \
                                                                        \
                      head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                      opahn = 1.0 + pow(alpha * head, n);               \
                      ahnm1 = pow(alpha * head, n - 1);                 \
                      prdat[ipr] = pow(1.0 - ahnm1 / (pow(opahn, m)), 2) \
                                   / pow(opahn, (m / 2));               \
                    }                                                   \
                  });                                                   \
                }                                                       \
              }    /* End if clause */                                  \
              else /* fcn = CALCDER */                                  \
              {                                                         \
                if (dummy1->lookup_tables[ir])                          \
                {                                                       \
                  switch (dummy1->lookup_tables[ir]->interpolation_method) \
                  {                                                     \
                                                                        \
                    case 0:                                             \
                    {                                                   \
                      GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
                      {                                                 \
                        /* Table Lookup */                              \
                        ipr = SubvectorEltIndex(pr_sub, i, j, k);       \
                        ipp = SubvectorEltIndex(pp_sub, i, j, k);       \
                        ipd = SubvectorEltIndex(pd_sub, i, j, k);       \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 0.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                                                                        \
                          prdat[ipr] = VanGLookupSpline(head,           \
                                                        dummy1->lookup_tables[ir], \
                                                        CALCDER);       \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                                                                        \
                    case 1:                                             \
                    {                                                   \
                      int pt = 0;                                       \
                      PhaseRelPermVanGTable *lookup_table = dummy1->lookup_tables[ir]; \
                      double interval = lookup_table->interval;         \
                      double min_pressure_head = lookup_table->min_pressure_head; \
                      int num_sample_points = lookup_table->num_sample_points; \
                      int max = num_sample_points + 1;                  \
                                                                        \
                      GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
                      {                                                 \
                        /* Table Lookup */                              \
                        ipr = SubvectorEltIndex(pr_sub, i, j, k);       \
                        ipp = SubvectorEltIndex(pp_sub, i, j, k);       \
                        ipd = SubvectorEltIndex(pd_sub, i, j, k);       \
                                                                        \
                        if (ppdat[ipp] >= 0.0)                          \
                          prdat[ipr] = 0.0;                             \
                        else                                            \
                        {                                               \
                          head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                                                                        \
                          if (head < fabs(min_pressure_head))           \
                          {                                             \
                            pt = (int)floor(head / interval);           \
                            assert(pt < max);                           \
                                                                        \
                            prdat[ipr] = lookup_table->a_der[pt] + lookup_table->slope_der[pt] * \
                                         (head - lookup_table->x[pt]);  \
                          }                                             \
                          else                                          \
                          {                                             \
                            prdat[ipr] = 0.0;                           \
                          }                                             \
                        }                                               \
                      });                                               \
                    }                                                   \
                    break;                                              \
                  }                                                     \
                }                                                       \
                else                                                    \
                {                                                       \
                  /* Compute VanG curve */                              \
                                                                        \
                  GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
                  {                                                     \
                    ipr = SubvectorEltIndex(pr_sub, i, j, k);           \
                    ipp = SubvectorEltIndex(pp_sub, i, j, k);           \
                    ipd = SubvectorEltIndex(pd_sub, i, j, k);           \
                                                                        \
                    if (ppdat[ipp] >= 0.0)                              \
                      prdat[ipr] = 0.0;                                 \
                    else                                                \
                    {                                                   \
                      alpha = alphas[ir];                               \
                      n = ns[ir];                                       \
                      m = 1.0e0 - (1.0e0 / n);                          \
                                                                        \
                      head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity); \
                      opahn = 1.0 + pow(alpha * head, n);               \
                      ahnm1 = pow(alpha * head, n - 1);                 \
                      coeff = 1.0 - ahnm1 * pow(opahn, -m);             \
                                                                        \
                      prdat[ipr] = 2.0 * (coeff / (pow(opahn, (m / 2)))) \
                                   * ((n - 1) * pow(alpha * head, n - 2) * alpha \
                                      * pow(opahn, -m)                  \
                                      - ahnm1 * m * pow(opahn, -(m + 1)) * n * alpha * ahnm1) \
                                   + pow(coeff, 2) * (m / 2) * pow(opahn, (-(m + 2) / 2)) \
                                   * n * alpha * ahnm1;                 \
                    }                                                   \
                  });                                                   \
                }                                                       \
              }   /* End else clause */                                 \
            }     /* End subgrid loop */                                \
          }       /* End subregion loop */                              \
        }         /* End if data not given by file */                   \
        else if (data_from_file == 1) /* alphas and ns given in pfb files */ \
        {                                                               \
          gr_solid = ProblemDataGrDomain(problem_data);                 \
          n_values = dummy1->n_values;                                  \
          alpha_values = dummy1->alpha_values;                          \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
                                                                        \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
            pp_sub = VectorSubvector(phase_pressure, sg);               \
            pd_sub = VectorSubvector(phase_density, sg);                \
                                                                        \
            n_values_sub = VectorSubvector(n_values, sg);               \
            alpha_values_sub = VectorSubvector(alpha_values, sg);       \
                                                                        \
            ix = SubgridIX(subgrid) - 1;                                \
            iy = SubgridIY(subgrid) - 1;                                \
            iz = SubgridIZ(subgrid) - 1;                                \
                                                                        \
            nx = SubgridNX(subgrid) + 2;                                \
            ny = SubgridNY(subgrid) + 2;                                \
            nz = SubgridNZ(subgrid) + 2;                                \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
            pddat = SubvectorData(pd_sub);                              \
                                                                        \
            n_values_dat = SubvectorData(n_values_sub);                 \
            alpha_values_dat = SubvectorData(alpha_values_sub);         \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                ipd = SubvectorEltIndex(pd_sub, i, j, k);               \
                                                                        \
                n_index = SubvectorEltIndex(n_values_sub, i, j, k);     \
                alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k); \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 1.0;                                     \
                else                                                    \
                {                                                       \
                  alpha = alpha_values_dat[alpha_index];                \
                  n = n_values_dat[n_index];                            \
                  m = 1.0e0 - (1.0e0 / n);                              \
                                                                        \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  opahn = 1.0 + pow(alpha * head, n);                   \
                  ahnm1 = pow(alpha * head, n - 1);                     \
                  prdat[ipr] = pow(1.0 - ahnm1 / (pow(opahn, m)), 2)    \
                               / pow(opahn, (m / 2));                   \
                }                                                       \
              });                                                       \
            }      /* End if clause */                                  \
            else   /* fcn = CALCDER */                                  \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                ipd = SubvectorEltIndex(pd_sub, i, j, k);               \
                                                                        \
                n_index = SubvectorEltIndex(n_values_sub, i, j, k);     \
                alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k); \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  alpha = alpha_values_dat[alpha_index];                \
                  n = n_values_dat[n_index];                            \
                  m = 1.0e0 - (1.0e0 / n);                              \
                                                                        \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  opahn = 1.0 + pow(alpha * head, n);                   \
                  ahnm1 = pow(alpha * head, n - 1);                     \
                  coeff = 1.0 - ahnm1 * pow(opahn, -m);                 \
                                                                        \
                  prdat[ipr] = 2.0 * (coeff / (pow(opahn, (m / 2))))    \
                               * ((n - 1) * pow(alpha * head, n - 2) * alpha \
                                  * pow(opahn, -m)                      \
                                  - ahnm1 * m * pow(opahn, -(m + 1)) * n * alpha * ahnm1) \
                               + pow(coeff, 2) * (m / 2) * pow(opahn, (-(m + 2) / 2)) \
                               * n * alpha * ahnm1;                     \
                }                                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End if data given by file */                       \
        break;                                                          \
      }        /* End case 1 */                                         \
                                                                        \
      case 2:  /* Haverkamp et.al. relative permeability */             \
      {                                                                 \
        double  *As, *gammas, head, tmp;                                \
                                                                        \
        GetDummyType(PhaseRelPerm, 2, (public_xtra->data), dummy2);     \
                                                                        \
        num_regions = (dummy2->num_regions);                            \
        region_indices = (dummy2->region_indices);                      \
        As = (dummy2->As);                                              \
        gammas = (dummy2->gammas);                                      \
                                                                        \
        /* Compute rel. perms. for Dirichlet BC's */                    \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
                                                                        \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
            pp_sub = VectorSubvector(phase_pressure, sg);               \
            pd_sub = VectorSubvector(phase_density, sg);                \
                                                                        \
            ix = SubgridIX(subgrid);                                    \
            iy = SubgridIY(subgrid);                                    \
            iz = SubgridIZ(subgrid);                                    \
                                                                        \
            nx = SubgridNX(subgrid);                                    \
            ny = SubgridNY(subgrid);                                    \
            nz = SubgridNZ(subgrid);                                    \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
            pddat = SubvectorData(pd_sub);                              \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipp = SubvectorEltIndex(pp_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipd = SubvectorEltIndex(pd_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 1.0;                                     \
                else                                                    \
                {                                                       \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  tmp = As[ir] + pow(head, gammas[ir]);                 \
                  prdat[ipr] = As[ir] / tmp;                            \
                }                                                       \
              });                                                       \
            }      /* End if clause */                                  \
            else   /* fcn = CALCDER */                                  \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipp = SubvectorEltIndex(pp_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipd = SubvectorEltIndex(pd_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  tmp = pow(head, gammas[ir]);                          \
                  prdat[ipr] = As[ir] * gammas[ir]                      \
                               * pow(head, gammas[ir] - 1) / pow(tmp, 2); \
                }                                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End subregion loop */                              \
                                                                        \
        /* Compute rel. perms. on interior */                           \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
                                                                        \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
            pp_sub = VectorSubvector(phase_pressure, sg);               \
            pd_sub = VectorSubvector(phase_density, sg);                \
                                                                        \
            ix = SubgridIX(subgrid) - 1;                                \
            iy = SubgridIY(subgrid) - 1;                                \
            iz = SubgridIZ(subgrid) - 1;                                \
                                                                        \
            nx = SubgridNX(subgrid) + 2;                                \
            ny = SubgridNY(subgrid) + 2;                                \
            nz = SubgridNZ(subgrid) + 2;                                \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
            pddat = SubvectorData(pd_sub);                              \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                ipd = SubvectorEltIndex(pd_sub, i, j, k);               \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 1.0;                                     \
                else                                                    \
                {                                                       \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  tmp = As[ir] + pow(head, gammas[ir]);                 \
                  prdat[ipr] = As[ir] / tmp;                            \
                }                                                       \
              });                                                       \
            }      /* End if clause */                                  \
            else   /* fcn = CALCDER */                                  \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                ipd = SubvectorEltIndex(pd_sub, i, j, k);               \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  prdat[ipr] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  tmp = pow(head, gammas[ir]);                          \
                  prdat[ipr] = As[ir] * gammas[ir]                      \
                               * pow(head, gammas[ir] - 1) / pow(tmp, 2); \
                }                                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End subregion loop */                              \
                                                                        \
        break;                                                          \
      }        /* End case 2 */                                         \
                                                                        \
      case 3:  /* Data relative permeability */                         \
      {                                                                 \
        GetDummyType(PhaseRelPerm, 3, (public_xtra->data), dummy3);     \
                                                                        \
        if (!amps_Rank(amps_CommWorld))                                 \
          printf("Data curves for rel perms not supported currently.\n"); \
        break;                                                          \
      }        /* End case 3 */                                         \
                                                                        \
      case 4:  /* Polynomial function of pressure relative permeability */ \
      {                                                                 \
        int     *degrees, dg;                                           \
        double **coefficients, *region_coeffs;                          \
                                                                        \
        GetDummyType(PhaseRelPerm, 4, (public_xtra->data), dummy4);     \
                                                                        \
        num_regions = (dummy4->num_regions);                            \
        region_indices = (dummy4->region_indices);                      \
        degrees = (dummy4->degrees);                                    \
        coefficients = (dummy4->coefficients);                          \
                                                                        \
        /* Compute rel. perms. for Dirichlet BC's */                    \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
          region_coeffs = coefficients[ir];                             \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
                                                                        \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
            pp_sub = VectorSubvector(phase_pressure, sg);               \
                                                                        \
            ix = SubgridIX(subgrid);                                    \
            iy = SubgridIY(subgrid);                                    \
            iz = SubgridIZ(subgrid);                                    \
                                                                        \
            nx = SubgridNX(subgrid);                                    \
            ny = SubgridNY(subgrid);                                    \
            nz = SubgridNZ(subgrid);                                    \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipp = SubvectorEltIndex(pp_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                if (ppdat[ipp] == 0.0)                                  \
                  prdat[ipr] = region_coeffs[0];                        \
                else                                                    \
                {                                                       \
                  prdat[ipr] = 0.0;                                     \
                  for (dg = 0; dg < degrees[ir] + 1; dg++)              \
                  {                                                     \
                    prdat[ipr] += region_coeffs[dg] * pow(ppdat[ipp], dg); \
                  }                                                     \
                }                                                       \
              });                                                       \
            }      /* End if clause */                                  \
            else   /* fcn = CALCDER */                                  \
            {                                                           \
              GrGeomSurfLoop(i, j, k, fdir, gr_solid, r, ix, iy, iz,    \
                             nx, ny, nz,                                \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                ipp = SubvectorEltIndex(pp_sub,                         \
                                        i + fdir[0], j + fdir[1], k + fdir[2]); \
                if (ppdat[ipp] == 0.0)                                  \
                  prdat[ipr] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  prdat[ipr] = 0.0;                                     \
                  for (dg = 1; dg < degrees[ir] + 1; dg++)              \
                  {                                                     \
                    prdat[ipr] += region_coeffs[dg] * dg                \
                                  * pow(ppdat[ipp], (dg - 1));          \
                  }                                                     \
                }                                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End subregion loop */                              \
                                                                        \
        /* Compute rel. perms. in interior */                           \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
          region_coeffs = coefficients[ir];                             \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
                                                                        \
            pr_sub = VectorSubvector(phase_rel_perm, sg);               \
            pp_sub = VectorSubvector(phase_pressure, sg);               \
                                                                        \
            ix = SubgridIX(subgrid) - 1;                                \
            iy = SubgridIY(subgrid) - 1;                                \
            iz = SubgridIZ(subgrid) - 1;                                \
                                                                        \
            nx = SubgridNX(subgrid) + 2;                                \
            ny = SubgridNY(subgrid) + 2;                                \
            nz = SubgridNZ(subgrid) + 2;                                \
                                                                        \
            r = SubgridRX(subgrid);                                     \
                                                                        \
            prdat = SubvectorData(pr_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                                                                        \
                if (ppdat[ipp] == 0.0)                                  \
                  prdat[ipr] = region_coeffs[0];                        \
                else                                                    \
                {                                                       \
                  prdat[ipr] = 0.0;                                     \
                  for (dg = 0; dg < degrees[ir] + 1; dg++)              \
                  {                                                     \
                    prdat[ipr] += region_coeffs[dg] * pow(ppdat[ipp], dg); \
                  }                                                     \
                }                                                       \
              });                                                       \
            }      /* End if clause */                                  \
            else   /* fcn = CALCDER */                                  \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ipr = SubvectorEltIndex(pr_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                                                                        \
                if (ppdat[ipp] == 0.0)                                  \
                  prdat[ipr] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  prdat[ipr] = 0.0;                                     \
                  for (dg = 1; dg < degrees[ir] + 1; dg++)              \
                  {                                                     \
                    prdat[ipr] += region_coeffs[dg] * dg                \
                                  * pow(ppdat[ipp], (dg - 1));          \
                  }                                                     \
                }                                                       \
              });                                                       \
            }     /* End else clause */                                 \
          }       /* End subgrid loop */                                \
        }         /* End subregion loop */                              \
                                                                        \
        break;                                                          \
      }        /* End case 4 */                                         \
    }          /* End switch */                                         \
                                                                        \
    IncFLOPCount(1);                                                    \
    EndTiming(public_xtra->time_index);                                 \
  }


#endif // _INLINE_PHASE_REL_PERM_H
