/* This definition _really_ screws with parsers, keep it in its own file */

#ifndef _INLINE_SATURATION_H
#define _INLINE_SATURATION_H

#define SATURATION_MODULE(module, saturation, pressure, density, ftype) \
  {                                                                     \
  GetModulePublicXtra(Saturation, module, public_xtra);                 \
  Vector *phase_saturation = saturation;                                \
  Vector *phase_pressure = pressure;                                    \
  Vector *phase_density = density;                                      \
  int fcn = ftype;                                                      \
                                                                        \
  GrGeomSolid   *gr_solid, *gr_domain;                                  \
  Subvector     *ps_sub;                                                \
  Subvector     *pp_sub;                                                \
  Subvector     *pd_sub;                                                \
  Subvector     *satRF_sub;                                             \
  Subvector     *n_values_sub;                                          \
  Subvector     *alpha_values_sub;                                      \
  Subvector     *s_res_values_sub;                                      \
  Subvector     *s_sat_values_sub;                                      \
                                                                        \
  double        *psdat, *ppdat, *pddat, *satRFdat;                      \
  double        *n_values_dat, *alpha_values_dat;                       \
  double        *s_res_values_dat, *s_sat_values_dat;                   \
                                                                        \
  int ix, iy, iz, r;                                                    \
  int nx, ny, nz;                                                       \
                                                                        \
  int sg;                                                               \
  int i, j, k, ips, ipp, ipd, ipRF;                                     \
                                                                        \
  int n_index, alpha_index, s_res_index, s_sat_index;                   \
                                                                        \
  int            *region_indices, num_regions, ir;                      \
                                                                        \
  Grid *grid = VectorGrid(phase_saturation);                            \
  SubgridArray *subgrids = GridSubgrids(grid);                          \
                                                                        \
  InitVectorAll(phase_saturation, -FLT_MAX);                            \
  switch (public_xtra->type) {                                             \
    case 0:                                                             \
    {                                                                   \
      SaturationType0 *dummy0 = (SaturationType0*)(public_xtra->data);     \
      num_regions = (dummy0->num_regions);                              \
      region_indices = (dummy0->region_indices);                        \
      double *values = (dummy0->values);                                \
                                                                        \
      for (ir = 0; ir < num_regions; ir++) {                            \
        GrGeomSolid *gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
        ForSubgridI(sg, subgrids)                                       \
        {                                                               \
          subgrid = SubgridArraySubgrid(subgrids, sg);                  \
          Subvector *ps_sub = VectorSubvector(phase_saturation, sg);    \
          ix = SubgridIX(subgrid);                                      \
          iy = SubgridIY(subgrid);                                      \
          iz = SubgridIZ(subgrid);                                      \
                                                                        \
          nx = SubgridNX(subgrid);                                      \
          ny = SubgridNY(subgrid);                                      \
          nz = SubgridNZ(subgrid);                                      \
                                                                        \
          r = SubgridRX(subgrid);                                       \
                                                                        \
          double *psdat = SubvectorData(ps_sub);                        \
          int ips;                                                      \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              psdat[ips] = values[ir];                                  \
            });                                                         \
          }                                                             \
          else   /* fcn = CALCDER */                                    \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              psdat[ips] = 0.0;                                         \
            });                                                         \
          }     /* End else clause */                                   \
        }                                                               \
      }                                                                 \
      break;                                                            \
    }                                                                   \
                                                                        \
    case 1: /* Van Genuchten saturation curve */                        \
    {                                                                   \
      int data_from_file;                                               \
      double *alphas, *ns, *s_ress, *s_difs;                            \
      double head, alpha, n, s_res, s_dif, s_sat, m;                    \
                                                                        \
      Vector *n_values, *alpha_values, *s_res_values, *s_sat_values;    \
                                                                        \
      GetDummyType(Saturation, 1, (public_xtra->data), dummy1);            \
                                                                        \
      num_regions = (dummy1->num_regions);                              \
      region_indices = (dummy1->region_indices);                        \
      alphas = (dummy1->alphas);                                        \
      ns = (dummy1->ns);                                                \
      s_ress = (dummy1->s_ress);                                        \
      s_difs = (dummy1->s_difs);                                        \
      data_from_file = (dummy1->data_from_file);                        \
                                                                        \
      if (data_from_file == 0) /* Soil parameters given by region */    \
      {                                                                 \
        for (ir = 0; ir < num_regions; ir++)                            \
        {                                                               \
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
          ForSubgridI(sg, subgrids)                                     \
          {                                                             \
            subgrid = SubgridArraySubgrid(subgrids, sg);                \
            ps_sub = VectorSubvector(phase_saturation, sg);             \
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
            psdat = SubvectorData(ps_sub);                              \
            ppdat = SubvectorData(pp_sub);                              \
            pddat = SubvectorData(pd_sub);                              \
                                                                        \
            if (fcn == CALCFCN)                                         \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                ipd = SubvectorEltIndex(pd_sub, i, j, k);               \
                                                                        \
                alpha = alphas[ir];                                     \
                n = ns[ir];                                             \
                m = 1.0e0 - (1.0e0 / n);                                \
                s_res = s_ress[ir];                                     \
                s_dif = s_difs[ir];                                     \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  psdat[ips] = s_dif + s_res;                           \
                else                                                    \
                {                                                       \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  psdat[ips] = s_dif / pow(1.0 + pow((alpha * head), n), m) \
                               + s_res;                                 \
                }                                                       \
              });                                                       \
            }    /* End if clause */                                    \
            else /* fcn = CALCDER */                                    \
            {                                                           \
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz, \
              {                                                         \
                ips = SubvectorEltIndex(ps_sub, i, j, k);               \
                ipp = SubvectorEltIndex(pp_sub, i, j, k);               \
                ipd = SubvectorEltIndex(pd_sub, i, j, k);               \
                                                                        \
                alpha = alphas[ir];                                     \
                n = ns[ir];                                             \
                m = 1.0e0 - (1.0e0 / n);                                \
                s_res = s_ress[ir];                                     \
                s_dif = s_difs[ir];                                     \
                                                                        \
                if (ppdat[ipp] >= 0.0)                                  \
                  psdat[ips] = 0.0;                                     \
                else                                                    \
                {                                                       \
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);     \
                  psdat[ips] = (m * n * alpha * pow(alpha * head, (n - 1))) * s_dif \
                               / (pow(1.0 + pow(alpha * head, n), m + 1)); \
                }                                                       \
              });                                                       \
            }   /* End else clause */                                   \
          }     /* End subgrid loop */                                  \
        }       /* End loop over regions */                             \
      }         /* End if data not from file */                         \
      else                                                              \
      {                                                                 \
        gr_solid = ProblemDataGrDomain(problem_data);                   \
        n_values = dummy1->n_values;                                    \
        alpha_values = dummy1->alpha_values;                            \
        s_res_values = dummy1->s_res_values;                            \
        s_sat_values = dummy1->s_sat_values;                            \
                                                                        \
        ForSubgridI(sg, subgrids)                                       \
        {                                                               \
          subgrid = SubgridArraySubgrid(subgrids, sg);                  \
          ps_sub = VectorSubvector(phase_saturation, sg);               \
          pp_sub = VectorSubvector(phase_pressure, sg);                 \
          pd_sub = VectorSubvector(phase_density, sg);                  \
                                                                        \
          n_values_sub = VectorSubvector(n_values, sg);                 \
          alpha_values_sub = VectorSubvector(alpha_values, sg);         \
          s_res_values_sub = VectorSubvector(s_res_values, sg);         \
          s_sat_values_sub = VectorSubvector(s_sat_values, sg);         \
                                                                        \
          ix = SubgridIX(subgrid);                                      \
          iy = SubgridIY(subgrid);                                      \
          iz = SubgridIZ(subgrid);                                      \
                                                                        \
          nx = SubgridNX(subgrid);                                      \
          ny = SubgridNY(subgrid);                                      \
          nz = SubgridNZ(subgrid);                                      \
                                                                        \
          r = SubgridRX(subgrid);                                       \
                                                                        \
          psdat = SubvectorData(ps_sub);                                \
          ppdat = SubvectorData(pp_sub);                                \
          pddat = SubvectorData(pd_sub);                                \
                                                                        \
          n_values_dat = SubvectorData(n_values_sub);                   \
          alpha_values_dat = SubvectorData(alpha_values_sub);           \
          s_res_values_dat = SubvectorData(s_res_values_sub);           \
          s_sat_values_dat = SubvectorData(s_sat_values_sub);           \
                                                                        \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              ipp = SubvectorEltIndex(pp_sub, i, j, k);                 \
              ipd = SubvectorEltIndex(pd_sub, i, j, k);                 \
                                                                        \
              n_index = SubvectorEltIndex(n_values_sub, i, j, k);       \
              alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k); \
              s_res_index = SubvectorEltIndex(s_res_values_sub, i, j, k); \
              s_sat_index = SubvectorEltIndex(s_sat_values_sub, i, j, k); \
                                                                        \
              alpha = alpha_values_dat[alpha_index];                    \
              n = n_values_dat[n_index];                                \
              m = 1.0e0 - (1.0e0 / n);                                  \
              s_res = s_res_values_dat[s_res_index];                    \
              s_sat = s_sat_values_dat[s_sat_index];                    \
                                                                        \
              if (ppdat[ipp] >= 0.0)                                    \
                psdat[ips] = s_sat;                                     \
              else                                                      \
              {                                                         \
                head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);       \
                psdat[ips] = (s_sat - s_res) /                          \
                             pow(1.0 + pow((alpha * head), n), m)       \
                             + s_res;                                   \
              }                                                         \
            });                                                         \
          }      /* End if clause */                                    \
          else   /* fcn = CALCDER */                                    \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              ipp = SubvectorEltIndex(pp_sub, i, j, k);                 \
              ipd = SubvectorEltIndex(pd_sub, i, j, k);                 \
                                                                        \
              n_index = SubvectorEltIndex(n_values_sub, i, j, k);       \
              alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k); \
              s_res_index = SubvectorEltIndex(s_res_values_sub, i, j, k); \
              s_sat_index = SubvectorEltIndex(s_sat_values_sub, i, j, k); \
                                                                        \
              alpha = alpha_values_dat[alpha_index];                    \
              n = n_values_dat[n_index];                                \
              m = 1.0e0 - (1.0e0 / n);                                  \
              s_res = s_res_values_dat[s_res_index];                    \
              s_sat = s_sat_values_dat[s_sat_index];                    \
              s_dif = s_sat - s_res;                                    \
                                                                        \
              if (ppdat[ipp] >= 0.0)                                    \
                psdat[ips] = 0.0;                                       \
              else                                                      \
              {                                                         \
                head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);       \
                psdat[ips] = (m * n * alpha * pow(alpha * head, (n - 1))) * s_dif \
                             / (pow(1.0 + pow(alpha * head, n), m + 1)); \
              }                                                         \
            });                                                         \
          }     /* End else clause */                                   \
        }       /* End subgrid loop */                                  \
      }         /* End if data_from_file */                             \
      break;                                                            \
    }        /* End case 1 */                                           \
                                                                        \
    case 2: /* Haverkamp et.al. saturation curve */                     \
    {                                                                   \
      double *alphas, *betas, *s_ress, *s_difs;                         \
      double head, alpha, beta, s_res, s_dif;                           \
                                                                        \
      GetDummyType(Saturation, 2, (public_xtra->data), dummy2);            \
                                                                        \
      num_regions = (dummy2->num_regions);                              \
      region_indices = (dummy2->region_indices);                        \
      alphas = (dummy2->alphas);                                        \
      betas = (dummy2->betas);                                          \
      s_ress = (dummy2->s_ress);                                        \
      s_difs = (dummy2->s_difs);                                        \
                                                                        \
      for (ir = 0; ir < num_regions; ir++)                              \
      {                                                                 \
        gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
                                                                        \
        ForSubgridI(sg, subgrids)                                       \
        {                                                               \
          subgrid = SubgridArraySubgrid(subgrids, sg);                  \
          ps_sub = VectorSubvector(phase_saturation, sg);               \
          pp_sub = VectorSubvector(phase_pressure, sg);                 \
          pd_sub = VectorSubvector(phase_density, sg);                  \
                                                                        \
          ix = SubgridIX(subgrid);                                      \
          iy = SubgridIY(subgrid);                                      \
          iz = SubgridIZ(subgrid);                                      \
                                                                        \
          nx = SubgridNX(subgrid);                                      \
          ny = SubgridNY(subgrid);                                      \
          nz = SubgridNZ(subgrid);                                      \
                                                                        \
          r = SubgridRX(subgrid);                                       \
                                                                        \
          psdat = SubvectorData(ps_sub);                                \
          ppdat = SubvectorData(pp_sub);                                \
          pddat = SubvectorData(pd_sub);                                \
                                                                        \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              ipp = SubvectorEltIndex(pp_sub, i, j, k);                 \
              ipd = SubvectorEltIndex(pd_sub, i, j, k);                 \
                                                                        \
              alpha = alphas[ir];                                       \
              beta = betas[ir];                                         \
              s_res = s_ress[ir];                                       \
              s_dif = s_difs[ir];                                       \
                                                                        \
              if (ppdat[ipp] >= 0.0)                                    \
                psdat[ips] = s_dif + s_res;                             \
              else                                                      \
              {                                                         \
                head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);       \
                psdat[ips] = alpha * s_dif / (alpha + pow(head, beta))  \
                             + s_res;                                   \
              }                                                         \
            });                                                         \
          }      /* End if clause */                                    \
          else   /* fcn = CALCDER */                                    \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              ipp = SubvectorEltIndex(pp_sub, i, j, k);                 \
              ipd = SubvectorEltIndex(pd_sub, i, j, k);                 \
                                                                        \
              alpha = alphas[ir];                                       \
              beta = betas[ir];                                         \
              s_res = s_ress[ir];                                       \
              s_dif = s_difs[ir];                                       \
                                                                        \
              if (ppdat[ipp] >= 0.0)                                    \
                psdat[ips] = 0.0;                                       \
              else                                                      \
              {                                                         \
                head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);       \
                psdat[ips] = alpha * s_dif * beta * pow(head, beta - 1) \
                             / pow((alpha + pow(head, beta)), 2);       \
              }                                                         \
            });                                                         \
          }     /* End else clause */                                   \
        }       /* End subgrid loop */                                  \
      }         /* End loop over regions */                             \
      break;                                                            \
    }        /* End case 2 */                                           \
                                                                        \
    case 3: /* Data points for saturation curve */                      \
    {                                                                   \
      GetDummyType(Saturation, 3, (public_xtra->data), dummy3);            \
                                                                        \
      num_regions = (dummy3->num_regions);                              \
      region_indices = (dummy3->region_indices);                        \
                                                                        \
      if (!amps_Rank(amps_CommWorld))                                   \
        printf("Data curves for sats not yet supported.\n");            \
                                                                        \
      break;                                                            \
    }        /* End case 3 */                                           \
                                                                        \
    case 4: /* Polynomial function of pressure saturation curve */      \
    {                                                                   \
      int     *degrees, dg;                                             \
      double **coefficients, *region_coeffs;                            \
                                                                        \
      GetDummyType(Saturation, 4, (public_xtra->data), dummy4);            \
                                                                        \
      num_regions = (dummy4->num_regions);                              \
      region_indices = (dummy4->region_indices);                        \
      degrees = (dummy4->degrees);                                      \
      coefficients = (dummy4->coefficients);                            \
                                                                        \
      for (ir = 0; ir < num_regions; ir++)                              \
      {                                                                 \
        gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]); \
        region_coeffs = coefficients[ir];                               \
                                                                        \
        ForSubgridI(sg, subgrids)                                       \
        {                                                               \
          subgrid = SubgridArraySubgrid(subgrids, sg);                  \
                                                                        \
          ps_sub = VectorSubvector(phase_saturation, sg);               \
          pp_sub = VectorSubvector(phase_pressure, sg);                 \
                                                                        \
                                                                        \
          ix = SubgridIX(subgrid);                                      \
          iy = SubgridIY(subgrid);                                      \
          iz = SubgridIZ(subgrid);                                      \
                                                                        \
          nx = SubgridNX(subgrid);                                      \
          ny = SubgridNY(subgrid);                                      \
          nz = SubgridNZ(subgrid);                                      \
                                                                        \
          r = SubgridRX(subgrid);                                       \
                                                                        \
          psdat = SubvectorData(ps_sub);                                \
          ppdat = SubvectorData(pp_sub);                                \
                                                                        \
          if (fcn == CALCFCN)                                           \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              ipp = SubvectorEltIndex(pp_sub, i, j, k);                 \
                                                                        \
              if (ppdat[ipp] == 0.0)                                    \
                psdat[ips] = region_coeffs[0];                          \
              else                                                      \
              {                                                         \
                psdat[ips] = 0.0;                                       \
                for (dg = 0; dg < degrees[ir] + 1; dg++)                \
                {                                                       \
                  psdat[ips] += region_coeffs[dg] * pow(ppdat[ipp], dg); \
                }                                                       \
              }                                                         \
            });                                                         \
          }      /* End if clause */                                    \
          else   /* fcn = CALCDER */                                    \
          {                                                             \
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,  \
            {                                                           \
              ips = SubvectorEltIndex(ps_sub, i, j, k);                 \
              ipp = SubvectorEltIndex(pp_sub, i, j, k);                 \
                                                                        \
              if (ppdat[ipp] == 0.0)                                    \
                psdat[ips] = 0.0;                                       \
              else                                                      \
              {                                                         \
                psdat[ips] = 0.0;                                       \
                for (dg = 0; dg < degrees[ir] + 1; dg++)                \
                {                                                       \
                  psdat[ips] += region_coeffs[dg] * dg                  \
                                * pow(ppdat[ipp], (dg - 1));            \
                }                                                       \
              }                                                         \
            });                                                         \
          }     /* End else clause */                                   \
        }       /* End subgrid loop */                                  \
      }         /* End loop over regions */                             \
      break;                                                            \
    }        /* End case 4 */                                           \
                                                                        \
    case 5: /* ParFlow binary file with spatially varying saturation values */ \
    {                                                                   \
      Vector *satRF;                                                    \
                                                                        \
      GetDummyType(Saturation, 5, (public_xtra->data), dummy5);            \
                                                                        \
      satRF = dummy5->satRF;                                            \
                                                                        \
      gr_domain = ProblemDataGrDomain(problem_data);                    \
                                                                        \
      ForSubgridI(sg, subgrids)                                         \
      {                                                                 \
        subgrid = SubgridArraySubgrid(subgrids, sg);                    \
        ps_sub = VectorSubvector(phase_saturation, sg);                 \
        satRF_sub = VectorSubvector(satRF, sg);                         \
                                                                        \
        ix = SubgridIX(subgrid);                                        \
        iy = SubgridIY(subgrid);                                        \
        iz = SubgridIZ(subgrid);                                        \
                                                                        \
        nx = SubgridNX(subgrid);                                        \
        ny = SubgridNY(subgrid);                                        \
        nz = SubgridNZ(subgrid);                                        \
                                                                        \
        /* RDF: assume resolution is the same in all 3 directions */    \
        r = SubgridRX(subgrid);                                         \
                                                                        \
        psdat = SubvectorData(ps_sub);                                  \
        satRFdat = SubvectorData(satRF_sub);                            \
                                                                        \
        if (fcn == CALCFCN)                                             \
        {                                                               \
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,   \
          {                                                             \
            ips = SubvectorEltIndex(ps_sub, i, j, k);                   \
            ipRF = SubvectorEltIndex(satRF_sub, i, j, k);               \
                                                                        \
            psdat[ips] = satRFdat[ipRF];                                \
          });                                                           \
        }     /* End if clause */                                       \
        else  /* fcn = CALCDER */                                       \
        {                                                               \
          GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,   \
          {                                                             \
            ips = SubvectorEltIndex(ps_sub, i, j, k);                   \
            ipRF = SubvectorEltIndex(satRF_sub, i, j, k);               \
                                                                        \
            psdat[ips] = 0.0;                                           \
          });                                                           \
        }    /* End else clause */                                      \
      }      /* End subgrid loop */                                     \
      break;                                                            \
    }        /* End case 5 */                                           \
  }                                                                     \
  }

#endif // _INLINE_SATURATION_H
