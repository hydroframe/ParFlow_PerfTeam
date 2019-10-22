#ifndef _INLINE_DEFS_H
#define _INLINE_DEFS_H

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


#define SATURATION_MODULE \
      switch (sat_xtra->type) {\
      case 0:\
      {\
        SaturationType0 *dummy0 = (SaturationType0*)(sat_xtra->data);\
        num_regions = (dummy0->num_regions);\
        region_indices = (dummy0->region_indices);\
        double *values = (dummy0->values);\
\
        for (ir = 0; ir < num_regions; ir++) {\
          GrGeomSolid *gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]);\
\
          ForSubgridI(sg, subgrids)\
          {\
            subgrid = SubgridArraySubgrid(subgrids, sg);\
            Subvector *ps_sub = VectorSubvector(phase_saturation, sg);\
            ix = SubgridIX(subgrid);\
            iy = SubgridIY(subgrid);\
            iz = SubgridIZ(subgrid);\
\
            nx = SubgridNX(subgrid);\
            ny = SubgridNY(subgrid);\
            nz = SubgridNZ(subgrid);\
\
            r = SubgridRX(subgrid);\
\
            double *psdat = SubvectorData(ps_sub);\
            int ips;\
            if (fcn == CALCFCN)\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                psdat[ips] = values[ir];\
              });\
            }\
            else   /* fcn = CALCDER */\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                psdat[ips] = 0.0;\
              });\
            }     /* End else clause */\
          }\
        }\
        break;\
      }\
\
          case 1: /* Van Genuchten saturation curve */\
    {\
      int data_from_file;\
      double *alphas, *ns, *s_ress, *s_difs;\
      double head, alpha, n, s_res, s_dif, s_sat, m;\
\
      Vector *n_values, *alpha_values, *s_res_values, *s_sat_values;\
\
      GetDummyType(Saturation, 1, (sat_xtra->data), dummy1);\
\
      num_regions = (dummy1->num_regions);\
      region_indices = (dummy1->region_indices);\
      alphas = (dummy1->alphas);\
      ns = (dummy1->ns);\
      s_ress = (dummy1->s_ress);\
      s_difs = (dummy1->s_difs);\
      data_from_file = (dummy1->data_from_file);\
\
      if (data_from_file == 0) /* Soil parameters given by region */\
      {\
        for (ir = 0; ir < num_regions; ir++)\
        {\
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]);\
\
          ForSubgridI(sg, subgrids)\
          {\
            subgrid = SubgridArraySubgrid(subgrids, sg);\
            ps_sub = VectorSubvector(phase_saturation, sg);\
            pp_sub = VectorSubvector(phase_pressure, sg);\
            pd_sub = VectorSubvector(phase_density, sg);\
\
            ix = SubgridIX(subgrid);\
            iy = SubgridIY(subgrid);\
            iz = SubgridIZ(subgrid);\
\
            nx = SubgridNX(subgrid);\
            ny = SubgridNY(subgrid);\
            nz = SubgridNZ(subgrid);\
\
            r = SubgridRX(subgrid);\
\
            psdat = SubvectorData(ps_sub);\
            ppdat = SubvectorData(pp_sub);\
            pddat = SubvectorData(pd_sub);\
\
            if (fcn == CALCFCN)\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                ipp = SubvectorEltIndex(pp_sub, i, j, k);\
                ipd = SubvectorEltIndex(pd_sub, i, j, k);\
\
                alpha = alphas[ir];\
                n = ns[ir];\
                m = 1.0e0 - (1.0e0 / n);\
                s_res = s_ress[ir];\
                s_dif = s_difs[ir];\
\
                if (ppdat[ipp] >= 0.0)\
                  psdat[ips] = s_dif + s_res;\
                else\
                {\
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);\
                  psdat[ips] = s_dif / pow(1.0 + pow((alpha * head), n), m)\
                               + s_res;\
                }\
              });\
            }    /* End if clause */\
            else /* fcn = CALCDER */\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                ipp = SubvectorEltIndex(pp_sub, i, j, k);\
                ipd = SubvectorEltIndex(pd_sub, i, j, k);\
\
                alpha = alphas[ir];\
                n = ns[ir];\
                m = 1.0e0 - (1.0e0 / n);\
                s_res = s_ress[ir];\
                s_dif = s_difs[ir];\
\
                if (ppdat[ipp] >= 0.0)\
                  psdat[ips] = 0.0;\
                else\
                {\
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);\
                  psdat[ips] = (m * n * alpha * pow(alpha * head, (n - 1))) * s_dif\
                               / (pow(1.0 + pow(alpha * head, n), m + 1));\
                }\
              });\
            }   /* End else clause */\
          }     /* End subgrid loop */\
        }       /* End loop over regions */\
      }         /* End if data not from file */\
      else\
      {\
        gr_solid = ProblemDataGrDomain(problem_data);\
        n_values = dummy1->n_values;\
        alpha_values = dummy1->alpha_values;\
        s_res_values = dummy1->s_res_values;\
        s_sat_values = dummy1->s_sat_values;\
\
        ForSubgridI(sg, subgrids)\
        {\
          subgrid = SubgridArraySubgrid(subgrids, sg);\
          ps_sub = VectorSubvector(phase_saturation, sg);\
          pp_sub = VectorSubvector(phase_pressure, sg);\
          pd_sub = VectorSubvector(phase_density, sg);\
\
          n_values_sub = VectorSubvector(n_values, sg);\
          alpha_values_sub = VectorSubvector(alpha_values, sg);\
          s_res_values_sub = VectorSubvector(s_res_values, sg);\
          s_sat_values_sub = VectorSubvector(s_sat_values, sg);\
\
          ix = SubgridIX(subgrid);\
          iy = SubgridIY(subgrid);\
          iz = SubgridIZ(subgrid);\
\
          nx = SubgridNX(subgrid);\
          ny = SubgridNY(subgrid);\
          nz = SubgridNZ(subgrid);\
\
          r = SubgridRX(subgrid);\
\
          psdat = SubvectorData(ps_sub);\
          ppdat = SubvectorData(pp_sub);\
          pddat = SubvectorData(pd_sub);\
\
          n_values_dat = SubvectorData(n_values_sub);\
          alpha_values_dat = SubvectorData(alpha_values_sub);\
          s_res_values_dat = SubvectorData(s_res_values_sub);\
          s_sat_values_dat = SubvectorData(s_sat_values_sub);\
\
          if (fcn == CALCFCN)\
          {\
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
            {\
              ips = SubvectorEltIndex(ps_sub, i, j, k);\
              ipp = SubvectorEltIndex(pp_sub, i, j, k);\
              ipd = SubvectorEltIndex(pd_sub, i, j, k);\
\
              n_index = SubvectorEltIndex(n_values_sub, i, j, k);\
              alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k);\
              s_res_index = SubvectorEltIndex(s_res_values_sub, i, j, k);\
              s_sat_index = SubvectorEltIndex(s_sat_values_sub, i, j, k);\
\
              alpha = alpha_values_dat[alpha_index];\
              n = n_values_dat[n_index];\
              m = 1.0e0 - (1.0e0 / n);\
              s_res = s_res_values_dat[s_res_index];\
              s_sat = s_sat_values_dat[s_sat_index];\
\
              if (ppdat[ipp] >= 0.0)\
                psdat[ips] = s_sat;\
              else\
              {\
                head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);\
                psdat[ips] = (s_sat - s_res) /\
                             pow(1.0 + pow((alpha * head), n), m)\
                             + s_res;\
              }\
            });\
          }      /* End if clause */\
          else   /* fcn = CALCDER */\
          {\
            GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
            {\
              ips = SubvectorEltIndex(ps_sub, i, j, k);\
              ipp = SubvectorEltIndex(pp_sub, i, j, k);\
              ipd = SubvectorEltIndex(pd_sub, i, j, k);\
\
              n_index = SubvectorEltIndex(n_values_sub, i, j, k);\
              alpha_index = SubvectorEltIndex(alpha_values_sub, i, j, k);\
              s_res_index = SubvectorEltIndex(s_res_values_sub, i, j, k);\
              s_sat_index = SubvectorEltIndex(s_sat_values_sub, i, j, k);\
\
              alpha = alpha_values_dat[alpha_index];\
              n = n_values_dat[n_index];\
              m = 1.0e0 - (1.0e0 / n);\
              s_res = s_res_values_dat[s_res_index];\
              s_sat = s_sat_values_dat[s_sat_index];\
              s_dif = s_sat - s_res;\
\
              if (ppdat[ipp] >= 0.0)\
                psdat[ips] = 0.0;\
              else\
              {\
                head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);\
                psdat[ips] = (m * n * alpha * pow(alpha * head, (n - 1))) * s_dif\
                             / (pow(1.0 + pow(alpha * head, n), m + 1));\
              }\
            });\
          }     /* End else clause */\
        }       /* End subgrid loop */\
      }         /* End if data_from_file */\
      break;\
    }        /* End case 1 */\
\
      case 2: /* Haverkamp et.al. saturation curve */\
      {\
        double *alphas, *betas, *s_ress, *s_difs;\
        double head, alpha, beta, s_res, s_dif;\
\
        GetDummyType(Saturation, 2, (sat_xtra->data), dummy2);\
\
        num_regions = (dummy2->num_regions);\
        region_indices = (dummy2->region_indices);\
        alphas = (dummy2->alphas);\
        betas = (dummy2->betas);\
        s_ress = (dummy2->s_ress);\
        s_difs = (dummy2->s_difs);\
\
        for (ir = 0; ir < num_regions; ir++)\
        {\
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]);\
\
          ForSubgridI(sg, subgrids)\
          {\
            subgrid = SubgridArraySubgrid(subgrids, sg);\
            ps_sub = VectorSubvector(phase_saturation, sg);\
            pp_sub = VectorSubvector(phase_pressure, sg);\
            pd_sub = VectorSubvector(phase_density, sg);\
\
            ix = SubgridIX(subgrid);\
            iy = SubgridIY(subgrid);\
            iz = SubgridIZ(subgrid);\
\
            nx = SubgridNX(subgrid);\
            ny = SubgridNY(subgrid);\
            nz = SubgridNZ(subgrid);\
\
            r = SubgridRX(subgrid);\
\
            psdat = SubvectorData(ps_sub);\
            ppdat = SubvectorData(pp_sub);\
            pddat = SubvectorData(pd_sub);\
\
            if (fcn == CALCFCN)\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                ipp = SubvectorEltIndex(pp_sub, i, j, k);\
                ipd = SubvectorEltIndex(pd_sub, i, j, k);\
\
                alpha = alphas[ir];\
                beta = betas[ir];\
                s_res = s_ress[ir];\
                s_dif = s_difs[ir];\
\
                if (ppdat[ipp] >= 0.0)\
                  psdat[ips] = s_dif + s_res;\
                else\
                {\
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);\
                  psdat[ips] = alpha * s_dif / (alpha + pow(head, beta))\
                               + s_res;\
                }\
              });\
            }      /* End if clause */\
            else   /* fcn = CALCDER */\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                ipp = SubvectorEltIndex(pp_sub, i, j, k);\
                ipd = SubvectorEltIndex(pd_sub, i, j, k);\
\
                alpha = alphas[ir];\
                beta = betas[ir];\
                s_res = s_ress[ir];\
                s_dif = s_difs[ir];\
\
                if (ppdat[ipp] >= 0.0)\
                  psdat[ips] = 0.0;\
                else\
                {\
                  head = fabs(ppdat[ipp]) / (pddat[ipd] * gravity);\
                  psdat[ips] = alpha * s_dif * beta * pow(head, beta - 1)\
                               / pow((alpha + pow(head, beta)), 2);\
                }\
              });\
            }     /* End else clause */\
          }       /* End subgrid loop */\
        }         /* End loop over regions */\
        break;\
      }        /* End case 2 */\
\
      case 3: /* Data points for saturation curve */\
      {\
        GetDummyType(Saturation, 3, (sat_xtra->data), dummy3);\
\
        num_regions = (dummy3->num_regions);\
        region_indices = (dummy3->region_indices);\
\
        if (!amps_Rank(amps_CommWorld))\
          printf("Data curves for sats not yet supported.\n");\
\
        break;\
      }        /* End case 3 */\
\
      case 4: /* Polynomial function of pressure saturation curve */\
      {\
        int     *degrees, dg;\
        double **coefficients, *region_coeffs;\
\
        GetDummyType(Saturation, 4, (sat_xtra->data), dummy4);\
\
        num_regions = (dummy4->num_regions);\
        region_indices = (dummy4->region_indices);\
        degrees = (dummy4->degrees);\
        coefficients = (dummy4->coefficients);\
\
        for (ir = 0; ir < num_regions; ir++)\
        {\
          gr_solid = ProblemDataGrSolid(problem_data, region_indices[ir]);\
          region_coeffs = coefficients[ir];\
\
          ForSubgridI(sg, subgrids)\
          {\
            subgrid = SubgridArraySubgrid(subgrids, sg);\
\
            ps_sub = VectorSubvector(phase_saturation, sg);\
            pp_sub = VectorSubvector(phase_pressure, sg);\
\
\
            ix = SubgridIX(subgrid);\
            iy = SubgridIY(subgrid);\
            iz = SubgridIZ(subgrid);\
\
            nx = SubgridNX(subgrid);\
            ny = SubgridNY(subgrid);\
            nz = SubgridNZ(subgrid);\
\
            r = SubgridRX(subgrid);\
\
            psdat = SubvectorData(ps_sub);\
            ppdat = SubvectorData(pp_sub);\
\
            if (fcn == CALCFCN)\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                ipp = SubvectorEltIndex(pp_sub, i, j, k);\
\
                if (ppdat[ipp] == 0.0)\
                  psdat[ips] = region_coeffs[0];\
                else\
                {\
                  psdat[ips] = 0.0;\
                  for (dg = 0; dg < degrees[ir] + 1; dg++)\
                  {\
                    psdat[ips] += region_coeffs[dg] * pow(ppdat[ipp], dg);\
                  }\
                }\
              });\
            }      /* End if clause */\
            else   /* fcn = CALCDER */\
            {\
              GrGeomInLoop(i, j, k, gr_solid, r, ix, iy, iz, nx, ny, nz,\
              {\
                ips = SubvectorEltIndex(ps_sub, i, j, k);\
                ipp = SubvectorEltIndex(pp_sub, i, j, k);\
\
                if (ppdat[ipp] == 0.0)\
                  psdat[ips] = 0.0;\
                else\
                {\
                  psdat[ips] = 0.0;\
                  for (dg = 0; dg < degrees[ir] + 1; dg++)\
                  {\
                    psdat[ips] += region_coeffs[dg] * dg\
                                  * pow(ppdat[ipp], (dg - 1));\
                  }\
                }\
              });\
            }     /* End else clause */\
          }       /* End subgrid loop */\
        }         /* End loop over regions */\
        break;\
      }        /* End case 4 */\
\
      case 5: /* ParFlow binary file with spatially varying saturation values */\
      {\
        Vector *satRF;\
\
        GetDummyType(Saturation, 5, (sat_xtra->data), dummy5);\
\
        satRF = dummy5->satRF;\
\
        gr_domain = ProblemDataGrDomain(problem_data);\
\
        ForSubgridI(sg, subgrids)\
        {\
          subgrid = SubgridArraySubgrid(subgrids, sg);\
          ps_sub = VectorSubvector(phase_saturation, sg);\
          satRF_sub = VectorSubvector(satRF, sg);\
\
          ix = SubgridIX(subgrid);\
          iy = SubgridIY(subgrid);\
          iz = SubgridIZ(subgrid);\
\
          nx = SubgridNX(subgrid);\
          ny = SubgridNY(subgrid);\
          nz = SubgridNZ(subgrid);\
\
          /* RDF: assume resolution is the same in all 3 directions */\
          r = SubgridRX(subgrid);\
\
          psdat = SubvectorData(ps_sub);\
          satRFdat = SubvectorData(satRF_sub);\
\
          if (fcn == CALCFCN)\
          {\
            GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,\
            {\
              ips = SubvectorEltIndex(ps_sub, i, j, k);\
              ipRF = SubvectorEltIndex(satRF_sub, i, j, k);\
\
              psdat[ips] = satRFdat[ipRF];\
            });\
          }     /* End if clause */\
          else  /* fcn = CALCDER */\
          {\
            GrGeomInLoop(i, j, k, gr_domain, r, ix, iy, iz, nx, ny, nz,\
            {\
              ips = SubvectorEltIndex(ps_sub, i, j, k);\
              ipRF = SubvectorEltIndex(satRF_sub, i, j, k);\
\
              psdat[ips] = 0.0;\
            });\
          }    /* End else clause */\
        }      /* End subgrid loop */\
        break;\
      }        /* End case 5 */\
    }


#endif // _INLINE_DEFS_H
