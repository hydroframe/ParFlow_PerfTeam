#ifndef _INLINE_OVERLAND_FLOW_DIFFUSIVE_H
#define _INLINE_OVERLAND_FLOW_DIFFUSIVE_H

#define OVERLAND_DIFFUSIVE_MODULE(module, grid, is, bc_struct, ipatch,\
                                  problem_data, pressure, old_pressure,\
                                  ke_der, kw_der, kn_der, ks_der,\
                                  kens_der, kwns_der, knns_der, ksns_der,\
                                  qx, qy, ftype)                        \
  {                                                                     \
    int sg = is;                                                        \
    int fcn = ftype;                                                    \
    double *ke_v = ke_der;                                              \
    double *kw_v = kw_der;                                              \
    double *kn_v = kn_der;                                              \
    double *ks_v = ks_der;                                              \
    double *ke_vns = kens_der;                                          \
    double *kw_vns = kwns_der;                                          \
    double *kn_vns = knns_der;                                          \
    double *ks_vns = ksns_der;                                          \
    double *qx_v = qx;                                                  \
    double *qy_v = qy;                                                  \
    Vector      *slope_x = ProblemDataTSlopeX(problem_data);            \
    Vector      *slope_y = ProblemDataTSlopeY(problem_data);            \
    Vector      *mannings = ProblemDataMannings(problem_data);          \
    Vector      *top = ProblemDataIndexOfDomainTop(problem_data);       \
                                                                        \
    Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub, *op_sub; \
                                                                        \
    Subgrid      *subgrid;                                              \
                                                                        \
    double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp, *opp;     \
                                                                        \
    double xdir, ydir;                                                  \
    double q_lo, q_mid, q_hi;                                           \
    double q_v[4], slope_fx_lo, slope_fx_hi, slope_fx_mid;              \
    double slope_fy_lo, slope_fy_hi, slope_fy_mid, dx, dy;              \
    double coeff, Pmean, P2, P3, Pdel, Pcen;                            \
    double slope_mean, manning, s1, s2, Sf_mag;                         \
    double Press_x, Press_y, Sf_x, Sf_y, Sf_xo, Sf_yo;                  \
    double Pupx, Pupy, Pupox, Pupoy, Pdown, Pdowno;                     \
    double ov_epsilon;                                                  \
                                                                        \
    int ival, sy_v, step;                                               \
    int            *fdir;                                               \
                                                                        \
    int i, ii, j, k, ip, ip2, ip3, ip4, ip0, io, itop,k1x, k1y, ipp1, ippsy; \
    int i1, j1, k1, k0x, k0y, iojm1, iojp1, ioip1, ioim1;               \
    /* @RMM get grid from global (assuming this is comp grid) to pass to CLM */ \
    int gnx = BackgroundNX(GlobalsBackground);                          \
    int gny = BackgroundNY(GlobalsBackground);                          \
                                                                        \
    p_sub = VectorSubvector(pressure, sg);                              \
    op_sub = VectorSubvector(old_pressure, sg);                         \
    sx_sub = VectorSubvector(slope_x, sg);                              \
    sy_sub = VectorSubvector(slope_y, sg);                              \
    mann_sub = VectorSubvector(mannings, sg);                           \
    top_sub = VectorSubvector(top, sg);                                 \
                                                                        \
    pp = SubvectorData(p_sub);                                          \
    opp = SubvectorData(op_sub);                                        \
                                                                        \
    sx_dat = SubvectorData(sx_sub);                                     \
    sy_dat = SubvectorData(sy_sub);                                     \
    mann_dat = SubvectorData(mann_sub);                                 \
    top_dat = SubvectorData(top_sub);                                   \
                                                                        \
    subgrid = GridSubgrid(grid, sg);                                    \
    dx = SubgridDX(subgrid);                                            \
    dy = SubgridDY(subgrid);                                            \
                                                                        \
    sy_v = SubvectorNX(top_sub);                                        \
                                                                        \
    ov_epsilon = GetDoubleDefault("Solver.OverlandDiffusive.Epsilon", 1.0e-5); \
                                                                        \
    if (fcn == CALCFCN)                                                 \
    {                                                                   \
      BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, sg, \
      {                                                                 \
        if (fdir[2] == 1)                                               \
        {                                                               \
          io = SubvectorEltIndex(sx_sub, i, j, 0);                      \
          itop = SubvectorEltIndex(top_sub, i, j, 0);                   \
                                                                        \
          k1 = (int)top_dat[itop];                                      \
          k0x = (int)top_dat[itop - 1];                                 \
          k0y = (int)top_dat[itop - sy_v];                              \
          k1x = (int)top_dat[itop + 1];                                 \
          k1y = (int)top_dat[itop + sy_v];                              \
                                                                        \
          if (k1 >= 0)                                                  \
          {                                                             \
            ip = SubvectorEltIndex(p_sub, i, j, k1);                    \
            ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);          \
            ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);         \
            Pupx = pfmax(pp[ipp1], 0.0);                                \
            Pupy = pfmax(pp[ippsy], 0.0);                               \
            Pupox = pfmax(opp[ipp1], 0.0);                              \
            Pupoy = pfmax(opp[ippsy], 0.0);                             \
            Pdown = pfmax(pp[ip], 0.0);                                 \
            Pdowno = pfmax(opp[ip], 0.0);                               \
                                                                        \
            Sf_x = sx_dat[io] + (Pupx - Pdown) / dx;                    \
            Sf_y = sy_dat[io] + (Pupy - Pdown) / dy;                    \
                                                                        \
            Sf_xo = sx_dat[io] + (Pupox - Pdowno) / dx;                 \
            Sf_yo = sy_dat[io] + (Pupoy - Pdowno) / dy;                 \
                                                                        \
            Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5);       \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
                                                                        \
            Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ipp1]), 0.0)); \
            Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ippsy]), 0.0)); \
                                                                        \
            qx_v[io] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0)); \
            qy_v[io] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0)); \
          }                                                             \
                                                                        \
          if (k0x < 0.0)                                                \
          {                                                             \
            Press_x = pfmax((pp[ip]), 0.0);                             \
            Sf_x = sx_dat[io] + (Press_x - 0.0) / dx;                   \
                                                                        \
            Pupox = pfmax(opp[ip], 0.0);                                \
            Sf_xo = sx_dat[io] + (Pupox - 0.0) / dx;                    \
                                                                        \
            double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
            if (Sf_x > 0.0)                                             \
            {                                                           \
              qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0)); \
            }                                                           \
          }                                                             \
                                                                        \
          if (k0y < 0.0)                                                \
          {                                                             \
            Press_y = pfmax((pp[ip]), 0.0);                             \
            Sf_y = sy_dat[io] + (Press_y - 0.0) / dx;                   \
                                                                        \
            Pupoy = pfmax(opp[ip], 0.0);                                \
            Sf_yo = sy_dat[io] + (Pupoy - 0.0) / dx;                    \
                                                                        \
            double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
                                                                        \
            if (Sf_y > 0.0)                                             \
            {                                                           \
              qy_v[io - sy_v] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0)); \
            }                                                           \
                                                                        \
            if (k0x < 0.0)                                              \
            {                                                           \
              if (Sf_x > 0.0)                                           \
              {                                                         \
                qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0)); \
              }                                                         \
            }                                                           \
          }                                                             \
        }                                                               \
      });                                                               \
                                                                        \
      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,     \
      {                                                                 \
        if (fdir[2] == 1)                                               \
        {                                                               \
          io = SubvectorEltIndex(sx_sub, i, j, 0);                      \
          ke_v[io] = qx_v[io];                                          \
          kw_v[io] = qx_v[io - 1];                                      \
          kn_v[io] = qy_v[io];                                          \
          ks_v[io] = qy_v[io - sy_v];                                   \
        }                                                               \
      });                                                               \
    }                                                                   \
    else                                                                \
    {                                                                   \
      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,     \
      {                                                                 \
        if (fdir[2] == 1)                                               \
        {                                                               \
          io = SubvectorEltIndex(sx_sub, i, j, 0);                      \
          itop = SubvectorEltIndex(top_sub, i, j, 0);                   \
                                                                        \
          k1 = (int)top_dat[itop];                                      \
          k0x = (int)top_dat[itop - 1];                                 \
          k0y = (int)top_dat[itop - sy_v];                              \
          k1x = (int)top_dat[itop + 1];                                 \
          k1y = (int)top_dat[itop + sy_v];                              \
                                                                        \
          if (k1 >= 0)                                                  \
          {                                                             \
            ip = SubvectorEltIndex(p_sub, i, j, k1);                    \
            ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);          \
            ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);         \
            Pupx = pfmax(pp[ipp1], 0.0);                                \
            Pupy = pfmax(pp[ippsy], 0.0);                               \
            Pupox = pfmax(opp[ipp1], 0.0);                              \
            Pupoy = pfmax(opp[ippsy], 0.0);                             \
            Pdown = pfmax(pp[ip], 0.0);                                 \
            Pdowno = pfmax(opp[ip], 0.0);                               \
                                                                        \
            Sf_x = sx_dat[io] + (Pupx - Pdown) / dx;                    \
            Sf_y = sy_dat[io] + (Pupy - Pdown) / dy;                    \
                                                                        \
            Sf_xo = sx_dat[io] + (Pupox - Pdowno) / dx;                 \
            Sf_yo = sy_dat[io] + (Pupoy - Pdowno) / dy;                 \
                                                                        \
            Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5);       \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
                                                                        \
            if (Sf_x < 0)                                               \
            {                                                           \
              ke_v[io] = (5.0 / 3.0) * (-sx_dat[io] - (Pupx / dx)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pdown, (2.0 / 3.0)) + \
                         (8.0 / 3.0) * RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
                                                                        \
              kw_v[io + 1] = -RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
                                                                        \
              ke_vns[io] = kw_v[io + 1];                                \
              kw_vns[io + 1] = ke_v[io];                                \
            }                                                           \
                                                                        \
            if (Sf_x >= 0)                                              \
            {                                                           \
              ke_v[io] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
                                                                        \
              kw_v[io + 1] = (5.0 / 3.0) * (-sx_dat[io] + (Pdown / dx)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) - \
                             (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
                                                                        \
              ke_vns[io] = kw_v[io + 1];                                \
              kw_vns[io + 1] = ke_v[io];                                \
            }                                                           \
                                                                        \
            if (Sf_y < 0)                                               \
            {                                                           \
              kn_v[io] = (5.0 / 3.0) * (-sy_dat[io] - (Pupy / dy)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pdown, (2.0 / 3.0)) + \
                         (8.0 / 3.0) * RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy); \
                                                                        \
              ks_v[io + sy_v] = -RPowerR(Pdown, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy); \
                                                                        \
              kn_vns[io] = ks_v[io + sy_v];                             \
              ks_vns[io + sy_v] = kn_v[io];                             \
            }                                                           \
                                                                        \
            if (Sf_y >= 0)                                              \
            {                                                           \
              kn_v[io] = RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy); \
                                                                        \
              ks_v[io + sy_v] = (5.0 / 3.0) * (-sy_dat[io] + (Pdown / dy)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupy, (2.0 / 3.0)) - \
                                (8.0 / 3.0) * RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy); \
                                                                        \
              kn_vns[io] = ks_v[io + sy_v];                             \
              ks_vns[io + sy_v] = kn_v[io];                             \
            }                                                           \
          }                                                             \
                                                                        \
          if (k0x < 0.0)                                                \
          {                                                             \
            Pupx = pfmax((pp[ip]), 0.0);                                \
            Sf_x = sx_dat[io] + (Pupx - 0.0) / dx;                      \
                                                                        \
            Pupox = pfmax(opp[ip], 0.0);                                \
            Sf_xo = sx_dat[io] + (Pupox - 0.0) / dx;                    \
                                                                        \
            double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
                                                                        \
            if (Sf_x < 0)                                               \
            {                                                           \
              ke_v[io - 1] = 0.0;                                       \
              kw_v[io] = 0.0;                                           \
              kw_vns[io] = 0.0;                                         \
              ke_vns[io - 1] = 0.0;                                     \
            }                                                           \
                                                                        \
            if (Sf_x >= 0)                                              \
            {                                                           \
              ke_v[io - 1] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
              kw_v[io] = (5.0 / 3.0) * (-sx_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) - \
                         (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
              ke_vns[io - 1] = kw_v[io];                                \
              kw_vns[io] = ke_v[io - 1];                                \
            }                                                           \
          }                                                             \
                                                                        \
          if (k0y < 0.0)                                                \
          {                                                             \
            Pupy = pfmax((pp[ip]), 0.0);                                \
            Sf_y = sy_dat[io] + (Pupy - 0.0) / dy;                      \
                                                                        \
            Pupoy = pfmax(opp[ip], 0.0);                                \
            Sf_yo = sy_dat[io] + (Pupoy - 0.0) / dy;                    \
                                                                        \
            double Sf_mag = RPowerR(Sf_xo * Sf_xo + Sf_yo * Sf_yo, 0.5); \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
                                                                        \
            if (Sf_y < 0)                                               \
            {                                                           \
              kn_v[io - sy_v] = 0.0;                                    \
              ks_v[io] = 0.0;                                           \
              ks_vns[io] = 0.0;                                         \
              kn_vns[io - sy_v] = 0.0;                                  \
            }                                                           \
                                                                        \
            if (Sf_y >= 0)                                              \
            {                                                           \
              kn_vns[io - sy_v] = RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy); \
              ks_v[io] = (5.0 / 3.0) * (-sy_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupy, (2.0 / 3.0)) - \
                         (8.0 / 3.0) * RPowerR(Pupy, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dy); \
              kn_vns[io - sy_v] = ks_v[io];                             \
              ks_vns[io] = kn_v[io - sy_v];                             \
            }                                                           \
                                                                        \
            if (k0x < 0.0)                                              \
            {                                                           \
              if (Sf_x < 0)                                             \
              {                                                         \
                kn_v[io - sy_v] = 0.0;                                  \
                ks_v[io] = 0.0;                                         \
                ks_vns[io] = 0.0;                                       \
                kn_vns[io - sy_v] = 0.0;                                \
              }                                                         \
                                                                        \
              if (Sf_x >= 0)                                            \
              {                                                         \
                ke_v[io - 1] = RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
                kw_v[io] = (5.0 / 3.0) * (-sx_dat[io] + 0.0) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io]) * RPowerR(Pupx, (2.0 / 3.0)) - \
                           (8.0 / 3.0) * RPowerR(Pupx, (5.0 / 3.0)) / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io] * dx); \
                ke_vns[io - 1] = kw_v[io];                              \
                kw_vns[io] = ke_v[io - 1];                              \
              }                                                         \
            }                                                           \
          }                                                             \
        }                                                               \
      });                                                               \
    }                                                                   \
  }



#endif // _INLINE_OVERLAND_FLOW_DIFFUSIVE_H
