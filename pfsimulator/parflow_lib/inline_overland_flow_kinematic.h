#ifndef _INLINE_OVERLAND_FLOW_KINEMATIC_H
#define _INLINE_OVERLAND_FLOW_KINEMATIC_H

#define OVERLAND_KINEMATIC_MODULE(module, grid, is, bc_struct, ipatch,\
                                  problem_data, pressure,\
                                  ke_der, kw_der, kn_der, ks_der,\
                                  kens_der, kwns_der, knns_der, ksns_der, \
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
    Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub;        \
                                                                        \
    Subgrid      *subgrid;                                              \
                                                                        \
    double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp;           \
                                                                        \
    double xdir, ydir;                                                  \
    double q_lo, q_mid, q_hi, qx_temp, qy_temp;                         \
    double q_v[4], slope_fx_lo, slope_fx_hi, slope_fx_mid;              \
    double slope_fy_lo, slope_fy_hi, slope_fy_mid, dx, dy;              \
    double coeff, Pmean, P2, P3, Pdel, Pcen;                            \
    double slope_mean, manning, s1, s2, Sf_mag;                         \
    double Press_x, Press_y, Sf_x, Sf_y, Sf_xo, Sf_yo;                  \
    double ov_epsilon;                                                  \
                                                                        \
    int ival, sy_v, step;                                               \
    int            *fdir;                                               \
                                                                        \
    int i, ii, j, k, ip, ip2, ip3, ip4, ip0, io, itop, k1x, k1y, ipp1, ippsy; \
    int i1, j1, k1, k0x, k0y, iojm1, iojp1, ioip1, ioim1;               \
    int gnx = BackgroundNX(GlobalsBackground);                          \
    int gny = BackgroundNY(GlobalsBackground);                          \
                                                                        \
    p_sub = VectorSubvector(pressure, sg);                              \
                                                                        \
    sx_sub = VectorSubvector(slope_x, sg);                              \
    sy_sub = VectorSubvector(slope_y, sg);                              \
    mann_sub = VectorSubvector(mannings, sg);                           \
    top_sub = VectorSubvector(top, sg);                                 \
                                                                        \
    pp = SubvectorData(p_sub);                                          \
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
    ov_epsilon = GetDoubleDefault("Solver.OverlandKinematic.Epsilon", 1.0e-5); \
                                                                        \
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
          k1x = pfmax((int)top_dat[itop + 1],0);                        \
          k1y = pfmax((int)top_dat[itop + sy_v],0);                     \
                                                                        \
          if (k1 >= 0)                                                  \
          {                                                             \
            ip = SubvectorEltIndex(p_sub, i, j, k1);                    \
            Sf_x = sx_dat[io];                                          \
            Sf_y = sy_dat[io];                                          \
            ipp1 = (int)SubvectorEltIndex(p_sub, i+1, j, k1x);          \
            ippsy = (int)SubvectorEltIndex(p_sub, i, j+1, k1y);         \
                                                                        \
            Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);           \
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
            if (k1 >= 0.0)                                              \
            {                                                           \
              Sf_x = sx_dat[io];                                        \
              Sf_y = sy_dat[io];                                        \
                                                                        \
              double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);  \
              if (Sf_mag < ov_epsilon)                                  \
                Sf_mag = ov_epsilon;                                    \
                                                                        \
              if (Sf_x > 0.0)                                           \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, i, j, k1);                \
                Press_x = pfmax((pp[ip]), 0.0);                         \
                qx_v[io - 1] = -(Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0)); \
              }                                                         \
            }                                                           \
          }                                                             \
                                                                        \
          if (k0y < 0.0)                                                \
          {                                                             \
            if (k1 >= 0.0)                                              \
            {                                                           \
              Sf_x = sx_dat[io];                                        \
              Sf_y = sy_dat[io];                                        \
                                                                        \
              double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);  \
              if (Sf_mag < ov_epsilon)                                  \
                Sf_mag = ov_epsilon;                                    \
                                                                        \
              if (Sf_y > 0.0)                                           \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, i, j, k1);                \
                Press_y = pfmax((pp[ip]), 0.0);                         \
                qy_v[io - sy_v] = -(Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0)); \
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
                                                                        \
            Sf_x = sx_dat[io];                                          \
            Sf_y = sy_dat[io];                                          \
                                                                        \
            Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);           \
            if (Sf_mag < ov_epsilon)                                    \
              Sf_mag = ov_epsilon;                                      \
                                                                        \
            Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ipp1]), 0.0)); \
            Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ippsy]), 0.0)); \
                                                                        \
            qx_temp = -(5.0 / 3.0) * (Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (2.0 / 3.0)); \
            qy_temp = -(5.0 / 3.0) * (Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (2.0 / 3.0)); \
                                                                        \
            ke_v[io] = pfmax(qx_temp, 0);                               \
            kw_v[io + 1] = -pfmax(-qx_temp, 0);                         \
            kn_v[io] = pfmax(qy_temp, 0);                               \
            ks_v[io + sy_v] = -pfmax(-qy_temp, 0);                      \
          }                                                             \
                                                                        \
          if (k0x < 0.0)                                                \
          {                                                             \
            if (k1 >= 0.0)                                              \
            {                                                           \
              Sf_x = sx_dat[io];                                        \
              Sf_y = sy_dat[io];                                        \
                                                                        \
              double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);  \
              if (Sf_mag < ov_epsilon)                                  \
                Sf_mag = ov_epsilon;                                    \
                                                                        \
              if (Sf_x > 0.0)                                           \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, i, j, k1);                \
                Press_x = pfmax((pp[ip]), 0.0);                         \
                qx_temp = -(5.0 / 3.0) * (Sf_x / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_x, (2.0 / 3.0)); \
                                                                        \
                kw_v[io] = qx_temp;                                     \
                ke_v[io - 1] = qx_temp;                                 \
              }                                                         \
            }                                                           \
          }                                                             \
                                                                        \
          if (k0y < 0.0)                                                \
          {                                                             \
            if (k1 >= 0.0)                                              \
            {                                                           \
              Sf_x = sx_dat[io];                                        \
              Sf_y = sy_dat[io];                                        \
                                                                        \
              double Sf_mag = RPowerR(Sf_x * Sf_x + Sf_y * Sf_y, 0.5);  \
              if (Sf_mag < ov_epsilon)                                  \
                Sf_mag = ov_epsilon;                                    \
                                                                        \
              if (Sf_y > 0.0)                                           \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, i, j, k1);                \
                Press_y = pfmax((pp[ip]), 0.0);                         \
                qy_temp = -(5.0 / 3.0) * (Sf_y / (RPowerR(fabs(Sf_mag), 0.5) * mann_dat[io])) * RPowerR(Press_y, (2.0 / 3.0)); \
                                                                        \
                ks_v[io] = qy_temp;                                     \
                kn_v[io - sy_v] = qy_temp;                              \
              }                                                         \
            }                                                           \
          }                                                             \
        }                                                               \
      });                                                               \
    }                                                                   \
  }

#endif // _INLINE_OVERLAND_FLOW_KINEMATIC_H
