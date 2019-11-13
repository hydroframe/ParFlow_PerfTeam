#ifndef _INLINE_OVERLAND_FLOW_H
#define _INLINE_OVERLAND_FLOW_H

#define OVERLAND_FLOW_MODULE(module, grid, is, bc_struct, ipatch,       \
                             problem_data, pressure, old_pressure,      \
                             ke_der, kw_der, kn_der, ks_der, qx, qy, ftype) \
  {                                                                     \
    /* This module does not make use of its PFModule pointer, don't bother setting that up */ \
    Vector      *slope_x = ProblemDataTSlopeX(problem_data);            \
    Vector      *slope_y = ProblemDataTSlopeY(problem_data);            \
    Vector      *mannings = ProblemDataMannings(problem_data);          \
    Vector      *top = ProblemDataIndexOfDomainTop(problem_data);       \
                                                                        \
    Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub;        \
                                                                        \
    double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp;           \
                                                                        \
    double xdir, ydir;                                                  \
    double q_lo, q_mid, q_hi;                                           \
    double q_v[3];                                                      \
                                                                        \
    int sg = is;                                                        \
    int fcn = ftype;                                                    \
    double *ke_v = ke_der;                                              \
    double *kw_v = kw_der;                                              \
    double *kn_v = kn_der;                                              \
    double *ks_v = ks_der;                                              \
    double *qx_v = qx;                                                  \
    double *qy_v = qy;                                                  \
                                                                        \
    int ival, sy_v, step;                                               \
    int            *fdir;                                               \
                                                                        \
    int i, ii, j, k, ip, io, itop;                                      \
    int i1, j1, k1;                                                     \
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
    sy_v = SubvectorNX(top_sub);                                        \
                                                                        \
    if (fcn == CALCFCN)                                                 \
    {                                                                   \
      if (qx_v == NULL || qy_v == NULL)  /* do not return velocity fluxes */ \
      {                                                                 \
        BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,   \
        {                                                               \
          if (fdir[2] == 1)                                             \
          {                                                             \
            io = SubvectorEltIndex(sx_sub, i, j, 0);                    \
            itop = SubvectorEltIndex(top_sub, i, j, 0);                 \
                                                                        \
            /* compute east and west faces */                           \
            /* First initialize velocities, q_v, for inactive region */ \
            q_v[0] = 0.0;                                               \
            q_v[1] = 0.0;                                               \
            q_v[2] = 0.0;                                               \
                                                                        \
            for (ii = -1; ii < 2; ii++)                                 \
            {                                                           \
              k1 = (int)top_dat[itop + ii];                             \
              if (k1 >= 0)                                              \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, (i + ii), j, k1);         \
                                                                        \
                if (sx_dat[io + ii] > 0.0)                              \
                  xdir = -1.0;                                          \
                else if (sx_dat[io + ii] < 0.0)                         \
                  xdir = 1.0;                                           \
                else                                                    \
                  xdir = 0.0;                                           \
                                                                        \
                q_v[ii + 1] = xdir * (RPowerR(fabs(sx_dat[io + ii]), 0.5) / mann_dat[io + ii]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0)); \
              }                                                         \
            }                                                           \
                                                                        \
            /* compute kw and ke - NOTE: io is for current cell */      \
            kw_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);        \
            ke_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);        \
                                                                        \
            /* compute north and south faces */                         \
            /* First initialize velocities, q_v, for inactive region */ \
            q_v[0] = 0.0;                                               \
            q_v[1] = 0.0;                                               \
            q_v[2] = 0.0;                                               \
                                                                        \
            for (ii = -1; ii < 2; ii++)                                 \
            {                                                           \
              step = ii * sy_v;                                         \
              k1 = (int)top_dat[itop + step];                           \
              if (k1 >= 0)                                              \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, i, (j + ii), k1);         \
                                                                        \
                if (sy_dat[io + step] > 0.0)                            \
                  ydir = -1.0;                                          \
                else if (sy_dat[io + step] < 0.0)                       \
                  ydir = 1.0;                                           \
                else                                                    \
                  ydir = 0.0;                                           \
                                                                        \
                q_v[ii + 1] = ydir * (RPowerR(fabs(sy_dat[io + step]), 0.5) / mann_dat[io + step]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0)); \
              }                                                         \
            }                                                           \
                                                                        \
            /* compute ks and kn - NOTE: io is for current cell */      \
            ks_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);        \
            kn_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);        \
          }                                                             \
        });                                                             \
      }                                                                 \
      else   /* return velocity fluxes */                               \
      {                                                                 \
        BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,   \
        {                                                               \
          if (fdir[2] == 1)                                             \
          {                                                             \
            io = SubvectorEltIndex(sx_sub, i, j, 0);                    \
            itop = SubvectorEltIndex(top_sub, i, j, 0);                 \
                                                                        \
            /* compute east and west faces */                           \
            /* First initialize velocities, q_v, for inactive region */ \
            q_v[0] = 0.0;                                               \
            q_v[1] = 0.0;                                               \
            q_v[2] = 0.0;                                               \
                                                                        \
            for (ii = -1; ii < 2; ii++)                                 \
            {                                                           \
              k1 = (int)top_dat[itop + ii];                             \
              if (k1 >= 0)                                              \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, (i + ii), j, k1);         \
                                                                        \
                if (sx_dat[io + ii] > 0.0)                              \
                  xdir = -1.0;                                          \
                else if (sx_dat[io + ii] < 0.0)                         \
                  xdir = 1.0;                                           \
                else                                                    \
                  xdir = 0.0;                                           \
                                                                        \
                q_v[ii + 1] = xdir * (RPowerR(fabs(sx_dat[io + ii]), 0.5) / mann_dat[io + ii]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0)); \
              }                                                         \
            }                                                           \
            qx_v[io] = q_v[1];                                          \
            /* compute kw and ke - NOTE: io is for current cell */      \
            kw_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);        \
            ke_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);        \
                                                                        \
            /* compute north and south faces */                         \
            /* First initialize velocities, q_v, for inactive region */ \
            q_v[0] = 0.0;                                               \
            q_v[1] = 0.0;                                               \
            q_v[2] = 0.0;                                               \
                                                                        \
            for (ii = -1; ii < 2; ii++)                                 \
            {                                                           \
              step = ii * sy_v;                                         \
              k1 = (int)top_dat[itop + step];                           \
              if (k1 >= 0)                                              \
              {                                                         \
                ip = SubvectorEltIndex(p_sub, i, (j + ii), k1);         \
                                                                        \
                if (sy_dat[io + step] > 0.0)                            \
                  ydir = -1.0;                                          \
                else if (sy_dat[io + step] < 0.0)                       \
                  ydir = 1.0;                                           \
                else                                                    \
                  ydir = 0.0;                                           \
                                                                        \
                q_v[ii + 1] = ydir * (RPowerR(fabs(sy_dat[io + step]), 0.5) / mann_dat[io + step]) * RPowerR(pfmax((pp[ip]), 0.0), (5.0 / 3.0)); \
              }                                                         \
            }                                                           \
            qy_v[io] = q_v[1];                                          \
            /* compute ks and kn - NOTE: io is for current cell */      \
            ks_v[io] = pfmax(q_v[0], 0.0) - pfmax(-q_v[1], 0.0);        \
            kn_v[io] = pfmax(q_v[1], 0.0) - pfmax(-q_v[2], 0.0);        \
          }                                                             \
        });                                                             \
      }                                                                 \
    }                                                                   \
    else  /* fcn == CALCDER: derivs of KE,KW,KN,KS w.r.t. current cell (i,j,k) */ \
    {                                                                   \
      if (qx_v == NULL || qy_v == NULL)  /* Do not return derivs of velocity fluxes */ \
      {                                                                 \
        BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,   \
        {                                                               \
          if (fdir[2] == 1)                                             \
          {                                                             \
            /* compute derivs for east and west faces */                \
                                                                        \
            /* current cell */                                          \
            io = SubvectorEltIndex(sx_sub, i, j, 0);                    \
            ip = SubvectorEltIndex(p_sub, i, j, k);                     \
                                                                        \
            if (sx_dat[io] > 0.0)                                       \
              xdir = -1.0;                                              \
            else if (sx_dat[io] < 0.0)                                  \
              xdir = 1.0;                                               \
            else                                                        \
              xdir = 0.0;                                               \
                                                                        \
            q_mid = xdir * (5.0 / 3.0) * (RPowerR(fabs(sx_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0)); \
            /* compute derivs of kw and ke - NOTE: io is for current cell */ \
            kw_v[io] = -pfmax(-q_mid, 0.0);                             \
            ke_v[io] = pfmax(q_mid, 0.0);                               \
                                                                        \
                                                                        \
            /* compute north and south faces */                         \
            if (sy_dat[io] > 0.0)                                       \
              ydir = -1.0;                                              \
            else if (sy_dat[io] < 0.0)                                  \
              ydir = 1.0;                                               \
            else                                                        \
              ydir = 0.0;                                               \
                                                                        \
            q_mid = ydir * (5.0 / 3.0) * (RPowerR(fabs(sy_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0)); \
            /* compute derivs of ks and kn - NOTE: io is for current cell */ \
            ks_v[io] = -pfmax(-q_mid, 0.0);                             \
            kn_v[io] = pfmax(q_mid, 0.0);                               \
          }                                                             \
        });                                                             \
      }                                                                 \
      else   /* return derivs of velocity fluxes */                     \
      {                                                                 \
        BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,   \
        {                                                               \
          if (fdir[2] == 1)                                             \
          {                                                             \
            /* compute derivs for east and west faces */                \
                                                                        \
            /* current cell */                                          \
            io = SubvectorEltIndex(sx_sub, i, j, 0);                    \
            ip = SubvectorEltIndex(p_sub, i, j, k);                     \
                                                                        \
            if (sx_dat[io] > 0.0)                                       \
              xdir = -1.0;                                              \
            else if (sx_dat[io] < 0.0)                                  \
              xdir = 1.0;                                               \
            else                                                        \
              xdir = 0.0;                                               \
                                                                        \
            q_mid = xdir * (5.0 / 3.0) * (RPowerR(fabs(sx_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0)); \
            qx_v[io] = q_mid;                                           \
            /* compute derivs of kw and ke - NOTE: io is for current cell */ \
            kw_v[io] = -pfmax(-q_mid, 0.0);                             \
            ke_v[io] = pfmax(q_mid, 0.0);                               \
                                                                        \
                                                                        \
            /* compute north and south faces */                         \
            if (sy_dat[io] > 0.0)                                       \
              ydir = -1.0;                                              \
            else if (sy_dat[io] < 0.0)                                  \
              ydir = 1.0;                                               \
            else                                                        \
              ydir = 0.0;                                               \
                                                                        \
            q_mid = ydir * (5.0 / 3.0) * (RPowerR(fabs(sy_dat[io]), 0.5) / mann_dat[io]) * RPowerR(pfmax((pp[ip]), 0.0), (2.0 / 3.0)); \
            qy_v[io] = q_mid;                                           \
            /* compute derivs of ks and kn - NOTE: io is for current cell */ \
            ks_v[io] = -pfmax(-q_mid, 0.0);                             \
            kn_v[io] = pfmax(q_mid, 0.0);                               \
          }                                                             \
        });                                                             \
      }                                                                 \
    }                                                                   \
  }


#endif // _INLINE_OVERLAND_FLOW_H
