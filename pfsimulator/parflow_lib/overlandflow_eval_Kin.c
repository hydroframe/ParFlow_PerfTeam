/*BHEADER*********************************************************************
 *
 *  Copyright (c) 1995-2009, Lawrence Livermore National Security,
 *  LLC. Produced at the Lawrence Livermore National Laboratory. Written
 *  by the Parflow Team (see the CONTRIBUTORS file)
 *  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.
 *
 *  This file is part of Parflow. For details, see
 *  http://www.llnl.gov/casc/parflow
 *
 *  Please read the COPYRIGHT file or Our Notice and the LICENSE file
 *  for the GNU Lesser General Public License.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (as published
 *  by the Free Software Foundation) version 2.1 dated February 1999.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
 *  and conditions of the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA
 **********************************************************************EHEADER*/
/*****************************************************************************
*
*  This module computes the contributions for the spatial discretization of the
*  diffusive wave approximation for the overland flow boundary condition:KE,KW,KN,KS.
*
*  It also computes the derivatives of these terms for inclusion in the Jacobian.
*
* Could add a switch statement to handle the Kinemative wave approx. also.
* -DOK
*****************************************************************************/
#include "parflow.h"
#include "llnlmath.h"
//#include "llnltyps.h"
/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef void PublicXtra;

typedef void InstanceXtra;

/*---------------------------------------------------------------------
 * Define macros for function evaluation
 *---------------------------------------------------------------------*/
#define RPMean(a, b, c, d)   UpstreamMean(a, b, c, d)

/*-------------------------------------------------------------------------
 * OverlandFlowEval
 *-------------------------------------------------------------------------*/

void    OverlandFlowEvalKin(
                             Grid *       grid, /* data struct for computational grid */
                             int          sg, /* current subgrid */
                             BCStruct *   bc_struct, /* data struct of boundary patch values */
                             int          ipatch, /* current boundary patch */
                             ProblemData *problem_data, /* Geometry data for problem */
                             Vector *     pressure, /* Vector of phase pressures at each block */
                             double *     ke_v, /* return array corresponding to the east face KE  */
                             double *     kw_v, /* return array corresponding to the west face KW */
                             double *     kn_v, /* return array corresponding to the north face KN */
                             double *     ks_v, /* return array corresponding to the south face KS */
                             double *     ke_vns, /* return array corresponding to the nonsymetric east face KE derivative  */
                             double *     kw_vns, /* return array corresponding to the nonsymetricwest face KW derivative */
                             double *     kn_vns, /* return array corresponding to the nonsymetricnorth face KN derivative */
                             double *     ks_vns, /* return array corresponding to the nonsymetricsouth face KS derivative*/
                             double *     qx_v, /* return array corresponding to the flux in x-dir */
                             double *     qy_v, /* return array corresponding to the flux in y-dir */
                             int          fcn) /* Flag determining what to calculate
                                                * fcn = CALCFCN => calculate the function value
                                                * fcn = CALCDER => calculate the function
                                                *                  derivative */
{
  PFModule      *this_module = ThisPFModule;

  Vector      *slope_x = ProblemDataTSlopeX(problem_data);
  Vector      *slope_y = ProblemDataTSlopeY(problem_data);
  Vector      *mannings = ProblemDataMannings(problem_data);
  Vector      *top = ProblemDataIndexOfDomainTop(problem_data);

  // printf("overland_eval_diffusive called\n");
  Subvector     *sx_sub, *sy_sub, *mann_sub, *top_sub, *p_sub;

  Subgrid      *subgrid;

  double        *sx_dat, *sy_dat, *mann_dat, *top_dat, *pp;

  double xdir, ydir;
  double q_lo, q_mid, q_hi;
  double q_v[4], slope_fx_lo, slope_fx_hi, slope_fx_mid;
  double slope_fy_lo, slope_fy_hi, slope_fy_mid, dx, dy;
  double coeff, Pmean, P2, P3, Pdel, Pcen;
  double slope_mean, manning, s1, s2;
  double Press_x, Press_y, Sf_x, Sf_y, Sf_xo, Sf_yo;

  int ival, sy_v, step;
  int            *fdir;

  int i, ii, j, k, ip, ip2, ip3, ip4, ip0, io, itop;
  int i1, j1, k1, k0x, k0y, iojm1, iojp1, ioip1, ioim1;
  /* @RMM get grid from global (assuming this is comp grid) to pass to CLM */
  int gnx = BackgroundNX(GlobalsBackground);
  int gny = BackgroundNY(GlobalsBackground);

  p_sub = VectorSubvector(pressure, sg);

  sx_sub = VectorSubvector(slope_x, sg);
  sy_sub = VectorSubvector(slope_y, sg);
  mann_sub = VectorSubvector(mannings, sg);
  top_sub = VectorSubvector(top, sg);

  pp = SubvectorData(p_sub);

  sx_dat = SubvectorData(sx_sub);
  sy_dat = SubvectorData(sy_sub);
  mann_dat = SubvectorData(mann_sub);
  top_dat = SubvectorData(top_sub);

  subgrid = GridSubgrid(grid, sg);
  dx = SubgridDX(subgrid);
  dy = SubgridDY(subgrid);

  sy_v = SubvectorNX(top_sub);



  if (fcn == CALCFCN)
  {
    // printf("Overland_Diffusive --CACLFCN\n");
    //    if(qx_v == NULL || qy_v == NULL) /* do not return velocity fluxes */
    //       {
    BCStructPatchLoopOvrlnd(i, j, k, fdir, ival, bc_struct, ipatch, sg,
    {
      if (fdir[2] == 1)
      {
        io = SubvectorEltIndex(sx_sub, i, j, 0);
        itop = SubvectorEltIndex(top_sub, i, j, 0);

        k1 = (int)top_dat[itop];
        k0x = (int)top_dat[itop - 1];
        k0y = (int)top_dat[itop - sy_v];
        double ov_epsilon= 1.0e-5;
        //printf("i=%d j=%d k=%d k1=%d k0x=%d k0y=%d\n",i,j,k,k1, k0x, k0y);


        if (k1 >= 0)
        {
          ip = SubvectorEltIndex(p_sub, i, j, k1);
          Sf_x = sx_dat[io];
          Sf_y = sy_dat[io];

          double Sf_mag = RPowerR(Sf_x*Sf_x+Sf_y*Sf_y,0.5); //+ov_epsilon;
          if (Sf_mag < ov_epsilon)
          Sf_mag = ov_epsilon;

          Press_x = RPMean(-Sf_x, 0.0, pfmax((pp[ip]), 0.0), pfmax((pp[ip+1]), 0.0));
          Press_y = RPMean(-Sf_y, 0.0, pfmax((pp[ip]), 0.0),pfmax((pp[ip+sy_v]), 0.0));

          qx_v[io] = -(Sf_x / (RPowerR(fabs(Sf_mag),0.5)*mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
          qy_v[io] = -(Sf_y / (RPowerR(fabs(Sf_mag),0.5)*mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));

        }

        //fix for lower x boundary
        if (k0x < 0.0) {
            if (k1 >= 0.0) {
              Sf_x = sx_dat[io];
              Sf_y = sy_dat[io];

              double Sf_mag = RPowerR(Sf_x*Sf_x+Sf_y*Sf_y,0.5); //+ov_epsilon;
              if (Sf_mag < ov_epsilon)
              Sf_mag = ov_epsilon;

              if (Sf_x > 0.0) {
                ip = SubvectorEltIndex(p_sub, i, j, k1);
                Press_x = pfmax((pp[ip]), 0.0);
                qx_v[io-1] = -(Sf_x / (RPowerR(fabs(Sf_mag),0.5)*mann_dat[io])) * RPowerR(Press_x, (5.0 / 3.0));
                //printf("i=%d j=%d k=%d k1=%d k0x=%d k0y=%d Sf_y=%f qx_v=%f\n",i,j,k,k1, k0x, k0y, Sf_x, qx_v[io-1]);
              }
            }
          }

          //fix for lower y boundary
          if (k0y < 0.0) {
              if (k1 >= 0.0) {

                Sf_x = sx_dat[io];
                Sf_y = sy_dat[io];

                double Sf_mag = RPowerR(Sf_x*Sf_x+Sf_y*Sf_y,0.5); //+ov_epsilon;
                if (Sf_mag < ov_epsilon)
                Sf_mag = ov_epsilon;

                if (Sf_y > 0.0) {
                  ip = SubvectorEltIndex(p_sub, i, j, k1);
                  Press_y = pfmax((pp[ip]), 0.0);
                  qy_v[io-sy_v] = -(Sf_y / (RPowerR(fabs(Sf_mag),0.5)*mann_dat[io])) * RPowerR(Press_y, (5.0 / 3.0));
                  //printf("i=%d j=%d k=%d k1=%d k0x=%d k0y=%d Sf_y=%f qy_v=%f\n",i,j,k,k1, k0x, k0y, Sf_y, qy_v[io-sy_v]);

                }
              }
            }
       }
    });

    BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
    {
      if (fdir[2] == 1)
      {
        io = SubvectorEltIndex(sx_sub, i, j, 0);
        ke_v[io] = qx_v[io];
        kw_v[io] = qx_v[io-1];
        kn_v[io] = qy_v[io];
        ks_v[io] = qy_v[io-sy_v];
        //printf("i=%d j=%d k=%d ke_v=%d kw_v=%d kn_v=%d ks_v=%f\n",i,j,k,ke_v[io],kw_v[io],kn_v[io],ks_v[io]);
      }
    });

  }
  else          //fcn = CALCDER calculates the derivs of KE KW KN KS wrt to current cell (i,j,k)
  {
    //     printf("ELSE... CALCDER");
    //  if(qx_v == NULL || qy_v == NULL) /* do not return velocity fluxes */
    // {
    BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
    {
      if (fdir[2] == 1)
      {
        //printf("Made it inside CALCDER");

        io = SubvectorEltIndex(sx_sub, i, j, 0);
        itop = SubvectorEltIndex(top_sub, i, j, 0);
        /* Current node */
        k1 = (int)top_dat[itop];
        ip0 = SubvectorEltIndex(p_sub, i, j, k1);
        manning = mann_dat[io];

        /*KW  - look at nodes i-1 and i*/
        k1 = (int)top_dat[itop - 1];
        ip = SubvectorEltIndex(p_sub, (i - 1), j, k1);


        /*Calcualte Friction Slope */
        if (i > 0)
        {
          slope_mean = sx_dat[io - 1];
          slope_fx_lo = slope_mean + (((pfmax((pp[ip0]), 0.0)) - (pfmax((pp[ip]), 0.0))) / dx);
        }
        else
        {
          slope_mean = sx_dat[io];
          slope_fx_lo = slope_mean;
        }

        manning = (1.0 + 0.0000001 / fabs(slope_fx_lo)) * mann_dat[io];
        Pcen = pfmax(pp[ip0], 0.0);                 //pressure of current cell
        Pdel = pfmax(pp[ip], 0.0);                  // pressure cell to the west

        /* Caluculate Derivative */
        if (fabs(slope_fx_lo) < 0.0000001)
        {
          kw_vns[io] = 0;
          kw_v[io] = 0;
        }

        else if (slope_fx_lo > 0.0)
        {
          xdir = -1.0;
          kw_vns[io] = xdir * (-1 / (2 * dx * manning)) * RPowerR(fabs(slope_fx_lo), -0.5) * RPowerR(Pcen, (5.0 / 3.0));

          kw_v[io] = xdir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fx_lo), 0.5) * RPowerR(Pcen, (2.0 / 3.0)) +
                             1 / (2 * dx * manning) * RPowerR(fabs(slope_fx_lo), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));
        }
        else if (slope_fx_lo < 0.0)
        {
          xdir = 1.0;
          /*This is dfi-1/di-1 East*/
          kw_vns[io] = xdir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fx_lo), 0.5) * RPowerR(Pdel, (2.0 / 3.0)) -
                               1 / (2 * dx * manning) * RPowerR(fabs(slope_fx_lo), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));
          /* This is dfi-1/di East */
          kw_v[io] = xdir * ((1 / (2 * dx * manning)) * RPowerR(fabs(slope_fx_lo), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));
          //kw_v[io]=0;
        }
        //else{
        //  kw_vns[io]=0;
        // kw_v[io]=0;
        //}

        //printf("WEST: i %d j %d %4.5e %4.5e %4.5e %4.5e %4.5e %4.5e \n", i, j, Pcen, Pdel, slope_fx_lo, slope_mean, kw_v[io], kw_vns[io]);


        /* KE - look at nodes i+1 and i */
        k1 = (int)top_dat[itop + 1];
        ip2 = SubvectorEltIndex(p_sub, (i + 1), j, k1);

        /*Calcualte Friction Slope */
        if (i < gnx - 1)
        {
          slope_mean = sx_dat[io];
          slope_fx_hi = slope_mean + (((pfmax((pp[ip2]), 0.0)) - (pfmax((pp[ip0]), 0.0))) / dx);
        }
        else
        {
          slope_mean = sx_dat[io];
          slope_fx_hi = slope_mean;
        }

        manning = (1.0 + 0.0000001 / fabs(slope_fx_hi)) * mann_dat[io];
        Pcen = pfmax(pp[ip0], 0.0);                 //pressure of current cel
        Pdel = pfmax(pp[ip2], 0.0);                  // pressure cell to the east

        /* Caluculate Derivative */
        if (fabs(slope_fx_hi) < 0.0000001)
        {
          ke_vns[io] = 0;
          ke_v[io] = 0;
        }

        else if (slope_fx_hi > 0.0)
        {
          xdir = -1.0;
          /*This is dfi+1/di+1 for kw */
          ke_vns[io] = xdir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fx_hi), 0.5) * RPowerR(Pdel, (2.0 / 3.0)) +
                               1 / (2 * dx * manning) * RPowerR(fabs(slope_fx_hi), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));

          /* This is dfi+1/di for kw */
          ke_v[io] = xdir * ((-1 / (2 * dx * manning)) * RPowerR(fabs(slope_fx_hi), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));
          //  ke_v[io]=0;
        }
        else if (slope_fx_hi < 0.0)
        {
          xdir = 1.0;
          ke_vns[io] = xdir * ((1 / (2 * dx * manning)) * RPowerR(fabs(slope_fx_hi), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));

          ke_v[io] = xdir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fx_hi), 0.5) * RPowerR(Pcen, (2.0 / 3.0)) -
                             1 / (2 * dx * manning) * RPowerR(fabs(slope_fx_hi), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));
        }
        // else{
        //   ke_vns[io]=0;
        //  ke_v[io]=0;
        //}

        // printf("i %d j %d %4.5f %4.5f %4.5f %4.5f %4.5e %4.5e %4.5e %4.5e \n", i, j, slope_mean, slope_fx_lo, Pcen, Pdel, kw_v[io], kw_vns[io], ke_v[io], ke_vns[io]);
        //printf("EAST: i %d j %d %4.5e %4.5e %4.5e %4.5e %4.5e %4.5e \n", i, j,  Pcen, Pdel, slope_fx_hi, slope_mean, ke_v[io], ke_vns[io]);


        /*KS  - look at nodes j-1 and j*/
        k1 = (int)top_dat[itop - sy_v];
        ip3 = SubvectorEltIndex(p_sub, i, (j - 1), k1);

        /*Calcualte Friction Slope */
        if (j > 0)
        {
          slope_mean = sy_dat[io - sy_v];
          slope_fy_lo = slope_mean + (((pfmax((pp[ip0]), 0.0)) - (pfmax((pp[ip3]), 0.0))) / dy);
        }
        else
        {
          slope_mean = sy_dat[io];
          slope_fy_lo = slope_mean;
        }
        manning = (1.0 + 0.0000001 / fabs(slope_fy_lo)) * mann_dat[io];
        Pcen = pfmax(pp[ip0], 0.0);                 //pressure of current cel
        Pdel = pfmax(pp[ip3], 0.0);                  // pressure cell to the south

        /* Caluculate Derivative */
        if (fabs(slope_fy_lo) < 0.0000001)
        {
          ks_vns[io] = 0;
          ks_v[io] = 0;
        }
        else if (slope_fy_lo > 0.0)
        {
          ydir = -1.0;
          ks_vns[io] = ydir * ((-1 / (2 * dy * manning)) * RPowerR(fabs(slope_fy_lo), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));

          ks_v[io] = ydir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fy_lo), 0.5) * RPowerR(Pcen, (2.0 / 3.0)) +
                             1 / (2 * dy * manning) * RPowerR(fabs(slope_fy_lo), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));
        }
        else if (slope_fy_lo < 0.0)
        {
          ydir = 1.0;
          ks_vns[io] = ydir * (5 / (3 * manning) * RPowerR(fabs(slope_fy_lo), 0.5) * RPowerR(Pdel, (2.0 / 3.0)) -
                               1 / (2 * dy * manning) * RPowerR(fabs(slope_fy_lo), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));

          //ks_v[io]=0.0;
          ks_v[io] = ydir * ((1 / (2 * dy * manning)) * RPowerR(fabs(slope_fy_lo), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));
        }
        //else{
        //  ks_vns[io]=0;
        // ks_v[io]=0;
        // }

        //printf("SOUTH: i %d j %d %4.5e %4.5e %4.5e %4.5e %4.5e %4.5e \n", i, j,  Pcen, Pdel, slope_fy_lo, slope_mean,  ks_v[io], ks_vns[io]);


        /* KN - look at nodes j+1 and j */
        k1 = (int)top_dat[itop + sy_v];
        ip4 = SubvectorEltIndex(p_sub, i, (j + 1), k1);

        /*Calcualte Friction Slope */
        if (j < gny - 1)
        {
          slope_mean = sy_dat[io];
          slope_fy_hi = slope_mean + (((pfmax((pp[ip4]), 0.0)) - (pfmax((pp[ip0]), 0.0))) / dy);
        }
        else
        {
          slope_mean = sy_dat[io];
          slope_fy_hi = slope_mean;
        }

        manning = (1.0 + 0.0000001 / fabs(slope_fy_hi)) * mann_dat[io];
        Pcen = pfmax(pp[ip0], 0.0);                 //pressure of current cel
        Pdel = pfmax(pp[ip4], 0.0);                  // pressure cell to the east

        /* Caluculate Derivative */
        if (fabs(slope_fy_hi) < 0.0000001)
        {
          kn_vns[io] = 0;
          kn_v[io] = 0;
        }

        else if (slope_fy_hi > 0.0)
        {
          ydir = -1.0;
          kn_vns[io] = ydir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fy_hi), 0.5) * RPowerR(Pdel, (2.0 / 3.0)) +
                               1 / (2 * dy * manning) * RPowerR(fabs(slope_fy_hi), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));

          //kn_v[io]=0.0;

          kn_v[io] = ydir * ((-1 / (2 * dy * manning)) * RPowerR(fabs(slope_fy_hi), -0.5) * RPowerR(Pdel, (5.0 / 3.0)));
        }
        else if (slope_fy_hi < 0.0)
        {
          ydir = 1.0;
          kn_vns[io] = ydir * ((1 / (2 * dy * manning)) * RPowerR(fabs(slope_fy_hi), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));

          kn_v[io] = ydir * ((5 / (3 * manning)) * RPowerR(fabs(slope_fy_hi), 0.5) * RPowerR(Pcen, (2.0 / 3.0)) -
                             1 / (2 * dy * manning) * RPowerR(fabs(slope_fy_hi), -0.5) * RPowerR(Pcen, (5.0 / 3.0)));
        }
        //else{
        //   kn_vns[io]=0;
        //   kn_v[io]=0;
        //}
        //printf("NORTH: i %d j %d %4.5e %4.5e %4.5e %4.5e %4.5e %4.5e \n", i, j,  Pcen, Pdel, slope_fy_hi, slope_mean,  kn_v[io], kn_vns[io]);
      }
    });
    //}
  }

  //  else /* return velocity fluxes */
  //  {
  //       printf("returning velocity fluxes \n");
  //  }

//	    BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
//	    {
//	       if (fdir[2] == 1)
//	       {
//	                   io = SubvectorEltIndex(sx_sub,i, j, 0);
//	                   itop = SubvectorEltIndex(top_sub, i, j, 0);
//
//	                   /* compute east and west faces */
//	                   /* First initialize velocities, q_v, for inactive region */
//	                   q_v[0] = 0.0;
//	                   q_v[1] = 0.0;
//	                   q_v[2] = 0.0;
//			   q_v[3] = 0.0;
//
//                          /* precompute some data for current node, to be reused */
//                         ip0 = SubvectorEltIndex(p_sub, i, j, k);
//                         //coeff = RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0));
//
//                        /* Dealing with Kw - look at nodes i-1 and i */
//                         k1 = (int)top_dat[itop-1];
//                         if(k1 >= 0)
//                         {                                //@RMM - moved down one line
//
//                                 ip = SubvectorEltIndex(p_sub,(i-1), j, k1);
//
//			  /* compute the friction slope for west node
//			   * this handles the pressure gradient correction
//			   * for the slope in Manning's equation
//			  */
//                           // @RMM  slope is PP(i)-PP(i-1)
//			      slope_fx_lo = sx_dat[io-1] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip]),0.0)))/dx);
//                           //@RMM upwind pressure based on slope direction
//			      if(slope_fx_lo > 0.0) {
//				     w = pp[ip];
//		                    xdir = -1.0;
//			      }
//			      else if(slope_fx_lo < 0.0) {
//                                  xdir = 1.0;
//			            w = pp[ip0];
//			      }
//                               else
//                                  xdir = 0.0;
//
//			      q_v[0] = xdir * (RPowerR(fabs(slope_fx_lo),0.5) / mann_dat[io-1]) * RPowerR(pfmax(w,0.0),(5.0/3.0));  //@RMM could use PMean here in place of w
//			   }
//
//			  /* compute the friction slope for center node */
//
//			    /*  slope_fx_mid = sx_dat[io] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip]),0.0)))/dx);
//			      if(slope_fx_mid > 0.0) {
//		                    xdir = -1.0;
//			      }
//			      else if(slope_fx_mid < 0.0) {
//	                            xdir = 1.0;
//			      }
//                            else
//                                  xdir = 0.0;
//
//			      q_v[1] = xdir * (RPowerR(fabs(slope_fx_mid),0.5) / mann_dat[io]) * RPowerR(pfmax(pp[ip0],0.0),(5.0/3.0));
//			      @RMM  */
//
//
//	                   /* compute kw - NOTE: io is for current cell */
//			   kw_v[io] = q_v[0];  //@RMM-- just i - (i-1) quantities
//
//			   /* initially set velocity flux */
//                         qx_v[io] = q_v[0];
//
//                        /* Dealing with Ke - look at nodes i and i+1 */
//                         /* Dealing with Ke - look at nodes i and i+1 */
//                         k1 = (int)top_dat[itop+1];
//			   ip2 = SubvectorEltIndex(p_sub,(i+1), j, k1);
//
//                         if(k1 >= 0)
//                         {
//
//			  /* compute the friction slope for east node
//			   * this handles the pressure gradient correction
//			   * for the slope in Manning's equation
//			        @RMM take PP(i+1)=PP(i)       */
//			      slope_fx_hi = sx_dat[io] - (((pfmax((pp[ip2]),0.0)) - (pfmax((pp[ip0]),0.0)))/dx);  //@RMM
//			      if(slope_fx_hi > 0.0) {
//		                    xdir = -1.0;
//				    w = pp[ip0];
//			      }
//			      else if(slope_fx_hi < 0.0) {
//                                  xdir = 1.0;
//				    w = pp[ip2];
//			      }
//                            else
//                                  xdir = 0.0;
//
//			      q_v[2] = xdir * (RPowerR(fabs(slope_fx_hi),0.5) / mann_dat[io+1]) * RPowerR(pfmax(w,0.0),(5.0/3.0));
//			   }
//
//			   /* compute the friction slope for center node */
//			   /*  @RMM
//			      slope_fx_mid = sx_dat[io] - (((pfmax((pp[ip2]),0.0)) - (pfmax((pp[ip0]),0.0)))/dx);
//
//			      if(slope_fx_mid > 0.0) {
//		                    xdir = -1.0;
//			      }
//			      else if(slope_fx_mid < 0.0) {
//	                            xdir = 1.0;
//			      }
//			      else
//                                  xdir = 0.0;
//
//			      q_v[3] = xdir * (RPowerR(fabs(slope_fx_mid),0.5) / mann_dat[io]) * RPowerR(pfmax(pp[ip0],0.0),(5.0/3.0));
//				*/
//
//	                   /* compute ke - NOTE: @RMM modified stencil to not upwind Q's but to upwind pressure-head */
//			   ke_v[io] = q_v[2];
//
//			   /* update velocity flux to take on upwinded value or else
//			    * the value associated with evaluating the pressure gradient
//			    * using nodes i and i+1
//			    */
//			   if(pfmax(q_v[2], 0.0))
//			      qx_v[io] = q_v[2];
//
//	               /* compute north and south faces */
//	               /* First initialize velocities, q_v, for inactive region */
//	                   q_v[4] = 0.0;
//		           q_v[5] = 0.0;
//	                   q_v[6] = 0.0;
//	                   q_v[7] = 0.0;
//
//                        /* Dealing with Ks - look at nodes j-1 and j */
//                         k1 = (int)top_dat[itop-sy_v];
//                         ip3 = SubvectorEltIndex(p_sub,i, (j-1), k1);
//
//			   if(k1 >= 0)
//                         {
//
//			  /* compute the friction slope for south node
//			   * this handles the pressure gradient correction
//			   * for the slope in Manning's equation
//			         now take slope based on PP(j)-PP(j-1)    */
//			      slope_fy_lo = sy_dat[io-sy_v] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip3]),0.0)))/dy);
//			      if(slope_fy_lo > 0.0) {
//		                    ydir = -1.0;
//				    w = pp[ip3];
//			      }
//			      else if(slope_fy_lo < 0.0) {
//                               ydir = 1.0;
//				 w = pp[ip0];
//			      }
//                            else
//                                  ydir = 0.0;
//
//			      q_v[4] = ydir * (RPowerR(fabs(slope_fy_lo),0.5) / mann_dat[io-sy_v]) * RPowerR(pfmax(w,0.0),(5.0/3.0));
//			        }
//			  /* compute the friction slope for center node */
//			     /* slope_fy_mid = sy_dat[io] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip3]),0.0)))/dy);
//			      if(slope_fy_mid > 0.0) {
//		                    ydir = -1.0;
//			      }
//			      else if(slope_fy_mid < 0.0) {
//                                  ydir = 1.0;
//			      }
//			      else
//                                  ydir = 0.0;
//
//			      q_v[5] = ydir * (RPowerR(fabs(slope_fy_mid),0.5) / mann_dat[io]) * RPowerR(pfmax(pp[ip0],0.0),(5.0/3.0));
//				*/
//
//	                   /* compute ks - NOTE: io is for current cell */
//			   ks_v[io] = q_v[4]; //@RMM
//
//			   /* initially set velocity flux */
//                         qy_v[io] = q_v[4];
//
//                        /* Dealing with Ke - look at nodes i and i+1 */
//                         k1 = (int)top_dat[itop+sy_v];
//                          ip4 = SubvectorEltIndex(p_sub,i, (j+1), k1);
//			   if(k1 >= 0)
//                         {
//
//			  /* compute the friction slope for east node
//			   * this handles the pressure gradient correction
//			   * for the slope in Manning's equation
//			  */
//			      slope_fy_hi = sy_dat[io+sy_v] - (((pfmax((pp[ip4]),0.0)) - (pfmax((pp[ip0]),0.0)))/dy);
//			      if(slope_fy_hi > 0.0) {
//		                    ydir = -1.0;
//				    w = pp[ip0];
//			      }
//			      else if(slope_fy_hi < 0.0) {
//	                            ydir = 1.0;
//				    w = pp[ip4];
//			      }
//			      else
//                                  ydir = 0.0;
//
//			      q_v[6] = ydir * (RPowerR(fabs(slope_fy_hi),0.5) / mann_dat[io+sy_v]) * RPowerR(pfmax(w,0.0),(5.0/3.0));
//			   }
//			  /* compute the friction slope for center node */
//			    /*  slope_fy_mid = sy_dat[io] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip4]),0.0)))/dy);
//			      if(slope_fy_mid > 0.0) {
//		                    ydir = -1.0;
//			      }
//			      else if(slope_fy_mid < 0.0) {
//	                            ydir = 1.0;
//			      }
//                               else
//                                  ydir = 0.0;
//
//			      q_v[7] = ydir * (RPowerR(fabs(slope_fy_mid),0.5) / mann_dat[io]) * RPowerR(pfmax(pp[ip0],0.0),(5.0/3.0));
//				*/
//
//	                   /* compute kn - NOTE: io is for current cell */
//			   kn_v[io] = q_v[6]; //@RMM
//
//			   printf(" %4.5f %4.5f %4.5f %4.5f \n", q_v[0], q_v[2],q_v[4], q_v[6]);
//
//			   /* update velocity flux to take on upwinded value or else
//			    * the value associated with evaluating the pressure gradient
//			    * using nodes i and i+1
//			    */
//			   if(pfmax(q_v[6], 0.0))
//			    qy_v[io] = q_v[6];
//
//	         }
//
//	     });
//
//	  }
//      }
/* ***********************************************************************
 * Still working on the jacobian!!!!
 ************************************************************************** */
//      else // fcn = CALCDER calculates the derivs of KE KW KN KS wrt to current cell (i,j,k)
//      {
//              if(qx_v == NULL || qy_v == NULL) /* do not return velocity fluxes */
//         {
//          BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
//	    {
//	       if (fdir[2] == 1)
//	       {
//	                   /* First initialize velocities, q_v, for inactive region */
//	                   q_v[0] = 0.0;
//	                   q_v[1] = 0.0;
//	                   q_v[2] = 0.0;
//	                   //printf("set velocities \n");
//	                   io = SubvectorEltIndex(sx_sub,i, j, 0);
//	                   ip = SubvectorEltIndex(top_sub, i, j, 0);
//
//                          /* precompute some data for current node, to be reused */
//                         ip0 = SubvectorEltIndex(p_sub, i, j, k);
//
//                         //   printf("initial data is recorded \n");
//		      //   ***********************************************************
//		      //
//	              /* compute east and  west face */
//	              //
//	              //   ***********************************************************
//
//                        /* Dealing with Kw - look at nodes i-1 and i */
//                         //printf("Why doesn't this work \n");
//                         //k1 = fabs((int)top_dat[itop-1]);
//                         //printf(" %4.2f \n", top_dat[itop-1]);
//                         //if(k1 >= 0)
//                         //{
//                            ip = SubvectorEltIndex(p_sub,(i-1), j, k);
//                        /* Compute the derivative wrt to node i,j,k  */
//                        // printf("calculating slope \n");
//                        slope_fx_lo = sx_dat[io] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip]),0.0)))/dx);
//		      if(slope_fx_lo > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_lo < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//                         q_v[0] = xdir * (1.0 / mann_dat[io]) * (((5.0/3.0)*RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) *(RPowerR(fabs(slope_fx_lo),0.5)))- ((1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fx_lo),-0.5))));
//
//
//			  /* compute the derivative wrt to node i-1,j,k */
//			      slope_fx_mid = sx_dat[io-1] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip]),0.0)))/dx);
//			      if(slope_fx_mid > 0.0)
//	                    xdir = -1.0;
//				 else if(slope_fx_mid < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//			      q_v[1] = xdir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fx_mid),-0.5));

  //}
  /* compute kw - NOTE: io is for current cell */

//			   kw_v[io] = pfmax(q_v[1],0.0) - pfmax(-q_v[0],0.0);

  /* Dealing with EAST NODE Ke - look at nodes i and i+1 */
//
//                         //k1 = fabs((int)top_dat[itop+1]);
//                         //if(k1 >= 0)
//                         //{
//                            ip2 = SubvectorEltIndex(p_sub,(i+1), j, k);
//
//                            /* compute the derivative wrt i,j,k*/
//			      slope_fx_hi = sx_dat[io] - (((pfmax((pp[ip2]),0.0)) - (pfmax((pp[ip0]),0.0)))/dx);
//			      if(slope_fx_hi > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_hi < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//			   q_v[2] = xdir * (1.0 / mann_dat[io]) * (((5.0/3.0)*RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) *(RPowerR(fabs(slope_fx_hi),0.5)))- ((1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fx_hi),-0.5))));
//
//			   /* Computer the derivative for the i+1,j,k */
//
//			   slope_fx_mid = sx_dat[io+1] - (((pfmax((pp[ip2]),0.0)) - (pfmax((pp[ip0]),0.0)))/dx);
//			      if(slope_fx_mid > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_mid < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//			   q_v[1] = xdir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fx_mid),-0.5));
//
//			   //}
//			   /* compute ke - NOTE: io is for current cell */
//
//			   ke_v[io] = pfmax(q_v[2],0.0) - pfmax(-q_v[1],0.0);
//
//
//		      //   ***********************************************************
//		      //
//	              /* compute north and south faces */
//	              //
//	              //   ***********************************************************
//
//	               /* First initialize velocities, q_v, for inactive region */
//	                   q_v[3] = 0.0;
//	                   q_v[4] = 0.0;
//	                   q_v[5] = 0.0;
//
//
//	                  /* Dealing with SOUTH NODE Ks - look at nodes j-1 and j */
//                         //k1 = fabs((int)top_dat[itop-sy_v]);
//                         //if(k1 >= 0)
//                         //{
//                            ip3 = SubvectorEltIndex(p_sub,i, (j-1), k);
//
//			  /* compute the derivative wrt to i,j,k */
//
//			      slope_fy_lo = sy_dat[io] - (((pfmax((pp[ip0]),0.0)) -(pfmax((pp[ip3]),0.0)))/dy);
//			      if(slope_fy_lo > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_lo < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//			      q_v[3] = ydir * (1.0 / mann_dat[io]) * (((5.0/3.0)*RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) *(RPowerR(fabs(slope_fy_lo),0.5)))- ((1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fy_lo),-0.5))));
//
//			     /* compute the derivative wrt to i, j-1 ,k */
//
//			      slope_fy_mid = sy_dat[io-sy_v] - (((pfmax((pp[ip0]),0.0)) -(pfmax((pp[ip3]),0.0)))/dy);
//			      if(slope_fy_mid > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_mid < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//			      q_v[4] = ydir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fy_mid),-0.5));
//
//			   //}
//
//	                   /* compute SOUTH NODE KS - NOTE: io is for current cell */
//			   ks_v[io] = pfmax(q_v[4],0.0)- pfmax(-q_v[3],0.0);
//
//                        /* Dealing with Ke - look at nodes i and i+1 */
//                         //k1 = fabs((int)top_dat[itop+sy_v]);
//                         //if(k1 >= 0)
//                         //{
//                            ip4 = SubvectorEltIndex(p_sub,i, (j+1), k);
//
//			  /* compute the derivative wrt to i j k */
//
//			      slope_fy_hi = sy_dat[io] - (((pfmax((pp[ip4]),0.0)) - (pfmax((pp[ip0]),0.0)))/dy);
//			      if(slope_fy_hi > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_hi < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//			      q_v[5] =ydir * (1.0 / mann_dat[io]) * (((5.0/3.0)*RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) *(RPowerR(fabs(slope_fy_hi),0.5)))- ((1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fy_hi),-0.5))));
//
//			    /* compute the derivative wrt to i j+1 k */
//
//			      slope_fy_mid = sy_dat[io+sy_v] - (((pfmax((pp[ip4]),0.0)) - (pfmax((pp[ip0]),0.0)))/dy);
//			      if(slope_fy_mid > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_mid < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//                          q_v[4] = ydir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fy_mid),-0.5));
//
//			   //}
//	                   /* compute kn - NOTE: io is for current cell */
//	                   printf(" %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f   %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f \n ",(pfmax(pp[ip0],0.0)), (pfmax(pp[ip],0.0)), (pfmax(pp[ip2],0.0)), (pfmax(pp[ip3],0.0)),(pfmax(pp[ip4],0.0)), slope_fy_lo, slope_fy_mid, slope_fy_hi, slope_fx_lo, slope_fx_mid, slope_fx_hi, q_v[0], q_v[1], q_v[2], q_v[3], q_v[4], q_v[5], kn_v[io], ke_v[io], ks_v[io],kw_v[io] );
//			   kn_v[io] = pfmax(q_v[5],0.0) - pfmax(-q_v[4],0.0) ;
//	         }
//
//	     });
//
//	 }
//	 else
//	 {
//              if(qx_v == NULL || qy_v == NULL) /* do not return velocity fluxes */
//         {
//            BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, sg,
//	    {
//	       if (fdir[2] == 1)
//	       {
//	                   /* First initialize velocities, q_v, for inactive region */
//	                   q_v[0] = 0.0;
//	                   q_v[1] = 0.0;
//	                   q_v[2] = 0.0;
//
//	                   io = SubvectorEltIndex(sx_sub,i, j, 0);
//	                   ip = SubvectorEltIndex(top_sub, i, j, 0);
//
//                          /* precompute some data for current node, to be reused */
//                         ip0 = SubvectorEltIndex(p_sub, i, j, k);
//                         coeff = RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0));
//
//
//		      //   ***********************************************************
//		      //
//	              /* compute east and west faces */
//	              //
//	              //   ***********************************************************
//
//
//                         /* Dealing with WEST NODE, Kw - look at nodes i-1 and i */
//                         //k1 = fabs((int)top_dat[itop-1]);
//                         //if(k1 >= 0)
//                         //{
//                            ip = SubvectorEltIndex(p_sub,(i-1), j, k);
//
//                            /* Compute the derivative wrt to node i,j,k  */
//
//                        slope_fx_lo = sx_dat[io] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip]),0.0)))/dx);
//			      if(slope_fx_lo > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_lo < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//                         q_v[0] = xdir * (1.0 / mann_dat[io]) * (((5.0/3.0) * RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) * (RPowerR(fabs(slope_fx_lo),0.5)))- ((1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) * (RPowerR(fabs(slope_fx_lo),-0.5))));
//
//
//			  /* compute the derivative wrt to node i-1,j,k */
//			      slope_fx_mid = sx_dat[io-1] - (((pfmax((pp[ip0]),0.0)) - (pfmax((pp[ip]),0.0)))/dx);
//			      if(slope_fx_mid > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_mid < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//			      q_v[1] = xdir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) * (RPowerR(fabs(slope_fx_mid),-0.5));
//
//			   //}
//	                   /* compute kw - NOTE: io is for current cell */
//			   kw_v[io] = pfmax(q_v[1],0.0) - pfmax(-q_v[0],0.0);
//			   qx_v[io] = q_v[0];					// IS THIS RIGHT????
//
//
//                        /* Dealing with EAST NODE, Ke - look at nodes i and i+1 */
//
//                         k1 = fabs((int)top_dat[itop+1]);
//                         if(k1 >= 0)
//                         {
//                            ip = SubvectorEltIndex(p_sub,(i+1), j, k);
//
//                            /* compute the derivative wrt i,j,k*/
//			      slope_fx_hi = sx_dat[io] - (((pfmax((pp[ip]),0.0)) - (pfmax((pp[ip0]),0.0)))/dx);
//			      if(slope_fx_hi > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_hi < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//			   q_v[2] = xdir * (1.0 / mann_dat[io]) * (((5.0/3.0) * RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) * (RPowerR(fabs(slope_fx_hi),0.5)))- ((1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fx_hi),-0.5))));
//
//			   /* Computer the derivative for the i+1,j,k */
//
//			   slope_fx_mid = sx_dat[io+1] - (((pfmax((pp[ip]),0.0)) - (pfmax((pp[ip0]),0.0)))/dx);
//			      if(slope_fx_mid > 0.0)
//		                    xdir = -1.0;
//				 else if(slope_fx_mid < 0.0)
//                                  xdir = 1.0;
//                               else
//                                  xdir = 0.0;
//
//			   q_v[1] = xdir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dx)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fx_mid),-0.5));
//
//			   }
//			   /* compute ke - NOTE: io is for current cell */
//			    ke_v[io] = pfmax(q_v[2],0.0) - pfmax(-q_v[1],0.0);
//
//			    if (pfmax(q_v[2],0.0))
//			            qx_v[io] = q_v[2];					//WHAT DO I PUT HERE??????????
//
//
//		      //   ***********************************************************
//		      //
//	              /* compute north and south faces */
//	              //
//	              //   ***********************************************************
//
//	               /* First initialize velocities, q_v, for inactive region */
//	                   q_v[3] = 0.0;
//	                   q_v[4] = 0.0;
//	                   q_v[5] = 0.0;
//
//
//	                  /* Dealing with SOUTH NODE Ks - look at nodes j-1 and j */
//                        // k1 = fabs((int)top_dat[itop-sy_v]);
//                        // if(k1 >= 0)
//                        // {
//                            ip = SubvectorEltIndex(p_sub,i, (j-1), k);
//
//			  /* compute the derivative wrt to i,j,k */
//
//			      slope_fy_lo = sy_dat[io] - (((pfmax((pp[ip0]),0.0)) -(pfmax((pp[ip]),0.0)))/dy);
//			      if(slope_fy_lo > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_lo < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//			      q_v[3] = ydir * (1.0 / mann_dat[io]) * (((5.0/3.0) * RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) * (RPowerR(fabs(slope_fy_lo),0.5))) - ((1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) * (RPowerR(fabs(slope_fy_lo),-0.5))));
//
//			     /* compute the derivative wrt to i, j-1 ,k */
//
//			      slope_fy_mid = sy_dat[io-sy_v] - (((pfmax((pp[ip0]),0.0)) -(pfmax((pp[ip]),0.0)))/dy);
//			      if(slope_fy_mid > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_mid < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//			      q_v[4] = ydir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fy_mid),-0.5));
//
//			  // }
//
//	                   /* compute SOUTH NODE KS - NOTE: io is for current cell */
//			   ks_v[io] = pfmax(q_v[4],0.0) - pfmax(-q_v[3],0.0);
//			   qy_v[io] = q_v[3];							//WHAT DO I PUT HERE??????????
//
//                        /* Dealing with Ke - look at nodes i and i+1 */
//                         //k1 = fabs((int)top_dat[itop+sy_v]);
//                         //if(k1 >= 0)
//                         //{
//                            ip = SubvectorEltIndex(p_sub,i, (j+1), k);
//
//			  /* compute the derivative wrt to i j k */
//
//			      slope_fy_hi = sy_dat[io] - (((pfmax((pp[ip]),0.0)) - (pfmax((pp[ip0]),0.0)))/dy);
//			      if(slope_fy_hi > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_hi < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//			      q_v[5] =ydir * (1.0 / mann_dat[io]) * (((5.0/3.0)*RPowerR(pfmax((pp[ip0]),0.0),(2.0/3.0)) * (RPowerR(fabs(slope_fy_hi),0.5)))- ((1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) * (RPowerR(fabs(slope_fy_hi),-0.5))));
//
//			    /* compute the derivative wrt to i j+1 k */
//
//			      slope_fy_mid = sy_dat[io+sy_v] - (((pfmax((pp[ip]),0.0)) - (pfmax((pp[ip0]),0.0)))/dy);
//			      if(slope_fy_mid > 0.0)
//		                    ydir = -1.0;
//				 else if(slope_fy_mid < 0.0)
//                                  ydir = 1.0;
//                               else
//                                  ydir = 0.0;
//
//                          q_v[4] = ydir * (1.0 / mann_dat[io]) * (1.0 / (2.0 * dy)) * RPowerR(pfmax((pp[ip0]),0.0),(5.0/3.0)) *(RPowerR(fabs(slope_fy_mid),-0.5));
//
//			   //}
//
//	                   /* compute kn - NOTE: io is for current cell */
//			   kn_v[io] = pfmax(q_v[5],0.0) - pfmax(-q_v[4],0.0);
//
//			   if (pfmax(q_v[4],0.0))
//			           qy_v[io] = q_v[4];
//	         }
//
//	     });
//
//	 }
//      }
}

//*/
/*--------------------------------------------------------------------------
 * OverlandFlowEvalKinInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *OverlandFlowEvalKinInitInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;

  instance_xtra = NULL;

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * OverlandFlowEvalKinFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  OverlandFlowEvalKinFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}

/*--------------------------------------------------------------------------
 * OverlandFlowEvalKinNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *OverlandFlowEvalKinNewPublicXtra()
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;

  public_xtra = NULL;

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}

/*-------------------------------------------------------------------------
 * OverlandFlowEvalKinFreePublicXtra
 *-------------------------------------------------------------------------*/

void  OverlandFlowEvalKinFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  if (public_xtra)
  {
    tfree(public_xtra);
  }
}

/*--------------------------------------------------------------------------
 * OverlandFlowEvalKinSizeOfTempData
 *--------------------------------------------------------------------------*/

int  OverlandFlowEvalKinSizeOfTempData()
{
  return 0;
}
