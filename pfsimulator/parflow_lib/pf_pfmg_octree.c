/*BHEADER**********************************************************************

  Copyright (c) 1995-2009, Lawrence Livermore National Security,
  LLC. Produced at the Lawrence Livermore National Laboratory. Written
  by the Parflow Team (see the CONTRIBUTORS file)
  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.

  This file is part of Parflow. For details, see
  http://www.llnl.gov/casc/parflow

  Please read the COPYRIGHT file or Our Notice and the LICENSE file
  for the GNU Lesser General Public License.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License (as published
  by the Free Software Foundation) version 2.1 dated February 1999.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
  and conditions of the GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
**********************************************************************EHEADER*/

#include "parflow.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

#ifdef HAVE_HYPRE
#include "hypre_dependences.h"

typedef struct
{
   int  max_iter;
   int  num_pre_relax;
   int  num_post_relax;
   int  smoother;
   
   int  box_size_power;

   int  time_index_pfmg;
   int  time_index_copy_hypre;

} PublicXtra;

typedef struct
{
   ProblemData         *problem_data;
   Grid                *grid;
   
   double               dxyz[3];

   HYPRE_StructGrid     hypre_grid;
   HYPRE_StructMatrix   hypre_mat;
   HYPRE_StructVector   hypre_b, hypre_x;
   HYPRE_StructStencil  hypre_stencil;
   
   HYPRE_StructSolver   hypre_pfmg_data;

} InstanceXtra;

/*--------------------------------------------------------------------------
 * PFMGOctree
 *--------------------------------------------------------------------------*/

void         PFMGOctree(
Vector      *soln,
Vector      *rhs,
double       tol,
int          zero)
{
   PFModule           *this_module    = ThisPFModule;
   InstanceXtra       *instance_xtra  = (InstanceXtra *)PFModuleInstanceXtra(this_module);
   PublicXtra         *public_xtra    = (PublicXtra *)PFModulePublicXtra(this_module);

   HYPRE_StructMatrix  hypre_mat      = instance_xtra -> hypre_mat;
   HYPRE_StructVector  hypre_b        = instance_xtra -> hypre_b;
   HYPRE_StructVector  hypre_x        = instance_xtra -> hypre_x;
   
   HYPRE_StructSolver  hypre_pfmg_data = instance_xtra -> hypre_pfmg_data;

   Grid               *grid           = VectorGrid(rhs);
   Subgrid            *subgrid;
   int                 sg;

   Subvector          *rhs_sub;
   Subvector          *soln_sub;

   double             *rhs_ptr;
   double             *soln_ptr;

   int                 ix,   iy,   iz;
   int                 nx,   ny,   nz;
   int                 nx_v, ny_v, nz_v;
   int                 i, j, k;
   int                 num_i, num_j, num_k;
   int                 iv;

   int                 num_iterations;
   double              rel_norm;

   int                 box_size_power = public_xtra -> box_size_power;

   GrGeomSolid        *gr_domain = ProblemDataGrDomain(instance_xtra -> problem_data);

   /* Copy rhs to hypre_b vector. */
   BeginTiming(public_xtra->time_index_copy_hypre);

   ForSubgridI(sg, GridSubgrids(grid))
   {
      int outside = 0;
      int boxnum  = -1;
      int action  = 0; // set values

      hypre_Box          *set_box;
      hypre_Box          *value_box;
      int                 ilo[3];
      int                 ihi[3];

      subgrid = SubgridArraySubgrid(GridSubgrids(grid), sg);
      rhs_sub = VectorSubvector(rhs, sg);

      rhs_ptr = SubvectorData(rhs_sub);

      ix = SubgridIX(subgrid);
      iy = SubgridIY(subgrid);
      iz = SubgridIZ(subgrid);

      nx = SubgridNX(subgrid);
      ny = SubgridNY(subgrid);
      nz = SubgridNZ(subgrid);

      nx_v = SubvectorNX(rhs_sub);
      ny_v = SubvectorNY(rhs_sub);
      nz_v = SubvectorNZ(rhs_sub);


      ilo[0] = SubvectorIX(rhs_sub);
      ilo[1] = SubvectorIY(rhs_sub);
      ilo[2] = SubvectorIZ(rhs_sub);
      ihi[0] = ilo[0] + nx_v - 1;
      ihi[1] = ilo[1] + ny_v - 1;
      ihi[2] = ilo[2] + nz_v - 1;

      value_box = hypre_BoxCreate();
      hypre_BoxSetExtents(value_box, ilo, ihi); 

      GrGeomInBoxLoop(i, j, k, 
		      num_i, num_j, num_k,
		      gr_domain, box_size_power,
		      ix, iy, iz, nx, ny, nz, 
		      {
			 ilo[0] = i;
			 ilo[1] = j;
			 ilo[2] = k;
			 ihi[0] = ilo[0] + num_i - 1;
			 ihi[1] = ilo[1] + num_j - 1;
			 ihi[2] = ilo[2] + num_k - 1;
			 
			 set_box = hypre_BoxCreate();
			 hypre_BoxSetExtents(set_box, ilo, ihi); 

			 hypre_StructVectorSetBoxValues ( hypre_b,
							  set_box, 
							  value_box, 
							  rhs_ptr, 
							  action, 
							  boxnum,
							  outside );
			 hypre_BoxDestroy(set_box);
			 
		      });

      hypre_BoxDestroy(value_box);
   }
   HYPRE_StructVectorAssemble(hypre_b);

   EndTiming(public_xtra->time_index_copy_hypre);

   if (tol > 0.0)
   {
      IfLogging(1)
      {
	 HYPRE_StructPFMGSetLogging(instance_xtra->hypre_pfmg_data, 1);
      }
   }

   /* Invoke the preconditioner using a zero initial guess */
   HYPRE_StructPFMGSetZeroGuess(hypre_pfmg_data);

   BeginTiming(public_xtra->time_index_pfmg);

   HYPRE_StructPFMGSolve(hypre_pfmg_data, hypre_mat, hypre_b, hypre_x);

   EndTiming(public_xtra->time_index_pfmg);

   if (tol > 0.0)
   {
      IfLogging(1)
      {
         FILE  *log_file;
	 
	 HYPRE_StructPFMGGetNumIterations(hypre_pfmg_data, &num_iterations);
	 HYPRE_StructPFMGGetFinalRelativeResidualNorm(hypre_pfmg_data, 
						     &rel_norm);

	 log_file = OpenLogFile("PFMG");
	 fprintf(log_file, "PFMGOctree num. its: %i  PFMGOctree Final norm: %12.4e\n", 
		 num_iterations, rel_norm);
	 CloseLogFile(log_file);
      }
   }

   /* Copy solution from hypre_x vector to the soln vector. */
   BeginTiming(public_xtra->time_index_copy_hypre);

   ForSubgridI(sg, GridSubgrids(grid))
   {
      int outside = 0;
      int boxnum  = -1;
      int action  = -1; // get values

      hypre_Box          *set_box;
      hypre_Box          *value_box;
      int                 ilo[3];
      int                 ihi[3];

      subgrid = SubgridArraySubgrid(GridSubgrids(grid), sg);
      soln_sub = VectorSubvector(soln, sg);

      soln_ptr = SubvectorData(soln_sub);

      ix = SubgridIX(subgrid);
      iy = SubgridIY(subgrid);
      iz = SubgridIZ(subgrid);

      nx = SubgridNX(subgrid);
      ny = SubgridNY(subgrid);
      nz = SubgridNZ(subgrid);

      nx_v = SubvectorNX(soln_sub);
      ny_v = SubvectorNY(soln_sub);
      nz_v = SubvectorNZ(soln_sub);

      iv  = SubvectorEltIndex(soln_sub, ix, iy, iz);

      ilo[0] = SubvectorIX(soln_sub);
      ilo[1] = SubvectorIY(soln_sub);
      ilo[2] = SubvectorIZ(soln_sub);
      ihi[0] = ilo[0] + nx_v - 1;
      ihi[1] = ilo[1] + ny_v - 1;
      ihi[2] = ilo[2] + nz_v - 1;

      value_box = hypre_BoxCreate();
      hypre_BoxSetExtents(value_box, ilo, ihi); 

      GrGeomInBoxLoop(i, j, k, 
		      num_i, num_j, num_k,
		      gr_domain, box_size_power,
		      ix, iy, iz, nx, ny, nz, 
		      {
			 ilo[0] = i;
			 ilo[1] = j;
			 ilo[2] = k;
			 ihi[0] = ilo[0] + num_i - 1;
			 ihi[1] = ilo[1] + num_j - 1;
			 ihi[2] = ilo[2] + num_k - 1;
			 
			 set_box = hypre_BoxCreate();
			 hypre_BoxSetExtents(set_box, ilo, ihi); 

			 hypre_StructVectorSetBoxValues ( hypre_x,
							  set_box, 
							  value_box, 
							  soln_ptr, 
							  action, 
							  boxnum,
							  outside );
			 hypre_BoxDestroy(set_box);
			 
		      });

      hypre_BoxDestroy(value_box);
   }
   EndTiming(public_xtra->time_index_copy_hypre);

}

/*--------------------------------------------------------------------------
 * PFMGOctreeInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule  *PFMGOctreeInitInstanceXtra(
Problem      *problem,
Grid         *grid,
ProblemData  *problem_data,
Matrix       *pf_matrix,
double       *temp_data)
{
   PFModule      *this_module        = ThisPFModule;
   PublicXtra    *public_xtra        = (PublicXtra *)PFModulePublicXtra(this_module);
   InstanceXtra  *instance_xtra;

   int                 max_iter      = public_xtra -> max_iter;
   int                 num_pre_relax = public_xtra -> num_pre_relax;
   int                 num_post_relax= public_xtra -> num_post_relax;
   int                 smoother      = public_xtra -> smoother;

   Grid               *mat_grid;
   Subgrid            *subgrid;
   int                 sg;

   Submatrix          *pf_sub;
   double             *cp, *wp, *ep, *sop, *np, *lp, *up;

   int                 i, j, k;
   int                 num_i, num_j, num_k;
   int                 ix, iy, iz;
   int                 nx, ny, nz;
   int                 nx_m, ny_m, nz_m;
   int                 im;
   int                 stencil_size;
   int                 symmetric;

   int                 full_ghosts[6]          = {1, 1, 1, 1, 1, 1};
   int                 no_ghosts[6]            = {0, 0, 0, 0, 0, 0};
   int                 stencil_indices[7]      = {0, 1, 2, 3, 4, 5, 6};
   int                 stencil_indices_symm[4] = {0, 1, 2, 3};
   int                 ilo[3];
   int                 ihi[3];
   
   int                 r;

   int                 box_size_power = public_xtra -> box_size_power;

   GrGeomSolid        *gr_domain = ProblemDataGrDomain(problem_data);

   if ( PFModuleInstanceXtra(this_module) == NULL )
      instance_xtra = ctalloc(InstanceXtra, 1);
   else
      instance_xtra = (InstanceXtra *)PFModuleInstanceXtra(this_module);

   instance_xtra -> problem_data = problem_data;

   if ( grid != NULL )
   {
      instance_xtra -> problem_data = problem_data;
      instance_xtra -> grid         = grid;
   }

   if(gr_domain != NULL) 
   {
      /* Free the HYPRE grid */
      if (instance_xtra -> hypre_grid) 
      {
	 HYPRE_StructGridDestroy(instance_xtra->hypre_grid);
	 instance_xtra->hypre_grid = NULL;
      }

      /* Set the HYPRE grid */
      HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &(instance_xtra->hypre_grid) );


      grid         = instance_xtra -> grid;
      
	 /* Set local grid extents as global grid values */
      ForSubgridI(sg, GridSubgrids(grid))
      {
	 subgrid = GridSubgrid(grid, sg);
	 	 
	 ix = SubgridIX(subgrid);
	 iy = SubgridIY(subgrid);
	 iz = SubgridIZ(subgrid);

	 nx = SubgridNX(subgrid);
	 ny = SubgridNY(subgrid);
	 nz = SubgridNZ(subgrid);

	 r = SubgridRX(subgrid);

	 instance_xtra->dxyz[0] = SubgridDX(subgrid);
	 instance_xtra->dxyz[1] = SubgridDY(subgrid);
	 instance_xtra->dxyz[2] = SubgridDZ(subgrid);

	 GrGeomInBoxLoop(i, j, k, 
			 num_i, num_j, num_k,
			 gr_domain, box_size_power,
			 ix, iy, iz, nx, ny, nz, 
			 {
			    ilo[0] = i;
			    ilo[1] = j;
			    ilo[2] = k;
			    ihi[0] = ilo[0] + num_i - 1;
			    ihi[1] = ilo[1] + num_j - 1;
			    ihi[2] = ilo[2] + num_k - 1;

			    HYPRE_StructGridSetExtents(instance_xtra->hypre_grid, ilo, ihi); 
			 });
      }

      HYPRE_StructGridAssemble(instance_xtra->hypre_grid);
   }

   /* Reset the HYPRE solver for each recompute of the PC matrix.  
      This reset will require a matrix copy from PF format to HYPRE format. */
   if ( pf_matrix != NULL )
   {
      /* Free old solver data because HYPRE requires a new solver if 
         matrix values change */
      if (instance_xtra->hypre_pfmg_data) {
	 HYPRE_StructPFMGDestroy(instance_xtra->hypre_pfmg_data);
	 instance_xtra->hypre_pfmg_data = NULL;
      }
      
      /* For remainder of routine, assume matrix is structured the same for
	 entire nonlinear solve process */
      /* Set stencil parameters */
      stencil_size = MatrixDataStencilSize(pf_matrix);
      if ( !(instance_xtra->hypre_stencil) )
      {
         HYPRE_StructStencilCreate(3, stencil_size, 
				&(instance_xtra->hypre_stencil) );

         for (i = 0; i < stencil_size; i++) 
         {
            HYPRE_StructStencilSetElement(instance_xtra->hypre_stencil, i,
                                         &(MatrixDataStencil(pf_matrix))[i*3]);
         }
      }

      /* Set up new matrix */
      symmetric = MatrixSymmetric(pf_matrix);
      if ( !(instance_xtra->hypre_mat) )
      {
         HYPRE_StructMatrixCreate(MPI_COMM_WORLD, instance_xtra->hypre_grid, 
			       instance_xtra->hypre_stencil,
			       &(instance_xtra->hypre_mat) );
	 HYPRE_StructMatrixSetNumGhost(instance_xtra->hypre_mat, full_ghosts);
         HYPRE_StructMatrixSetSymmetric(instance_xtra->hypre_mat, symmetric);
         HYPRE_StructMatrixInitialize(instance_xtra->hypre_mat);
      }

      /* Set up new right-hand-side vector */
      if ( !(instance_xtra->hypre_b) )
      {
         HYPRE_StructVectorCreate(MPI_COMM_WORLD, 
			       instance_xtra->hypre_grid, 
			       &(instance_xtra->hypre_b) );
	 HYPRE_StructVectorSetNumGhost(instance_xtra->hypre_b, no_ghosts);
	 HYPRE_StructVectorInitialize(instance_xtra->hypre_b);
      }

      /* Set up new solution vector */
      if ( !(instance_xtra->hypre_x) )
      {
         HYPRE_StructVectorCreate(MPI_COMM_WORLD, 
			       instance_xtra->hypre_grid, 
			       &(instance_xtra->hypre_x) );
	 HYPRE_StructVectorSetNumGhost(instance_xtra->hypre_x, full_ghosts);
	 HYPRE_StructVectorInitialize(instance_xtra->hypre_x);
      }
      HYPRE_StructVectorSetConstantValues(instance_xtra->hypre_x, 0.0e0);
      HYPRE_StructVectorAssemble(instance_xtra->hypre_x);

      /* Copy the matrix entries */
      BeginTiming(public_xtra->time_index_copy_hypre);

      mat_grid = MatrixGrid(pf_matrix);
      ForSubgridI(sg, GridSubgrids(mat_grid))
      {
	 subgrid = GridSubgrid(mat_grid, sg);

	 pf_sub  = MatrixSubmatrix(pf_matrix, sg);

	 if (symmetric)
	 {
	    /* Pull off upper diagonal coeffs here for symmetric part */
	    cp      = SubmatrixStencilData(pf_sub, 0);
	    ep      = SubmatrixStencilData(pf_sub, 2);
	    np      = SubmatrixStencilData(pf_sub, 4);
	    up      = SubmatrixStencilData(pf_sub, 6);
	 }
	 else
	 {
	    cp      = SubmatrixStencilData(pf_sub, 0);
	    wp      = SubmatrixStencilData(pf_sub, 1);
	    ep      = SubmatrixStencilData(pf_sub, 2);
	    sop     = SubmatrixStencilData(pf_sub, 3);
	    np      = SubmatrixStencilData(pf_sub, 4);
	    lp      = SubmatrixStencilData(pf_sub, 5);
	    up      = SubmatrixStencilData(pf_sub, 6);
	 }

	 ix = SubgridIX(subgrid);
	 iy = SubgridIY(subgrid);
	 iz = SubgridIZ(subgrid);
	 
	 nx = SubgridNX(subgrid);
	 ny = SubgridNY(subgrid);
	 nz = SubgridNZ(subgrid);
	 
	 nx_m  = SubmatrixNX(pf_sub);
	 ny_m  = SubmatrixNY(pf_sub);
	 nz_m  = SubmatrixNZ(pf_sub);

	 im  = SubmatrixEltIndex(pf_sub,  ix, iy, iz);

	 if (symmetric)
	 {
	    int outside = 0;
	    int boxnum  = -1;
	    int action  = 0; // set values
	    int stencil;

	    hypre_Box          *set_box;
	    hypre_Box          *value_box;

	    ilo[0] = SubmatrixIX(pf_sub);
	    ilo[1] = SubmatrixIY(pf_sub);
	    ilo[2] = SubmatrixIZ(pf_sub);
	    ihi[0] = ilo[0] + nx_m - 1;
	    ihi[1] = ilo[1] + ny_m - 1;
	    ihi[2] = ilo[2] + nz_m - 1;

	    value_box = hypre_BoxCreate();
	    hypre_BoxSetExtents(value_box, ilo, ihi); 

	    GrGeomInBoxLoop(i, j, k, 
			    num_i, num_j, num_k,
			    gr_domain, box_size_power,
			    ix, iy, iz, nx, ny, nz, 
			    {
			       
			       ilo[0] = i;
			       ilo[1] = j;
			       ilo[2] = k;
			       ihi[0] = ilo[0] + num_i - 1;
			       ihi[1] = ilo[1] + num_j - 1;
			       ihi[2] = ilo[2] + num_k - 1;
			       
			       set_box = hypre_BoxCreate();
			       hypre_BoxSetExtents(set_box, ilo, ihi); 
			       
                               /* IMF: commented print statement
			       amps_Printf("hypre symm matrix value box : %d (%d, %d, %d) (%d, %d, %d)\n", PV_l,
					   ilo[0], ilo[1], ilo[2], ihi[0], ihi[1], ihi[2]);
			       */
                               
			       /*
				 Note that loop over stencil's is necessary due to hypre
				 interface wanting stencil values to be contiguous.
				 FORTRAN ordering of (stencil, i, j, k).  PF stores as
				 (i, j, k, stencil)
			       */
			       for(stencil = 0; stencil < stencil_size; ++stencil) {
				  
				  /* 
				     symmetric stencil values are at 0, 2, 4, 6
				  */
				  double *values = SubmatrixStencilData(pf_sub, stencil*2);
				  
				  hypre_StructMatrixSetBoxValues( instance_xtra->hypre_mat,
								  set_box,
								  value_box,
								  1, 
								  &stencil_indices_symm[stencil], 
								  values,
								  action,
								  boxnum,
								  outside );
			       }
			       
			       hypre_BoxDestroy(set_box);
			    });
	    
	    hypre_BoxDestroy(value_box);

	 }
	 else
	 {
	    int outside = 0;
	    int boxnum  = -1;
	    int action  = 0; // set values
	    int stencil;

	    hypre_Box          *set_box;
	    hypre_Box          *value_box;

	    ilo[0] = ix;
	    ilo[1] = iy;
	    ilo[2] = iz;
	    ihi[0] = ilo[0] + nx - 1;
	    ihi[1] = ilo[1] + ny - 1;
	    ihi[2] = ilo[2] + nz - 1;

	    value_box = hypre_BoxCreate();
	    hypre_BoxSetExtents(value_box, ilo, ihi); 

	    GrGeomInBoxLoop(i, j, k, 
			    num_i, num_j, num_k,
			    gr_domain, box_size_power,
			    ix, iy, iz, nx, ny, nz, 
			    {
			       ilo[0] = i;
			       ilo[1] = j;
			       ilo[2] = k;
			       ihi[0] = ilo[0] + num_i - 1;
			       ihi[1] = ilo[1] + num_j - 1;
			       ihi[2] = ilo[2] + num_k - 1;

			       set_box = hypre_BoxCreate();
			       hypre_BoxSetExtents(set_box, ilo, ihi); 

			       /*
				 Note that loop over stencil's is necessary due to hypre
				 interface wanting stencil values to be contiguous.
				 FORTRAN ordering of (stencil, i, j, k).  PF stores as
				 (i, j, k, stencil)
			       */
			       for(stencil = 0; stencil < stencil_size; ++stencil) {
				  
				  double *values = SubmatrixStencilData(pf_sub, stencil);
				  
				  hypre_StructMatrixSetBoxValues( instance_xtra->hypre_mat,
								  set_box,
								  value_box,
								  1, 
								  &stencil_indices[stencil], 
								  values,
								  action,
								  boxnum,
								  outside );
			       }

			       hypre_BoxDestroy(set_box);
			    });

	    hypre_BoxDestroy(value_box);
	 }
      }   /* End subgrid loop */
      HYPRE_StructMatrixAssemble(instance_xtra->hypre_mat);

      EndTiming(public_xtra->time_index_copy_hypre);

      /* Set up the PFMG preconditioner */
      HYPRE_StructPFMGCreate(MPI_COMM_WORLD,
				&(instance_xtra->hypre_pfmg_data) );

      HYPRE_StructPFMGSetTol(instance_xtra->hypre_pfmg_data, 1.0e-30);
      /* Set user parameters for PFMG */
      HYPRE_StructPFMGSetMaxIter(instance_xtra->hypre_pfmg_data, max_iter);
      HYPRE_StructPFMGSetNumPreRelax(instance_xtra->hypre_pfmg_data, 
				     num_pre_relax);
      HYPRE_StructPFMGSetNumPostRelax(instance_xtra->hypre_pfmg_data, 
				      num_post_relax);
      /* Jacobi = 0; weighted Jacobi = 1; red-black GS symmetric = 2; red-black GS non-symmetric = 3 */
      HYPRE_StructPFMGSetRelaxType(instance_xtra->hypre_pfmg_data, smoother);

      /* Use non-Galkerkin option */
      HYPRE_StructPFMGSetRAPType(instance_xtra->hypre_pfmg_data, 1);

      HYPRE_StructPFMGSetSkipRelax(instance_xtra->hypre_pfmg_data, 1);

      HYPRE_StructPFMGSetDxyz(instance_xtra->hypre_pfmg_data, 	
			      instance_xtra->dxyz);

      HYPRE_StructPFMGSetup(instance_xtra->hypre_pfmg_data, 
			    instance_xtra->hypre_mat, 
			    instance_xtra->hypre_b, instance_xtra->hypre_x);
   }

   PFModuleInstanceXtra(this_module) = instance_xtra;
   return this_module;
}


/*--------------------------------------------------------------------------
 * PFMGOctreeFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  PFMGOctreeFreeInstanceXtra()
{
   PFModule      *this_module   = ThisPFModule;
   InstanceXtra  *instance_xtra = (InstanceXtra *)PFModuleInstanceXtra(this_module);

   if (instance_xtra)
   {
      if (instance_xtra->hypre_pfmg_data)
	 HYPRE_StructPFMGDestroy(instance_xtra->hypre_pfmg_data);
      if (instance_xtra->hypre_mat)
	 HYPRE_StructMatrixDestroy(instance_xtra->hypre_mat);
      if (instance_xtra->hypre_b)
	 HYPRE_StructVectorDestroy(instance_xtra->hypre_b);
      if (instance_xtra->hypre_x)
	 HYPRE_StructVectorDestroy(instance_xtra->hypre_x);
      if (instance_xtra->hypre_stencil)
	 HYPRE_StructStencilDestroy(instance_xtra->hypre_stencil);
      if (instance_xtra->hypre_grid)
	 HYPRE_StructGridDestroy(instance_xtra->hypre_grid);

      tfree(instance_xtra);
   }
}

/*--------------------------------------------------------------------------
 * PFMGOctreeNewPublicXtra
 *--------------------------------------------------------------------------*/

PFModule  *PFMGOctreeNewPublicXtra(char *name)
{
   PFModule      *this_module   = ThisPFModule;
   PublicXtra    *public_xtra;

   char           key[IDB_MAX_KEY_LEN];
   char          *smoother_name;

   NameArray      smoother_switch_na;

   int            smoother;

   public_xtra = ctalloc(PublicXtra, 1);

   sprintf(key, "%s.MaxIter", name);
   public_xtra -> max_iter = GetIntDefault(key, 1);

   sprintf(key, "%s.NumPreRelax", name);
   public_xtra -> num_pre_relax = GetIntDefault(key, 1);

   sprintf(key, "%s.NumPostRelax", name);
   public_xtra -> num_post_relax = GetIntDefault(key, 1);

   sprintf(key, "%s.BoxSizePowerOf2", name);
   public_xtra -> box_size_power = GetIntDefault(key, 4);

   /* Use a dummy place holder so that cardinalities match 
      with what HYPRE expects */
   smoother_switch_na = NA_NewNameArray("Dummy Jacobi WJacobi RBGaussSeidelSymmetric RBGaussSeidelNonSymmetric");
   sprintf(key, "%s.Smoother", name);
   smoother_name = GetStringDefault(key, "RBGaussSeidelNonSymmetric");
   smoother = NA_NameToIndex(smoother_switch_na, smoother_name);
   if (smoother != 0)
   {
      public_xtra->smoother = NA_NameToIndex(smoother_switch_na, 
					     smoother_name) - 1;
   }
   else
   {
      InputError("Error: Invalid value <%s> for key <%s>.\n", 
		 smoother_name, key);
   }
   NA_FreeNameArray(smoother_switch_na);

   public_xtra -> time_index_pfmg = RegisterTiming("PFMGOctree");
   public_xtra -> time_index_copy_hypre = RegisterTiming("HYPRE_Copies");

   PFModulePublicXtra(this_module) = public_xtra;

   return this_module;
}

/*-------------------------------------------------------------------------
 * PFMGOctreeFreePublicXtra
 *-------------------------------------------------------------------------*/

void  PFMGOctreeFreePublicXtra()
{
   PFModule    *this_module   = ThisPFModule;
   PublicXtra  *public_xtra   = (PublicXtra *)PFModulePublicXtra(this_module);

   if ( public_xtra )
   {
      tfree(public_xtra);
   }
}

/*--------------------------------------------------------------------------
 * PFMGOctreeSizeOfTempData
 *--------------------------------------------------------------------------*/

int  PFMGOctreeSizeOfTempData()
{
   return 0;
}

#else 

/*
 Hyper is not available.   
*/

void         PFMGOctree(soln, rhs, tol, zero)
Vector      *soln;
Vector      *rhs;
double       tol;
int          zero;
{
   amps_Printf("Error: Parflow not compiled with hypre, can't use pfmg\n");   
   exit(1);
}

PFModule  *PFMGOctreeInitInstanceXtra(problem, grid, problem_data,  
				      pf_matrix, temp_data)
Problem      *problem;
Grid         *grid;
ProblemData  *problem_data;
Matrix       *pf_matrix;
double       *temp_data;
{
   amps_Printf("Error: Parflow not compiled with hypre, can't use pfmg\n");   
   exit(1);
   return NULL;
}

void  PFMGOctreeFreeInstanceXtra()
{
   amps_Printf("Error: Parflow not compiled with hypre, can't use pfmg\n");   
   exit(1);
}

PFModule  *PFMGOctreeNewPublicXtra(char *name)
{
   amps_Printf("Error: Parflow not compiled with hypre, can't use pfmg\n");   
   exit(1);
   return NULL;
}

void  PFMGOctreeFreePublicXtra()
{
   amps_Printf("Error: Parflow not compiled with hypre, can't use pfmg\n");   
   exit(1);
}

int  PFMGOctreeSizeOfTempData()
{
   amps_Printf("Error: Parflow not compiled with hypre, can't use pfmg\n");   
   exit(1);
   return 0;
}

#endif

