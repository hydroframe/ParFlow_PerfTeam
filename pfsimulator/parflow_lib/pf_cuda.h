/*------------------------------------------------------------------------
  Header file for all CUDA specific functions for linking
  ----------------------------------------------------------------------*/

#ifndef _PF_CUDA_H
#define _PF_CUDA_H

/* cu_axpy.c */
void CU_Axpy(double, Vector*, Vector*);

/* cu_copy.c */
void CU_Copy(Vector*, Vector*);

/* cu_innerprod.c */
double CU_InnerProd(Vector *x, Vector *y);

/* cu_matrix.c */
Stencil *CU_NewStencil(int shape [][3], int sz);
CommPkg *CU_NewMatrixUpdatePkg(Matrix *matrix, Stencil *ghost);
CommHandle *CU_InitMatrixUpdate(Matrix *matrix);
void CU_FinalizeMatrixUpdate(CommHandle *handle);
Matrix *CU_NewMatrixType(Grid *grid, SubregionArray *range, Stencil *stencil, int symmetry, Stencil *ghost, enum matrix_type type);
void CU_FreeStencil(Stencil *stencil);
void CU_FreeMatrix(Matrix *matrix);
void CU_InitMatrix(Matrix *A, double value);

/* cu_matvec.c */
void CU_Matvec(double alpha, Matrix *A, Vector *x, double beta, Vector *y);
void CU_MatvecSubMat(void *  current_state,
                     double  alpha,
                     Matrix *JB,
                     Matrix *JC,
                     Vector *x,
                     double  beta,
                     Vector *y);
void CU_MatvecJacF(ProblemData *problem_data,
                   double       alpha,
                   Matrix *     JF,
                   Vector *     x,
                   double       beta,
                   Vector *     y);
void CU_MatvecJacE(ProblemData *problem_data,
                   double       alpha,
                   Matrix *     JE,
                   Vector *     x,
                   double       beta,
                   Vector *     y);


/* cu_vector.c */
CommPkg *CU_NewVectorCommPkg(Vector *vector, ComputePkg *compute_pkg);
VectorUpdateCommHandle  *CU_InitVectorUpdate(Vector *vector,
                                             int update_mode);
void CU_FinalizeVectorUpdate(VectorUpdateCommHandle *handle);
Vector  *CU_NewVector(Grid *grid,
                   int   nc,
                   int   num_ghost);
Vector  *CU_NewVectorType(Grid *           grid,
                          int              nc,
                          int              num_ghost,
                          enum vector_type type);
void CU_FreeVector(Vector *vector);
void CU_InitVector(Vector *v, double value);
void CU_InitVectorAll(Vector *v, double value);
void CU_InitVectorInc(Vector *v, double value, double inc);
void CU_InitVectorRandom(Vector *v, long seed);

/* cu_vector_utilities.c */
void CU_PFVLinearSum(double a, Vector *x, double b, Vector *y, Vector *z);
void CU_PFVConstInit(double c, Vector *z);
void CU_PFVProd(Vector *x, Vector *y, Vector *z);
void CU_PFVDiv(Vector *x, Vector *y, Vector *z);
void CU_PFVScale(double c, Vector *x, Vector *z);
void CU_PFVAbs(Vector *x, Vector *z);
void CU_PFVInv(Vector *x, Vector *z);
void CU_PFVAddConst(Vector *x, double b, Vector *z);
double CU_PFVDotProd(Vector *x, Vector *y);
double CU_PFVMaxNorm(Vector *x);
double CU_PFVWrmsNorm(Vector *x, Vector *w);
double CU_PFVWL2Norm(Vector *x, Vector *w);
double CU_PFVL1Norm(Vector *x);
double CU_PFVMin(Vector *x);
double CU_PFVMax(Vector *x);
int CU_PFVConstrProdPos(Vector *c, Vector *x);
void CU_PFVCompare(double c, Vector *x, Vector *z);
int CU_PFVInvTest(Vector *x, Vector *z);
void CU_PFVCopy(Vector *x, Vector *y);
void CU_PFVSum(Vector *x, Vector *y, Vector *z);
void CU_PFVDiff(Vector *x, Vector *y, Vector *z);
void CU_PFVNeg(Vector *x, Vector *z);
void CU_PFVScaleSum(double c, Vector *x, Vector *y, Vector *z);
void CU_PFVScaleDiff(double c, Vector *x, Vector *y, Vector *z);
void CU_PFVLin1(double a, Vector *x, Vector *y, Vector *z);
void CU_PFVLin2(double a, Vector *x, Vector *y, Vector *z);
void CU_PFVAxpy(double a, Vector *x, Vector *y);
void CU_PFVScaleBy(double a, Vector *x);
void CU_PFVLayerCopy(int a, int b, Vector *x, Vector *y);

/* cu_mg_semi.c */
void CU_MGSemi(Vector *x, Vector *b, double tol, int zero);
void CU_SetupCoarseOps(Matrix **A_l, Matrix **P_l, int num_levels, SubregionArray **f_sra_l, SubregionArray **c_sra_l);
PFModule *CU_MGSemiInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void CU_MGSemiFreeInstanceXtra(void);
PFModule *CU_MGSemiNewPublicXtra(char *name);
void CU_MGSemiFreePublicXtra(void);
int CU_MGSemiSizeOfTempData(void);

/* cu_mg_semi_prolong.c */
void CU_MGSemiProlong(Matrix *A_f, Vector *e_f, Vector *e_c, Matrix *P, SubregionArray *f_sr_array, SubregionArray *c_sr_array, ComputePkg *compute_pkg, CommPkg *e_f_comm_pkg);
ComputePkg *CU_NewMGSemiProlongComputePkg(Grid *grid, Stencil *stencil, int sx, int sy, int sz, int c_index, int f_index);

/* cu_mg_semi_restrict.c */
void CU_MGSemiRestrict(Matrix *A_f, Vector *r_f, Vector *r_c, Matrix *P, SubregionArray *f_sr_array, SubregionArray *c_sr_array, ComputePkg *compute_pkg, CommPkg *r_f_comm_pkg);
ComputePkg *CU_NewMGSemiRestrictComputePkg(Grid *grid, Stencil *stencil, int sx, int sy, int sz, int c_index, int f_index);

/* cu_rb_GS_point.c */
void CU_RedBlackGSPoint(Vector *x, Vector *b, double tol, int zero);
PFModule *CU_RedBlackGSPointInitInstanceXtra(Problem *problem, Grid *grid, ProblemData *problem_data, Matrix *A, double *temp_data);
void CU_RedBlackGSPointFreeInstanceXtra(void);
PFModule *CU_RedBlackGSPointNewPublicXtra(char *name);
void CU_RedBlackGSPointFreePublicXtra(void);
int CU_RedBlackGSPointSizeOfTempData(void);


/* cu_problem_saturation.c */
void CU_Saturation(Vector *phase_saturation, Vector *phase_pressure, Vector *phase_density, double gravity, ProblemData *problem_data, int fcn);
PFModule *CU_SaturationInitInstanceXtra(Grid *grid, double *temp_data);
void CU_SaturationFreeInstanceXtra(void);
PFModule *CU_SaturationNewPublicXtra(void);
void CU_SaturationFreePublicXtra(void);
int CU_SaturationSizeOfTempData(void);

/* cu_problem_phase_source.c */
void CU_PhaseSource(Vector *phase_source, int phase, Problem *problem, ProblemData *problem_data, double time);
PFModule *CU_PhaseSourceInitInstanceXtra(void);
void CU_PhaseSourceFreeInstanceXtra(void);
PFModule *CU_PhaseSourceNewPublicXtra(int num_phases);
void CU_PhaseSourceFreePublicXtra(void);
int CU_PhaseSourceSizeOfTempData(void);

/* cu_problem_phase_rel_perm.c */
void CU_PhaseRelPerm(Vector *phase_rel_perm, Vector *phase_pressure, Vector *phase_density, double gravity, ProblemData *problem_data, int fcn);
PFModule *CU_PhaseRelPermInitInstanceXtra(Grid *grid, double *temp_data);
void CU_PhaseRelPermFreeInstanceXtra(void);
PFModule *CU_PhaseRelPermNewPublicXtra(void);
void CU_PhaseRelPermFreePublicXtra(void);
int CU_PhaseRelPermSizeOfTempData(void);

/* cu_problem_phase_density.c */
void CU_PhaseDensity(int phase, Vector *phase_pressure, Vector *density_v, double *pressure_d, double *density_d, int fcn);
PFModule *CU_PhaseDensityInitInstanceXtra(void);
void CU_PhaseDensityFreeInstanceXtra(void);
PFModule *CU_PhaseDensityNewPublicXtra(int num_phases);
void CU_PhaseDensityFreePublicXtra(void);
int CU_PhaseDensitySizeOfTempData(void);


/* cu_vector.c */
void CU_FreeSubvector(Subvector *);
void CU_FreeTempVector(Vector *);
void CU_FreEVector(Vector *);


#endif // _PF_CUDA_H
