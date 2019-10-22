#ifndef _GET_MODULES_XTRA_H
#define _GET_MODULES_XTRA_H

#define GET_PUBLICXTRA_FUNC(module)                               \
  module ## PublicXtra *Get ## module ## Xtra()                   \
  {                                                               \
    PFModule *this_module = ThisPFModule;                         \
    return (module##PublicXtra*)PFModulePublicXtra(this_module);  \
  }

#define GET_INSTANCEXTRA_FUNC(module)                                 \
  module ## InstanceXtra *Get ## module ## InstanceXtra()             \
  {                                                                   \
    PFModule *this_module = ThisPFModule;                             \
    return (module##InstanceXtra*)PFModuleInstanceXtra(this_module);  \
  }

#define _GetInstanceXtra(module, pf_module)     \
  (ThisPFModule = pf_module, Get##module##InstanceXtra())

#define _GetPublicXtra(module, pf_module)       \
  (ThisPFModule = pf_module, Get ## module ## Xtra())

#define GetModulePublicXtra(module, pf_module, varname)                       \
  module ## PublicXtra* varname = _GetPublicXtra(module, pf_module);

#define GetModuleInstanceXtra(module, pf_module, varname) \
  module ## InstanceXtra* varname = _GetInstanceXtra(module, pf_module);

#define GetDummyType(module, type, xtra, varname)    \
  module ## Type ## type* varname = (module ## Type ## type*)(xtra);

// Phase Density Struct Mirrors
typedef struct {
  int num_phases;
  int *type;
  void **data;
} PhaseDensityPublicXtra;

typedef struct {
  double constant;
} PhaseDensityType0;

typedef struct {
  double reference_density;
  double compressibility_constant;
} PhaseDensityType1;

PhaseDensityPublicXtra *GetPhaseDensityXtra(void);

// Saturation Mirrors
typedef struct {
  int type;     /* input type */
  void  *data;  /* pointer to Type structure */

  NameArray regions;
} SaturationPublicXtra;

typedef struct {
  Grid    *grid;

  double  *temp_data;
} SaturationInstanceXtra;

typedef struct {
  int num_regions;
  int    *region_indices;
  double *values;
} SaturationType0;

typedef struct {
  int num_regions;
  int    *region_indices;
  int data_from_file;
  char   *alpha_file;
  char   *n_file;
  char   *s_sat_file;
  char   *s_res_file;
  double *alphas;
  double *ns;
  double *s_ress;
  double *s_difs;
  Vector *alpha_values;
  Vector *n_values;
  Vector *s_res_values;
  Vector *s_sat_values;
} SaturationType1;                      /* Van Genuchten Saturation Curve */

typedef struct {
  int num_regions;
  int    *region_indices;
  double *alphas;
  double *betas;
  double *s_ress;
  double *s_difs;
} SaturationType2;                      /* Haverkamp et.al. Saturation Curve */

typedef struct {
  int num_regions;
  int    *region_indices;
} SaturationType3;                      /* Data points for Saturation Curve */

typedef struct {
  int num_regions;
  int     *region_indices;
  int     *degrees;
  double **coefficients;
} SaturationType4;                      /* Polynomial function for Saturation Curve */

typedef struct {
  char    *filename;

  Vector  *satRF;
} SaturationType5;                      /* Spatially varying field over entire domain
                               * read from a file */
SaturationPublicXtra *GetSaturationXtra(void);

// BC_Pressure xtras
typedef struct {
  int num_phases;
  //int     iflag;   //@RMM
} BCPressurePublicXtra;

typedef struct {
  Problem  *problem;
  PFModule *phase_density;

  double     ***elevations;
  ProblemData  *problem_data;
  Grid         *grid;
} BCPressureInstanceXtra;

BCPressurePublicXtra *GetBCPressureXtra(void);
BCPressureInstanceXtra *GetBCPressureInstanceXtra(void);

// PhaseRelPerm Structs
typedef struct {
  NameArray regions;

  int type;     /* input type */
  void  *data;  /* pointer to Type structure */

  int time_index;
} PhaseRelPermPublicXtra;

typedef struct {
  Grid   *grid;
  double *temp_data;
} PhaseRelPermInstanceXtra;

typedef struct {
  int num_regions;
  int    *region_indices;
  double *values;
} PhaseRelPermType0;


typedef struct {
  double min_pressure_head;
  int num_sample_points;

  double *x;
  double *a;
  double *d;
  double *a_der;
  double *d_der;

  /* used by linear interpolation method */
  double *slope;
  double *slope_der;

  int interpolation_method;

  double interval;
} PhaseRelPermVanGTable;

typedef struct {
  int num_regions;
  int    *region_indices;
  int data_from_file;
  double *alphas;
  double *ns;
  char   *alpha_file;
  char   *n_file;
  Vector *alpha_values;
  Vector *n_values;

  PhaseRelPermVanGTable **lookup_tables;

#ifdef PF_PRINT_VG_TABLE
  int     *print_table;
#endif
} PhaseRelPermType1;                      /* Van Genuchten Rel. Perm. */

typedef struct {
  int num_regions;
  int    *region_indices;
  double *As;
  double *gammas;
} PhaseRelPermType2;                      /* Haverkamp, et.al. Rel. Perm. */

typedef struct {
  int num_regions;
  int    *region_indices;
} PhaseRelPermType3;                      /* Data points for Rel. Perm. */

typedef struct {
  int num_regions;
  int    *region_indices;
  int    *degrees;
  double **coefficients;
} PhaseRelPermType4;                      /* Polynomial Function for Rel. Perm. */

PhaseRelPermPublicXtra *GetPhaseRelPermXtra(void);
PhaseRelPermInstanceXtra *GetPhaseRelPermInstanceXtra(void);

// Richards BC Internal Structs

typedef struct {
  NameArray internal_bc_names;
  int num_conditions;

  int     *type;
  void   **data;
} RichardsBCInternalPublicXtra;

typedef struct {
  double xlocation;
  double ylocation;
  double zlocation;
  double value;
} RichardsBCInternalType0;                      /* basic point condition */

RichardsBCInternalPublicXtra *GetRichardsBCInternalXtra(void);

#endif // _GET_MODULES_XTRA_H
