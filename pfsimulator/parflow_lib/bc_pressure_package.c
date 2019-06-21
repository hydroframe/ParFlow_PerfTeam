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

//#include "parflow.h"

#include <string.h>

#include "bc_pressure_package.h"

#define BC_TYPE(type, values) typedef struct values Type ## type;
BC_TYPE_TABLE
#undef BC_TYPE


/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  int num_phases;

  /* cycling info */
  int num_cycles;

  int    *interval_divisions;
  int   **intervals;
  int    *repeat_counts;

  /* patch info */
  NameArray patches;
  int num_patches;

  int    *input_types;    /* num_patches input types */
  int    *patch_indexes;  /* num_patches patch indexes */
  int    *cycle_numbers;  /* num_patches cycle numbers */
  void  **data;           /* num_patches pointers to Type structures */
} PublicXtra;

typedef struct {
  Problem *problem;
} InstanceXtra;

/*--------------------------------------------------------------------------
 * BCPressurePackage
 *--------------------------------------------------------------------------*/

void         BCPressurePackage(
                               ProblemData *problem_data)
{
  PFModule         *this_module = ThisPFModule;
  PublicXtra       *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  BCPressureData   *bc_pressure_data
    = ProblemDataBCPressureData(problem_data);

  TimeCycleData    *time_cycle_data;

  int num_patches;
  int i;
  int cycle_length, cycle_number, interval_division,
    interval_number;

  /* Allocate the bc data */
  BCPressureDataNumPhases(bc_pressure_data) = (public_xtra->num_phases);

  BCPressureDataNumPatches(bc_pressure_data) = (public_xtra->num_patches);

  if ((public_xtra->num_patches) > 0)
  {
    /* Load the time cycle data */
    time_cycle_data = NewTimeCycleData((public_xtra->num_cycles),
                                       (public_xtra->interval_divisions));

    for (cycle_number = 0; cycle_number < (public_xtra->num_cycles);
         cycle_number++)
    {
      TimeCycleDataIntervalDivision(time_cycle_data, cycle_number)
        = (public_xtra->interval_divisions[cycle_number]);
      cycle_length = 0;
      for (interval_number = 0;
           interval_number < (public_xtra->
                              interval_divisions[cycle_number]);
           interval_number++)
      {
        cycle_length += (public_xtra->intervals[cycle_number])[interval_number];
        TimeCycleDataInterval(time_cycle_data, cycle_number, interval_number) = (public_xtra->intervals[cycle_number])[interval_number];
      }
      TimeCycleDataRepeatCount(time_cycle_data, cycle_number) = (public_xtra->repeat_counts[cycle_number]);
      TimeCycleDataCycleLength(time_cycle_data, cycle_number) = cycle_length;
    }

    BCPressureDataTimeCycleData(bc_pressure_data) = time_cycle_data;

    /* Load the Boundary Condition Data */

    num_patches = BCPressureDataNumPatches(bc_pressure_data);
    BCPressureDataTypes(bc_pressure_data) = ctalloc(int, num_patches);
    BCPressureDataPatchIndexes(bc_pressure_data) = ctalloc(int, num_patches);
    BCPressureDataCycleNumbers(bc_pressure_data) = ctalloc(int, num_patches);
    BCPressureDataBCTypes(bc_pressure_data) = ctalloc(int, num_patches);
    BCPressureDataValues(bc_pressure_data) = ctalloc(void **,
                                                     num_patches);

    ForEachPatch(num_patches, i)
    {
      BCPressureDataType(bc_pressure_data, i) = (public_xtra->input_types[i]);
      BCPressureDataPatchIndex(bc_pressure_data, i) = (public_xtra->patch_indexes[i]);
      BCPressureDataCycleNumber(bc_pressure_data, i) = (public_xtra->cycle_numbers[i]);

      interval_division = TimeCycleDataIntervalDivision(time_cycle_data, BCPressureDataCycleNumber(bc_pressure_data, i));
      BCPressureDataIntervalValues(bc_pressure_data, i) = ctalloc(void *, interval_division);
      ForEachInterval(interval_division, interval_number)
      {
        Do_SetupPatchIntervals(public_xtra, interval, i,
        {
          /* Setup a fixed pressure condition structure */
          SetupPatchInterval(DirEquilRefPatch,
          {
            BCPressureType0 *bc_pressure_type0;

            int phase;
            int num_phases = BCPressureDataNumPhases(bc_pressure_data);

            BCPressureDataBCType(bc_pressure_data, i) = DirichletBC;

            GetTypeStruct(DirEquilRefPatch, data, public_xtra, i);

            bc_pressure_type0 = ctalloc(BCPressureType0, 1);

            BCPressureType0RefSolid(bc_pressure_type0) =
              (data->reference_solid);

            BCPressureType0RefPatch(bc_pressure_type0) =
              (data->reference_patch);

            BCPressureType0Value(bc_pressure_type0) = (data->values[interval_number]);

            if (num_phases > 1)
            {
              BCPressureType0ValueAtInterfaces(bc_pressure_type0) = ctalloc(double, (num_phases - 1));
              for (phase = 1; phase < num_phases; phase++)
              {
                BCPressureType0ValueAtInterface(bc_pressure_type0, phase)
                  = ((data->value_at_interface[interval_number])[phase - 1]);
              }
            }
            else
            {
              BCPressureType0ValueAtInterfaces(bc_pressure_type0) = NULL;
            }

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number) = (void*)bc_pressure_type0;
          });

          /* Setup a piecewise linear pressure condition structure */
          SetupPatchInterval(DirEquilPLinear,
          {
            BCPressureType1   *bc_pressure_type1;
            int point, phase;
            int num_points;
            int num_phases = BCPressureDataNumPhases(bc_pressure_data);

            BCPressureDataBCType(bc_pressure_data, i) = DirichletBC;

            GetTypeStruct(DirEquilPLinear, data, public_xtra, i);

            num_points = (data->num_points[interval_number]);

            bc_pressure_type1 = ctalloc(BCPressureType1, 1);

            BCPressureType1XLower(bc_pressure_type1) = (data->xlower[interval_number]);
            BCPressureType1YLower(bc_pressure_type1) = (data->ylower[interval_number]);
            BCPressureType1XUpper(bc_pressure_type1) = (data->xupper[interval_number]);
            BCPressureType1YUpper(bc_pressure_type1) = (data->yupper[interval_number]);
            BCPressureType1NumPoints(bc_pressure_type1) = (data->num_points[interval_number]);

            BCPressureType1Points(bc_pressure_type1) = ctalloc(double, num_points);
            BCPressureType1Values(bc_pressure_type1) = ctalloc(double, num_points);

            for (point = 0; point < num_points; point++)
            {
              BCPressureType1Point(bc_pressure_type1, point) = ((data->points[interval_number])[point]);
              BCPressureType1Value(bc_pressure_type1, point) = ((data->values[interval_number])[point]);
            }

            if (num_phases > 1)
            {
              BCPressureType1ValueAtInterfaces(bc_pressure_type1) = ctalloc(double, (num_phases - 1));

              for (phase = 1; phase < num_phases; phase++)
              {
                BCPressureType1ValueAtInterface(bc_pressure_type1, phase) = ((data->value_at_interface[interval_number])[phase - 1]);
              }
            }
            else
            {
              BCPressureType1ValueAtInterfaces(bc_pressure_type1) = NULL;
            }

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number) = (void*)bc_pressure_type1;
          });

          /* Setup a constant flux condition structure */
          SetupPatchInterval(FluxConst,
          {
            BCPressureType2   *bc_pressure_type2;

            BCPressureDataBCType(bc_pressure_data, i) = FluxBC;

            GetTypeStruct(FluxConst, data, public_xtra, i);

            bc_pressure_type2 = ctalloc(BCPressureType2, 1);

            BCPressureType2Value(bc_pressure_type2)
              = (data->values[interval_number]);

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type2;
          });

          /* Setup a volumetric flux condition structure */
          SetupPatchInterval(FluxVolumetric,
          {
            BCPressureType3   *bc_pressure_type3;

            BCPressureDataBCType(bc_pressure_data, i) = FluxBC;

            GetTypeStruct(FluxVolumetric, data, public_xtra, i);

            bc_pressure_type3 = ctalloc(BCPressureType3, 1);

            BCPressureType3Value(bc_pressure_type3)
              = (data->values[interval_number]);

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type3;
          });

          /* Setup a file defined pressure condition structure */
          SetupPatchInterval(PressureFile,
          {
            BCPressureType4   *bc_pressure_type4;

            BCPressureDataBCType(bc_pressure_data, i) = DirichletBC;

            GetTypeStruct(PressureFile, data, public_xtra, i);

            bc_pressure_type4 = ctalloc(BCPressureType4, 1);

            BCPressureType4FileName(bc_pressure_type4)
              = ctalloc(char, strlen((data->filenames)[interval_number]) + 1);

            strcpy(BCPressureType4FileName(bc_pressure_type4),
                   ((data->filenames)[interval_number]));

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type4;
          });

          /* Setup a file defined flux condition structure */
          SetupPatchInterval(FluxFile,
          {
            BCPressureType5   *bc_pressure_type5;

            BCPressureDataBCType(bc_pressure_data, i) = FluxBC;

            GetTypeStruct(FluxFile, data, public_xtra, i);

            bc_pressure_type5 = ctalloc(BCPressureType5, 1);

            BCPressureType5FileName(bc_pressure_type5)
              = ctalloc(char, strlen((data->filenames)[interval_number]) + 1);

            strcpy(BCPressureType5FileName(bc_pressure_type5),
                   ((data->filenames)[interval_number]));

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type5;
          });

          /* Setup a Dir. pressure MATH problem condition */
          SetupPatchInterval(ExactSolution,
          {
            BCPressureType6   *bc_pressure_type6;

            BCPressureDataBCType(bc_pressure_data, i) = DirichletBC;

            GetTypeStruct(ExactSolution, data, public_xtra, i);

            bc_pressure_type6 = ctalloc(BCPressureType6, 1);

            BCPressureType6FunctionType(bc_pressure_type6)
              = (data->function_type);

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type6;
          });

          /*//sk  Setup a overland flow condition structure */
          SetupPatchInterval(OverlandFlow,
          {
            BCPressureType7   *bc_pressure_type7;

            BCPressureDataBCType(bc_pressure_data, i) = OverlandBC;

            GetTypeStruct(OverlandFlow, data, public_xtra, i);

            bc_pressure_type7 = ctalloc(BCPressureType7, 1);

            BCPressureType7Value(bc_pressure_type7)
              = (data->values[interval_number]);

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type7;
          });

          /* Setup a file defined flux condition structure for overland flow BC*/
          SetupPatchInterval(OverlandFlowPFB,
          {
            BCPressureType8   *bc_pressure_type8;

            BCPressureDataBCType(bc_pressure_data, i) = OverlandBC;

            GetTypeStruct(OverlandFlowPFB, data, public_xtra, i);

            bc_pressure_type8 = ctalloc(BCPressureType8, 1);

            BCPressureType8FileName(bc_pressure_type8)
              = ctalloc(char, strlen((data->filenames)[interval_number]) + 1);

            strcpy(BCPressureType8FileName(bc_pressure_type8),
                   ((data->filenames)[interval_number]));

            BCPressureDataIntervalValue(bc_pressure_data, i, interval_number)
              = (void*)bc_pressure_type8;

            break;
          });
        }); /* End Do_SetupPatchIntervals */
      } /* End ForEachInterval */
    } /* End ForEachPatch */
  }
}


/*--------------------------------------------------------------------------
 * BCPressurePackageInitInstanceXtra
 *--------------------------------------------------------------------------*/

PFModule *BCPressurePackageInitInstanceXtra(
                                            Problem *problem)
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra;

  if (PFModuleInstanceXtra(this_module) == NULL)
    instance_xtra = ctalloc(InstanceXtra, 1);
  else
    instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (problem != NULL)
  {
    (instance_xtra->problem) = problem;
  }

  PFModuleInstanceXtra(this_module) = instance_xtra;
  return this_module;
}


/*--------------------------------------------------------------------------
 * BCPressurePackageFreeInstanceXtra
 *--------------------------------------------------------------------------*/

void  BCPressurePackageFreeInstanceXtra()
{
  PFModule      *this_module = ThisPFModule;
  InstanceXtra  *instance_xtra = (InstanceXtra*)PFModuleInstanceXtra(this_module);

  if (instance_xtra)
  {
    tfree(instance_xtra);
  }
}


/*--------------------------------------------------------------------------
 * BCPressurePackageNewPublicXtra
 *--------------------------------------------------------------------------*/
PFModule  *BCPressurePackageNewPublicXtra(
                                          int num_phases)
{
  PFModule      *this_module = ThisPFModule;
  PublicXtra    *public_xtra;

  char *patch_names;

  char *patch_name;

  char key[IDB_MAX_KEY_LEN];

  char *switch_name;

  char *cycle_name;

  int global_cycle;

  int domain_index;

  int phase;
  char *interval_name;

  int num_patches;
  int num_cycles;

  int i, interval_number, interval_division;

  NameArray type_na;
  NameArray function_na;

  //sk
  //type_na = NA_NewNameArray("DirEquilRefPatch DirEquilPLinear FluxConst FluxVolumetric PressureFile FluxFile ExactSolution OverlandFlow OverlandFlowPFB");
  type_na = NA_NewNameArray(TOSTR(BC_TYPE_NAMES));

  function_na = NA_NewNameArray("dum0 X XPlusYPlusZ X3Y2PlusSinXYPlus1 X3Y4PlusX2PlusSinXYCosYPlus1 XYZTPlus1 XYZTPlus1PermTensor");

  /* allocate space for the public_xtra structure */
  public_xtra = ctalloc(PublicXtra, 1);

  (public_xtra->num_phases) = num_phases;

  patch_names = GetString("BCPressure.PatchNames");

  public_xtra->patches = NA_NewNameArray(patch_names);
  num_patches = NA_Sizeof(public_xtra->patches);

  (public_xtra->num_patches) = num_patches;
  (public_xtra->num_cycles) = num_patches;

  if (num_patches > 0)
  {
    (public_xtra->num_cycles) = num_cycles = num_patches;

    (public_xtra->interval_divisions) = ctalloc(int, num_cycles);
    (public_xtra->intervals) = ctalloc(int *, num_cycles);
    (public_xtra->repeat_counts) = ctalloc(int, num_cycles);

    (public_xtra->input_types) = ctalloc(int, num_patches);
    (public_xtra->patch_indexes) = ctalloc(int, num_patches);
    (public_xtra->cycle_numbers) = ctalloc(int, num_patches);
    (public_xtra->data) = ctalloc(void *, num_patches);

    /* Determine the domain geom index from domain name */
    switch_name = GetString("Domain.GeomName");
    domain_index = NA_NameToIndex(GlobalsGeomNames, switch_name);
    if (domain_index < 0)
    {
      InputError("Error: invalid geometry name <%s> for key <%s>\n",
                 switch_name, "Domain.GeomName");
    }

    ForEachPatch(num_patches, i)
    {
      patch_name = NA_IndexToName(public_xtra->patches, i);

      public_xtra->patch_indexes[i] =
        NA_NameToIndex(GlobalsGeometries[domain_index]->patches,
                       patch_name);

      sprintf(key, "Patch.%s.BCPressure.Type", patch_name);
      switch_name = GetString(key);

      public_xtra->input_types[i] = NA_NameToIndex(type_na, switch_name);

      if (public_xtra->input_types[i] < 0)
      {
        InputError("Error: invalid type name <%s> for <%s>\n",
                   switch_name, key);
      }

      sprintf(key, "Patch.%s.BCPressure.Cycle", patch_name);
      cycle_name = GetString(key);

      public_xtra->cycle_numbers[i] = i;

      global_cycle = NA_NameToIndex(GlobalsCycleNames, cycle_name);

      if (global_cycle < 0)
      {
        InputError("Error: invalid cycle name <%s> for <%s>\n",
                   cycle_name, key);
      }


      interval_division = public_xtra->interval_divisions[i] =
        GlobalsIntervalDivisions[global_cycle];

      public_xtra->repeat_counts[i] =
        GlobalsRepeatCounts[global_cycle];

      (public_xtra->intervals[i]) = ctalloc(int, interval_division);

      ForEachInterval(interval_division, interval_number)
      {
        public_xtra->intervals[i][interval_number] =
          GlobalsIntervals[global_cycle][interval_number];
      }

      Do_SetupPatchTypes(public_xtra, interval, i,
      {
        SetupPatchType(DirEquilRefPatch,
        {
          int size;

          NewTypeStruct(DirEquilRefPatch, data);

          (data->values) = ctalloc(double,
                                     interval_division);
          (data->value_at_interface) = ctalloc(double *,
                                                 interval_division);

          sprintf(key, "Patch.%s.BCPressure.RefGeom", patch_name);
          switch_name = GetString(key);

          data->reference_solid = NA_NameToIndex(GlobalsGeomNames,
                                                   switch_name);

          if (data->reference_solid < 0)
          {
            InputError("Error: invalid geometry name <%s> for reference solid <%s>\n", switch_name, key);
          }

          sprintf(key, "Patch.%s.BCPressure.RefPatch", patch_name);
          switch_name = GetString(key);

          data->reference_patch =
            NA_NameToIndex(GeomSolidPatches(
                                            GlobalsGeometries[data->reference_solid]),
                           switch_name);

          if (data->reference_patch < 0)
          {
            InputError("Error: invalid reference patch name <%s> for key <%s>\n", switch_name, key);
          }

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.Value",
                    patch_name,
                    NA_IndexToName(
                                   GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->values[interval_number] = GetDouble(key);

            if (num_phases > 1)
            {
              size = (num_phases - 1);

              (data->value_at_interface[interval_number]) =
                ctalloc(double, size);

              for (phase = 1; phase < num_phases; phase++)
              {
                sprintf(key, "Patch.%s.BCPressure.%s.%s.IntValue",
                        patch_name,
                        NA_IndexToName(
                                       GlobalsIntervalNames[global_cycle],
                                       interval_number),
                        NA_IndexToName(GlobalsPhaseNames, phase));

                data->value_at_interface[interval_number][phase] =
                  GetDouble(key);
              }
            }
            else
            {
              (data->value_at_interface[interval_number]) = NULL;
            }
          }

          (public_xtra->data[i]) = (void*)data;
        });

        SetupPatchType(DirEquilPLinear,
        {
          int k, num_points, size;

          NewTypeStruct(DirEquilPLinear, data);

          (data->xlower) = ctalloc(double, interval_division);
          (data->ylower) = ctalloc(double, interval_division);
          (data->xupper) = ctalloc(double, interval_division);
          (data->yupper) = ctalloc(double, interval_division);

          (data->num_points) = ctalloc(int, interval_division);

          (data->points) = ctalloc(double *, interval_division);
          (data->values) = ctalloc(double *, interval_division);
          (data->value_at_interface) = ctalloc(double *, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            interval_name =
              NA_IndexToName(GlobalsIntervalNames[global_cycle],
                             interval_number);

            /* read in the xy-line */
            sprintf(key, "Patch.%s.BCPressure.%s.XLower", patch_name,
                    interval_name);
            data->xlower[interval_number] = GetDouble(key);

            sprintf(key, "Patch.%s.BCPressure.%s.YLower", patch_name,
                    interval_name);
            data->ylower[interval_number] = GetDouble(key);

            sprintf(key, "Patch.%s.BCPressure.%s.XUpper", patch_name,
                    interval_name);
            data->xupper[interval_number] = GetDouble(key);

            sprintf(key, "Patch.%s.BCPressure.%s.YUpper", patch_name,
                    interval_name);
            data->yupper[interval_number] = GetDouble(key);

            /* read num_points */
            sprintf(key, "Patch.%s.BCPressure.%s.NumPoints", patch_name,
                    interval_name);
            num_points = GetInt(key);

            (data->num_points[interval_number]) = num_points;

            (data->points[interval_number]) = ctalloc(double, num_points);
            (data->values[interval_number]) = ctalloc(double, num_points);

            for (k = 0; k < num_points; k++)
            {
              sprintf(key, "Patch.%s.BCPressure.%s.%d.Location",
                      patch_name, interval_name, k);
              data->points[interval_number][k] = GetDouble(key);

              sprintf(key, "Patch.%s.BCPressure.%s.%d.Value",
                      patch_name, interval_name, k);
              data->values[interval_number][k] = GetDouble(key);
            }

            if (num_phases > 1)
            {
              size = (num_phases - 1);

              (data->value_at_interface[interval_number]) =
                ctalloc(double, size);

              for (phase = 1; phase < num_phases; phase++)
              {
                sprintf(key, "Patch.%s.BCPressure.%s.%s.IntValue",
                        patch_name,
                        NA_IndexToName(
                                       GlobalsIntervalNames[global_cycle],
                                       interval_number),
                        NA_IndexToName(GlobalsPhaseNames, phase));

                data->value_at_interface[interval_number][phase] =
                  GetDouble(key);
              }
            }
            else
            {
              (data->value_at_interface[interval_number]) = NULL;
            }
          }

          StoreTypeStruct(public_xtra, data, i);
        });

        SetupPatchType(FluxConst,
        {
          NewTypeStruct(FluxConst, data);

          (data->values) = ctalloc(double, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.Value",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->values[interval_number] = GetDouble(key);
          }

          StoreTypeStruct(public_xtra, data, i);
        });

        SetupPatchType(FluxVolumetric,
        {
          NewTypeStruct(FluxVolumetric, data);

          (data->values) = ctalloc(double, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.Value",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->values[interval_number] = GetDouble(key);
          }

          StoreTypeStruct(public_xtra, data, i);
        });

        SetupPatchType(PressureFile,
        {
          NewTypeStruct(PressureFile, data);

          (data->filenames) = ctalloc(char *, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.FileName",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->filenames[interval_number] = GetString(key);
          }

          StoreTypeStruct(public_xtra, data, i);
        });

        SetupPatchType(FluxFile,
        {
          NewTypeStruct(FluxFile, data);

          (data->filenames) = ctalloc(char *, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.FileName",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->filenames[interval_number] = GetString(key);
          }

          StoreTypeStruct(public_xtra, data, i);
        });

        SetupPatchType(ExactSolution,
        {
          NewTypeStruct(ExactSolution, data);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.PredefinedFunction",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));
            switch_name = GetString(key);

            data->function_type = NA_NameToIndex(function_na,
                                                   switch_name);

            if (data->function_type < 0)
            {
              InputError("Error: invalid function type <%s> for key <%s>\n",
                         switch_name, key);
            }

            // MCB: This is overwriting the data struct inside the for loop
            // Also structured this way in master branch without changes
            // Bug? Intentional?
            StoreTypeStruct(public_xtra, data, i);
          }
        });


        SetupPatchType(OverlandFlow,
        {
          NewTypeStruct(OverlandFlow, data);

          (data->values) = ctalloc(double, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.Value",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->values[interval_number] = GetDouble(key);
          }

          StoreTypeStruct(public_xtra, data, i);
        });

        SetupPatchType(OverlandFlowPfB,
        {
          NewTypeStruct(OverlandFlowPFB, data);

          (data->filenames) = ctalloc(char *, interval_division);

          ForEachInterval(interval_division, interval_number)
          {
            sprintf(key, "Patch.%s.BCPressure.%s.FileName",
                    patch_name,
                    NA_IndexToName(GlobalsIntervalNames[global_cycle],
                                   interval_number));

            data->filenames[interval_number] = GetString(key);
          }

          StoreTypeStruct(public_xtra, data, i);
        });
      });      /* End Do_SetupPatchTypes */
    }
  }

  NA_FreeNameArray(type_na);
  NA_FreeNameArray(function_na);

  PFModulePublicXtra(this_module) = public_xtra;
  return this_module;
}


/*-------------------------------------------------------------------------
 * BCPressurePackageFreePublicXtra
 *-------------------------------------------------------------------------*/

void  BCPressurePackageFreePublicXtra()
{
  PFModule    *this_module = ThisPFModule;
  PublicXtra  *public_xtra = (PublicXtra*)PFModulePublicXtra(this_module);

  Type0         *dummy0;
  Type1         *dummy1;
  Type2         *dummy2;
  Type3         *dummy3;
  Type4         *dummy4;
  Type5         *dummy5;
  Type6         *dummy6;
  Type7         *dummy7;
  Type8         *dummy8;


  int num_patches, num_cycles;
  int i, interval_number, interval_division;

  if (public_xtra)
  {
    /* Free the well information */
    num_patches = (public_xtra->num_patches);
    NA_FreeNameArray(public_xtra->patches);

    if (num_patches > 0)
    {
      ForEachPatch(num_patches, i)
      {
        interval_division = (public_xtra->interval_divisions[(public_xtra->cycle_numbers[i])]);
        Do_FreePatches(public_xtra, i,
        {
          FreePatch(DirEquilRefPatch,
          {
            GetTypeStruct(DirEquilRefPatch, data, public_xtra, i);

            ForEachInterval(interval_division, interval)
            {
              tfree((data->value_at_interface[interval_number]));
            }

            tfree((data->value_at_interface));
            tfree((data->values));

            tfree(data);
          });

          FreePatch(DirEquilPLinear,
          {
            int interval_number;

            GetTypeStruct(DirEquilPLinear, data, public_xtra, i);

            ForEachInterval(interval_division, interval)
            {
              tfree((data->value_at_interface[interval_number]));
              tfree((data->values[interval_number]));
              tfree((data->points[interval_number]));
            }

            tfree((data->value_at_interface));
            tfree((data->points));
            tfree((data->values));

            tfree((data->num_points));

            tfree((data->yupper));
            tfree((data->xupper));
            tfree((data->ylower));
            tfree((data->xlower));

            tfree(data);
          });

          FreePatch(FluxConst,
          {
            GetTypeStruct(FluxConst, data, public_xtra, i);
            tfree((data->values));
            tfree(data);
          });

          FreePatch(FluxVolumetric,
          {
            GetTypeStruct(FluxVolumetric, data, public_xtra, i);
            tfree((data->values));
            tfree(data);
          });

          FreePatch(PressureFile,
          {
            int interval_number;
            GetTypeStruct(PressureFile, data, public_xtra, i);

            ForEachInterval(interval_division, interval)
            {
              tfree(((data->filenames)[interval_number]));
            }
            // @RMM had to remove to not error our
            tfree((data->filenames));
            tfree(data);
          });

          FreePatch(FluxFile,
          {
            int interval_number;
            GetTypeStruct(FluxFile, data, public_xtra, i);
            ForEachInterval(interval_division, interval)
            {
              tfree(((data->filenames)[interval_number]));
            }

            tfree((data->filenames));
            tfree(data);
          });

          FreePatch(ExactSolution,
          {
            GetTypeStruct(ExactSolution, data, public_xtra, i);
            tfree(data);
          });

          //sk
          FreePatch(OverlandFlow,
          {
            GetTypeStruct(OverlandFlow, data, public_xtra, i);
            tfree((data->values));
            tfree(data);
          });

          //RMM
          FreePatch(OverlandFlowPFB,
          {
            int interval_number;

            GetTypeStruct(OverlandFlowPFB, data, public_xtra, i);

            /* for(interval_number = 0; interval_number < interval_division; interval_number++)
             * {
             *   tfree(((dummy8 -> filenames)[interval_number]));
             * }  */
            // @RMM had to remove to not error our

            tfree((data->filenames));
            tfree(data);
          });
        }); /* End Do_FreePatches */
      } /* End ForEachPatch */

      tfree((public_xtra->data));
      tfree((public_xtra->cycle_numbers));
      tfree((public_xtra->patch_indexes));
      tfree((public_xtra->input_types));

      /* Free the time cycling information */
      num_cycles = (public_xtra->num_cycles);

      tfree((public_xtra->repeat_counts));

      for (i = 0; i < num_cycles; i++)
      {
        tfree((public_xtra->intervals[i]));
      }
      tfree((public_xtra->intervals));

      tfree((public_xtra->interval_divisions));
    }

    tfree(public_xtra);
  }
}


/*--------------------------------------------------------------------------
 * BCPressurePackageSizeOfTempData
 *--------------------------------------------------------------------------*/

int  BCPressurePackageSizeOfTempData()
{
  return 0;
}
