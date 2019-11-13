#ifndef _PROBLEM_PHASE_REL_PERM_H
#define _PROBLEM_PHASE_REL_PERM_H

static inline double VanGLookupLinear(
                                      double     pressure_head,
                                      VanGTable *lookup_table,
                                      int        fcn)
{
  double rel_perm = 0.0;
  int pt = 0;
  int num_sample_points = lookup_table->num_sample_points;
  double min_pressure_head = lookup_table->min_pressure_head;
  int max = num_sample_points + 1;

  // This table goes from 0 to fabs(min_pressure_head)
  assert(pressure_head >= 0);

  // SGS TODO add warning in output?
  // Make sure values are in the table range, if lower then set to the 0.0 which is limit of VG curve
  if (pressure_head < fabs(min_pressure_head))
  {
    double interval = lookup_table->interval;

    // Use direct table lookup to avoid using this binary search since
    // we have uniformly spaced points.

    pt = (int)floor(pressure_head / interval);
    assert(pt < max);

    // using cubic Hermite interpolation


    if (fcn == CALCFCN)
    {
      rel_perm = lookup_table->a[pt] + lookup_table->slope[pt] * (pressure_head - lookup_table->x[pt]);
    }
    else
    {
      rel_perm = lookup_table->a_der[pt] + lookup_table->slope_der[pt] * (pressure_head - lookup_table->x[pt]);
    }
  }

  return rel_perm;
}

static inline double VanGLookupSpline(
                                      double     pressure_head,
                                      VanGTable *lookup_table,
                                      int        fcn)
{
  double rel_perm, t;
  int pt = 0;
  int num_sample_points = lookup_table->num_sample_points;
  double min_pressure_head = lookup_table->min_pressure_head;
  int max = num_sample_points + 1;

  // This table goes from 0 to fabs(min_pressure_head)
  assert(pressure_head >= 0);

  // SGS TODO add warning in output?
  // Make sure values are in the table range, if lower then set to the 0.0 which is limit of VG curve
  if (pressure_head >= fabs(min_pressure_head))
  {
    return 0.0;
  }
  else
  {
    // Use direct table lookup to avoid using this binary search since
    // we have uniformly spaced points.
    double interval = lookup_table->interval;
    pt = (int)floor(pressure_head / interval);
    if (pt > max)
    {
      pt = max - 1;
    }

#if 0
    // When using variably spaced interpolation points, use binary
    // search to find the interval
    {
      int min = 0;
      int mid;

      while (max != min + 1)
      {
        mid = min + floor((max - min) / 2);
        if (pressure_head == lookup_table->x[mid])
        {
          min = mid;
          max = min + 1;
        }
        if (pressure_head < lookup_table->x[mid])
        {
          max = mid;
        }
        else
        {
          min = mid;
        }
      }
      pt = min;
    }
#endif
  }

  double x = lookup_table->x[pt];
  double a = lookup_table->a[pt];
  double d = lookup_table->d[pt];
  double a_der = lookup_table->a_der[pt];
  double d_der = lookup_table->d_der[pt];

  // using cubic Hermite interpolation
  if (fcn == CALCFCN)
  {
    t = (pressure_head - x) / (lookup_table->x[pt + 1] - x);
    rel_perm = (2.0 * pow(t, 3) - 3.0 * pow(t, 2) + 1.0) * a
               + (pow(t, 3) - 2.0 * pow(t, 2) + t)
               * (lookup_table->x[pt + 1] - x) * d + (-2.0 * pow(t, 3)
                                                      + 3.0 * pow(t, 2)) * (lookup_table->a[pt + 1])
               + (pow(t, 3) - pow(t, 2)) * (lookup_table->x[pt + 1] - x)
               * (lookup_table->d[pt + 1]);
  }
  else
  {
    t = (pressure_head - x) / (lookup_table->x[pt + 1] - x);
    rel_perm = (2.0 * pow(t, 3) - 3.0 * pow(t, 2) + 1.0) * a_der
               + (pow(t, 3) - 2.0 * pow(t, 2) + t)
               * (lookup_table->x[pt + 1] - x) * d_der + (-2.0 * pow(t, 3)
                                                          + 3.0 * pow(t, 2)) * (lookup_table->a_der[pt + 1])
               + (pow(t, 3) - pow(t, 2)) * (lookup_table->x[pt + 1] - x)
               * (lookup_table->d_der[pt + 1]);
  }

  return rel_perm;
}

#endif // _PROBLEM_PHASE_REL_PERM_H
