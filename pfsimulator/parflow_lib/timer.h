#ifndef _TIMER_H
#define _TIMER_H

#define GETTIME(tval)                           \
  (double)tval.tv_sec + (double)tval.tv_usec * 1.e-6;


#define START_TIMER(start_time)                 \
  struct timeval tp;                            \
  gettimeofday(&tp, NULL);                      \
  start_time = GETTIME(tp);

#define STOP_TIMER(start_time, fp)              \
  gettimeofday(&tp, NULL);                      \
  double end_time = GETTIME(tp);                \
  fprintf(fp, "%lf\n", end_time - start_time);

#endif // _TIMER_H
