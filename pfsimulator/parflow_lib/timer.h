#ifndef _TIMER_H
#define _TIMER_H

#define GETTIME(tval)                                 \
  (double)tval.tv_sec + (double)tval.tv_usec * 1.e-6;

#define START_TIMER()                           \
  double start_time = 0.0;                      \
  double end_time = 0.0;                        \
  double total_time = 0.0;                      \
  struct timeval tp;                            \
  gettimeofday(&tp, NULL);                      \
  start_time = GETTIME(tp);

#define END_TIMER(fp) STOP_TIMER(fp)

#define STOP_TIMER(fp)                          \
  gettimeofday(&tp, NULL);                      \
  end_time = GETTIME(tp);                       \
  total_time += (end_time - start_time);        \
  fprintf(fp, "%lf\n", total_time);

#endif // _TIMER_H
