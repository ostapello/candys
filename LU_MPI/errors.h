#ifndef ERRORS_H
#define ERRORS_H

#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <fenv.h>
#include <unistd.h>
#include <memory>
#include <stdarg.h>
#include <iostream>
#include <math.h>

#define WRONG_ARG -1

enum LU_errors
{
  OPEN_FILE = -2,
  READ_ERR  = -3,
  BAD_MATRIX = -4,
  BM_FIND_LR = -5,
  FIND_X = -6
};



#define LEN 1234
#define MAX_PRINT 20
#define MAIN_PROCESS 0

#define PRINT_PROCESS 2
#endif // ERRORS_H
