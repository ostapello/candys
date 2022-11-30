#ifndef PROLE_H
#define PROLE_H

#include "errors.h"

int get_arguments (int argc, char *argv[], int &n, int &m, int &r, int &s, char **file_name);
void print_err (int err, char *argv[]);

void printf_main_process (const char *format, ...);
void printf_i_process (int i, const char *format, ...);
void print_i_array_in_i_proc (int *A, int size_A, int proc_number, const char *text = nullptr);
void print_d_array_in_i_proc (double *A, int size_A, int proc_number, int max_print = MAX_PRINT);

#endif // PROLE_H
