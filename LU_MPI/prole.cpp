#include "prole.h"

void print_err (int err, char *argv[])
{
  int my_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (my_rank == MAIN_PROCESS)
    {
      switch (err)
        {
        case WRONG_ARG:
          printf ("USAGE   : ./a.out n m r s [filename.txt]\n");
          printf ("\t  n - size of matrix (n > 0)\n"
                  "\t  m - size of block  (m > 0 && m <= n && m %% 3 == 0)\n"
                  "\t  r - count of elements to display (r > 0)\n"
                  "\t  s - number of formula for input matrix (s >= 0)\n"
                  "\t  [filename.txt] - name of file with marix if s = 0\n"
                  "                         (._.)                      \n");
          break;
        case OPEN_FILE:
          printf ("ERROR: cannot open file: %s   (._.)\n", argv[5]);
          break;
        case READ_ERR:
          printf ("ERROR: some problem with reading elements in file: %s     (._.)\n", argv[5]);
          break;
        case BAD_MATRIX:
          printf ("ERROR: cannot make LU factorization     (._.)\n");
          break;
        case BM_FIND_LR:
          printf ("ERROR: some problem with LU factorization     (._.)\n");
          break;
        case FIND_X:
          printf ("ERROR: some problem with finding solution     (._.)\n");
          break;
        }
    }
}


int get_arguments (int argc, char *argv[], int &n, int &m, int &r, int &s, char **file_name)
{
  if (!((argc == 5) || (argc == 6)) ||
      sscanf (argv[1], "%d", &n) != 1 || n <= 0 ||
      sscanf (argv[2], "%d", &m) != 1 || m <= 0 || m > n || (m % 3) ||
      sscanf (argv[3], "%d", &r) != 1 || r <= 0 ||
      sscanf (argv[4], "%d", &s) != 1 || s < 0)
    {
      return WRONG_ARG;
    }
  if (argc == 5 && s == 0)
    return WRONG_ARG;
  if (argc == 6)
    {
      if (s != 0)
        return WRONG_ARG;
      *file_name = argv[5];
    }

  return 0;
}


void printf_main_process (const char *format, ...)
{
  int my_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  va_list args;
  va_start (args, format);

  if (my_rank == MAIN_PROCESS)
    vprintf (format, args);

  va_end (args);
}

void printf_i_process (int i, const char *format, ...)
{
  int my_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  va_list args;
  va_start (args, format);

  if (my_rank == i)
    vprintf (format, args);

  va_end (args);
}

void print_i_array_in_i_proc (int *A, int size_A, int proc_number, const char *text)
{
  int my_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (my_rank == proc_number)
    {
      if (text)
        {
          printf ("print_i_arr %s; size = %d; proc %d\n", text, size_A, proc_number);
        }
      else
        {
          printf ("print_i_arr; size = %d; proc %d\n", size_A, proc_number);
        }
      for (int i = 0; i < size_A; i++)
        {
          printf ("%3d", A[i]);
        }
      printf ("\n");
    }
}

void print_d_array_in_i_proc (double *A, int size_A, int proc_number, int max_print)
{
  int my_rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (my_rank == proc_number)
    {
      /*
      if (text)
        {
          printf ("print_d_arr %s; size = %d; proc %d\n", text, size_A, proc_number);
        }
      else
        {
          printf ("print_d_arr; size = %d; proc %d\n", size_A, proc_number);
        }
        */
      for (int i = 0; i < size_A && i < max_print; i++)
        {
          printf ("%10.3e ", A[i]);
        }
      printf ("\n");
    }
}

