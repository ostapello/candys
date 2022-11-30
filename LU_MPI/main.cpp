#include "lu_matrix.h"
#include "prole.h"

int MPI_main (int argc, char *argv[]);

int main (int argc, char *argv[])
{
  feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW| FE_UNDERFLOW);
  int err = MPI_main (argc, argv);
  if (err < 0)
    {
      print_err (err, argv);
    }
  MPI_Finalize ();
  return 0;
}

int MPI_main (int argc, char *argv[])
{
  MPI_Init (&argc, &argv);

  int k = 0, p = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size (comm, &p); //взять число процессов
  MPI_Comm_rank (comm, &k); //взять свой номер

  // get main parametrs
  int n, m, r, s;
  char *file_name = nullptr;
  int err = 0;

  err = get_arguments (argc, argv, n, m, r, s, &file_name);
  if (err < 0)
    {
      return err;
    }

  LU_matrix A (n, m);
  std::unique_ptr <double []> B_ptr (new double [n]);
  std::unique_ptr <double []> X_ptr (new double [n]);
  std::unique_ptr <double []> B_buf_ptr (new double [n]);

  printf_main_process ("... initiliazation of matrix\n");
  err = A.init_matrix (s, file_name);
  if (err < 0)
    {
      return err;
    }
  printf_main_process ("... success! norm = %10.3e\n", A.get_norm ());

  A.print_matrix ("ORIGIN MATRIX", r);

  printf_main_process ("... building of vector B\n");
  A.build_rhs (B_ptr.get (), X_ptr.get ());
  printf_main_process ("... success! vector B:\n");
  print_d_array_in_i_proc (B_ptr.get (), n, MAIN_PROCESS, r);


  printf_main_process ("... making factorization\n");
  double t = MPI_Wtime ();
  err = A.make_factorization ();
  if (err < 0)
    {
      return err;
    }
  t = MPI_Wtime () - t;
  printf_main_process ("... success! time = %.2f\n", t);

  //printf_main_process ("TEST: LU_matrix\n");
  //A.print_matrix ("CHANGED MATRIX", r, "out.txt");
  //A.print_matrix ("CHANGED MATRIX", r);


  printf_main_process ("... finding vector X:\n");
  err = A.find_x (B_ptr.get (), X_ptr.get ());
  if (err < 0)
    {
      return err;
    }
  printf_main_process ("... success! vector X:\n");
  print_d_array_in_i_proc (X_ptr.get (), n, MAIN_PROCESS, r);


  err = A.init_matrix (s, file_name);
  if (err < 0)
    {
      return err;
    }

  A.build_rhs (B_ptr.get (), B_buf_ptr.get ());

  double residual = A.find_res (X_ptr.get (), B_ptr.get ());

  if (k == MAIN_PROCESS)
    printf("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n", argv[0], residual, t, s, n, m, p);

  return 0;
}
