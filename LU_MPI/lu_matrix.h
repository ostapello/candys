#ifndef LU_MATRIX_H
#define LU_MATRIX_H

#include "errors.h"
#include "prole.h"

double f (int s, int n, int i, int j);

class LU_matrix
{
private:
  // matrix's parametrs
  const int size_of_system = 0;
  const int size_of_block = 0;
  const int size_of_tail = 0;
  const int number_of_blocks_row = 0;
  double norm = 0;
  static constexpr double eps = 1e-16;

  // utility parametrs
  int max_number_of_offdiag_block = 0;
  int number_of_proc = 0;
  int proc_number = 0;
  int respon_proc_for_tail = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  int size_of_box_of_blocks = 0;
  int tail = 0; // tail = (size_of_tail? 1 : 0);
  int mm_size = 0;
  int ml_size = 0;
  int ll_size = 0;
  std::unique_ptr <int []> L_offsets_offdiag_block_ptr;
  std::unique_ptr <int []> U_offsets_offdiag_block_ptr;

  // value's arrays
  std::unique_ptr <double []> L_off_diag_full_block_ptr;
  std::unique_ptr <double []> U_off_diag_full_block_ptr;
  std::unique_ptr <double []> diag_full_block_ptr;

  std::unique_ptr <double []> L_tail_bloc_ptr;
  std::unique_ptr <double []> U_tail_bloc_ptr;
  std::unique_ptr <double []> diag_tail_ptr;

  // buffer's arrays
  std::unique_ptr <double []> buf_block_row_ptr_1;
  std::unique_ptr <double []> buf_block_row_ptr_2;
  std::unique_ptr <double []> buf_max_block_row_for_proc_ptr;
  //std::unique_ptr <double []> buf_block_row_ptr_3;
  std::unique_ptr <int []> pointers_array_ptr;
  std::unique_ptr <double []> buf_5_blocks;

public:
  LU_matrix (const int size_of_system_arg, const int size_of_block_arg) :
    size_of_system (size_of_system_arg),
    size_of_block (size_of_block_arg),
    size_of_tail (size_of_system_arg % size_of_block_arg),
    number_of_blocks_row (size_of_system_arg / size_of_block_arg),
    tail ((size_of_tail ? 1 : 0)),
    mm_size (size_of_block * size_of_block),
    ml_size (size_of_block * size_of_tail),
    ll_size (size_of_tail * size_of_tail)
  {
    MPI_Comm_size (comm, &number_of_proc); //take number of procces
    MPI_Comm_rank (comm, &proc_number); //take procces's number

    init_offsets_offdiag_block (); // also determine respon_proc_for_tail
    allocate_memory ();
  }

  int init_matrix (int formula_number, char *file_name = nullptr);
  void print_matrix (const char *text = "MATRIX", int max_print = MAX_PRINT, const char *file_name = nullptr);
  void build_rhs (double *B, double *X);
  int make_factorization ();
  int find_x (double *B, double *X);
  double find_res (double *X, double *B);
private:
  void allocate_memory ();
  void init_offsets_offdiag_block ();

  void fill_first_blocks_for_L (int count_of_block, int formula_number);
  void unpack_first_blocks_for_L (int count_of_block);
  void place_buf_pointers (int start_block);
  void fill_blocks_for_diag_and_U (int i_block_row, int formula_number);
  void unpack_blocks_for_diag_and_U (int i_block_row);
  void fill_L_tail (int formula_number);
  void fill_matrix (int formula_number);
  int read_first_blocks_for_L (int count_of_block, FILE *fp);
  int read_blocks_for_diag_and_U (int i_block_row, FILE *fp);
  int read_L_tail (FILE *fp);
  void unpack_L_tail ();
  int read_matrix (char *file_name);
  void collet_first_L_and_diag_blocks (int i_block_row);
  void print_block_row_in_i_proc (double *A, int row_size, int col_size, int count_of_blocks, int i_proc,
                                  int max_print = MAX_PRINT);
  void print_block_i_row (double *block_row, int i_block_row, int row_size, int max_print = MAX_PRINT, FILE *fp = stdout);
  int next_proc (int cur_proc);
  void sum_for_vector_B (double *block, double *vector,int size_of_row, int size_of_col, int start);
  void work_with_L (int i_block_row, double *B);
  int get_global_number_of_block_column (int i_block_row);
  void work_with_U (int i_block_row, double *B);
  void work_with_L_tail (double *B);

  int get_owner_of_big_part (int i_block_row);
  void copy_offsets ();
  void move_l2r (double *left_buf, int num_row, int num_col, double *right_buf);

  int make_prepare_mult_for_diag_block (int i_block_row, int start_block);
  int mult_and_factor_block (int i_block_row, int start_block);
  void fill_column_for_send (int i_block_row);
  void fill_row_for_send (int i_block_row);
  void find_L_column (int i_block_row);
  void find_U_row (int i_block_row);
  void find_L_column_tail (int i_block_row);
  void find_U_row_tail (int i_block_row);
  int find_last_diag_tail_block ();

  void mult_blocks_block (double *a, double *b, double *c, int I, int S, int J);
  void mult_blocks (double *a, double *b, double *c, int i, int s, int j);
  void mult_block_U (double *a, double *u, double *c, int I, int J);
  void mult_block_L (double *l, double *a, double *c, int I, int J);
  void add_blocks (double *a, double *b, double *c, int M, int L);
  void subst_blocks (double *a, double *b, double *c, int M, int L);
  void find_Ur (double *a, double *b, int l);
  int find_Lr (double *a, double *b, int l);
  int factor_block (double *a, int m);

  int find_block_y (int i_block_row, double *X, double *B);
  void find_block_x (int i_block_row, double *X, double *B);
  void L_update_b (int i_block_row, double *X, double *B);
  void U_update_b (int i_block_row, double *X, double *B);
  void mult_LX (double *L, double *X, double *B, int n_row, int n_col);
  int calc_tail_L (double *X, double *B);
  void calc_tail_U (double *X, double *B);

  void mult_block_i_row (double *recv_buf, int i_block_row, int size_of_block, double *X, double *B, double &max_num, double &max_det);
public:
  int    * get_offsets_offdiag_block   ()       {return U_offsets_offdiag_block_ptr.get ();}
  int      get_number_of_blocks_row    () const {return number_of_blocks_row;}
  int      get_number_of_offdiag_block () const {return U_offsets_offdiag_block_ptr[number_of_blocks_row];}
  double   get_norm                    () const {return norm;}
  double * get_prepare_buf             ()       {return buf_max_block_row_for_proc_ptr.get ();}
  double * get_block_1                 ()       {return buf_5_blocks.get ();}
  double * get_rev_LU_block            ()       {return buf_5_blocks.get () + mm_size;}
  double * get_block_2                 ()       {return buf_5_blocks.get () + 2 * mm_size;}
  double * get_diag_full_block (int i_row = 0) {return diag_full_block_ptr.get () + (i_row / number_of_proc)
                                                                                  * mm_size;}
  double * get_L_off_diag (int i_row = 0) {return L_off_diag_full_block_ptr.get () +
                                                  U_offsets_offdiag_block_ptr[i_row] * mm_size;}
  double * get_U_off_diag (int i_row = 0) {return U_off_diag_full_block_ptr.get ()
                                                + U_offsets_offdiag_block_ptr[i_row] * mm_size;}
  double * get_L_tail (int i_col = 0) {if (respon_proc_for_tail) return L_tail_bloc_ptr.get ()
                                                                      + i_col * ml_size;
                                       return nullptr;}
  double * get_U_tail (int i_row = 0) {if (respon_proc_for_tail) return U_tail_bloc_ptr.get ()
                                                                      + i_row * ml_size;
                                       return nullptr;}
  double * get_diag_tail () {if (respon_proc_for_tail) return diag_tail_ptr.get ();
                             return nullptr;}
  int get_number_of_offdiag_block_in_i_row (int i_block_row) const {return U_offsets_offdiag_block_ptr[i_block_row + 1]
                                                                         - U_offsets_offdiag_block_ptr[i_block_row];}
};


#endif // LU_MATRIX_H
