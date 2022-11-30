#include "lu_matrix.h"

double f (int s, int n, int i, int j)
{
  switch (s)
    {
    case 1:
      return n - (i > j ? i : j) + 1;
    case 2:
      return (i > j ? i : j);
    case 3:
      return (i > j ? i - j : j - i);
    case 4:
      return ((double) 1) / (i + j - 1);
    default:
      printf ("Something WRONG!!!!!\n");
      return 0;
    }
  return 0;
}

void LU_matrix::mult_blocks_block (double *a, double *b, double *c, int I, int S, int J)
{
  int i, j, s;
  double *ai, *ci, *bj;
  double c00, c01, c02, c10, c11, c12, c20, c21, c22;
  int m = size_of_block;

  for (i = 0; i < I; i += 3)
    {
      ai = a + i * m;
      ci = c + i * m;
      for (j = 0; j < J; j += 3)
        {
          bj = b + j;
          c00 = 0; c01 = 0; c02 = 0;
          c10 = 0; c11 = 0; c12 = 0;
          c20 = 0; c21 = 0; c22 = 0;
          for (s = 0; s < S; s++)
            {
              c00 += ai[s] * bj[s * m];
              c01 += ai[s] * bj[s * m + 1];
              c02 += ai[s] * bj[s * m + 2];

              c10 += ai[s + m] * bj[s * m];
              c11 += ai[s + m] * bj[s * m + 1];
              c12 += ai[s + m] * bj[s * m + 2];

              c20 += ai[s + m * 2] * bj[s * m];
              c21 += ai[s + m * 2] * bj[s * m + 1];
              c22 += ai[s + m * 2] * bj[s * m + 2];
            }
          ci[j] = c00;
          ci[j + 1] = c01;
          ci[j + 2] = c02;

          ci[j + m] = c10;
          ci[j + 1 + m] = c11;
          ci[j + 2 + m] = c12;

          ci[j + m * 2] = c20;
          ci[j + 1 + m * 2] = c21;
          ci[j + 2 + m * 2] = c22;
        }
    }
}

void LU_matrix::mult_blocks (double *a, double *b, double *c, int I, int S, int J)
{
  int i, j, s;
  double *ai, *bj, *ci, t;

  //printf ("I = %d, S = %d, J = %d\n", I, S, J);
  for (i = 0; i < I; i++)
    {
      ai = a + i * S;
      ci = c + i * J;
      for (j = 0; j < J; j++)
        {
          bj = b + j;
          t = 0;
          for (s = 0; s < S; s++)
            {
              t += ai[s] * bj[s * J];
            }
          ci[j] = t;
        }
    }
}

void LU_matrix::mult_block_U (double *a, double *u, double *c, int I, int J)
{
  int i, j, s;
  double *ai, *uj, *ci, t;
  int m = J;
  //посчитаем частями
  //сначала для всех элементов выше диагонали j > i;
  if (I == J)
    {
      for (i = 0; i < I - 1; i++)
        {
          ai = a + i * m;
          ci = c + i * m;
          for (j = i + 1; j < J; j++)
            {
              t = 0;
              uj = u + j;
              for (s = 0; s < j; s++)
                {
                  t += ai[s] * uj[s * m];
                  //printf ("i = %d, j = %d, s = %d, t = %lf\n", i, j, s, t);
                }
              ci[j] = t + ai[s];
            }
        }
    }
  else
    {
      //I < J
      //те, кто в квадратике
      for (i = 0; i < I; i++)
        {
          ai = a + i * m;
          ci = c + i * m;
          for (j = i + 1; j < J; j++)
            {
              t = 0;
              uj = u + j;
              for (s = 0; s < j; s++)
                {
                  t += ai[s] * uj[s * m];
                  //printf ("i = %d, j = %d, s = %d, t = %lf\n", i, j, s, t);
                }
              ci[j] = t + ai[s];
            }
        }
    }
  //теперь те, которые ниже диагонали
  for (i = 1; i < I; i++)
    {
      ai = a + i * m;
      ci = c + i * m;
      for (j = 0; j < i; j++)
        {
          t = 0;
          uj = u + j;
          for (s = 0; s < j; s++)
            {
              t += ai[s] * uj[s * m];
            }
          ci[j] = ai[j] + t;
        }
    }
  //диагональные
  for (i = 0; i < I; i++)
    {
      t = 0;
      ai = a + i * m;
      uj = u + i;
      for (s = 0; s < i; s++)
        {
          t += ai[s] * uj[s * m];
        }
      c[i * m + i] = t + ai[i];
    }
}

void LU_matrix::mult_block_L (double *l, double *a, double *c, int I, int J)
{
  int i, j, s;
  double *li, *aj, *ci, t;
  //посчитаем частями
  //сначала для всех элементов выше диагонали или на диагонали j >= i;
  for (i = 0; i < I; i++)
    {
      li = l + i * I;
      ci = c + i * J;
      for (j = i; j < J; j++)
        {
          t = 0;
          aj = a + j;
          for (s = 0; s < i + 1; s++)
            {
              t += li[s] * aj[s * J];
            }
          ci[j] = t;
        }
    }
  //теперь те, которые ниже диагонали i > j
  if (I == J)
    {
      for (i = 1; i < I; i++)
        {
          li = l + i * I;
          ci = c + i * J;
          for (j = 0; j < i; j++)
            {
              t = 0;
              aj = a + j;
              for (s = 0; s <= i; s++)
                {
                  t += li[s] * aj[s * J];
                }
              ci[j] = t;
            }
        }
    }
  else
    {
      //I > J
      for (i = 1; i < J; i++)
        {
          li = l + i * I;
          ci = c + i * J;
          for (j = 0; j < i; j++)
            {
              t = 0;
              aj = a + j;
              for (s = 0; s <= i; s++)
                {
                  t += li[s] * aj[s * J];
                }
              ci[j] = t;
            }
        }
      for (; i < I; i++)
        {
          li = l + i * I;
          ci = c + i * J;
          for (j = 0; j < J; j++)
            {
              t = 0;
              aj = a + j;
              for (s = 0; s <= i; s++)
                {
                  t += li[s] * aj[s * J];
                }
              ci[j] = t;
            }
        }
    }
}

void LU_matrix::add_blocks (double *a, double *b, double *c, int M, int L)
{
  int i, j;
  double *ai, *bi, *ci;
  int m = L;

  for (i = 0; i < M; i++)
    {
      ai = a + i * m;
      bi = b + i * m;
      ci = c + i * m;
      for (j = 0; j < L; j++)
        {
          ci[j] = ai[j] + bi[j];
        }
    }
}

void LU_matrix::subst_blocks (double *a, double *b, double *c, int M, int L)
{
  int i, j;
  double *ai, *bi, *ci;
  int m = L;

  for (i = 0; i < M; i++)
    {
      ai = a + i * m;
      bi = b + i * m;
      ci = c + i * m;
      for (j = 0; j < L; j++)
        {
          ci[j] = ai[j] - bi[j];
        }
    }
}

int LU_matrix::factor_block (double *a, int m)
{
  int i, k, l;
  double t, aii;
  double *lk, *ui, *uk, *ai;
  int col = m;

  if (m == 1)
    {
      //printf ("fabs(a[0]) = %e, norm = %e, eps = %e\n", fabs (a[0]), norm, eps);
      if (fabs(a[0]) < norm * eps)
        return BAD_MATRIX;
    }
  for (i = 1; i < m; i++)
    {
      ai = a + (i - 1) * col;
      aii = ai[i - 1];

      if (fabs(aii) < norm * eps)
        {
          return BAD_MATRIX;
        }

      for (k = i; k < m; k++)
        {
          t = 0;
          uk = a + k;

          for (l = 0; l < i - 1; l++)
            {
              t += ai[l] * uk[l * col];
            }

          ai[k] = (ai[k] - t) / aii;
        }
      ui = a + i;

      for (k = i; k < m; k++)
        {
          t = 0;
          lk = a + k * col;

          for (l = 0; l < i; l++)
            {
              t += lk[l] * ui[l * col];
            }

          ui[k * col] = ui[k * col] - t;
        }
    }
  return 0;
}

void LU_matrix::find_Ur (double *a, double *b, int l)
{
  int j, i, s;
  double *aj, *bj, bjj, *bi;
  int m = l;
  for (j = l - 1; j > 0; j--)
    {
      aj = a + j;
      bj = b + j;
      for (i = 0; i < j; i++)
        {
          bj[i * m] = -1 * aj[i * m];
        }
      /*
      printf ("1)b matrix:\n");
      print_matrix (b, m, m, m, m);
      */
      for (i = j + 1; i < l; i++)
        {
          bi = b + i;
          bjj = bi[j * m];
          for (s = 0; s < j; s++)
            {
              //printf ("j = %d, i = %d, s = %d, bjj = %lf\n", j, i, s, bjj);
              bi[s * m] += bj[s * m] * bjj;
            }
        }
      /*
      printf ("2)b matrix:\n");
      print_matrix (b, m, m, m, m);
      */
    }
}

int LU_matrix::find_Lr (double *a, double *b, int l)
{
  int j, i, s;
  double *aj, *bj, aii, bii, *bi, bjj;
  int m = l;

  for (j = 0; j < l; j++)
    {
      aj = a + j;
      bj = b + j;
      aii = aj[j * m];
      //printf ("aii = %lf j = %d\n", aii, j);
      if (fabs (aii) < norm * eps)
        return BM_FIND_LR;

      bii = 1 / aii;
      bj[j * m] = bii;
      bi = b + j * m;
      for (i = 0; i < j; i++)
        {
          //printf ("1)j = %d, i = %d, bii = %lf\n",j, i, bii);
          bi[i] *= bii;
        }
      /*
      printf ("1)b matrix:\n");
      print_matrix (b, m, m, m, m);
      */
      for (i = j + 1; i < l; i++)
        {
          bj[i * m] = -1 * aj[i * m] * bii;
        }
      /*
      printf ("2)b matrix:\n");
      print_matrix (b, m, m, m, m);
      */
      for (i = 0; i < j; i++)
        {
          bi = b + i;
          bjj = bi[j * m];
          for (s = j + 1; s < l; s++)
            {
              //printf ("j = %d, i = %d, s = %d, bjj = %lf\n", j, i, s, bjj);
              bi[s * m] -= bjj * aj[s * m];
            }
        }
      /*
      printf ("3)b matrix:\n");
      print_matrix (b, m, m, m, m);
      */
    }
  return 0;
}

void LU_matrix::copy_offsets ()
{
  memcpy (L_offsets_offdiag_block_ptr.get (), U_offsets_offdiag_block_ptr.get (), (number_of_blocks_row + 1) * sizeof (int));
}

void LU_matrix::move_l2r (double *left_buf, int num_row, int num_col, double *right_buf)
{
  memcpy (right_buf, left_buf, num_row * num_col * sizeof (double));
}

int LU_matrix::get_owner_of_big_part (int i_block_row)
{
  if (i_block_row == 0)
    return 0;
  if (i_block_row < number_of_proc)
    return i_block_row;
  return i_block_row % number_of_proc;
}

void LU_matrix::init_offsets_offdiag_block ()
{
  L_offsets_offdiag_block_ptr = std::make_unique <int []> (number_of_blocks_row + 1);
  U_offsets_offdiag_block_ptr = std::make_unique <int []> (number_of_blocks_row + 1);

  int *offsets_offdiag_block = U_offsets_offdiag_block_ptr.get ();
  int blocks_tail = number_of_blocks_row % number_of_proc;
  int max_num_of_offdiag_block_for_proc = number_of_blocks_row / number_of_proc
                                                 + (blocks_tail > proc_number ? 1 : 0);


  if (size_of_tail)
    {
      respon_proc_for_tail = (blocks_tail == proc_number ? 1 : 0);
    }
  //printf_i_process (0, "max_num_of_offdiag_block_for_proc = %d\n", max_num_of_offdiag_block_for_proc);
  int flag = proc_number;

  offsets_offdiag_block[0] = 0;
  for (int i_row = 0; i_row < number_of_blocks_row; i_row++)
    {
      if (i_row == flag)
        {
          max_num_of_offdiag_block_for_proc--;
          flag += number_of_proc;
        }
      offsets_offdiag_block[i_row + 1] = max_num_of_offdiag_block_for_proc;
    }
  //print_i_array_in_i_proc (offsets_offdiag_block, number_of_blocks_row + 1, PRINT_PROCESS, "in begining L");

  for (int i = 1; i < number_of_blocks_row + 1; i++)
    {
      offsets_offdiag_block[i] += offsets_offdiag_block[i - 1];
    }
  memcpy (L_offsets_offdiag_block_ptr.get (), offsets_offdiag_block, (number_of_blocks_row + 1) * sizeof (int));
}

void LU_matrix::allocate_memory ()
{
  int number_of_offdiag_block = get_number_of_offdiag_block ();
  MPI_Allreduce (&number_of_offdiag_block, &max_number_of_offdiag_block, 1, MPI_INT, MPI_MAX, comm);
  size_of_box_of_blocks = (number_of_blocks_row / number_of_proc) + 1;
  //printf_i_process (0, "boxs size = %d\n", size_of_box_of_blocks);


  int LU_offdiag_size = get_number_of_offdiag_block () * size_of_block * size_of_block;
  L_off_diag_full_block_ptr = std::make_unique <double []> (LU_offdiag_size);
  U_off_diag_full_block_ptr = std::make_unique <double []> (LU_offdiag_size);

  int diag_size = size_of_box_of_blocks * size_of_block * size_of_block;
  diag_full_block_ptr = std::make_unique <double []> (diag_size);


  if (respon_proc_for_tail)
    {
      int LU_tail_size = ml_size * number_of_blocks_row;
      L_tail_bloc_ptr = std::make_unique <double []> (LU_tail_size);
      U_tail_bloc_ptr = std::make_unique <double []> (LU_tail_size);
      diag_tail_ptr = std::make_unique <double []> (ll_size);

      memset (L_tail_bloc_ptr.get (), 0, sizeof (double) * LU_tail_size);
      memset (U_tail_bloc_ptr.get (), 0, sizeof (double) * LU_tail_size);
      memset (diag_tail_ptr.get (), 0, sizeof (double) * ll_size);
    }

  int buf_size = size_of_box_of_blocks * number_of_proc * mm_size;
  buf_block_row_ptr_1 = std::make_unique <double []> (buf_size);
  buf_block_row_ptr_2 = std::make_unique <double []> (buf_size);

  int buf_size_2 = U_offsets_offdiag_block_ptr.get ()[1] * mm_size;
  buf_max_block_row_for_proc_ptr = std::make_unique <double []> (buf_size_2);
  buf_5_blocks = std::make_unique <double []> (mm_size * 5);
  //buf_block_row_ptr_3 = std::make_unique <double []> (buf_size);

  pointers_array_ptr = std::make_unique <int []> (number_of_blocks_row + 1);
  memset (pointers_array_ptr.get (), 0, sizeof (int) * (number_of_blocks_row + 1));

  memset (L_off_diag_full_block_ptr.get (), 0, sizeof (double) * LU_offdiag_size);
  memset (U_off_diag_full_block_ptr.get (), 0, sizeof (double) * LU_offdiag_size);
  memset (diag_full_block_ptr.get (), 0, sizeof (double) * diag_size);


  memset (buf_block_row_ptr_1.get (), 0, sizeof (double) * buf_size);
  memset (buf_block_row_ptr_2.get (), 0, sizeof (double) * buf_size);

  memset (buf_max_block_row_for_proc_ptr.get (), 0, sizeof (double) * buf_size_2);
  //memset (buf_block_row_ptr_3.get (), 0, sizeof (double) * buf_size);
  memset (buf_5_blocks.get (), 0, sizeof (double) * 5 * mm_size);
}

void LU_matrix::fill_first_blocks_for_L (int count_of_block, int formula_number)
{
  double *buf = buf_block_row_ptr_1.get ();
  int i_global = count_of_block * size_of_block;
  //printf ("count_of_block = %d\n", count_of_block);
  double max = norm;
  int *p = pointers_array_ptr.get ();
  for (int i = 0; i < count_of_block; i++)
    {
      p[i]= size_of_block * size_of_block * i;
    }

  for (int i_loc = 0; i_loc < size_of_block; i_loc++)
    {
      int j_global = 0;
      double sum = 0;

      for (int j_bloc_col = 0; j_bloc_col < count_of_block; j_bloc_col++)
        {
          double *buf_p = buf + p[j_bloc_col];
          for (int j_loc = 0; j_loc < size_of_block; j_loc++)
            {
              //printf ("i_global = %d, j_global = %d, ", i_global, j_global);
              double t = f (formula_number, size_of_system, i_global + 1, ++j_global);
              //printf ("1 - t = %10.3e\n", t);
              buf_p[j_loc] = t;
              sum += fabs (t);
            }
          p[j_bloc_col] += size_of_block;
        }
      max = (sum > max ? sum : max);
      i_global++;
    }
  norm = max;
}

void LU_matrix::unpack_first_blocks_for_L (int count_of_block)
{
  //printf_i_process (PRINT_PROCCES, "tut? proc = %d, count_of_block = %d\n", proc_number, count_of_block);
  double * recv_buf = buf_block_row_ptr_2.get ();
  //print_block_row_in_i_proc (recv_buf, size_of_block, size_of_block, count_of_block, PRINT_PROCCES);
  double * L = get_L_off_diag ();
  int * L_offsets = L_offsets_offdiag_block_ptr.get ();
  //print_i_array_in_i_proc (L_offsets, number_of_blocks_row + 1, PRINT_PROCESS);
  int sq_size = size_of_block * size_of_block;

  //print_d_array_in_i_proc (recv_buf, count_of_block * size_of_block * size_of_block, PRINT_PROCCES, "recv_buf in unpack");
  for (int i_block = 0; i_block < count_of_block; i_block++)
    {
      double *p2begin_i_block = recv_buf + i_block * sq_size;
      double *p2begin_L_block = L + size_of_block * size_of_block * L_offsets[i_block];
      //printf_i_process (PRINT_PROCESS, "L_offsets[i_block] = %d\n", L_offsets[i_block]);
      move_l2r (p2begin_i_block, size_of_block, size_of_block, p2begin_L_block);
      //printf_i_process (PRINT_PROCESS, "DONE!\n");
      L_offsets[i_block]++;
    }
  //print_d_array_in_i_proc (L, size_of_block * size_of_block * get_number_of_offdiag_block (), PRINT_PROCCES, "L");
}

void LU_matrix::place_buf_pointers (int start_block)
{
  int sq_size = size_of_block * size_of_block;
  int *p = pointers_array_ptr.get ();

  memset (p, 0, sizeof (int) * (number_of_blocks_row + 1));

  int j = start_block % number_of_proc;
  int start_j = j;
  //printf ("j = %d\n", j);
  for (int i = start_block, diff = 0; i < number_of_blocks_row + tail; i++, j++)
    {
      if (j == number_of_proc)
        {
          j = 0;
        }
      if (j == start_j && i != start_block)
        diff += sq_size;
      p[i] = size_of_box_of_blocks * sq_size * j + diff;
    }
}

void LU_matrix::fill_blocks_for_diag_and_U (int i_block_row, int formula_number)
{
  place_buf_pointers (i_block_row);
  //print_i_array_in_i_proc (pointers_array_ptr.get (), number_of_blocks_row + 1, PRINT_PROCESS, "pointers");
  int *p = pointers_array_ptr.get ();
  double *buf = buf_block_row_ptr_1.get ();

  //print_i_array_in_i_proc (p, number_of_blocks_row + tail, 0);

  int i_global = i_block_row * size_of_block;
  double max = norm;
  for (int i_loc = 0; i_loc < size_of_block; i_loc++)
    {
      int j_global = i_block_row * size_of_block;
      double sum = 0;
      for (int j_bloc_col = i_block_row; j_bloc_col < number_of_blocks_row; j_bloc_col++)
        {
          double * buf_p = buf + p[j_bloc_col];
          for (int j_loc = 0; j_loc < size_of_block; j_loc++)
            {
              //printf ("i_global = %d, j_global = %d, ", i_global, j_global);
              double t = f (formula_number, size_of_system, i_global + 1, ++j_global);
              //printf ("1 - t = %10.3e\n", t);
              buf_p[j_loc] = t;
              sum += fabs (t);
            }
          p[j_bloc_col] += size_of_block;
        }
      if (tail)
        {
          double * buf_p = buf + p[number_of_blocks_row];
          for (int j_loc = 0; j_loc < size_of_tail; j_loc++)
            {
              double t = f (formula_number, size_of_system, i_global + 1, ++j_global);
              //printf ("2 - t = %10.3e\n", t);
              buf_p[j_loc] = t;
              sum += fabs (t);
            }
          p[number_of_blocks_row] += size_of_tail;
        }
      max = (sum > max ? sum : max);
      i_global++;
      //print_d_array_in_i_proc (buf, size_of_system, 0);
      //print_i_array_in_i_proc (p, number_of_blocks_row + tail, 0);
    }
  norm = max;
}

void LU_matrix::unpack_blocks_for_diag_and_U (int i_block_row)
{
  int count_of_block_for_proc = get_number_of_offdiag_block_in_i_row (i_block_row);
  double *recv_buf = buf_block_row_ptr_2.get ();
  int d_block = 0;
  int sq_size = size_of_block * size_of_block;

  //printf_i_process (PRINT_PROCESS, "START\n");
  if (proc_number == get_owner_of_big_part (i_block_row))
    {
      int needed_shift = i_block_row / number_of_proc;
      double *diag_block = diag_full_block_ptr.get () + needed_shift * sq_size;
      move_l2r (recv_buf, size_of_block, size_of_block, diag_block);
      d_block++;
    }
  double *U_off_diag_blocks = U_off_diag_full_block_ptr.get () + U_offsets_offdiag_block_ptr[i_block_row] * sq_size;
  double *blocks_for_move = recv_buf + d_block * sq_size;
  memcpy (U_off_diag_blocks, blocks_for_move, count_of_block_for_proc * sq_size * sizeof (double));

  if (respon_proc_for_tail)
    {
      double *U_tail_block = U_tail_bloc_ptr.get () + i_block_row * size_of_block * size_of_tail;
      double *block_for_move = recv_buf + (count_of_block_for_proc + d_block) * sq_size;
      move_l2r (block_for_move, size_of_block, size_of_tail, U_tail_block);
    }

  //printf_i_process (PRINT_PROCESS, "END\n");
}

void LU_matrix::fill_L_tail (int formula_number)
{
  double *L_tail = get_L_tail ();
  double *diag_tail = get_diag_tail ();
  int i_global = number_of_blocks_row * size_of_block;
  double max = norm;
  int sq_size = size_of_block * size_of_tail;
  //print_block_row_in_i_proc (get_L_tail (), size_of_tail, size_of_block, number_of_blocks_row, PRINT_PROCESS, 30);

  for (int i_loc = 0; i_loc < size_of_tail; i_loc++)
    {
      int j_global = 0;
      double sum = 0;
      double *L_tail_i = L_tail + i_loc * size_of_block;

      for (int j_bloc_col = 0; j_bloc_col < number_of_blocks_row; j_bloc_col++)
        {
          double * L_tail_ij = L_tail_i + sq_size * j_bloc_col;
          for (int j_loc = 0; j_loc < size_of_block; j_loc++)
            {
              //printf ("i_global = %d, j_global = %d, ", i_global, j_global);
              double t = f (formula_number, size_of_system, i_global + 1, ++j_global);
              //printf ("1 - t = %10.3e\n", t);
              L_tail_ij[j_loc] = t;
              sum += fabs (t);
            }
          //print_block_row_in_i_proc (get_L_tail (), size_of_tail, size_of_block, number_of_blocks_row, PRINT_PROCESS, 30);
        }
      double * diag_tail_i = diag_tail + i_loc * size_of_tail;
      for (int j_loc = 0; j_loc < size_of_tail; j_loc++)
        {
          //printf ("i_global = %d, j_global = %d, ", i_global, j_global);
          double t = f (formula_number, size_of_system, i_global + 1, ++j_global);
          //printf ("1 - t = %10.3e\n", t);
          diag_tail_i[j_loc] = t;
          sum += fabs (t);
        }
      max = (sum > max ? sum : max);
      i_global++;
    }
  norm = max;
}

void LU_matrix::fill_matrix (int formula_number)
{
  double *send_buf = buf_block_row_ptr_1.get ();
  double *recv_buf = buf_block_row_ptr_2.get ();
  int shared_parts_size = size_of_box_of_blocks * size_of_block * size_of_block;
  copy_offsets ();
  //print_i_array_in_i_proc (L_offsets_offdiag_block_ptr.get (), number_of_blocks_row + 1, PRINT_PROCESS, "L_array");

  //print_block_row_in_i_proc (get_L_tail (), size_of_tail, size_of_block, number_of_blocks_row, PRINT_PROCESS, 30);
  for (int i_bloc_row = 0; i_bloc_row < number_of_blocks_row; i_bloc_row++)
    {
      int size_of_first_blocks = i_bloc_row * size_of_block * size_of_block;
      int owner = get_owner_of_big_part (i_bloc_row);
      //printf_i_process (MAIN_PROCESS, "i_block_row = %d, owner = %d\n", i_bloc_row, owner);
      if (proc_number == MAIN_PROCESS)
        {
          fill_first_blocks_for_L (i_bloc_row, formula_number);
          //print_d_array_in_i_proc (send_buf, size_of_first_blocks, MAIN_PROCESS, "send_buf");

          // send
          if (owner == MAIN_PROCESS)
            {
              memcpy (recv_buf, send_buf, size_of_first_blocks * sizeof (double)); // may be better
            }
          else
            {
              //printf_i_process (owner, "0 tut?\n");
              MPI_Send (send_buf, size_of_first_blocks, MPI_DOUBLE, owner, 0 , comm);
            }
        }
      else
        {
          if (proc_number == owner)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_of_first_blocks, MPI_DOUBLE, MAIN_PROCESS, 0, comm, &st);
            }
        }

      //MPI_Barrier (comm);
      //printf_i_process (owner, "tut?\n");
      //print_d_array_in_i_proc (recv_buf, size_of_first_blocks, owner, "recv_buf");

      //printf_i_process (MAIN_PROCESS, "onwer = %d\n", owner);

      if (proc_number == owner)
        {
          // recv
          unpack_first_blocks_for_L (i_bloc_row);
        }


      //printf_i_process (PRINT_PROCESS, "==== IT = %d\n", i_bloc_row);
      //printf_i_process (PRINT_PROCESS, "size_of_box_of_blocks = %d\n", size_of_box_of_blocks);
      //for test
      memset (send_buf, 0, size_of_system * size_of_block * sizeof (double));
      if (proc_number == MAIN_PROCESS)
        {
          fill_blocks_for_diag_and_U (i_bloc_row, formula_number);
          //print_d_array_in_i_proc (buf_block_row_ptr.get (), size_of_block * size_of_system, 0);
          //print_block_row_in_i_proc (send_buf, size_of_block, size_of_block, size_of_box_of_blocks * number_of_proc, PRINT_PROCESS, 100);
        }


      MPI_Scatter (send_buf, shared_parts_size, MPI_DOUBLE, recv_buf, shared_parts_size, MPI_DOUBLE,
                   MAIN_PROCESS, comm);

      //print_block_row_in_i_proc (recv_buf, size_of_block, size_of_block, size_of_box_of_blocks, PRINT_PROCESS, 27);
      unpack_blocks_for_diag_and_U (i_bloc_row);


  }

  if (respon_proc_for_tail)
    {
      fill_L_tail (formula_number);
    }

/*
  printf_i_process (PRINT_PROCESS, "DIAG ");
  print_block_row_in_i_proc (get_diag_full_block (), size_of_block, size_of_block, size_of_box_of_blocks, PRINT_PROCESS, 30);

  printf_i_process (PRINT_PROCESS, "L full ");
  print_block_row_in_i_proc (get_L_off_diag (), size_of_block, size_of_block, get_number_of_offdiag_block (), PRINT_PROCESS, 30);
  if (respon_proc_for_tail)
    {
      printf_i_process (PRINT_PROCESS, "L tail");
      print_block_row_in_i_proc (get_L_tail (), size_of_tail, size_of_block, number_of_blocks_row, PRINT_PROCESS, 30);
    }

  printf_i_process (PRINT_PROCESS, "U offdiag ");
  print_block_row_in_i_proc (get_U_off_diag (), size_of_block, size_of_block, get_number_of_offdiag_block (), PRINT_PROCESS, 30);
  if (respon_proc_for_tail)
    {
      printf_i_process (PRINT_PROCESS, "U tail ");
      print_block_row_in_i_proc (get_U_tail (), size_of_block, size_of_tail, number_of_blocks_row, PRINT_PROCESS, 30);
      printf_i_process (PRINT_PROCESS, "diag tail ");
      print_block_row_in_i_proc (get_diag_tail (), size_of_tail, size_of_tail, 1, PRINT_PROCESS, 30);
    }
*/
}

int LU_matrix::read_first_blocks_for_L (int count_of_block, FILE *fp)
{
  double *buf = buf_block_row_ptr_1.get ();
  double max = norm;
  int *p = pointers_array_ptr.get ();
  double t;

  for (int i = 0; i < count_of_block; i++)
    {
      p[i]= size_of_block * size_of_block * i;
    }

  for (int i_loc = 0; i_loc < size_of_block; i_loc++)
    {
      double sum = 0;
      for (int j_bloc_col = 0; j_bloc_col < count_of_block; j_bloc_col++)
        {
          double *buf_p = buf + p[j_bloc_col];
          for (int j_loc = 0; j_loc < size_of_block; j_loc++)
            {
              if (fscanf (fp, "%lf", &t) != 1)
                {
                  return READ_ERR;
                }
              //printf ("1 - t = %10.3e\n", t);
              buf_p[j_loc] = t;
              sum += fabs (t);
            }
          p[j_bloc_col] += size_of_block;
        }
      max = (sum > max ? sum : max);
    }
  norm = max;
  return 0;
}

int LU_matrix::read_blocks_for_diag_and_U (int i_block_row, FILE *fp)
{
  place_buf_pointers (i_block_row);
  int *p = pointers_array_ptr.get ();
  double *buf = buf_block_row_ptr_1.get ();
  double t;
  double max = norm;

  for (int i_loc = 0; i_loc < size_of_block; i_loc++)
    {
      double sum = 0;
      for (int j_bloc_col = i_block_row; j_bloc_col < number_of_blocks_row; j_bloc_col++)
        {
          double * buf_p = buf + p[j_bloc_col];
          for (int j_loc = 0; j_loc < size_of_block; j_loc++)
            {
              if (fscanf (fp, "%lf", &t) != 1)
                {
                  return READ_ERR;
                }
              buf_p[j_loc] = t;
              sum += fabs (t);
            }
          p[j_bloc_col] += size_of_block;
        }
      if (tail)
        {
          double * buf_p = buf + p[number_of_blocks_row];
          for (int j_loc = 0; j_loc < size_of_tail; j_loc++)
            {
              if (fscanf (fp, "%lf", &t) != 1)
                {
                  return READ_ERR;
                }
              buf_p[j_loc] = t;
              sum += fabs (t);
            }
          p[number_of_blocks_row] += size_of_tail;
        }
      max = (sum > max ? sum : max);
    }
  norm = max;
  return 0;
}

int LU_matrix::read_L_tail (FILE *fp)
{
  double *L_tail = buf_block_row_ptr_1.get ();
  double *diag_tail = L_tail + size_of_tail * size_of_block * number_of_blocks_row;
  double max = norm;
  int sq_size = size_of_block * size_of_tail;
  double t;

  for (int i_loc = 0; i_loc < size_of_tail; i_loc++)
    {
      double sum = 0;
      double *L_tail_i = L_tail + i_loc * size_of_block;

      for (int j_bloc_col = 0; j_bloc_col < number_of_blocks_row; j_bloc_col++)
        {
          double * L_tail_ij = L_tail_i + sq_size * j_bloc_col;
          for (int j_loc = 0; j_loc < size_of_block; j_loc++)
            {
              if (fscanf (fp, "%lf", &t) != 1)
                {
                  return READ_ERR;
                }
              L_tail_ij[j_loc] = t;
              sum += fabs (t);
            }
        }
      double * diag_tail_i = diag_tail + i_loc * size_of_tail;
      for (int j_loc = 0; j_loc < size_of_tail; j_loc++)
        {
          if (fscanf (fp, "%lf", &t) != 1)
            {
              return READ_ERR;
            }
          diag_tail_i[j_loc] = t;
          sum += fabs (t);
        }
      max = (sum > max ? sum : max);
    }
  norm = max;
  return 0;
}

void LU_matrix::unpack_L_tail ()
{
  int size_L_tail = size_of_tail * size_of_block * number_of_blocks_row;
  double *buf_L_tail = buf_block_row_ptr_2.get ();
  double *buf_diag_tail = buf_block_row_ptr_2.get () + size_L_tail;
  double *L_tail = get_L_tail ();
  double *diag_tail = get_diag_tail ();

  memcpy (L_tail, buf_L_tail, size_L_tail * sizeof (double));
  memcpy (diag_tail, buf_diag_tail, size_of_tail * size_of_tail * sizeof (double));
}

int LU_matrix::read_matrix (char *file_name)
{
  double *send_buf = buf_block_row_ptr_1.get ();
  double *recv_buf = buf_block_row_ptr_2.get ();
  int shared_parts_size = size_of_box_of_blocks * size_of_block * size_of_block;

  copy_offsets ();
  int err = 0;

  FILE *fp = nullptr;

  if (proc_number == MAIN_PROCESS)
    {
      fp = fopen (file_name, "r");
      if (!fp)
        err = OPEN_FILE;
    }

  MPI_Bcast (&err, 1, MPI_INT, MAIN_PROCESS, comm);
  if (err < 0)
    return err;

  for (int i_bloc_row = 0; i_bloc_row < number_of_blocks_row; i_bloc_row++)
    {
      int size_of_first_blocks = i_bloc_row * size_of_block * size_of_block;
      int owner = get_owner_of_big_part (i_bloc_row);
      if (proc_number == MAIN_PROCESS)
        {
          err = read_first_blocks_for_L (i_bloc_row, fp);

          if (owner == MAIN_PROCESS)
            {
              memcpy (recv_buf, send_buf, size_of_first_blocks * sizeof (double)); // may be better
            }
          else
            {
              MPI_Send (send_buf, size_of_first_blocks, MPI_DOUBLE, owner, 0 , comm);
            }
        }
      else
        {
          if (proc_number == owner)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_of_first_blocks, MPI_DOUBLE, MAIN_PROCESS, 0, comm, &st);
            }
        }

      MPI_Bcast (&err, 1, MPI_INT, MAIN_PROCESS, comm);
      if (err < 0)
        {
          if (proc_number == MAIN_PROCESS)
            fclose (fp);
          return err;
        }

      if (proc_number == owner)
        {
          unpack_first_blocks_for_L (i_bloc_row);
        }

      if (proc_number == MAIN_PROCESS)
        {
          err = read_blocks_for_diag_and_U (i_bloc_row, fp);
        }

      MPI_Bcast (&err, 1, MPI_INT, MAIN_PROCESS, comm);
      if (err < 0)
        {
          if (proc_number == MAIN_PROCESS)
            fclose (fp);
          return err;
        }

      MPI_Scatter (send_buf, shared_parts_size, MPI_DOUBLE, recv_buf, shared_parts_size, MPI_DOUBLE,
                   MAIN_PROCESS, comm);

      unpack_blocks_for_diag_and_U (i_bloc_row);
  }

  if (tail)
    {
      int size_L_tail = size_of_tail * size_of_block * number_of_blocks_row + size_of_tail * size_of_tail;
      int owner = get_owner_of_big_part (number_of_blocks_row);
      if (proc_number == MAIN_PROCESS)
        {
          err = read_L_tail (fp);
          if (owner == MAIN_PROCESS)
            {
              memcpy (recv_buf, send_buf, size_L_tail * sizeof (double)); // may be better
            }
          else
            {
              MPI_Send (send_buf, size_L_tail, MPI_DOUBLE, owner, 0 , comm);
            }
        }
      else
        {
          if (proc_number == owner)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_L_tail, MPI_DOUBLE, MAIN_PROCESS, 0, comm, &st);
            }
        }
      MPI_Bcast (&err, 1, MPI_INT, MAIN_PROCESS, comm);
      if (err < 0)
        {
          if (proc_number == MAIN_PROCESS)
            fclose (fp);
          return err;
        }

      if (proc_number == owner)
        {
          unpack_L_tail ();
        }
    }

  if (proc_number == MAIN_PROCESS)
    fclose (fp);
  return 0;
}

int LU_matrix::init_matrix (int formula_number, char *file_name)
{
  int err = 0;
  if (formula_number)
    {
      fill_matrix (formula_number);
    }
  else
    {
      err = read_matrix (file_name);
    }
  MPI_Bcast (&norm, 1, MPI_DOUBLE, MAIN_PROCESS, comm);
  return err;
}

void LU_matrix::print_block_row_in_i_proc (double *A, int row_size, int col_size, int count_of_blocks, int i_proc,
                                           int max_print)
{
  if (proc_number == i_proc)
    {
      int j = 0;
      printf ("BLOCK PRINT TEST\n");
      int sq_size = row_size * col_size;
      row_size = (row_size > max_print ? max_print : row_size);
      for (int i_row = 0; i_row < row_size; i_row++)
        {
          j = 0;
          double *Ai = A + col_size * i_row;
          for (int i_block_col = 0; i_block_col < count_of_blocks && j < max_print; i_block_col++)
            {
              double *Ai_bloc = Ai + i_block_col * sq_size;
              for (int j_loc = 0; j_loc < col_size && j < max_print; j_loc++)
                {
                  printf ("%6.3lf ", Ai_bloc[j_loc]);
                  //printf ("%10.3e ", Ai_bloc[j_loc]);
                  j++;
                }
              printf (" ");
            }
          printf ("\n");
        }
      printf ("END TEST\n");
    }
}

void LU_matrix::collet_first_L_and_diag_blocks (int i_block_row)
{
  double *L = L_off_diag_full_block_ptr.get ();
  double *send_buf = buf_block_row_ptr_1.get ();
  int *L_offsets = L_offsets_offdiag_block_ptr.get ();
  int sq_size = size_of_block * size_of_block;

  for (int i_block = 0; i_block < i_block_row; i_block++)
    {
      move_l2r (L + L_offsets[i_block] * sq_size, size_of_block, size_of_block, send_buf + i_block * sq_size);
      L_offsets[i_block]++;
    }
  move_l2r (get_diag_full_block (i_block_row), size_of_block, size_of_block, send_buf + i_block_row * sq_size);
}

void LU_matrix::print_block_i_row (double *block_row, int i_block_row, int row_size, int max_print, FILE *fp)
{
  if (proc_number == MAIN_PROCESS)
    {
      int j = 0;
      int sq_size = size_of_block * row_size;
      row_size = (row_size > max_print - i_block_row * size_of_block ? max_print - i_block_row * size_of_block : row_size);
      for (int i_row = 0; i_row < row_size; i_row++)
        {
          j = 0;
          double *Ai = block_row + size_of_block * i_row;
          for (int i_block_col = 0; i_block_col < number_of_blocks_row && j < max_print; i_block_col++)
            {
              double *Ai_bloc = Ai + i_block_col * sq_size;
              for (int j_loc = 0; j_loc < size_of_block && j < max_print; j_loc++)
                {
                  //fprintf (fp, "%8.3lf ", Ai_bloc[j_loc]);
                  //fprintf (fp, "%3.lf ", Ai_bloc[j_loc]);
                  fprintf (fp, "%10.3e", Ai_bloc[j_loc]);
                  j++;
                }
              fprintf (fp,"   ");
            }
          if (tail)
            {
              double *Ai_bloc = block_row + number_of_blocks_row * sq_size + i_row * size_of_tail;
              for (int j_loc = 0; j_loc < size_of_tail && j < max_print; j_loc++)
                {
                  //fprintf (fp, "%8.3lf ", Ai_bloc[j_loc]);
                  //fprintf (fp, "%3.lf ", Ai_bloc[j_loc]);
                  fprintf (fp, "%10.3e", Ai_bloc[j_loc]);
                  j++;
                }
            }
          fprintf (fp, "\n");
        }
      fprintf (fp, "\n");
    }
}

int LU_matrix::next_proc (int cur_proc)
{
  return (cur_proc + 1 == number_of_proc ? 0 : cur_proc + 1); //may be it is equal (cur_proc + 1) % number_of_proc
}

void LU_matrix::print_matrix (const char *text, int max_print, const char *file_name)
{
  copy_offsets ();
  double *send_buf = buf_block_row_ptr_1.get ();
  double *recv_buf = buf_block_row_ptr_2.get ();
  int sq_size = size_of_block * size_of_block;

  FILE *fp = stdout;
  if (proc_number == MAIN_PROCESS)
    {
      if (file_name)
        {
          printf ("PRINT MATRIX IN %s\n", file_name);
          fp = fopen (file_name, "w");
          if (!fp)
            return;
        }
      else
        {
          printf_main_process ("%s: \n", text);
        }
    }



  for (int i_block_row = 0; i_block_row < number_of_blocks_row && (i_block_row * size_of_block) < max_print; i_block_row++)
    {
      /// only for test
      //memset (send_buf, 0, sizeof (double) * size_of_box_of_blocks * number_of_proc * size_of_block * size_of_block);
      //memset (recv_buf, 0, sizeof (double) * size_of_box_of_blocks * number_of_proc * size_of_block * size_of_block);
      int owner = get_owner_of_big_part (i_block_row);

      int size_of_first_part = sq_size * (i_block_row + 1);

      if (proc_number == owner)
        {
          collet_first_L_and_diag_blocks (i_block_row);
          if (proc_number == MAIN_PROCESS)
            {
              memcpy (recv_buf, send_buf, sizeof (double) * size_of_first_part);
            }
          else
            {
              MPI_Send (send_buf, size_of_first_part, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
            }
        }
      else
        {
          if (proc_number == MAIN_PROCESS)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_of_first_part, MPI_DOUBLE, owner, 0 /* tag */, comm, &st);
            }
        }

      int send_proc = owner;
      int shift_of_block = 0;

      for (int i_block = i_block_row + 1; i_block < number_of_blocks_row && (i_block * size_of_block < max_print); i_block++)
        {
          double *U_i = get_U_off_diag (i_block_row) + shift_of_block * sq_size;
          double *recv_buf_i = recv_buf + i_block * sq_size;
          send_proc = next_proc (send_proc);
          if (proc_number == send_proc)
            {
              if (proc_number == MAIN_PROCESS)
                {
                  memcpy (recv_buf_i, U_i, sizeof (double) * sq_size);
                }
              else
                {
                  MPI_Send (U_i, sq_size, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
                }
              shift_of_block++;
            }
          else
            {
              if (proc_number == MAIN_PROCESS)
                {
                  MPI_Status st;
                  MPI_Recv (recv_buf_i, sq_size, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
                }
            }
        }
      if (tail && (number_of_blocks_row * size_of_block) < max_print)
        {
          double *U_tail = get_U_tail (i_block_row);
          double *recv_buf_i = recv_buf + number_of_blocks_row * sq_size;
          send_proc = next_proc (send_proc);
          if (proc_number == send_proc)
            {
              if (proc_number == MAIN_PROCESS)
                {
                  memcpy (recv_buf_i, U_tail, sizeof (double) * size_of_block * size_of_tail);
                }
              else
                {
                  MPI_Send (U_tail, size_of_block * size_of_tail, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
                }
            }
          else
            {
              if (proc_number == MAIN_PROCESS)
                {
                  MPI_Status st;
                  MPI_Recv (recv_buf_i, size_of_block * size_of_tail, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
                }
            }
        }

      print_block_i_row (recv_buf, i_block_row, size_of_block, max_print, fp);
    }

  if (tail && (number_of_blocks_row * size_of_block) < max_print)
    {
      double *L_tail = get_L_tail ();
      int send_proc = get_owner_of_big_part (number_of_blocks_row);

      if (respon_proc_for_tail)
        {

          if (proc_number == MAIN_PROCESS)
            {
              memcpy (recv_buf, L_tail, sizeof (double) * size_of_block * size_of_tail * number_of_blocks_row);
            }
          else
            {
              MPI_Send (L_tail, size_of_block * size_of_tail * number_of_blocks_row, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
            }
        }
      else
        {
          if (proc_number == MAIN_PROCESS)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_of_block * size_of_tail * number_of_blocks_row, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
            }
        }

      double * diag_tail = get_diag_tail ();
      double * recv_tail = recv_buf + size_of_block * size_of_tail * number_of_blocks_row;

      if (respon_proc_for_tail)
        {
          if (proc_number == MAIN_PROCESS)
            {
              memcpy (recv_tail, diag_tail, sizeof (double) * size_of_tail * size_of_tail);
            }
          else
            {
              MPI_Send (diag_tail, size_of_tail * size_of_tail, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
            }
        }
      else
        {
          if (proc_number == MAIN_PROCESS)
            {
              MPI_Status st;
              MPI_Recv (recv_tail, size_of_tail * size_of_tail, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
            }
        }
    }
  print_block_i_row (recv_buf, number_of_blocks_row, size_of_tail, max_print, fp);
  if (proc_number == MAIN_PROCESS)
    {
      if (file_name)
        {
          fclose (fp);
        }
    }
}

void LU_matrix::sum_for_vector_B (double *block, double *vector,int size_of_row, int size_of_col, int start)
{
  for (int i_row = 0; i_row < size_of_row; i_row++)
    {
      double sum = 0;
      double *A_i = block + i_row * size_of_col;
      for (int j_col = start; j_col < size_of_col; j_col += 2)
        {
          //printf_i_process (1, "A_i[%d] = %10.3e\n", j_col, A_i[j_col]);
          sum += A_i[j_col];
        }
      //printf_i_process (1, "i = %d, sum = %10.3e\n", i_row, sum);
      vector[i_row] += sum;
    }
}

void LU_matrix::work_with_L (int i_block_row, double *B)
{
  double *L = get_L_off_diag ();
  int *L_offsets = L_offsets_offdiag_block_ptr.get ();
  for (int i_block = 0; i_block < i_block_row; i_block++)
    {
      double *L_i_block = L + L_offsets[i_block] * size_of_block * size_of_block;
      sum_for_vector_B (L_i_block, B, size_of_block, size_of_block, (i_block * size_of_block) & 1);
      L_offsets[i_block]++;
    }
}

int LU_matrix::get_global_number_of_block_column (int i_block_row)
{
  int m = proc_number;
  int k = get_owner_of_big_part (i_block_row);
  int p = number_of_proc;
  return (m <= k ? p - k + m : m - k);
}

void LU_matrix::work_with_U (int i_block_row, double *B)
{
  int count_of_block = get_number_of_offdiag_block_in_i_row (i_block_row);
  double *U = get_U_off_diag (i_block_row);
  int st = (get_global_number_of_block_column (i_block_row) + i_block_row);
  for (int i_block = 0; i_block < count_of_block; i_block++)
    {
      double *U_i = U + i_block * size_of_block * size_of_block;
      int start = (st + number_of_proc * i_block) * size_of_block;
      sum_for_vector_B (U_i, B, size_of_block, size_of_block, start & 1);
    }
}

void LU_matrix::work_with_L_tail (double *B)
{
  double *L_tail = get_L_tail ();
  double *diag_tail = get_diag_tail ();

  for (int i_block = 0; i_block < number_of_blocks_row; i_block++)
    {
      double *L_tail_i = L_tail + i_block * size_of_tail * size_of_block;
      sum_for_vector_B (L_tail_i, B, size_of_tail, size_of_block, (i_block * size_of_block) & 1);
    }
  sum_for_vector_B (diag_tail, B, size_of_tail, size_of_tail, (number_of_blocks_row * size_of_block) & 1);
}

void LU_matrix::build_rhs (double *B, double *X)
{
  memset (X, 0, sizeof (double) * size_of_system);
  memset (B, 0, sizeof (double) * size_of_system);
  copy_offsets ();

  //print_i_array_in_i_proc (L_offsets_offdiag_block_ptr.get (), number_of_blocks_row + 1, MAIN_PROCESS);

  for (int i_block_row = 0; i_block_row < number_of_blocks_row; i_block_row++)
    {
      double *B_i = X + i_block_row * size_of_block;
      //work with L's blocks and with diagonal block
      int owner = get_owner_of_big_part (i_block_row);
      if (proc_number == owner)
        {
          work_with_L (i_block_row, B_i);
          //print_d_array_in_i_proc (B_i, 3, MAIN_PROCESS);
          sum_for_vector_B (get_diag_full_block (i_block_row), B_i, size_of_block, size_of_block,
                            1 & (i_block_row * size_of_block));
          //print_d_array_in_i_proc (B_i, 3, MAIN_PROCESS);
        }
      //work with U's blocks
      work_with_U (i_block_row, B_i);
      if (respon_proc_for_tail)
        {
          sum_for_vector_B (get_U_tail (i_block_row), B_i, size_of_block, size_of_tail,
                            1 & (number_of_blocks_row * size_of_block));
        }
    }
  if (respon_proc_for_tail)
    {
      work_with_L_tail (X + number_of_blocks_row * size_of_block);
    }

  MPI_Reduce (X, B, size_of_system, MPI_DOUBLE, MPI_SUM, MAIN_PROCESS, comm);
}


int LU_matrix::make_prepare_mult_for_diag_block (int i_block_row, int start_block)
{
  if (!i_block_row)
    return 0;

  //printf_i_process (PRINT_PROCESS, "i_block_row = %d, start_block = %d\n", i_block_row, start_block);
  double *result = get_prepare_buf () + ((i_block_row / number_of_proc) - 1
                                         + (proc_number ? 1 : 0)) * mm_size;
  int number_of_block_in_row = get_number_of_offdiag_block_in_i_row (i_block_row) + (!proc_number ? 1 : 0);
  double *buf_block = get_block_1 ();

  for (int i_row = start_block; i_row < i_block_row; i_row++)
    {
      int shift = (proc_number > get_owner_of_big_part (i_row) ? 1 : 0);
      double *L_i = get_L_off_diag (i_row) + shift * mm_size;
      double *U_i = get_U_off_diag (i_row) + shift * mm_size;
      for (int i_block = 0; i_block < number_of_block_in_row; i_block++)
        {
          int i_shift = i_block * mm_size;
          //printf_i_process (PRINT_PROCESS, "i_block = %d, number_of_block_in_row = %d\n", i_block, number_of_block_in_row);
          //printf_i_process (PRINT_PROCESS, "L_block");
          //print_block_row_in_i_proc (L_i + i_shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
          //printf_i_process (PRINT_PROCESS, "U_block");
          //print_block_row_in_i_proc (U_i + i_shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
          mult_blocks_block (L_i + i_shift, U_i + i_shift, buf_block, size_of_block, size_of_block, size_of_block);
          add_blocks (buf_block, result + i_shift, result + i_shift, size_of_block, size_of_block);
        }
    }
  return i_block_row;
}

int LU_matrix::mult_and_factor_block (int i_block_row, int start_block)
{
  double *buf = get_block_1 ();
  double *prev_res_of_mult = get_prepare_buf () + (i_block_row / number_of_proc + (proc_number ? 0 : -1)) * mm_size;
  double *diag_block = get_diag_full_block (i_block_row);
  for (int i_block = start_block; i_block < i_block_row; i_block++)
    {
      double *L_i = get_L_off_diag (i_block);
      double *U_i = get_U_off_diag (i_block);
      mult_blocks_block (L_i, U_i, buf, size_of_block, size_of_block, size_of_block);
      add_blocks (buf, prev_res_of_mult, prev_res_of_mult, size_of_block, size_of_block);
    }
  subst_blocks (diag_block, prev_res_of_mult, diag_block, size_of_block, size_of_block);
  return factor_block (diag_block, size_of_block);
}

void LU_matrix::fill_column_for_send (int i_block_row)
{
  double *send_buf = buf_block_row_ptr_1.get ();
  double *U = get_U_off_diag ();
  int *offset = U_offsets_offdiag_block_ptr.get ();

  int shift = i_block_row / number_of_proc;

  //printf_i_process (PRINT_PROCESS, "Start with i_bloc_row = %d\n", i_block_row);
  for (int i_block = 0; i_block < i_block_row; i_block++)
    {
      if (i_block % number_of_proc == proc_number)
        shift--;
      //printf_i_process (PRINT_PROCESS, "shift = %d\n", shift);
      memcpy (send_buf + i_block * mm_size, U + (offset[i_block] + shift) * mm_size, sizeof (double) * mm_size);
    }
  memcpy (send_buf + i_block_row * mm_size, get_rev_LU_block (), sizeof (double) * mm_size);
}

void LU_matrix::find_L_column (int i_block_row)
{
  int *shifted_offset = L_offsets_offdiag_block_ptr.get ();
  int number_of_rows = U_offsets_offdiag_block_ptr[i_block_row + 1] - shifted_offset[i_block_row];


  double *L = get_L_off_diag ();
  double *U_col = buf_block_row_ptr_1.get();
  double *buf_row = buf_block_row_ptr_2.get ();
  double *buf_block = get_block_1 ();
  double *U_r = U_col + i_block_row * mm_size;

  memset (buf_row, 0, sizeof (double) * number_of_rows * mm_size);
//  print_block_row_in_i_proc (buf_row,  size_of_block, size_of_block, number_of_rows, PRINT_PROCESS, 30);
//  print_i_array_in_i_proc (shifted_offset, number_of_blocks_row + 1, PRINT_PROCESS, "shifted offsets");
  for (int j_col = 0; j_col < i_block_row; j_col++)
    {
      double *L_col = L + shifted_offset[j_col] * mm_size;
      double *U_row = U_col + j_col * mm_size;
      for (int i_row = 0; i_row < number_of_rows; i_row++)
        {
          int shift = i_row * mm_size;
/*
          printf_i_process (PRINT_PROCESS, "L_col + shift");
          print_block_row_in_i_proc (L_col + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
          printf_i_process (PRINT_PROCESS, "U_row");
          print_block_row_in_i_proc (U_row, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
*/
          mult_blocks_block (L_col + shift, U_row, buf_block, size_of_block, size_of_block, size_of_block);
          add_blocks (buf_block, buf_row + shift, buf_row + shift, size_of_block, size_of_block);
        }
    }

  double *L_col = get_L_off_diag (i_block_row);
/*
  printf_i_process (PRINT_PROCESS, "i_block = %d: number_of_rows = %d\n", i_block_row, number_of_rows);
  print_block_row_in_i_proc (buf_row,  size_of_block, size_of_block, number_of_rows, PRINT_PROCESS, 30);
*/
  for (int i_row = 0; i_row < number_of_rows; i_row++)
    {
      int shift = i_row * mm_size;
/*
      printf_i_process (PRINT_PROCESS, "L_block");
      print_block_row_in_i_proc (L_col + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
      printf_i_process (PRINT_PROCESS, "buf_row + shift");
      print_block_row_in_i_proc (buf_row + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
      printf_i_process (PRINT_PROCESS, "RES_block");
      print_block_row_in_i_proc (get_block_1 (), size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
*/


      if (i_block_row)
        {
          /*
          printf_i_process (PRINT_PROCESS, "L_col + shift");
          print_block_row_in_i_proc (L_col + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
          printf_i_process (PRINT_PROCESS, "buf_row + shift");
          print_block_row_in_i_proc (buf_row + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
          */
          subst_blocks (L_col + shift, buf_row + shift, L_col + shift, size_of_block, size_of_block);
        }
/*
      printf_i_process (PRINT_PROCESS, "L_col + shift");
      print_block_row_in_i_proc (L_col + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
      printf_i_process (PRINT_PROCESS, "U_r");
      print_block_row_in_i_proc (U_r, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
*/
      mult_block_U (L_col + shift, U_r, buf_block, size_of_block, size_of_block);
      move_l2r (buf_block, size_of_block, size_of_block, L_col + shift);
    }
}

void LU_matrix::fill_row_for_send (int i_block_row)
{
  double *send_buf = buf_block_row_ptr_1.get ();
  double *L = get_L_off_diag ();
  int *offset = U_offsets_offdiag_block_ptr.get ();

  int shift = i_block_row / number_of_proc;

  //printf_i_process (PRINT_PROCESS, "Start with i_bloc_row = %d\n", i_block_row);
  for (int i_block = 0; i_block < i_block_row; i_block++)
    {
      if (i_block % number_of_proc == proc_number)
        shift--;
      //printf_i_process (PRINT_PROCESS, "shift = %d\n", shift);
      memcpy (send_buf + i_block * mm_size, L + (offset[i_block] + shift) * mm_size, sizeof (double) * mm_size);
    }
}

void LU_matrix::find_U_row (int i_block_row)
{
  int *shifted_offset = L_offsets_offdiag_block_ptr.get ();
  int number_of_cols = U_offsets_offdiag_block_ptr[i_block_row + 1] - shifted_offset[i_block_row];

  double *U = get_U_off_diag ();
  double *L_row = buf_block_row_ptr_1.get();
  double *buf_row = buf_block_row_ptr_2.get ();
  double *buf_block = get_block_1 ();
  double *L_r = L_row + i_block_row * mm_size;

  memset (buf_row, 0, sizeof (double) * number_of_cols * mm_size);
  for (int i_row = 0; i_row < i_block_row; i_row++)
    {
      double *U_row = U + shifted_offset[i_row] * mm_size;
      double *L_col = L_row + i_row * mm_size;
      for (int j_col = 0; j_col < number_of_cols; j_col++)
        {
          int shift = j_col * mm_size;
          mult_blocks_block (L_col, U_row + shift, buf_block, size_of_block, size_of_block, size_of_block);
          add_blocks (buf_block, buf_row + shift, buf_row + shift, size_of_block, size_of_block);
        }
    }

  double *U_row = get_U_off_diag (i_block_row);

  //printf_i_process (PRINT_PROCESS, "i_block = %d: number_of_rows = %d\n", i_block_row, number_of_rows);

  for (int i_col = 0; i_col < number_of_cols; i_col++)
    {
      int shift = i_col * mm_size;
/*
      printf_i_process (PRINT_PROCESS, "L_r");
      print_block_row_in_i_proc (L_r, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
      printf_i_process (PRINT_PROCESS, "U_row + shift");
      print_block_row_in_i_proc (U_row + shift, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
*/

      if (i_block_row)
        {
          subst_blocks (U_row + shift, buf_row + shift, U_row + shift, size_of_block, size_of_block);
        }

      mult_block_L (L_r, U_row + shift, buf_block, size_of_block, size_of_block);
      move_l2r (buf_block, size_of_block, size_of_block, U_row + shift);
    }
  if (proc_number == get_owner_of_big_part (i_block_row + 1))
    {
      for (int i_block = 0; i_block <= i_block_row; i_block++)
        {
          shifted_offset[i_block]++;
        }
    }
}

void LU_matrix::find_L_column_tail (int i_block_row)
{
  double *L_tail = get_L_tail ();
  double *U_col = buf_block_row_ptr_1.get ();
  double *buf_block_1 = get_block_1 ();
  double *buf_block_2 = get_block_2 ();
  double *U_r = buf_block_row_ptr_1.get() + i_block_row * mm_size;

  memset (buf_block_1, 0, sizeof (double) * ml_size);
  for (int i_block = 0; i_block < i_block_row; i_block++)
    {
      mult_blocks (L_tail + i_block * ml_size, U_col + i_block * mm_size, buf_block_2, size_of_tail, size_of_block, size_of_block);
      add_blocks (buf_block_2, buf_block_1, buf_block_1, size_of_tail, size_of_block);
    }
/*
  printf_i_process (PRINT_PROCESS, "L_tail + i_block_row * ml_size");
  print_block_row_in_i_proc (L_tail + i_block_row * ml_size, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
  printf_i_process (PRINT_PROCESS, "U_r");
  print_block_row_in_i_proc (U_r, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
*/
  subst_blocks (L_tail + i_block_row * ml_size, buf_block_1, L_tail + i_block_row * ml_size, size_of_tail, size_of_block);
  mult_block_U (L_tail + i_block_row * ml_size, U_r, buf_block_2, size_of_tail, size_of_block);
  move_l2r (buf_block_2, size_of_tail, size_of_block, L_tail + i_block_row * ml_size);
}

void LU_matrix::find_U_row_tail (int i_block_row)
{
  double *U_tail = get_U_tail ();
  double *L_row = buf_block_row_ptr_1.get ();
  double *buf_block_1 = get_block_1 ();
  double *buf_block_2 = get_block_2 ();
  double *L_r = buf_block_row_ptr_1.get() + i_block_row * mm_size;

  memset (buf_block_1, 0, sizeof (double) * ml_size);
  for (int i_block = 0; i_block < i_block_row; i_block++)
    {
      mult_blocks (L_row + i_block * mm_size, U_tail + i_block * ml_size, buf_block_2, size_of_block, size_of_block, size_of_tail);
      add_blocks (buf_block_2, buf_block_1, buf_block_1, size_of_block, size_of_tail);
    }
/*
  printf_i_process (PRINT_PROCESS, "U_tail + i_block_row * ml_size");
  print_block_row_in_i_proc (U_tail + i_block_row * ml_size, size_of_block, size_of_tail, 1, PRINT_PROCESS, 30);
  printf_i_process (PRINT_PROCESS, "L_r");
  print_block_row_in_i_proc (L_r, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
*/
  subst_blocks (U_tail + i_block_row * ml_size, buf_block_1, U_tail + i_block_row * ml_size, size_of_block, size_of_tail);
  mult_block_L (L_r, U_tail + i_block_row * ml_size, buf_block_2, size_of_block, size_of_tail);
  move_l2r (buf_block_2, size_of_tail, size_of_block, U_tail + i_block_row * ml_size);
}

int LU_matrix::find_last_diag_tail_block ()
{
  double *L_tail = get_L_tail ();
  double *U_tail = get_U_tail ();
  double *buf_block_1 = get_block_1 ();
  double *buf_block_2 = get_block_2 ();
  double *diag_block = get_diag_tail ();

  memset (buf_block_2, 0, sizeof (double) * ll_size);
  for (int i_block = 0; i_block < number_of_blocks_row; i_block++)
    {
      int shift = i_block * ml_size;
      mult_blocks (L_tail + shift, U_tail + shift, buf_block_1, size_of_tail, size_of_block, size_of_tail);
      add_blocks (buf_block_1, buf_block_2, buf_block_2, size_of_tail, size_of_tail);
    }
  subst_blocks (diag_block, buf_block_2, diag_block, size_of_tail, size_of_tail);
  return factor_block (diag_block, size_of_tail);
}

int LU_matrix::make_factorization ()
{
  int start_block = 0;
  int err = 0;
  double *send_buf = buf_block_row_ptr_1.get ();
  copy_offsets ();
  //double *recv_buf = buf_block_row_ptr_2.get ();



  /*
  mult_blocks_block (get_L_off_diag (), get_U_off_diag (), get_block_1 (), size_of_block, size_of_block, size_of_block);
  printf_i_process (PRINT_PROCESS, "L_block");
  print_block_row_in_i_proc (get_L_off_diag (), size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
  printf_i_process (PRINT_PROCESS, "U_block");
  print_block_row_in_i_proc (get_U_off_diag (), size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
  printf_i_process (PRINT_PROCESS, "RES_block");
  print_block_row_in_i_proc (get_block_1 (), size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
  return -1;
  */

  for (int i_block_row = 0; i_block_row < number_of_blocks_row; i_block_row++)
    {
      //printf_i_process (PRINT_PROCESS, "====== IT = %d =======\n", i_block_row);
      // only for testing
      //memset (send_buf, 0, sizeof (double) * number_of_blocks_row * mm_size);
      // some speed thing
      if (!(i_block_row % number_of_proc))
        {
          start_block = make_prepare_mult_for_diag_block (i_block_row, start_block);
        }
      // find diagonal
      int owner = get_owner_of_big_part (i_block_row);
      if (proc_number == owner)
        {
          if (i_block_row)
            {
              err = mult_and_factor_block (i_block_row, start_block);
            }
          else
            {
              err = factor_block (get_diag_full_block (), size_of_block);
            }
          if (err >= 0)
            {
              err = find_Lr (get_diag_full_block (i_block_row), get_rev_LU_block (), size_of_block);
            }
          if (err >= 0)
            {
              find_Ur (get_diag_full_block (i_block_row), get_rev_LU_block (), size_of_block);
            }
        }
      MPI_Bcast (&err, 1, MPI_INT, MAIN_PROCESS, comm);
      if (err < 0)
        {
          return err;
        }

      // find column
      if (proc_number == owner)
        {
          fill_column_for_send (i_block_row);
        }
      MPI_Bcast (send_buf, (i_block_row + 1) * mm_size, MPI_DOUBLE, owner, comm);
      //print_block_row_in_i_proc (send_buf, size_of_block, size_of_block, number_of_blocks_row, PRINT_PROCESS);

      find_L_column (i_block_row);

      if (respon_proc_for_tail)
        {
          find_L_column_tail (i_block_row);
        }

      // find row
      if (proc_number == owner)
        {
          fill_row_for_send (i_block_row);
        }
      MPI_Bcast (send_buf, i_block_row * mm_size, MPI_DOUBLE, owner, comm);

      //print_block_row_in_i_proc (send_buf, size_of_block, size_of_block, number_of_blocks_row, PRINT_PROCESS);

      find_U_row (i_block_row);
      if (respon_proc_for_tail)
        {
          find_U_row_tail (i_block_row);
        }
    }
  int proc = get_owner_of_big_part (number_of_blocks_row);
  if (respon_proc_for_tail)
    {
      err = find_last_diag_tail_block ();
    }
  MPI_Bcast (&err, 1, MPI_INT, proc, comm);
  if (err < 0)
    {
      return err;
    }

  return 0;
}

int LU_matrix::find_block_y (int i_block_row, double *Y, double *B)
{
  double *a = get_diag_full_block (i_block_row);

  for (int i_row = 0; i_row < size_of_block; i_row++)
    {
      double aii = a[i_row * size_of_block + i_row];
      if (fabs(aii) < norm * eps)
        return FIND_X;
      double t = 0;
      double *li = a + i_row * size_of_block;
      for (int j_col = 0; j_col < i_row; j_col++)
        {
          t += li[j_col] * Y[j_col];
        }
      Y[i_row] = (B[i_row] - t) / aii;
    }
  return 0;
}

void LU_matrix::mult_LX (double *L, double *X, double *B, int n_row, int n_col)
{
  // B - LX
  for (int i_row = 0; i_row < n_row; i_row++)
    {
      double t = 0;
      double *Li = L + i_row * n_col;
      for (int j_col = 0; j_col < n_col; j_col++)
        {
          t += Li[j_col] * X[j_col];
        }
      B[i_row] -= t;
    }
}

void LU_matrix::L_update_b (int i_block_row, double *X, double *B)
{
  if (i_block_row == number_of_blocks_row - 1)
    return;
  double *L_i = get_L_off_diag (i_block_row);
  int number_of_block_in_col = get_number_of_offdiag_block_in_i_row (i_block_row);
  //printf ("number = %d\n", number_of_block_in_col);
  for (int i_block = 0; i_block < number_of_block_in_col; i_block++)
    {
      double *B_i = B + size_of_block * (i_block_row + get_global_number_of_block_column (i_block_row) + i_block * number_of_proc);
      /*
      printf_i_process (PRINT_PROCESS, "B_i before, shift = %d\n", (i_block_row + get_global_number_of_block_column (i_block_row) + i_block * number_of_proc));
      print_d_array_in_i_proc (B_i, size_of_block, PRINT_PROCESS);
      print_block_row_in_i_proc (L_i + i_block * mm_size, size_of_block, size_of_block, 1, PRINT_PROCESS, 30);
      */
      mult_LX (L_i + i_block * mm_size, X, B_i, size_of_block, size_of_block);
      /*
      printf_i_process (PRINT_PROCESS, "B_i after\n");
      print_d_array_in_i_proc (B_i, size_of_block, PRINT_PROCESS);
      */
    }
}


void LU_matrix::find_block_x (int i_block_row, double *X, double *B)
{
  double *a = get_diag_full_block (i_block_row);

  //print_block_row_in_i_proc (a, size_of_block, size_of_block, 1, PRINT_PROCESS);
  //print_d_array_in_i_proc (B, size_of_block, PRINT_PROCESS);

  for (int i_row = size_of_block - 1; i_row >= 0; i_row--)
    {
      double t = 0;
      double *ui = a + i_row * size_of_block;
      for (int j_col = size_of_block - 1; j_col > i_row; j_col--)
        {
          t += ui[j_col] * X[j_col];
        }
      X[i_row] = B[i_row] - t;
    }
  //print_d_array_in_i_proc (X, size_of_block, PRINT_PROCESS);

}

void LU_matrix::U_update_b (int i_block_row, double *X, double *B)
{
  if (!i_block_row)
    return;

  int *shifted_offset = L_offsets_offdiag_block_ptr.get ();


  for (int i_block = i_block_row - 1; i_block >= 0; i_block--)
    {
      double *U = get_U_off_diag () + (shifted_offset[i_block + 1] - 1) * mm_size;
      mult_LX (U, X, B + i_block * size_of_block, size_of_block, size_of_block);
      shifted_offset[i_block + 1]--;
    }
}

int LU_matrix::calc_tail_L (double *X, double *B)
{
  double *L_tail = get_L_tail ();
  double *B_tail = B + number_of_blocks_row * size_of_block;
  for (int i_block = 0; i_block < number_of_blocks_row; i_block++)
    {
      mult_LX (L_tail + i_block * ml_size, X + i_block * size_of_block, B_tail, size_of_tail, size_of_block);
    }

  double *a = get_diag_tail ();
  double *Y_tail = X + number_of_blocks_row * size_of_block;

  for (int i_row = 0; i_row < size_of_tail; i_row++)
    {
      double aii = a[i_row * size_of_tail + i_row];
      if (fabs(aii) < norm * eps)
        return FIND_X;
      double t = 0;
      double *li = a + i_row * size_of_tail;
      for (int j_col = 0; j_col < i_row; j_col++)
        {
          t += li[j_col] * Y_tail[j_col];
        }
      Y_tail[i_row] = (B_tail[i_row] - t) / aii;
    }
  return 0;
}

void LU_matrix::calc_tail_U (double *X, double *B)
{
  double *a = get_diag_tail ();

  //print_block_row_in_i_proc (a, size_of_tail, size_of_tail, 1, PRINT_PROCESS);
  //print_d_array_in_i_proc (B, size_of_system, PRINT_PROCESS);
  //print_d_array_in_i_proc (X, size_of_system, PRINT_PROCESS);

  double *B_tail = B + number_of_blocks_row * size_of_block;
  double *X_tail = X + number_of_blocks_row * size_of_block;

  for (int i_row = size_of_tail - 1; i_row >= 0; i_row--)
    {
      double t = 0;
      double *ui = a + i_row * size_of_tail;
      for (int j_col = size_of_tail - 1; j_col > i_row; j_col--)
        {
          t += ui[j_col] * X_tail[j_col];
        }
      X_tail[i_row] = B_tail[i_row] - t;
    }
  //print_d_array_in_i_proc (X, size_of_system, PRINT_PROCESS);

  for (int i_block = number_of_blocks_row - 1; i_block >= 0; i_block--)
    {
      double *U = get_U_tail () + i_block * ml_size;
      mult_LX (U, X_tail, B + i_block * size_of_block, size_of_block, size_of_tail);
    }
  //print_d_array_in_i_proc (B, size_of_system, PRINT_PROCESS);
  //print_d_array_in_i_proc (X, size_of_system, PRINT_PROCESS);
}

int LU_matrix::find_x (double *B, double *X)
{
  int err = 0;

  //printf ("syze_of_system = %d\n", size_of_system);
  MPI_Bcast (B, size_of_system, MPI_DOUBLE, MAIN_PROCESS, comm);

  // Ly = b
  for (int i_block = 0; i_block < number_of_blocks_row; i_block++)
    {
      int owner = get_owner_of_big_part (i_block);
      int shift = i_block * size_of_block;
      if (proc_number == owner)
        {
          err = find_block_y (i_block, X + shift, B + shift);
        }
      MPI_Bcast (&err, 1, MPI_INT, owner, comm);
      if (err < 0)
        {
          return err;
        }
      MPI_Bcast (X + shift, size_of_block, MPI_DOUBLE, owner, comm);
      L_update_b (i_block, X + shift, B);
    }


  int proc = get_owner_of_big_part (number_of_blocks_row);


  if (respon_proc_for_tail)
    {
      err = calc_tail_L (X, B);
    }

  MPI_Bcast (&err, 1, MPI_INT, proc, comm);
  if (err < 0)
    {
      return err;
    }

  MPI_Bcast (X + number_of_blocks_row * size_of_block, size_of_tail, MPI_DOUBLE, proc, comm);

  // Ux = y
  //printf_main_process ("... vector Y:\n");
  //print_d_array_in_i_proc (X, size_of_system, MAIN_PROCESS, 30);


  if (respon_proc_for_tail)
    {
      calc_tail_U (B, X);
    }

  MPI_Bcast (X, size_of_system, MPI_DOUBLE, proc, comm);
  MPI_Bcast (B, size_of_system, MPI_DOUBLE, proc, comm);



  copy_offsets ();
  for (int i_block = number_of_blocks_row - 1; i_block >= 0; i_block--)
    {
      int owner = get_owner_of_big_part (i_block);
      int shift = i_block * size_of_block;
      if (proc_number == owner)
        {
          find_block_x (i_block, B + shift, X + shift);
          U_update_b (i_block, B + shift, X);
        }
      MPI_Bcast (X, size_of_system, MPI_DOUBLE, owner, comm);
      MPI_Bcast (B + shift, size_of_block, MPI_DOUBLE, owner, comm);
    }

  if (proc_number == MAIN_PROCESS)
    {
      memcpy (X, B, size_of_system * sizeof (double));
    }
  return 0;
}

void LU_matrix::mult_block_i_row (double *recv_buf, int i_block_row, int row_size, double *X, double *B, double &max_num, double &max_det)
{
  if (proc_number == MAIN_PROCESS)
    {
      int sq_size = size_of_block * row_size;
      double *Bi = B + i_block_row * size_of_block;

      for (int i_row = 0; i_row < row_size; i_row++)
        {
          double *Ai = recv_buf + size_of_block * i_row;

          double t = 0;
          for (int i_block_col = 0; i_block_col < number_of_blocks_row; i_block_col++)
            {
              double *Ai_bloc = Ai + i_block_col * sq_size;
              double *X_loc = X + i_block_col * size_of_block;

              for (int j_loc = 0; j_loc < size_of_block; j_loc++)
                {
                  t += Ai_bloc[j_loc] * X_loc[j_loc];
                }
            }
          if (tail)
            {
              double *Ai_bloc = recv_buf + number_of_blocks_row * sq_size + i_row * size_of_tail;
              double *X_loc = X + number_of_blocks_row * size_of_block;
              for (int j_loc = 0; j_loc < size_of_tail; j_loc++)
                {
                  t += Ai_bloc[j_loc] * X_loc[j_loc];
                }
            }
          double cand = Bi[i_row];
          t -= cand;

          max_num = (max_num > fabs(t) ? max_num : fabs(t));
          max_det = (max_det > fabs(cand) ? max_det : fabs(cand));
        }
    }
}

double LU_matrix::find_res (double *X, double *B)
{
  copy_offsets ();
  double *send_buf = buf_block_row_ptr_1.get ();
  double *recv_buf = buf_block_row_ptr_2.get ();
  int sq_size = size_of_block * size_of_block;

  double max_num = 0;
  double max_den = 0;

  for (int i_block_row = 0; i_block_row < number_of_blocks_row; i_block_row++)
    {
      /// only for test
      //memset (send_buf, 0, sizeof (double) * size_of_box_of_blocks * number_of_proc * size_of_block * size_of_block);
      //memset (recv_buf, 0, sizeof (double) * size_of_box_of_blocks * number_of_proc * size_of_block * size_of_block);
      int owner = get_owner_of_big_part (i_block_row);

      int size_of_first_part = sq_size * (i_block_row + 1);

      if (proc_number == owner)
        {
          collet_first_L_and_diag_blocks (i_block_row);
          if (proc_number == MAIN_PROCESS)
            {
              memcpy (recv_buf, send_buf, sizeof (double) * size_of_first_part);
            }
          else
            {
              MPI_Send (send_buf, size_of_first_part, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
            }
        }
      else
        {
          if (proc_number == MAIN_PROCESS)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_of_first_part, MPI_DOUBLE, owner, 0 /* tag */, comm, &st);
            }
        }

      int send_proc = owner;
      int shift_of_block = 0;

      for (int i_block = i_block_row + 1; i_block < number_of_blocks_row; i_block++)
        {
          double *U_i = get_U_off_diag (i_block_row) + shift_of_block * sq_size;
          double *recv_buf_i = recv_buf + i_block * sq_size;
          send_proc = next_proc (send_proc);
          if (proc_number == send_proc)
            {
              if (proc_number == MAIN_PROCESS)
                {
                  memcpy (recv_buf_i, U_i, sizeof (double) * sq_size);
                }
              else
                {
                  MPI_Send (U_i, sq_size, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
                }
              shift_of_block++;
            }
          else
            {
              if (proc_number == MAIN_PROCESS)
                {
                  MPI_Status st;
                  MPI_Recv (recv_buf_i, sq_size, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
                }
            }
        }
      if (tail)
        {
          double *U_tail = get_U_tail (i_block_row);
          double *recv_buf_i = recv_buf + number_of_blocks_row * sq_size;
          send_proc = next_proc (send_proc);
          if (proc_number == send_proc)
            {
              if (proc_number == MAIN_PROCESS)
                {
                  memcpy (recv_buf_i, U_tail, sizeof (double) * size_of_block * size_of_tail);
                }
              else
                {
                  MPI_Send (U_tail, size_of_block * size_of_tail, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
                }
            }
          else
            {
              if (proc_number == MAIN_PROCESS)
                {
                  MPI_Status st;
                  MPI_Recv (recv_buf_i, size_of_block * size_of_tail, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
                }
            }
        }

      mult_block_i_row (recv_buf, i_block_row, size_of_block, X, B, max_num, max_den);
    }

  if (tail)
    {
      double *L_tail = get_L_tail ();
      int send_proc = get_owner_of_big_part (number_of_blocks_row);

      if (respon_proc_for_tail)
        {

          if (proc_number == MAIN_PROCESS)
            {
              memcpy (recv_buf, L_tail, sizeof (double) * size_of_block * size_of_tail * number_of_blocks_row);
            }
          else
            {
              MPI_Send (L_tail, size_of_block * size_of_tail * number_of_blocks_row, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
            }
        }
      else
        {
          if (proc_number == MAIN_PROCESS)
            {
              MPI_Status st;
              MPI_Recv (recv_buf, size_of_block * size_of_tail * number_of_blocks_row, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
            }
        }

      double * diag_tail = get_diag_tail ();
      double * recv_tail = recv_buf + size_of_block * size_of_tail * number_of_blocks_row;

      if (respon_proc_for_tail)
        {
          if (proc_number == MAIN_PROCESS)
            {
              memcpy (recv_tail, diag_tail, sizeof (double) * size_of_tail * size_of_tail);
            }
          else
            {
              MPI_Send (diag_tail, size_of_tail * size_of_tail, MPI_DOUBLE, MAIN_PROCESS, 0 /* tag */, comm);
            }
        }
      else
        {
          if (proc_number == MAIN_PROCESS)
            {
              MPI_Status st;
              MPI_Recv (recv_tail, size_of_tail * size_of_tail, MPI_DOUBLE, send_proc, 0 /* tag */, comm, &st);
            }
        }
    }
  mult_block_i_row (recv_buf, number_of_blocks_row, size_of_tail, X, B, max_num, max_den);
  if (proc_number == MAIN_PROCESS)
    return max_num / max_den;
  return 0;
}





















