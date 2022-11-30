#include "functions.h"

unsigned int get_next_func_id (unsigned int k)
{
  return (k >= NUM_OF_FUNC - 1 ? 0 : k + 1);
}

unsigned int get_next_graph_id (unsigned int k)
{
  return (k >= NUM_OF_GRAPH - 1 ? 0 : k + 1);
}


int eq (double a, double b)
{
  if (a < b)
    return 0;
  if (b < a)
    return 0;
  return 1;
}

double f_0 (double /*x*/)
{
  return 1;
}
double f_1 (double x)
{
  return x;
}
double f_2 (double x)
{
  return x * x;
}
double f_3 (double x)
{
  return x * x * x;
}
double f_4 (double x)
{
  return x * x * x * x;
}
double f_5 (double x)
{
  return exp (x);
}
double f_6 (double x)
{
  return 1 / (25 * x * x + 1);
}


double df_0 (double )
{
  return 0;
}
double df_1 (double )
{
  return 1;
}
double df_2 (double x)
{
  return 2 * x;
}
double df_3 (double x)
{
  return 3 * x * x;
}
double df_4 (double x)
{
  return 4 * x * x * x;
}
double df_5 (double x)
{
  return exp (x);
}
double df_6 (double x)
{
  return -50 * x / ((25 * x * x + 1) * (25 * x * x + 1));
}


Information::Information (double a_arg, double b_arg, unsigned int n_arg, int func_id_arg)
{
  a = a_arg;
  b = b_arg;
  n = n_arg;
  func_id = func_id_arg;
}

void Information::set_graph_mode (unsigned int graph_id)
{
  switch (graph_id)
    {
    case 0:
      paint_f = 1;
      if (n <= MAX_NUM_OF_NODES)
        {
          paint_Pf1 = 1;
        }
      else
        {
          printf ("I can't build Pf1 for more than 50 nodes\n");
          paint_Pf1 = 0;
        }
      paint_Pf2 = 0;
      paint_Dis1 = 0;
      paint_Dis2 = 0;
      break;
    case 1:
      paint_f = 1;
      paint_Pf1 = 0;
      paint_Pf2 = 1;
      paint_Dis1 = 0;
      paint_Dis2 = 0;
      break;
    case 2:
      paint_f = 1;
      if (n <= MAX_NUM_OF_NODES)
        {
          paint_Pf1 = 1;
        }
      else
        {
          printf ("I can't build Pf1 for more than 50 nodes\n");
          paint_Pf1 = 0;
        }
      paint_Pf2 = 1;
      paint_Dis1 = 0;
      paint_Dis2 = 0;
      break;
    case 3:
      paint_f = 0;
      paint_Pf1 = 0;
      paint_Pf2 = 0;
      if (n <= MAX_NUM_OF_NODES)
        {
          paint_Dis1 = 1;
        }
      else
        {
          printf ("I can't build Pf1 for more than 50 nodes\n");
          paint_Dis1 = 0;
        }
      paint_Dis2 = 1;
      break;
    }
}

double ApproxFunx::F (double x)
{
  return f(x);
}

double ApproxFunx::Pf1 (double value)
{
  double *a = A.data ();
  double *x = X.data ();
  double t = a[2 * number_of_nodes_Pf1 - 1];
  for (int i = 2 * number_of_nodes_Pf1 - 2; i >= 0; i--)
    {
      t *= (value - x[i]);
      t += a[i];
    }
  return t;
}

void ApproxFunx::fill_X_by_x ()
{
  if (X.size () < number_of_nodes_Pf1 * 2)
    X.resize (number_of_nodes_Pf1 * 2);
  const double step = (b - a) / (number_of_nodes_Pf1 - 1);
  double *x = X.data ();
  //printf ("step = %f, number_of_nodes = %u\n", step, number_of_nodes);
  for (unsigned int i = 0; i < (number_of_nodes_Pf1 - 1); i++)
    {
      //printf ("i = %u, a[i] = %f\n", i, a + step * i);
      x[2 * i] = a + step * i;
      x[2 * i + 1] = a + step * i;
    }
  x[2 * number_of_nodes_Pf1 - 2] = b;
  x[2 * number_of_nodes_Pf1 - 1] = b;
}


void ApproxFunx::fill_A_by_fx (double eps)
{
  if (A.size () < number_of_nodes_Pf1 * 2)
    A.resize (number_of_nodes_Pf1 * 2);
  const double *x = X.data ();
  double *a = A.data ();
  for (unsigned int i = 0; i < 2 * number_of_nodes_Pf1; i++)
    {
      a[i] = f(x[i]);
    }
  if (number_of_nodes_Pf1 & 1)
    {
      a[number_of_nodes_Pf1] += eps;
      a[number_of_nodes_Pf1 - 1] += eps;
    }
  else
    {
      a[number_of_nodes_Pf1] += eps;
      a[number_of_nodes_Pf1 + 1] += eps;
    }
}

void ApproxFunx::first_step_of_filling_table ()
{
  double *a = A.data ();
  double *x = X.data ();
  for (int i = 2 * number_of_nodes_Pf1 - 1; i >= 2 ; i -= 2)
    {
      a[i] = df (x[i]);
      a[i - 1] = (a[i - 1] - a[i - 2]) / (x[i - 1] - x[i - 2]);
    }
  a[1] = df (x[1]);
}

void ApproxFunx::calc_A ()
{
  first_step_of_filling_table ();
  //print_vector (A);

  double *a = A.data ();
  const double *x = X.data ();
  //print_vector (A);
  //print_vector (X);
  for (unsigned int i = 1; i < 2 * number_of_nodes_Pf1; i++)
    {
      //printf ("====== i = %u\n", i);
      for (unsigned int j = 2 * number_of_nodes_Pf1 - 1; j > i; j--)
        {
          //printf ("a[%u] = %7.2f; a[%u] = %7.2f; x[%u] = %7.2f; x[%u] = %7.2f\n", j, a[j], j - 1, a[j - 1], j, x[j], j - i - 1, x[j - i - 1]);
          a[j] = (a[j] - a[j - 1]) / (x[j] - x[j - i - 1]);
        }
      //print_vector (A);
    }
}

int ApproxFunx::build_approx_func_Pf1 (unsigned int number_of_nodes_arg, double eps)
{
  if (number_of_nodes_arg > MAX_NUM_OF_NODES)
    return CAN_NOT_DO;
  number_of_nodes_Pf1 = number_of_nodes_arg;


  fill_X_by_x ();
  //print_vector (A);
  fill_A_by_fx (eps);
  //print_vector (A);
  calc_A ();
  return DONE;
}

void ApproxFunx::fill_FD (double eps)
{
  if (FD.size () < number_of_nodes_Pf2)
    FD.resize (number_of_nodes_Pf2);

  Eps = eps;

  double *F = FD.data ();
  double x = a;

  // fill FD as f(x)
  for (unsigned int i = 0; i < number_of_nodes_Pf2; i++)
    {
      F[i] = f(x + h * i);
    }
  F[number_of_nodes_Pf2 / 2] += eps;

  //print_vector (FD);

  double f_1 = F[0];
  for (unsigned int i = 0; i < number_of_nodes_Pf2 - 1; i++)
    {
      double f_2 = F[i + 1];
      F[i] = (f_2 - f_1) / h;
      f_1 = f_2;
    }
}

void ApproxFunx::fill_border_coef (double x0, double x_n1)
{
  if (B.size () < number_of_nodes_Pf2)
    B.resize (number_of_nodes_Pf2);

  double FD0 = FD[0];
  double l = x0 - a;
  double l2 = l * l;
  double l3 = l2 * l;
  A1 = h * h * l - 2 * h * l2 + l3;
  A2 = l3 - h * l2;
  B[0] = h * h * (f(x0) - f(a)) - 3 * h * FD0 * l2 + 2 * FD0 * l3;

  double FDN = FD[number_of_nodes_Pf2 - 2];
  double r = x_n1 - (b - h);
  double r2 = r * r;
  double r3 = r2 * r;
  A3 = h * h * r - 2 * h * r2 + r3;
  A4 = r3  - h * r2;
  // bad variation of + eps
  if (number_of_nodes_Pf2 == 3)
    {
      B[number_of_nodes_Pf2 - 1] = h * h * (f(x_n1) - f(b - h) - Eps) - 3 * h * FDN * r2 + 2 * FDN * r3;
    }
  else
    {
      B[number_of_nodes_Pf2 - 1] = h * h * (f(x_n1) - f(b - h)) - 3 * h * FDN * r2 + 2 * FDN * r3;
    }

  //printf ("A1 = %e\nA2 = %e\nA3 = %e\nA4 = %e, B[0] = %e, B[n] = %e\n", A1, A2, A3, A4, B[0], B[number_of_nodes - 1]);
}

void ApproxFunx::LU_fact ()
{
  if (L.size () < number_of_nodes_Pf2)
    L.resize (number_of_nodes_Pf2);
  if (U.size () < number_of_nodes_Pf2)
    U.resize (number_of_nodes_Pf2);
  double *l = L.data ();
  double *u = U.data ();

  l[0] = A1;
  u[0] = A2 / A1;

  for (unsigned int i = 1; i < number_of_nodes_Pf2 - 1; i++)
    {
      double t = h * (4 - u[i - 1]);
      l[i] = t;
      u[i] = h / t;
    }

  u[number_of_nodes_Pf2 - 1] = A3; // actualy it is l_nn-1
  l[number_of_nodes_Pf2 - 1] = A4 - A3 * u[number_of_nodes_Pf2 - 2];
  //printf ("A3 = %e\n", A3);

}

// B[0] and B[n-1] we already know
void ApproxFunx::fill_rhs ()
{
  double *b = B.data ();
  double *fd = FD.data ();
  double fd0 = fd[0];
  for (unsigned int i = 1; i < number_of_nodes_Pf2 - 1; i++)
    {
      double fd1 = fd[i];
      b[i] = h * (3 * fd0 + 3 * fd1);
      fd0 = fd1;
    }
}

int ApproxFunx::Ly_B ()
{
  if (D.size () < number_of_nodes_Pf2)
    D.resize (number_of_nodes_Pf2);

  double *y = D.data ();
  double *b = B.data ();
  double *l = L.data ();

  double t = l[0];
  if (eq (t, 0))
    return -1;
  double y0 = b[0] / t;
  y[0] = y0;
  for (unsigned int i = 1; i < number_of_nodes_Pf2 - 1; i++)
    {
      double t = l[i];
      if (eq (t, 0))
        return -1;
      double y1 = (b[i] - h * y0) / t;
      y[i] = y1;
      y0 = y1;
    }

  t = l[number_of_nodes_Pf2 - 1];
  if (eq (t, 0))
    return -1;

  y[number_of_nodes_Pf2 - 1] = (b[number_of_nodes_Pf2 - 1] - U[number_of_nodes_Pf2 - 1] * y0) / t;
  return 0;
}

void ApproxFunx::Ux_Y ()
{
  double *y = D.data ();
  double *x = B.data ();
  double *u = U.data ();

  double xn = y[number_of_nodes_Pf2 - 1];
  x[number_of_nodes_Pf2 - 1] = xn;

  for (int i = number_of_nodes_Pf2 - 2; i >= 0; i--)
    {
      double x0 = y[i] - u[i] * xn;
      //printf ("y[i] = %e\nu[i] = %e\nxn = %e\n", y[i], u[i], xn);
      //printf ("u[i] * xn = %e\n", u[i] * xn);
      x[i] = x0;
      xn = x0;
    }
}

int ApproxFunx::find_d ()
{
  /// Ly = B
  /// unswer in D
  if (int res = Ly_B (); res < 0)
    return res;

  //print_vector (D, "D_");

  /// Ux = y
  /// unswer in B
  Ux_Y ();

  std::copy_n (B.data (), number_of_nodes_Pf2, D.data ());
  //print_vector (D, "D");
  return 0;
}

int ApproxFunx::build_approx_func_Pf2 (unsigned int number_of_nodes_arg, double eps)
{
  number_of_nodes_Pf2 = number_of_nodes_arg;

  h = (b - a) / (number_of_nodes_Pf2 - 1);

  fill_FD (eps);
  //print_vector (FD);
  fill_border_coef (a - h * 0.5, b + h * 0.5);


  LU_fact ();
  //print_vector (L, "L");
  //print_vector (U, "U");
  //print_vector (B, "B");
  fill_rhs ();
  //print_vector (B, "B");

  if (int res = find_d (); res < 0)
    {
      printf ("I can't do LU factorization for Pf2 with number of nodes = %u\n", number_of_nodes_Pf2);
      return -1;
    }

  return 0;
}

double ApproxFunx::Pf2 (double x)
{

  double l = x - a;
  int k = l / h; // spline's number
  if (k >= (int) number_of_nodes_Pf2 - 1)
    {
      k = number_of_nodes_Pf2 - 2;
    }
  if (k < 0)
    {
      k = 0;
    }

  double x_0 = a + k * h;
  double f_x0 = f (x_0);
  if (k == ((int) number_of_nodes_Pf2) / 2)
    f_x0 += Eps;
  double x_1 = x_0 + h;
  double di = D[k];
  double Fi = FD[k];
  return f_x0 + (x - x_0) * di + (x - x_0) * (x - x_0) * (Fi - di) / h + (x - x_0) * (x - x_0) * (x - x_1) *
      (di + D[k + 1] - 2 * Fi) / (h * h);
}

// min_y and max_y can has starting values
void ApproxFunx::get_min_max_on_ab (double &min_y, double &max_y, double &global_fabs, double a, double b, int n_x_pieces,
                                    double (ApproxFunx::*FF) (double))
{
  const double delta_x = (b - a) / n_x_pieces;
  double x1 = a;

  global_fabs = 0;

  for (int i = 0; i < n_x_pieces; i++)
    {
      double y1 = (this->*FF) (x1 + i * delta_x);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
      double abs_y1 = fabs (y1);
      if (abs_y1 > global_fabs)
        global_fabs = abs_y1;
    }
}


void ApproxFunx::set_func_param (unsigned int func_id)
{
  switch (func_id)
    {
    case 0:
      f_name = "f(x) = 1";
      f = f_0;
      df = df_0;
      break;
    case 1:
      f_name = "f(x) = x";
      f = f_1;
      df = df_1;
      break;
    case 2:
      f_name = "f(x) = x^2";
      f = f_2;
      df = df_2;
      break;
    case 3:
      f_name = "f(x) = x^3";
      f = f_3;
      df = df_3;
      break;
    case 4:
      f_name = "f(x) = x^4";
      f = f_4;
      df = df_4;
      break;
    case 5:
      f_name = "f(x) = e^x";
      f = f_5;
      df = df_5;
      break;
    case 6:
      f_name = "f(x) = 1 / (25x^2 + 1)";
      f = f_6;
      df = df_6;
      break;
    default:
      printf ("(._.) Something wrong in set_func_param (._.)\n");
      f_name = "(._.) UNKNOWN ERROR (._.)";
      f = f_0;
      df = df_0;
      break;
    }
}

double ApproxFunx::DIS1 (double x)
{
  return fabs(f(x) - Pf1(x));
}

double ApproxFunx::DIS2 (double x)
{
  return fabs(f(x) - Pf2(x));
}

ApproxFunx::ApproxFunx (double a_arg, double b_arg, unsigned int number_of_nodes_arg, int func_id_arg)
  : a (a_arg), b(b_arg)
{
  set_func_param (func_id_arg);
  if (number_of_nodes_arg > MAX_NUM_OF_NODES)
    {
      printf ("I build first approximate function for 50 nodes\n");
      build_approx_func_Pf1 (MAX_NUM_OF_NODES);
    }
  else
    {
      build_approx_func_Pf1 (number_of_nodes_arg);
    }
  build_approx_func_Pf2 (number_of_nodes_arg);
}


void ApproxFunx::print_vector (std::vector <double> &T, const char *text)
{
  if (text)
    {
      printf ("%s: ", text);
    }
  for (unsigned int i = 0; i < T.size (); i++)
    {
      if (fabs (T[i]) < 1.e-2)
        printf ("%8.2e ", T[i]);
      else
        printf ("%8.16e ", T[i]);
    }
  printf ("\n");
}
