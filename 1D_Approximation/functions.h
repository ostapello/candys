#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <math.h>
#include <stdio.h>
#include <algorithm>


#define EPS 1e-1
#define MAX_PLUS_EPS 10
#define MAX_DELTA_BA 1e-8
#define MAX_NUM_OF_NODES 50
#define MAX_ZOOM 20
#define NUM_OF_FUNC 7
#define NUM_OF_GRAPH 4
#define DEF_A -10
#define DEF_B 10
#define DEF_N 3

enum report {DONE = 0, CAN_NOT_DO = 1, ALREADY_DONE = 2};

int eq (double a, double b);

double f_0 (double x);
double f_1 (double x);
double f_2 (double x);
double f_3 (double x);
double f_4 (double x);
double f_5 (double x);
double f_6 (double x);

double df_0 (double x);
double df_1 (double x);
double df_2 (double x);
double df_3 (double x);
double df_4 (double x);
double df_5 (double x);
double df_6 (double x);



const char * get_f_name (unsigned int k);

unsigned int get_next_func_id (unsigned int k);
unsigned int get_next_graph_id (unsigned int k);

class Information
{
public:
  int func_id = 0;
  double a = 0;
  double b = 0;
  double min_y = 0;
  double max_y = 0;
  int zoom = 1;
  unsigned int n = 0;

  double MAX = std::numeric_limits <double>::max ();
  double MIN = -MAX;


  double global_abs_max_f_on_ab = 0;
  unsigned int p = 0;

  double res_f = 0;
  double res_Pf = 0;
  double res_Pf2 = 0;
  double res_Dis1 = 0;
  double res_Dis2 = 0;

  unsigned int paint_f = 1;
  unsigned int paint_Pf1 = 1;
  unsigned int paint_Pf2 = 1;
  unsigned int paint_Dis1 = 1;
  unsigned int paint_Dis2 = 1;

  int graph_id = 0;


  Information (double a_arg, double b_arg, unsigned int n_arg, int func_id_arg);
  void set_graph_mode (unsigned int graph_id);
};

class ApproxFunx
{
private:
  double a;
  double b;
  double h = 0;
  const char *f_name;
  double (*f) (double);
  double (*df) (double);
  unsigned int number_of_nodes_Pf1;
  unsigned int number_of_nodes_Pf2;
  // Pf
  std::vector <double> X; // a = x1 < ... < xn = b
  std::vector <double> A; // coefficients of app_func
  // Pf2
  std::vector <double> D; // coefficients for spline
  std::vector <double> FD; // coefficients for 1 divided differences
  std::vector <double> L; // for diag of L
  std::vector <double> U; // for diag of U
  std::vector <double> B; // for rhs
  double Eps = 0;

  // border's coefficients

  double A1;
  double A2;
  double A3;
  double A4;

public:
  ApproxFunx (double a = DEF_A, double b = DEF_B, unsigned int number_of_nodes_arg = DEF_N, int func_id_arg = 0);
  ~ApproxFunx () {};

  void set_func_param (unsigned int func_id);

  int build_approx_func_Pf1 (unsigned int number_of_nodes_arg, double eps = 0);
  void get_min_max_on_ab (double &min_y, double &max_y, double &global_fabs, double a, double b, int n_x_pieces,
                          double (ApproxFunx::*FF) (double));

  int build_approx_func_Pf2 (unsigned int number_of_nodes_arg, double eps = 0);

  double F (double x);
  double Pf1 (double x);
  double Pf2 (double x);
  double DIS1 (double x);
  double DIS2 (double x);

  /// getters
  double (*get_f ()) (double) {return f;}
  double get_a () {return a;}
  double get_b () {return b;}
  const char * get_f_name () {return f_name;}

private:
  // for Pf1
  void fill_A_by_fx (double eps = 0);
  void fill_X_by_x ();
  void calc_A ();
  void first_step_of_filling_table ();

  // for Pf2
  void fill_FD (double eps = 0);
  void fill_border_coef (double x0, double x_n1);
  void LU_fact ();
  void fill_rhs ();
  int find_d ();
  int Ly_B ();
  void Ux_Y ();

  // for test
  void print_vector (std::vector <double> &T, const char *text = nullptr);
};



#endif // FUNCTIONS_H
