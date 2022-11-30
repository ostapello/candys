#include "window.h"

void MainWindow::set_residual ()
{
  double min, max;
  graph->set_start_max_min (min, max);
  max = 0;
  if (info.n <= MAX_NUM_OF_NODES)
    Pf.get_min_max_on_ab (min, max, info.res_Dis1, info.a, info.b, graph->width (), &ApproxFunx::DIS1);
  graph->set_start_max_min (min, max);
  max = 0;
  Pf.get_min_max_on_ab (min, max, info.res_Dis2, info.a, info.b, graph->width (), &ApproxFunx::DIS2);
  printf ("== residual 1 = %e; residual 2 = %e\n", info.res_Dis1, info.res_Dis2);
}

void MainWindow::set_global_max ()
{
  graph->set_start_max_min (info.min_y, info.max_y);
  Pf.get_min_max_on_ab (info.min_y, info.max_y, info.res_f, info.a, info.b, graph->width (), &ApproxFunx::F);
  info.global_abs_max_f_on_ab = info.res_f;
}

void MainWindow::set_local_max ()
{
  graph->set_start_max_min (info.min_y, info.max_y);
  if (info.paint_f)
    Pf.get_min_max_on_ab (info.min_y, info.max_y, info.res_f, info.a, info.b, graph->width (), &ApproxFunx::F);
  if (info.paint_Pf1)
    Pf.get_min_max_on_ab (info.min_y, info.max_y, info.res_Pf, info.a, info.b, graph->width (), &ApproxFunx::Pf1);
  if (info.paint_Pf2)
    Pf.get_min_max_on_ab (info.min_y, info.max_y, info.res_Pf2, info.a, info.b, graph->width (), &ApproxFunx::Pf2);
  double global;
  if (info.paint_Dis1)
      Pf.get_min_max_on_ab (info.min_y, info.max_y, global, info.a, info.b, graph->width (), &ApproxFunx::DIS1);
  if (info.paint_Dis2)
    Pf.get_min_max_on_ab (info.min_y, info.max_y, global, info.a, info.b, graph->width (), &ApproxFunx::DIS2);
  graph->correcting_min_max (info.min_y, info.max_y);
}

void MainWindow::plus_eps ()
{
  unsigned int p = info.p;
  if (p + 1 > MAX_PLUS_EPS)
    {
      printf ("Epsilon is |f| already. I assume it's too much (You can change #define MAX_PLUS_EPS 10 in functions.h)\n");
      return;
    }
  info.p++;

  Pf.build_approx_func_Pf1 (info.n, EPS * info.p * info.global_abs_max_f_on_ab);
  Pf.build_approx_func_Pf2 (info.n, EPS * info.p * info.global_abs_max_f_on_ab);
  //printf ("Epsilon is %25.18f\n", EPS * info.p * info.global_abs_max_f_on_ab);
  set_residual ();
  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::minus_eps ()
{
  unsigned int p = info.p;
  if (!p)
    {
      return;
    }
  info.p--;
  Pf.build_approx_func_Pf1 (info.n, EPS * info.p * info.global_abs_max_f_on_ab);
  Pf.build_approx_func_Pf2 (info.n, EPS * info.p * info.global_abs_max_f_on_ab);

  set_residual ();
  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::plus_nodes ()
{
  unsigned int n = info.n * 2;
  info.p = 0;
  if (n > MAX_NUM_OF_NODES)
    {
      printf ("n is equal %u already. I can't build first approximate function for more then 50 nodes\n", info.n);
    }
  else
    {
      Pf.build_approx_func_Pf1 (n);
    }
  info.n = n;
  Pf.build_approx_func_Pf2 (info.n);

  set_residual ();
  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::minus_nodes ()
{
  if (info.n * 0.5 < 2)
    {
      info.n = 2;
    }
  else
    {
      info.n *= 0.5;
    }
  info.p = 0;

  Pf.build_approx_func_Pf1 (info.n);
  Pf.build_approx_func_Pf2 (info.n);

  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::plus_size ()
{
  if (info.zoom == MAX_ZOOM)
    {
      printf ("Zoom is already x20. I assume it's too much (You can change #define MAX_ZOOM 20 in functions.h\n");
      return;
    }
  double b = info.b;
  double a = info.a;
  double shift = (b - a) * 0.25;
  a += shift;
  b -= shift;
  if (b - a < MAX_DELTA_BA)
    {
      printf ("So big zoom. I can't do it\n");
      return;
    }
  info.zoom++;
  info.a += shift;
  info.b -= shift;

  set_residual ();
  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::minus_size ()
{
  if (info.zoom == 1)
    return;
  info.zoom--;
  double shift = (info.b - info.a) * 0.5;
  info.a -= shift;
  info.b += shift;

  set_residual ();
  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::change_graph ()
{
  info.graph_id = get_next_graph_id (info.graph_id);
  info.set_graph_mode (info.graph_id);

  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::change_func ()
{
  info.func_id = get_next_func_id (info.func_id);
  Pf.set_func_param (info.func_id);
  info.n = start_n; //
  Pf.build_approx_func_Pf1 (info.n);
  Pf.build_approx_func_Pf2 (info.n);
  info.zoom = 1;
  info.a = Pf.get_a ();
  info.b = Pf.get_b ();

  set_residual ();
  set_global_max ();
  set_local_max ();
  graph->update ();
  buttons->update ();
}

void MainWindow::connect_buttons (ButtonsArea *buttons)
{
  QAction action;
  // chane_func
  connect (buttons->change_func, &QPushButton::clicked, this, &MainWindow::change_func);
  QShortcut *cf_0 = new QShortcut (QKeySequence (tr ("0")), this);
  connect (cf_0, &QShortcut::activated, this, &MainWindow::change_func);

  connect (buttons->change_graph, &QPushButton::clicked, this, &MainWindow::change_graph);
  QShortcut *cg_1 = new QShortcut (QKeySequence (tr ("1")), this);
  connect (cg_1, &QShortcut::activated, this, &MainWindow::change_graph);

  // zooming
  connect (buttons->plus_size, &QPushButton::clicked, this, &MainWindow::plus_size);
  QShortcut *ps_2 = new QShortcut (QKeySequence (tr ("2")), this);
  connect (ps_2, &QShortcut::activated, this, &MainWindow::plus_size);
  connect (buttons->minus_size, &QPushButton::clicked, this, &MainWindow::minus_size);
  QShortcut *ps_3 = new QShortcut (QKeySequence (tr ("3")), this);
  connect (ps_3, &QShortcut::activated, this, &MainWindow::minus_size);

  // nodes
  connect (buttons->plus_nodes, &QPushButton::clicked, this, &MainWindow::plus_nodes);
  QShortcut *pn_4 = new QShortcut (QKeySequence (tr ("4")), this);
  connect (pn_4, &QShortcut::activated, this, &MainWindow::plus_nodes);
  connect (buttons->minus_nodes, &QPushButton::clicked, this, &MainWindow::minus_nodes);
  QShortcut *mn_5 = new QShortcut (QKeySequence (tr ("5")), this);
  connect (mn_5, &QShortcut::activated, this, &MainWindow::minus_nodes);

  // eps
  connect (buttons->plus_eps, &QPushButton::clicked, this, &MainWindow::plus_eps);
  QShortcut *pe_6 = new QShortcut (QKeySequence (tr ("6")), this);
  connect (pe_6, &QShortcut::activated, this, &MainWindow::plus_eps);
  connect (buttons->minus_eps, &QPushButton::clicked, this, &MainWindow::minus_eps);
  QShortcut *me_7 = new QShortcut (QKeySequence (tr ("7")), this);
  connect (me_7, &QShortcut::activated, this, &MainWindow::minus_eps);
}

void MainWindow::set_layout (QWidget *central_box)
{
  QHBoxLayout *hbox = new QHBoxLayout (central_box);

  graph = new GraphArea (&Pf, &info, central_box);
  buttons = new ButtonsArea (&Pf, &info, central_box);

  hbox->addWidget (buttons, 0, Qt::AlignLeft);
  hbox->addSpacing (1);
  hbox->addWidget (graph);
}

MainWindow::MainWindow(double a_arg, double b_arg, unsigned int n_arg, unsigned int k_arg, QWidget *parent)
  : QWidget (parent),
    info (a_arg, b_arg, n_arg, k_arg),
    Pf(a_arg, b_arg, n_arg, k_arg),
    start_n (n_arg)
{
  resize (1000, 800);
  set_layout (this);
  connect_buttons (buttons);
  setWindowTitle ("1D by Ostap Karkovskiy");
}
