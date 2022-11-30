#include "grapharea.h"

void GraphArea::set_start_max_min (double &min_y, double &max_y)
{
  min_y = info->MAX;
  max_y = info->MIN;
}

void GraphArea::correcting_min_max (double &min_y, double &max_y)
{
  // if you want to see axis X always
  if (min_y > 0)
    min_y = 0;
  if (max_y < 0)
    max_y = 0;
  if (eq (min_y, max_y) && eq(min_y, 0))
    {
      min_y = -5;
      max_y = 5;
    }
  // if you wan't to see axis X always
  /*
  if (eq (min_y, max_y))
    {
      min_y -= 5;
      max_y += 5;
    }
  */
  const double delta = (max_y - min_y) * 0.1;
  max_y += delta;
  min_y -= delta;
}

GraphArea::GraphArea (ApproxFunx *F_arg, Information *info_arg, QWidget *parent)
  : QWidget (parent),
    F (F_arg),
    info (info_arg)
{
  setMinimumSize (QSize (500, 500));
  setAttribute (Qt::WA_StyledBackground, true);
  setStyleSheet ("background-color : black");
  setSizePolicy (QSizePolicy::Expanding, QSizePolicy::Expanding);
  info->a = F->get_a ();
  info->b = F->get_b ();

  info->paint_f = 1;
  if (info->n <= MAX_NUM_OF_NODES)
    {
      info->paint_Pf1 = 1;
    }
  else
    {
      printf ("I can't build Pf1 for more than 50 nodes\n");
      info->paint_Pf1 = 0;
    }
  info->paint_Pf2 = 0;
  info->paint_Dis1 = 0;
  info->paint_Dis2 = 0;

  double min, max;
  set_start_max_min (min, max);
  F->get_min_max_on_ab (min, max, info->res_Dis1, info->a, info->b, width (), &ApproxFunx::DIS1);
  set_start_max_min (min, max);
  F->get_min_max_on_ab (min, max, info->res_Dis2, info->a, info->b, width (), &ApproxFunx::DIS2);

  set_start_max_min (info->min_y, info->max_y);
  F->get_min_max_on_ab (info->min_y, info->max_y, info->res_f, info->a, info->b, width (), &ApproxFunx::F);
  info->global_abs_max_f_on_ab = info->res_f;
  if (info->paint_Pf1) F->get_min_max_on_ab (info->min_y, info->max_y, info->res_Pf,info->a, info->b, width (), &ApproxFunx::Pf1);
  correcting_min_max (info->min_y, info->max_y);
}

void GraphArea::transform_coordinate_system (QPainter &painter)
{
  double b = info->b, a = info->a;
  double min_y = info->min_y, max_y = info->max_y;
  painter.translate (0.5 * width (), 0.5 * height ());
  painter.scale (width () / (b - a), -height () / (max_y - min_y));
  painter.translate (-0.5 * (a + b), -0.5 * (min_y + max_y));
}

void GraphArea::paint_graph (QPainter &painter, double (ApproxFunx::*Pf) (double), const QColor &color)
{
  QPen orig_func_pen (color, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
  orig_func_pen.setCosmetic (true);
  painter.setPen (orig_func_pen);

  double b = info->b, a = info->a;

  const double delta_x = (b - a) / width ();

  //if (color == Qt::darkBlue) printf ("starte DARKBLUE\n");
  //if (color == Qt::darkYellow) printf ("starte DARKYellow\n");
  double x1 = a;
  double y1 = (F->*Pf) (x1);
  for (int i = 1; i < width (); i++)
    {
      double x2 = a + i * delta_x;
      double y2 = (F->*Pf) (x2);
      //printf ("x1 = %e, y1 = %e, x2 = %e, y2 = %e\n", x1, y1, x2, y2);
      //printf ("tut\n");
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
      //printf ("?tut\n");
      x1 = x2, y1 = y2;
    }
  double x2 = b;
  double y2 = (F->*Pf) (x2);
  painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
}

void GraphArea::paint_axis (QPainter &painter)
{
  QPen axis_pen (Qt::red, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (axis_pen);
  painter.drawLine (QPointF (info->a, 0), QPointF (info->b, 0));
  painter.drawLine (QPointF (0, info->max_y), QPointF (0, info->min_y));
}

void GraphArea::paintEvent (QPaintEvent */*event*/)
{
  QPainter painter (this);

  painter.save ();
  transform_coordinate_system (painter);

  if (info->paint_f) paint_graph (painter, &ApproxFunx::F, Qt::green);
  if (info->paint_Pf1) paint_graph (painter, &ApproxFunx::Pf1, (Qt::blue));
  if (info->paint_Pf2) paint_graph (painter, &ApproxFunx::Pf2, (Qt::yellow));
  if (info->paint_Dis1) paint_graph (painter, &ApproxFunx::DIS1, (Qt::darkBlue));
  if (info->paint_Dis2) paint_graph (painter, &ApproxFunx::DIS2, (Qt::darkYellow));
  paint_axis (painter);

  painter.restore ();
}

