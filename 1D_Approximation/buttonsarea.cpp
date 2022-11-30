#include "buttonsarea.h"

InfoWindow::InfoWindow (ApproxFunx *F_arg, Information *info_arg, QWidget *parent)
  : QWidget (parent),
    info (info_arg),
    F (F_arg)
{
  setAttribute (Qt::WA_StyledBackground, true);
  setStyleSheet ("background-color : black");
  setSizePolicy (QSizePolicy::Minimum, QSizePolicy::Expanding);
}

void InfoWindow::paintEvent (QPaintEvent */*event*/)
{
  QPainter painter (this);
  //QFont info_font ("Times", 10, QFont::Bold);
  QPen info_pen (Qt::white, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (info_pen);
  //painter.setFont (info_font);
  //printf ("here: f_name = %s, width = %u\n", F->get_f_name (), width ());
  QString str;
  unsigned int middle = width () * 0.5;
  painter.drawText (0, 15, F->get_f_name ());
  painter.drawLine (middle, 18, middle, 137 * 2);
  painter.drawLine (0, 18, width (), 18);

  unsigned int step = 30;
  str.setNum (info->a, 'e', 1);
  painter.drawText (0, 15 + 1 * step, "a"); painter.drawText (middle + 5, 15 + 1 * step, str);
  painter.drawLine (0, 17 + 1 * step, width (), 17 + 1 * step);

  str.setNum (info->b, 'e', 1);
  painter.drawText (0, 15 + 2 * step, "b"); painter.drawText (middle + 5, 15 + 2 * step, str);
  painter.drawLine (0, 17 + step * 2, width (), 17 + step * 2);
  str.setNum (info->min_y, 'e', 1);
  painter.drawText (0, 15 + 3 * step, "min_y"); painter.drawText (middle + 5, 14 + 3 * step, str);
  painter.drawLine (0, 17 + step * 3, width (), 17 + step * 3);
  str.setNum (info->max_y, 'e', 1);
  painter.drawText (0, 15 + 4 * step, "max_y"); painter.drawText (middle + 5, 15 + 4 * step, str);
  painter.drawLine (0, 17 + step * 4, width (), 17 + step * 4);
  str.setNum (info->zoom);
  painter.drawText (0, 15 + 5 * step, "zoom"); painter.drawText (middle + 5, 15 + 5 * step, str);
  painter.drawLine (0, 17 + step * 5, width (), 17 + step * 5);
  str.setNum (info->n);
  painter.drawText (0, 15 + 6 * step, "n_nodes"); painter.drawText (middle + 5, 15 + 6 * step, str);
  painter.drawLine (0, 17 + step * 6, width (), 17 + step * 6);
  str.setNum (info->p);
  painter.drawText (0, 15 + 7 * step, "p_eps"); painter.drawText (middle + 5, 15 + 7 * step, str);
  painter.drawLine (0, 17 + step * 7, width (), 17 + step * 7);
  str.setNum (info->graph_id);
  painter.drawText (0, 15 + 8 * step , "grap_id"); painter.drawText (middle + 5, 15 + 8 * step, str);
  painter.drawLine (0, 17 + step * 8, width (), 17 + step * 8);

  QPen f_pen (Qt::green, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (f_pen);
  painter.drawText (0, height () - 10, "||f||=");
  painter.setPen (info_pen);
  str.setNum (info->res_f, 'e', 3);
  painter.drawText (middle - 5, height () - 10, str);

  unsigned int step2 = 34;

  QPen Pf1_pen (Qt::blue, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (Pf1_pen);
  painter.drawText (0, height () - 10 - step2 * 1, "||Pf||="); //painter.drawText (middle + 5, 120, str);
  painter.setPen (info_pen);
  str.setNum (info->res_Pf, 'e', 3);
  painter.drawText (middle - 5, height () - 10 - step2 * 1, str);

  QPen Pf2_pen (Qt::yellow, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (Pf2_pen);
  painter.drawText (0, height () - 10 - step2 * 2, "||Pf2||=");
  painter.setPen (info_pen);
  str.setNum (info->res_Pf2, 'e', 3);
  painter.drawText (middle - 5, height () - 10 - step2 * 2, str);

  QPen Dis1_pen (Qt::darkBlue, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (Dis1_pen);
  painter.drawText (0, height () - 10 - step2 * 3, "||Dis1||="); //painter.drawText (middle + 5, 120, str);
  painter.setPen (info_pen);
  str.setNum (info->res_Dis1, 'e', 3);
  painter.drawText (middle - 5, height () - 10 - step2 * 3, str);

  QPen Dis2_pen (Qt::darkYellow, 0, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
  painter.setPen (Dis2_pen);
  painter.drawText (0, height () - 10 - step2 * 4, "||Dis2||=");
  painter.setPen (info_pen);
  str.setNum (info->res_Dis2, 'e', 3);
  painter.drawText (middle - 5, height () - 10 - step2 * 4, str);
}
void ButtonsArea::set_button (QPushButton **button, QWidget *parent, const char *buttons_name, QVBoxLayout *vbox)
{
  *button = new QPushButton (buttons_name, parent);
  (*button)->setSizePolicy (QSizePolicy::Minimum, QSizePolicy::Minimum);
  vbox->addWidget (*button);
}

ButtonsArea::ButtonsArea (ApproxFunx *F_arg, Information *info, QWidget *parent)
  : QWidget (parent)
{
  QVBoxLayout *vbox = new QVBoxLayout (this);
  vbox->setSpacing (1);

  set_button (&change_func , this, "0. Change function", vbox);
  set_button (&change_graph, this, "1. Change graphic" , vbox);
  set_button (&plus_size   , this, "2. x2 size"        , vbox);
  set_button (&minus_size  , this, "3. /2 size"        , vbox);
  set_button (&plus_nodes  , this, "4. x2 nodes"       , vbox);
  set_button (&minus_nodes , this, "5. /2 nodes"       , vbox);
  set_button (&plus_eps    , this, "6. + epsilon"      , vbox);
  set_button (&minus_eps   , this, "7. - epsilon"      , vbox);


  text_info = new InfoWindow (F_arg, info, this);
  vbox->addWidget (text_info);
}
