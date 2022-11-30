#ifndef GRAPHAREA_H
#define GRAPHAREA_H

#include <QWidget>
#include <QPainter>
#include <QMouseEvent>
#include <QLabel>

#include "functions.h"

class GraphArea : public QWidget
{
  Q_OBJECT
private:
  ApproxFunx *F = nullptr;
  Information *info = nullptr;


public:
  GraphArea (ApproxFunx *F_arg, Information *info, QWidget *parent = nullptr);
  ~GraphArea () {};

protected:
  void paintEvent (QPaintEvent *event);

private:
  void find_minY_max_Y (double (*f) (double));
  void transform_coordinate_system (QPainter &painter);
  void paint_graph (QPainter &painter, double (ApproxFunx::*Pf) (double), const QColor &color);
  void paint_axis (QPainter &painter);

public:
  void set_start_max_min (double &min_y, double &max_y);
  void correcting_min_max (double &min_y, double &max_y);
};

#endif // GRAPHAREA_H
