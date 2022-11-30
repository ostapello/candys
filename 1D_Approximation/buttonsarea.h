#ifndef BUTTONSAREA_H
#define BUTTONSAREA_H

#include <QWidget>
#include <QPushButton>
#include <QLabel>
#include <QBoxLayout>
#include <QPainter>
#include <QString>

#include "functions.h"

class MainWindow;

class InfoWindow : public QWidget
{
  Q_OBJECT

private:
  Information *info;
  ApproxFunx *F;

public:
  InfoWindow (ApproxFunx *F_arg = nullptr, Information *info = nullptr, QWidget *parent = nullptr);
  ~InfoWindow () {}

protected:
  void paintEvent (QPaintEvent *event);

};


class ButtonsArea : public QWidget
{
  Q_OBJECT

private:
  QPushButton *change_func = nullptr;
  QPushButton *change_graph = nullptr;
  QPushButton *plus_size = nullptr;
  QPushButton *minus_size = nullptr;
  QPushButton *plus_nodes = nullptr;
  QPushButton *minus_nodes = nullptr;
  QPushButton *plus_eps = nullptr;
  QPushButton *minus_eps = nullptr;
  InfoWindow *text_info = nullptr;

public:
  ButtonsArea (ApproxFunx *F = nullptr, Information *info = nullptr, QWidget *parent = nullptr);
  ~ButtonsArea () {};

  friend class MainWindow;

private:
  void set_button (QPushButton **button, QWidget *parent, const char *buttons_name, QVBoxLayout *vbox);

};

#endif // BUTTONSAREA_H
