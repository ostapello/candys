#ifndef WINDOW_H
#define WINDOW_H

#include <memory>
#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QPainter>
#include <QGridLayout>
#include <QAction>
#include <QShortcut>
#include <QKeySequence>

#include "functions.h"
#include "buttonsarea.h"
#include "grapharea.h"

class MainWindow : public QWidget
{
  Q_OBJECT

private:
  GraphArea *graph = nullptr;
  ButtonsArea *buttons = nullptr;
  Information info;
  ApproxFunx Pf;
  unsigned int start_n = 0;


public:
  MainWindow (double a_arg = DEF_A, double b_arg = DEF_B, unsigned int n_arg = DEF_N, unsigned int k_arg = 0,
              QWidget *parent = nullptr);
  ~MainWindow () {};
private:
  void set_layout (QWidget *central_box);
  void connect_buttons (ButtonsArea *buttons);

  void set_global_max ();
  void set_local_max ();
  void set_residual ();

private slots:
  void change_func ();
  void change_graph ();
  void plus_size ();
  void minus_size ();
  void plus_nodes ();
  void minus_nodes();
  void plus_eps ();
  void minus_eps();
};

#endif // WINDOW_H
