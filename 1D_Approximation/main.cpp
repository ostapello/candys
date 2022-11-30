#include "window.h"
//#include "functions.h"
#include <fenv.h>

#include <QApplication>

int parse_arg (int argc, char *argv[], double &a, double &b, unsigned int &n, int &k)
{
  if (argc != 5)
    {
      printf ("Usage : %s a b n k\n"
              "        a - left border\n"
              "        b - right border\n"
              "        n - starting number of nodes\n"
              "        k - starting function's id\n", argv[0]);
      return -1;
    }
  if (sscanf (argv[1], "%lf", &a) != 1)
    {
      printf ("INPUT ERROR   : a should be double\n"
              "your input is : %s\n", argv[1]);
      return -2;
    }
  if (sscanf (argv[2], "%lf", &b) != 1)
    {
      printf ("INPUT ERROR   : b should be double\n"
              "your input is : %s\n", argv[2]);
      return -3;
    }
  if (b - a <= MAX_DELTA_BA)
    {
      printf ("INPUT ERROR   : b - a should be more then %e\n"
              "your input is : a = %e, b = %e\n", MAX_DELTA_BA, a, b);
      return -4;
    }
  if (sscanf (argv[3], "%u", &n) != 1 || n < 2)
    {
      printf ("INPUT ERROR   : n is integer number more then 1\n"
              "your input is : n = %s\n", argv[3]);
      return -5;
    }
  if (sscanf (argv[4], "%d", &k) != 1 || k < 0 || k > 6)
    {
      printf ("INPUT ERROR   : k is nonegative integer number less then 7\n"
              "your input is : k = %s\n", argv[4]);
      return -5;
    }
  //printf ("k = %u\n", k);
  return 0;
}

int main(int argc, char *argv[])
{
  double a, b;
  unsigned int n;
  int k;

  feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW| FE_UNDERFLOW);

  if (int res = parse_arg (argc, argv, a, b, n, k); res < 0)
    {
      return 0;
    }

  QApplication app(argc, argv);
  MainWindow Mwindow (a, b, n, k);
  Mwindow.show();
  return app.exec();
}
