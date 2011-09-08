#include "TGraph.h"
///////////////////////////////////////////////////////////////////
class BenharDrawer {
public :
  BenharDrawer();
  ~BenharDrawer();
  void draw();
  TGraph graph_fg_ma10_;
  TGraph graph_fg_ma12_;
  TGraph graph_sf_ma10_;
  TGraph graph_sf_ma12_;

private :
  void adjust_graph(TGraph& graph);
  double x0_ps_;
  double x1_ps_;

  double y1_ps_;
  double y0_ps_;

  double x0_ph_;
  double x1_ph_;

  double y0_ph_;
  double y1_ph_;
};
