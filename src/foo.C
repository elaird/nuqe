#include <iostream>
#include "Rtypes.h"
#include "TT_params.h"
#include "TT_nucleus.h"
#include "TT_generator.h"

////////////////////////////////////////////////////////////////////////
Int_t main(Int_t argc,char **argv) {
  Bool_t help=0;

  if (help) {
    std::cout << "--USAGE: " << argv[0] << " CONFIG_FILE" << std::endl;
    std::cout << "--for further help, see readme" << std::endl;
    return 1;
  }

  Int_t* init_styles=0;
  Bool_t* init_styles_on=0;
  Int_t NRuns=0;
  if (argc==1) {
    NRuns=1;
    init_styles=new Int_t[NRuns];
    init_styles_on=new Bool_t[NRuns];
    init_styles[0]=-1;
    init_styles_on[0]=kTRUE;
  }
  if (argc>1) {
    NRuns=10;
    init_styles=new Int_t[NRuns];
    init_styles_on=new Bool_t[NRuns];
    for (Int_t iRun=0;iRun<NRuns;iRun++) {
      init_styles[iRun]=iRun;
      init_styles_on[iRun]=kTRUE;
    }
    init_styles_on[9]=kTRUE;
  }

  TT_params params(argv[1]);
  params.Set_rand_type_and_seed(3,0);
  //TT_nucleus nucleus(8,8);
  TT_nucleus nucleus(6,6);
  //nucleus.n_plot();

  for (Int_t iRun=0;iRun<NRuns;iRun++) {
    if (!init_styles_on[iRun]) continue;
  
    params.Init(init_styles[iRun]);
    Bool_t params_good=params.Check_for_problems();
    if (!params_good) return 1;
  
    TT_generator gen(&params,&nucleus);
    gen.Setup_processes();
    Bool_t status=gen.Generate_events();
    gen.Make_graphs();
    gen.Write_shuffled_tree();
    if (!status) return 1;
  }

  delete [] init_styles;
  delete [] init_styles_on;
  return 0;
}
