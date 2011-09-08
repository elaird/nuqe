#define TT_generator_cxx
#include "TT_generator.h"
#include "TGraph.h"

////////////////////////////////////////////////////////////////////////
void TT_generator::Setup_processes()
{
  f_files[1]->cd();
  //TDirectory init_histos("init_histos","init_histos");
  //init_histos.cd();

  for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
    if (!f_params->f_processes_on[iProcess]) continue;
    for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
      //if (!f_flux_bins_on[iBin]) continue;
      f_rate_est[iProcess][iBin]=1.0/(f_N_Processes*f_N_rate_bins);
      f_drawer  [iProcess][iBin]=new TT_drawer(f_params,f_nucleus,f_process[iProcess],iBin);
      f_drawer  [iProcess][iBin]->Init_randomly();
      f_drawer  [iProcess][iBin]->Compute_integral();
      f_drawer  [iProcess][iBin]->f_init_histo->Write();
    }
  }
  f_files[1]->cd();
  std::cout << "Generator initialized." << std::endl;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_generator::Generate_events()
{
  f_keep_going=1;
  f_event_display=1;
  Double_t EVENT_LIMIT=pow(2,53);

  while (f_keep_going) {
    Double_t total_accepted_old=f_total_accepted;

    for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
      for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
	//if (!f_flux_bins_on[iBin]) continue;
	if ((!f_params->f_processes_on[iProcess]) || (!f_keep_going) || Reject_process(iProcess,iBin)) continue;

	Bool_t got_one=0;
	while (!got_one) {
	  f_total_points[iProcess][iBin]++;
	  if (f_total_points[iProcess][iBin]==EVENT_LIMIT) {
	    std::cerr << Form("EVENT LIMIT reached for iProcess=%d,iBin=%d",iProcess,iBin) << std::endl;
	    return kFALSE;
	  }
	  got_one=f_drawer[iProcess][iBin]->Draw_point();
	  if (f_drawer[iProcess][iBin]->f_DIE_NOW) return kFALSE;
	}
	
	f_accepted_points[iProcess][iBin]++;
	f_total_accepted++;
	if (f_total_accepted==EVENT_LIMIT) {
	  std::cerr << "EVENT LIMIT reached overall" << std::endl;
	  return kFALSE;
	}
	
	f_keep_going=Keep_going();
	Report_progress();
	Fill_tree(f_drawer[iProcess][iBin]->f_event);
	Update_rates(iProcess,iBin);    
      } //end bin loop
    } //end process loop
    if (total_accepted_old==f_total_accepted) {
      std::cout << "caught in an infinite loop" << std::endl;
      return kFALSE;
    }

  } //end accepted_event loop

  f_files[0]->cd();
  f_trees[0]->Write();
  
  return kTRUE;
}


////////////////////////////////////////////////////////////////////////
Bool_t TT_generator::Keep_going()
{
  Bool_t ans=f_total_accepted<f_total_to_accept;
  for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
    for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
      //if (!f_flux_bins_on[iBin]) continue;
      ans = ans || (f_params->f_processes_on[iProcess] && f_accepted_points[iProcess][iBin]<2);
    }
  }
  return ans;
}

////////////////////////////////////////////////////////////////////////
void TT_generator::Report_progress()
{
  if (f_total_accepted==f_event_display || !f_keep_going) {
    TDatime dt;
    std::cout << "(" << dt.AsString() << ") total accepted: " << f_total_accepted << std::endl;
    f_event_display*=2;
  }
}   

////////////////////////////////////////////////////////////////////////
Bool_t TT_generator::Reject_process(Int_t process,Int_t bin)
{
  Double_t accepted_total=0;
  Double_t rate_total=0.0;
  for (Int_t jProcess=0;jProcess<f_N_Processes;jProcess++) {
    for (Int_t jBin=0;jBin<f_N_rate_bins;jBin++) {
      //if (!f_flux_bins_on[jBin]) continue;
      accepted_total+=f_accepted_points[jProcess][jBin];
      rate_total+=f_rate_est[jProcess][jBin];
    }
  }
  Double_t frac_desired=f_rate_est[process][bin]/rate_total;
  Double_t frac_actual =(f_accepted_points[process][bin]+0.0)/f_total_accepted;
  return frac_actual>frac_desired;
}

////////////////////////////////////////////////////////////////////////
void TT_generator::Fill_tree(TT_event *event)
{
  for (Int_t iComp=0;iComp<4;iComp++) {
    f_tree_k     [iComp]=event->f_k_lower     [iComp];
    f_tree_kprime[iComp]=event->f_kprime_lower[iComp];
    f_tree_q     [iComp]=event->f_q_lower     [iComp];
    f_tree_p     [iComp]=event->f_p_lower     [iComp];
    f_tree_pprime[iComp]=event->f_pprime_lower[iComp];
  }
  f_tree_xs=event->f_dsigma_dall;
  for (Int_t iMA=0;iMA<f_tree_N_MA;iMA++) {
    f_tree_MA_weights[iMA]=event->f_MA_weights[iMA];
  }
  f_tree_mag_p=event->f_mag_p;
  f_tree_cth_pq=event->f_cos_theta_pq;
  f_tree_phi_p=event->f_phi_p;
  f_tree_process=event->f_process;
  f_tree_enuqe=event->f_enuqe;
  f_tree_Q2qe=event->f_Q2qe;
  f_tree_kappa=event->f_kappa;
  f_tree_lambda=event->f_lambda;
  f_tree_tau=event->f_tau;
  f_tree_psi=event->f_psi;
  f_trees[0]->Fill();
}

////////////////////////////////////////////////////////////////////////
void TT_generator::Update_rates(Int_t process,Int_t bin)
{
  Double_t eff=(f_accepted_points[process][bin]+0.0)/f_total_points[process][bin];
  f_rate_est[process][bin]=eff*f_drawer[process][bin]->f_integral;
  f_rate_err[process][bin]=sqrt(eff*(1.0-eff)/f_total_points[process][bin])*f_drawer[process][bin]->f_integral;
}

////////////////////////////////////////////////////////////////////////
void TT_generator::Make_graphs()
{
  Double_t processes[f_N_Processes];
  Double_t rates    [f_N_Processes];
  Double_t rate_errs[f_N_Processes];

  for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
    processes[iProcess]=iProcess;
    rates    [iProcess]=0.0;
    rate_errs[iProcess]=0.0;
    if (!f_params->f_processes_on[iProcess]) continue;

    for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
      //if (!f_flux_bins_on[iBin]) continue;
      rates    [iProcess]+=f_rate_est[iProcess][iBin];
      rate_errs[iProcess]+=pow(f_rate_err[iProcess][iBin],2);
      Double_t eff=(f_accepted_points[iProcess][iBin]+0.0)/f_total_points[iProcess][iBin];
      printf("process: %d, bin: %d, eff.: %6.4g, rate: %6.4g, rate_error: %6.4g\n",f_process[iProcess],iBin,eff,f_rate_est[iProcess][iBin],f_rate_err[iProcess][iBin]);
    }
    rate_errs[iProcess]=sqrt(rate_errs[iProcess]);
  }
  
  TGraph graph_rate      (f_N_Processes,processes,rates);
  TGraph graph_rate_error(f_N_Processes,processes,rate_errs);
  
  graph_rate.SetName("ccqe_rate");
  graph_rate_error.SetName("ccqe_rate_error");
  
  f_files[1]->cd();
  graph_rate.Write();
  graph_rate_error.Write();
}

////////////////////////////////////////////////////////////////////////
void TT_generator::Write_shuffled_tree()
{
  Long64_t Long64_Entries=f_trees[0]->GetEntries();
  UInt_t UInt_Max_Entries=-1;
  Long64_t Long64_Max_Entries=UInt_Max_Entries;

  if (Long64_Entries>Long64_Max_Entries) {
    std::cerr << "The tree has " << Long64_Entries << " entries.  The maximum allowed is " << Long64_Max_Entries << "." << std::endl;
  }
  UInt_t NEntries=Long64_Entries;
  Double_t mem_usage=NEntries*4.0/1.0e6;
  printf("Shuffling will take %6.3g MB of memory.\n ",mem_usage);

  UInt_t *shuffled_indices=new UInt_t[NEntries];

  for (UInt_t iEntry=0;iEntry<NEntries;iEntry++) {
    shuffled_indices[iEntry]=iEntry;
  }

  ////shuffle indices
  //for (UInt_t iEntry=0;iEntry<NEntries;iEntry++) {
  //  UInt_t offset=f_params->f_rand.Integer(NEntries-iEntry);
  //  UInt_t temp=shuffled_indices[iEntry];
  //  shuffled_indices[iEntry]=shuffled_indices[iEntry+offset];
  //  shuffled_indices[iEntry+offset]=temp;
  //}

  //restrict to the number desired
  for (UInt_t iEntry=0; iEntry<NEntries;iEntry++) {
    f_trees[0]->GetEntry(shuffled_indices[iEntry]);
    f_trees[1]->Fill();
  }

  delete [] shuffled_indices;

  f_files[1]->cd();
  f_trees[1]->Write();
  f_files[1]->Close();
  f_files[0]->Close();
}
