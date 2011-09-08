This is a program to generate CCQE neutrino events.  It ran in 2007,
but by 2011 was unmaintaned and had many uncommitted changes.  The
changes found were committed and the code was moved to github in
Sept. 2011.

////list of ambiguities/issues
//--write some docs (implementing nucl-th/0512004v4)
//--ksi in F2?
//--Mpi+-0?
//--poles in Smith-Moniz form factors?
//--M value fixed?
//--use delta function to eliminate the integral over E?
//--overall sign of epsilon
//--epsilon upper/lower -1 proof
//--what if contraction is complex?
//--args, options, seed
//--form factors
//--nuclear magneton
//--Double_t ok?
//--negative E/limits of integration for E?
//----has implications for norm. of SM spectral function
//----E lower threshold for binding e near 0
//--virtual functions?/fix up OO stuff
//--single nucleon lower energy thershold
//--Llewellyn-Smith does not make sense at high Q2
//--test small values of eb, pf
//--W boson mass?
//--how to find p_max for spectral functions
//--5 GeV comparison
//--cmake compiler setting
//--w<e_b ok?
//--muon phi/direction of qbold
//--factor of pi between AS and their reference?
//--LS comparison
//--try nu_e events
//--fix height finding
//--check sign of E_b in TT_event::compute_enuqe_and_q2qe
//--check enuqe and Q2qe formulas
//--process==3 in setup_kinematics
//--run valgrind
//--in particular, see Init_valgrind
//--efficiency at high neutrino energy
//--e scattering
//--types in write_shuffled_tree
//--make more events than we want?
//--shuffle tree
//--put date/random label on tmp file
//--directory for init_histos
//--zero event bins (check boundaries)
//--flux quantiles
//--flux normalization
//--delete tree in TT_generator.h?
//--somehow ensure that flux histo is a TH1F
//--SM isobar?
//--anti-nu mode residual nucleus

//--sigma, omega
//http://www.slac.stanford.edu/spires/find/hep/www?rawcmd=FIND+K+SIGMA+OMEGA+MODEL&FORMAT=www&SEQUENCE=ds
//http://www.slac.stanford.edu/spires/find/hep/www?irn=1301462
//http://www.slac.stanford.edu/spires/find/hep/www?irn=1333011
//http://www.slac.stanford.edu/spires/find/hep/www?irn=1498711
