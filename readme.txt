Brent R. Hickman and Timothy P. Hubbard, "Replacing Sample Trimming with 
Boundary Correction in Nonparametric Estimation of First-Price 
Auctions", Journal of Applied Econometrics, Vol. 30, No. 5, 2015, pp. 
739-762.

All data we used are already in the Journal of Applied Econometrics
Data Archive, because we revisited the application considered in

  Sandra Campo, Isabelle M. Perrigne, and Quang H. Vuong, "Asymmetry in
  First-Price Auctions With Affiliated Private Values", Journal of
  Applied Econometrics, Vol. 18, No. 2, 2003, pp. 179-207.

See:  http://qed.econ.queensu.ca/jae/2003-v18.2/campo-perrigne-vuong/ 

In addition, we have provided three sets of replication files coded in
Matlab. These are stored in three zip files. Below we describe the
files in each zip file.  Briefly, the files in "MC-related.zip" can be
used to replicate the Monte Carlo experiments presented in our paper;
the files in "Application-related.zip" can be used to replicate the
application of the proposed estimator to the OCS data noted above; the
files in "Extensions.zip" are the core files we think a user
interested in applying boundary correction to a new environment and/or
application might want to start with.

Note that the programs can be used to execute the boundary-corrected
approach we propose as well as a trimming-based approach. We ask
researchers who use our code, in whatever capacity, to please cite
this paper.

All the Matlab files are ASCII files in DOS format. Unix/Linux users
should use "unzip -a".


Description of files in "MC-related.zip"

hh_mcstudy.m: main driver file to run and replicate Monte Carlo
  experiments

computeise.m: compute integrated squared error between kernel-smoothed
  density estimate and true density

computeise_stg1.m: compute integrated squared error between
  kernel-smoothed density estimate and true density using a grid of
  points for truth

computeise_stg1_postrun.m: variation of computeise.m called by
  hh_mcstudy.m

computeise_stg2_postrun.m: variation of computeise.m called by
  hh_mcstudy.m

evalkspdf.m: evaluate specified kernel at a grid of points

evalskpdf_num.m: evaluate function of a specified kernel at a grid of
  points

evalkspdf_denom.m: evaluate function of a specified kernel at a grid
  of points

kscdf.m: compute empirical distribution function

kspdf_bc.m: compute boundary-corrected kernel-smoothed estimate of
  density function

kspdf_no_bc.m: compute kernel-smoothed estimate of a density function

mixbetarnd.m: function to help generate (pseudo) random valuations
  from a mixture of Beta distributions

processoutput.m: takes as input the beta_mc_ex#.mat (where # =
  {1,2,3,4}) files which are output from hh_mcstudy.m to replicate
  results presented in the paper

risamint_mixbetas.m: function for computing the integral portion of
  the symmetric IPVP bid function

trim.m: trim bids within a bandwidth of the minimum and maximum
  observed bids


Description of files in "Application-related.zip"

Remember to download the OCS data from Campo, Perrigne, and Vuong as
described above before running these programs which read in a file
called Ocs702.dat.

bcappliedtoCPV.m: boundary correction applied to the CPV data

bcappliedtoCPVfcn: a function called by compare_figs.m which
  implements boundary correction on the CPV data

compare_figs.m: a file which compares the boundary correction and
  trimming-based approaches

evalkspdf.m: evaluate specified kernel at a grid of points

evalskpdf_num.m: evaluate function of a specified kernel at a grid of
  points

evalkspdf_denom.m: evaluate function of a specified kernel at a grid
  of points

jointDistbcfcn.m: compute boundary-corrected kernel-smoothed estimate
  of joint density function

kscdf.m: compute empirical distribution function

kscdf2dCPV.m: compute two-dimensional empirical distribution function
  using standard kernel

kscdf2dCPVbc.m: compute two-dimensional empirical distribution
  function using boundary-corrected kernel

kspdf_bc.m: compute boundary-corrected kernel smoothed estimate of
  density function

kspdf_no_bc_CPV.m: compute kernel-smoothed estimate of a density
  function

kspdf2dCPV.m: compute two-dimensional kernel-smoothed estimate of a
  density function

kspdf2dCPV.m: compute two-dimensional kernel-smoothed estimate of a
  density function using boundary correction

kspdf2dCPVFAST_joint.m: compute two-dimensional kernel-smoothed
  estimate of a density function using boundary correction in a
  faster way

kzroutine.m: KZ routine for boundary corrected kernel

kzker_joint.m: another KZ routine for boundary corrected kernel

replicateCPV.m: a replication file of the trimming-based approach of
  Campo, Perrigne, and Vuong

replicateCPVfcn.m: a function called by compare_figs.m which
  implements trimming on the CPV data

trim_comparison_scatter.m: a program to compare the different ways of
  trimming

trimasymmetricCPV.m: trim auctions for asymmetric APV case

trimsymmetricCPV.m: trim auctions for symmetric APV case


Description of files in "Extensions.zip"

bcgpv_driver.m: generates uniformly distributed data and applies
  boundary correction

evalkspdf.m: evaluate specified kernel at a grid of points

evalskpdf_num.m: evaluate function of a specified kernel at a grid of
  points

evalkspdf_denom.m: evaluate function of a specified kernel at a grid
  of points

kscdf.m: compute empirical distribution function

kspdf_bc.m: compute boundary-corrected kernel-smoothed estimate of
  density function


 
