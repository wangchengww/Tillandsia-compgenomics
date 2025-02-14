seqfile = sample.txt * sequence data filename
treefile = samplephy.phy     * tree structure file name
outfile = mlc           * main result file name

 noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
verbose = 0  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
			 * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 0  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
  clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
 aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
			 * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
  model = 2
			 * models for codons:
				 * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
			 * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
			 * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
			 * 13:3normal>0

  icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
  kappa = 3  * initial or fixed kappa
fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
  omega = 1  * initial or fixed omega, for codons or codon-based AAs

fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
  alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
 Malpha = 0  * different alphas for genes
  ncatG = 8  * # of categories in dG of NSsites models

  getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

Small_Diff = .5e-6
cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
 method = 1  * Optimization method 0: simultaneous; 1: one branch a time
