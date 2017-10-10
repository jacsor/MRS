plot1DSigWindows = function(ans_mrs,fdr=0.1,precision = 0.001,...) { # Plot windows significant at given Bayesian FDR level
  ans_mrs2 = ans_mrs
  c_vec = seq(0,1,by=precision)
  FDR = double(length(c_vec))
  for (i in 1:length(c_vec)) {
    FDR[i] = bFDR(ans_mrs2,c_vec[i])
  }
  FDR[length(c_vec)] = 0
  c = c_vec[min(which(FDR<fdr))]
  print(c)
  sig_window_ind = which(ans_mrs2$RepresentativeTree$AltProbs > c)
  print(sig_window_ind)
  n_windows = length(ans_mrs2$RepresentativeTree$AltProbs)
  ans_mrs2$RepresentativeTree$AltProbs[setdiff(1:n_windows,sig_window_ind)] = 0
  ans_mrs2$RepresentativeTree$EffectSizes[setdiff(1:n_windows,sig_window_ind)] = 0
  
  plot1D(ans_mrs2,...)
}
