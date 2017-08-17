bFDR = function(ans_mrs,c) {
  PMAPs = ans_mrs$RepresentativeTree$AltProbs
  return(1-sum(PMAPs[PMAPs>c])/sum(PMAPs>c))
}
