calciumcorrection=function(rawtrace){
background=rawtrace[,grepl("BG",names(rawtrace))]
background$average=(background$BG_1+background$BG_2+background$BG_3)/3
bgcorrectedtraces=rawtrace[,grepl("Cell",names(rawtrace))]-background$average
normalizationvalues=vector(length = ncol(bgcorrectedtraces))
for (i in 1:ncol(bgcorrectedtraces)){
  normalizationvalues[i]=mean(bgcorrectedtraces[,i])
}
normalizedtraces=data.frame(matrix(ncol = ncol(bgcorrectedtraces), nrow = nrow(bgcorrectedtraces)))
for (i in 1:ncol(bgcorrectedtraces)){
  normalizedtraces[,i]=bgcorrectedtraces[,i]/normalizationvalues[i]
}
colnames(normalizedtraces)=colnames(bgcorrectedtraces)
loessmatrix=data.frame(matrix(NA, nrow = nrow(bgcorrectedtraces), ncol = ncol(bgcorrectedtraces)-1))
Time=(1:nrow(bgcorrectedtraces))
for (i in 1:ncol(bgcorrectedtraces))
{
  Testintermediate=as.matrix(bgcorrectedtraces[,i])
  loessintermediate=loess(Testintermediate ~ Time, span = 0.45)
  smoothed=predict(loessintermediate)
  loessmatrix[,i]=smoothed-2
}
correctedmatrix=bgcorrectedtraces[,1:ncol(bgcorrectedtraces)]-loessmatrix
correctedmatrix$Time=Time
return(correctedmatrix)
}
