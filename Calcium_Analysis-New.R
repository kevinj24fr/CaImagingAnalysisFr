
############################################################################################## 
#       Welcome to the calcium analysis script provided by the MILO laboratory Freiburg      # 
#       Author: Kevin Joseph                                                                 #
#       1. Segment your cells in FIJI, and use the multimeasure function in FIJI to obtain   #
#       raw calcium traces.                                                                  #
#       2. Make sure to include 3 background spots from spatially nonoverlapping locations   #
#       3. Write the data into an excel file, with the trace columns named Cell_ and        #
#          background columns as BG_                                                        #
##############################################################################################
#
#
#

calciumcorrection=function(rawtrace){

#Background corrected Traces
background=rawtrace[,grepl("BG",names(rawtrace))]
background$average=(background$BG_1+background$BG_2+background$BG_3)/3
bgcorrectedtraces=rawtrace[,grepl("Cell",names(rawtrace))]-background$average

#Normalized Traces
normalizationvalues=vector(length = ncol(bgcorrectedtraces))
for (i in 1:ncol(bgcorrectedtraces)){
normalizationvalues[i]=mean(bgcorrectedtraces[,i])
}
normalizedtraces=data.frame(matrix(ncol = ncol(bgcorrectedtraces), nrow = nrow(bgcorrectedtraces)))
for (i in 1:ncol(bgcorrectedtraces)){
  normalizedtraces[,i]=bgcorrectedtraces[,i]/normalizationvalues[i]
}
colnames(normalizedtraces)=colnames(bgcorrectedtraces)
# Baseline Drift corrected Traces using a loess regression


loessmatrix=data.frame(matrix(NA, nrow = nrow(bgcorrectedtraces), ncol = ncol(bgcorrectedtraces)-1))
#loessmatrix=data.frame(matrix(NA, nrow = nrow(normalizedtraces), ncol = ncol(normalizedtraces)-1))
Time=(1:nrow(bgcorrectedtraces))
for (i in 1:ncol(bgcorrectedtraces))
{
  Testintermediate=as.matrix(bgcorrectedtraces[,i])
  #plot(x=Time, y=Testintermediate, type="l")
  loessintermediate=loess(Testintermediate ~ Time, span = 0.45)
  smoothed=predict(loessintermediate)
  #points(y=smoothed-2, x=Time, type="l", col="red")
  loessmatrix[,i]=smoothed-2
}
correctedmatrix=bgcorrectedtraces[,1:ncol(bgcorrectedtraces)]-loessmatrix
correctedmatrix$Time=Time
return(correctedmatrix)
}

correctedtraces=calciumcorrection(rawtrace = readxl::read_excel("~/Documents/Documents – Kevin’s MacBook Pro/Microglia Project/Data/Microglia Project/Kevin/233+MG/B3.xlsx"))

#write.csv(correctedtraces,"A2_corrected.csv")

#plot(correctedmatrix$Cell_14,type="l")
#points(correctedmatrix$Cell_45,type="l", col="red")


#plot(rawtrace$Cell_56,type="l")



#correctedggplotdf=reshape2::melt(correctedtraces, id.vars="Time")

#ggplot(correctedggplotdf, aes(x=Time,y=value))+
#  geom_line(aes(color=variable))+
#  theme_bw()+theme(legend.position = "none")


findthreshold=function(inputdf, thresholdmultiplier){
  thresholdingmatrix=inputdf[,!grepl("Time",names(inputdf))]
  thresholdvalues=vector(length = ncol(thresholdingmatrix))
    for (i in 1: ncol(thresholdingmatrix)){
    thresholdvalues[i]=thresholdmultiplier*sd(thresholdingmatrix[,i])
    }
  return(thresholdvalues)
}

thresholdvector=findthreshold(correctedtraces, thresholdmultiplier = 2)

findpeaks=function(inputdf, thresholdvector){
  peakmatrix=inputdf[,!grepl("Time",names(inputdf))]
  detectedspikes<-purrr::map(.x=1:ncol(peakmatrix), .f=function(i){
    out<-pracma::findpeaks(peakmatrix[,i],nups = 1, 
                                          ndowns = 1, 
                                          peakpat = NULL, 
                                          minpeakheight = -Inf, 
                                          minpeakdistance = 3, 
                                          threshold = thresholdvector[i], 
                                          sortstr = FALSE)
    if(!is.null(out)){
      out<-as.data.frame(out)
      names(out) <- c("Amp", "Peak_time", "Peak_start", "Peak_end")
    }
    
    return(out)
    
    
  })
  names(detectedspikes)=colnames(peakmatrix)
  
  return(detectedspikes)
  }

peaks=findpeaks(correctedtraces, thresholdvector)

#list<-peaks[[111]]

peaks<-purrr::map(.x=1:length(peaks), .f=function(i){
  list<-peaks[[i]]
  if(!is.null(list)){
    
    comp<-purrr::map_df(.x=1:nrow(list), .f=function(ii){
      
      #isolate the raw signal
      peak.val<-correctedtraces[,i][list$Peak_start[ii]:list$Peak_end[ii]]
      val<-data.frame(
        peakwidth=list$Peak_end[ii]-list$Peak_start[ii],
        skewness=moments::skewness(peak.val),
        kurtosis=moments::kurtosis(peak.val),
        risetime=list$Peak_time[ii]-list$Peak_start[ii],
        falltime=list$Peak_end[ii]-list$Peak_time[ii]
      
      )
      
    })
    
    list<-cbind(list, comp)
    
    list$Peak_ID<-paste0(names(peaks)[i], "_P", 1:nrow(list))
    list$Cell<-names(peaks)[i]
    
    return(list)
  }
})

spike.df_B3<-purrr::map_df(.x=peaks, .f=function(i){if(!is.null(i)){i}})

spike.df_B3$Experiment=paste0("MG+233")
spike.df_B3$Trial=paste0("Trial3")

spike.df_Control=rbind(spike.df_A1,spike.df_A2)
spike.df_MGcoculture=rbind(spike.df_B1,spike.df_B2,spike.df_B3)

spike.df_all=rbind(spike.df_Control,spike.df_MGcoculture)
spike.df_all[spike.df_all == "MG+233"] <- "MG_Coculture"



ggplot(spike.df_all)+geom_point(aes(x=skewness, y=kurtosis, color=Experiment))
ggplot(spike.df_all)+geom_point(aes(x=peakwidth, y=Amp, color=Experiment))

ggplot(spike.df_all)+geom_boxplot(aes(x=Experiment, y=Amp))

ggpubr::stat_compare_means(comparisons = list(c("Control","MG_Coculture")))



mat<-spike.df_all[1:9]
for(i in 1:ncol(mat)){mat[,i]<-scales::rescale(mat[,i], c(0,1))}

ggplot(mat)+geom_point(aes(x=log(skewness), y=log(kurtosis), color=spike.df_all$Experiment))

ggplot(spike.df_all,aes(x=Experiment, y=Amp))+geom_boxplot(outlier.shape = NA)+ggpubr::stat_compare_means(comparisons = list(c("Control","MG_Coculture")))

umap<-umap::umap(mat)
colnames(umap$layout)=c("UMAP_1","UMAP_2")

ggplot(data.frame(umap$layout))+
  geom_point(aes(x=UMAP_1, y=UMAP_2, color=mat$Amp))+
  theme_classic()

pca<-pcaMethods::pca(mat, nPcs=3)

ggplot()+
  geom_point(data=pca@scores %>% as.data.frame(), mapping=aes(x=PC1, y=PC2+PC3, color=spike.df_all$Experiment))+
  theme_classic()+
  Seurat::NoLegend()

a=table(spike.df_Control$Cell)
b=table(spike.df_MGcoculture$Cell)
table(spike.df_MGcoculture$Cell)

plot(density(a),col="blue", type="l")
points(density(b), col="red", type="l")


prcomp(na.omit(mat), scale=TRUE)$x %>%
  as_tibble() %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(aes(color=spike.df_all$Amp,shape=spike.df_all$Experiment))

correctedggplotdf=reshape2::melt(correctedtraces, id.vars="Time")

ggplot(correctedggplotdf, aes(x=Time,y=value))+
  geom_line(aes(color=variable))+
  theme_bw()+theme(legend.position = "none")


testdfspike.all=reshape2::melt(spike.df_all, id.vars=c("Experiment","Trial"), measure.vars="Cell")
