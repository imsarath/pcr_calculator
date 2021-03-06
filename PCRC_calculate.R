###
# Rscript for prostate cancer risk calculation
###

# calDCRwithMRI - calculation for detectable cancer risk patients with MRI and prior biopsy
calDCRwithMRI <- function(psa, dre_outcome, volume, age, pirads) {
  preBiopsy = 1.0
  psaLog = log2(psa)
  volumeLog = log2(volume)
  pirad = as.character(pirads-1)
  fiPirad = ifelse(pirad == '2', 0.08258, 
              ifelse(pirad == '3', 1.40371, 
                ifelse(pirad == '4', 1.91967,
                      0.0)))
  val_prior = -2.62636 + (-0.745) + ((((-1.470 - 0.677 * preBiopsy + (0.576 - 0.423 * preBiopsy) * (psaLog -2)-1.043 * (volumeLog - 5.5) + 0.68 * dre_outcome))*1.4857+2.3426)) * 0.78749+age*0.02467+ fiPirad
  
  PCRC_phi = 1/(1+exp(-(val_prior)))
  
  return( round(PCRC_phi * 100))
}

# calSCRwithMRI - calculation for significant cancer risk patients with MRI and prior biopsy
calSCRwithMRI <- function(psa, dre_outcome, volume, age, pirads) {
  preBiopsy = 1.0
  psaLog = log2(psa)
  volumeLog = log2(volume)
  pirad = as.character(pirads-1)
  fiPirad = ifelse(pirad == '2', 0.94909,
                   ifelse(pirad == '3', 2.67536, 
                          ifelse(pirad == '4', 3.05243,
                                 0.0)))
  
  valPriorSig = -5.19208 + (-0.67) +((((-3.489 - 1.136 *  preBiopsy + (1.075 - 0.434 * preBiopsy) *( psaLog-2)- 1.501 * (volumeLog - 5.5) + 1.311 * dre_outcome)) * 0.9420+ 2.1452))* 0.74027+ age * 0.04286 + fiPirad;
  ProRC_SCR = 1/(1+exp(-(valPriorSig)));
  
  return( round(ProRC_SCR*100))
}

# calDCRwithMRInoPrebio - calculation for detectable cancer risk score for patients with MRI and no prior biopsy
calDCRwithMRInoPrebio <- function(psa, dre_outcome, volume, age, pirads ) {
  psaLog = log2(psa)
  volumeLog = log2(volume)
  pirad = as.character(pirads-1)
  fiPirad = ifelse(pirad == '2', -0.38169, 
                   ifelse(pirad == '3', 0.61290, 
                          ifelse(pirad == '4',  2.01363,
                                 0.0)))

  aux =  ((-1.826 + 1.024 * (psaLog - 2.0)  - 1.50 * (volumeLog - 5.4) + 0.992 * dre_outcome) * 1.139 ) + 1.325
  valPrior = -3.52647 + (-1.165) + aux * 0.79560 + age * 0.04738 + fiPirad
  
  PRC_DRE_noPreBio = 1/(1+exp(-(valPrior)))
  
  return( round(PRC_DRE_noPreBio * 100))
}

# calSCRwithMRInoPrebio - calculation for significant cancer risk score for patients with MRI and no prior biopsy
calSCRwithMRInoPrebio <- function(psa, dre_outcome, volume, age, pirads ) {
  psaLog = log2(psa)
  volumeLog = log2(volume)
  pirad = as.character(pirads-1)
  fiPirad = ifelse(pirad == '2', 0.44510, 
                   ifelse(pirad == '3', 1.47778, 
                          ifelse(pirad == '4',  2.61551,
                                 0.0)))

  aux =  ((-3.457 + 1.177 * ((psaLog - 2.0) - 1.526 * (volumeLog - 5.4) + 1.813 * dre_outcome)) * 0.718)+ 1.069
  
  valPriorSig = -4.05287 + (-1.61) + aux * 0.60973 + age * 0.04026 + fiPirad;
  
  PRC_DRE_noPreBio_Sig = 1/(1+exp(-(valPriorSig)))
  
  return( round(PRC_DRE_noPreBio_Sig * 100))
}

calDCRwithoutMRI <- function(psa, dre_outcome, volume) {
  preBiopsy = 1.0
  psaLog = log2(psa)
  volumeLog = log2(volume)
  
  valPrior = -1.470 - 0.677 * preBiopsy + (0.576 - 0.423 * preBiopsy) * (psaLog - 2.0) - 1.043 * (volumeLog - 5.5) + (0.68 * dre_outcome)
  PCRS_DRE = 1 / (1 + exp( -(valPrior)))
  
  return(round(PCRS_DRE * 100))
}

calSCRwithoutMRI <- function(psa, dre_outcome, volume) {
  preBiopsy = 1.0
  psaLog = log2(psa)
  volumeLog = log2(volume)
  
  valPriorSig = -3.489 - 1.136 * preBiopsy + (1.075 - 0.434 * preBiopsy) * (psaLog - 2.0) - 1.501 * (volumeLog - 5.5) + (1.311 * dre_outcome)
  
  PCRS_DRE_Sig = 1.0 / (1.0 + exp(-valPriorSig))
  
  return(round(PCRS_DRE_Sig * 100))
}

calDCRwithoutMRInoPreBio <- function(psa, dre_outcome, volume) {
  psaLog = log2(psa)
  volumeLog = log2(volume)
  
  valPrior  = -1.826 + 1.024 * (psaLog - 2.0) - 1.50 * (volumeLog - 5.4) + 0.992 * (dre_outcome)
  
  PCRS_DRE_noPrebio = 1 / (1 + exp( -(valPrior)))
  
  return(round(PCRS_DRE_noPrebio * 100))

}

calSCRwithoutMRInoPreBio <- function(psa, dre_outcome, volume) {
  psaLog = log2(psa)
  volumeLog = log2(volume)
  
  valPrior = -3.457 + 1.177 * (psaLog - 2.0) - 1.526 * (volumeLog - 5.4) + 1.813 * (dre_outcome)
  PCRS_DRE_noPrebio_ser = 1 / (1 + exp( -valPrior))
  
  return(round(PCRS_DRE_noPrebio_ser * 100))
}

pat_df <- read.csv('/home/chimera/Downloads/ashkan_mri.tsv', header = T, sep = '\t')

pat_df$df.DCRSwithMRI <- ifelse(pat_df$df.prevBiop == 0, calDCRwithMRInoPrebio(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol, pat_df$df.age.x, pat_df$df.Tumor1_PIRAD), calDCRwithMRI(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol, pat_df$df.age.x, pat_df$df.Tumor1_PIRAD))

pat_df$df.SCRSwithMRI <- ifelse(pat_df$df.prevBiop == 0, calSCRwithMRInoPrebio(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol, pat_df$df.age.x, pat_df$df.Tumor1_PIRAD), calSCRwithMRI(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol, pat_df$df.age.x, pat_df$df.Tumor1_PIRAD))

pat_df$df.DCRSnoMRI <- ifelse(pat_df$df.prevBiop == 0, calDCRwithoutMRInoPreBio(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol), calDCRwithoutMRI(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol))

pat_df$df.SCRSnoMRI <- ifelse(pat_df$df.prevBiop == 0, calSCRwithoutMRInoPreBio(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol), calSCRwithoutMRI(pat_df$df.psa, pat_df$df.dre, pat_df$df.dreVol))

