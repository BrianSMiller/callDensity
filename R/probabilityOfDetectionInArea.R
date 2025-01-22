# library(MASS, quietly = T, warn.conflicts = F)     # for the mvrnorm

#'  Simulation to predict average probability of detection
#' The simulation has several steps:
#' (1) SNR regression and TL data now passed in via function call
#' (a)Thin original TL data so have a 100m range step, not 20m
#'    The resolution given in the propagation model output may be too fine a
#'    resolution, and will increase the code run time.
#' (2) GAM for the SNR regression passed in via function call -
#'    initial model does not need to be in the loop
#'
#' (3) The outer loop (to be run 1000 times): (a) Generate a possible mean and
#'    standard deviations for source level using a parametric bootstrap approach (b)
#'    Generate a possible mean and standard deviations for noise level using a
#'    parametric bootstrap approach (c) Generate a set of coefficients for the
#'    detector characterisation curve
#'
#' (4) For each combination of outer loop parameters: For each column of TL
#'    lookup table : (a) draw ~10000 source level values from a distribution
#'    generated in outer loop and apply to all of the virtual calls along the
#'    transect (b) draw ~10000 noise level values from a distribution generated
#'    in outer loop and apply to all of the range steps along the transect (c)
#'    using the TL values, calculate the RL given the SL for each virtual call
#'    (d) calculate the SNR of each virtual call, given the assigned noise
#'    levels (e) predict prob(detect) given the SNR for all SNR values (~10000
#'    per profile)
#'
#' (5) for all p(det), calculate a weighted average (weighting by
#'    distance) Weight each prob(detect) by the range (i.e. multiply the p(det) by
#'    the range) Then divide all weighted p(dets) for a transect by the sum of ALL
#'    ranges across each transect. Having run the outer loop 1000 times, have 1000
#'    weighted-average p(det) per transect
#' (6) Final p(det) calculations
#'  (a) Calculate a transect-specific mean p(det) for each transect
#'  (b) Calculate an overall p(det) for the whole simulation
#'  (c) Save results in a results table
#' (7) Variance calculations
#'  (a) Calculate the standard error of transect-specific p(det) values
#'  (b) Calculate the standard deviation of transect-specific p(det) values
#'  (c) Save results in a table and print the table for future calculations
#'
#' @param res.1 SNR-Detection function from fitSNRdetectionFunc. Must be an
#'   object of class GLM, GAM, SCAM, or VGLM that can be passed to functions
#'   coef() and vcov()
#' @param SL Distribution of source levels (data.frame with columns mean,sd)
#' @param TLlookup Transmission loss lookup table with range in the first column
#'     and TL in dB along an azimuth for columns 2-n
#' @param NL Noise level distribution (data.frame with columns mean,sd)
#' @param transectFile File name base for writing transect outputs
#' @param simResultsFile File name base for writing simulation results
#' @param paFile File name for writing probability of detection in the area, pa
#' @param output.resolution.m Spatial resolution of density/transect output in m
#' @param outerloop Number of iterations in outerloop (default 1000)
#' @param truncationDistance scalar or matrix of truncation distances.
#'    If a matrix is provided, then the dimensions should be 1xN with N being
#'    the number of transects
#' @param snrTruncationThreshold scalar SNR in dB below which the probability of
#'    detection will be forcibly set to zero.
#'
#' @export
pDetInArea <-
  function(res.1, SL, TLlookup,  NL, # Sonar equation inputs
           transectFile, simResultsFile, paFile, # file output names
           output.resolution.m = 100, outerloop = 1000,
           truncationDistance=(Inf), snrTruncationThreshold=-Inf) {
    ### A5.2 Code used for the Monte Carlo Simulation

  # Written (started) 04 Jul 10 by Danielle Harris
  # This version edited 01 Sept 11
  # This version edited by FRC up to apr21
  # Latest version adapted by BSM on 12 Oct 2022
  #   INPUT/OUTPUT changes
  #   -Bring all user-selectable parameters
  #   -Change names of output files based on siteCode and season
  #   -Replace all hard-coded "magic" numbers with correct, and full dimensions
  #     from input files
  #   -Load NL from a separate csv that contains datetime, t, and noise level,
  #     with nl being in the same units & frequency band as SL (dB re 1 uPa RMS).
  #     Previously calibrated NL and RL were loaded from the capture-history file.
  #     This provides the option to use a NL distribution that is representative
  #     throughout the study period, without being coupled to detections.
  #   COSMETIC changes
  #   -Lower 'verbosity' by removing most calls to str()/head()/print()/etc.
  #   -Layout 8 transects in 3x3 grid with respect to true north
  # Further modified by BSM on 12 Dec 2022:
  #   -Move all user-adjustable parameters into parameter data.frame, p.
  #     This includes all input and output file names, as well as numeric params.
  # Further modified by BSM in Apr 2023 to make into a function
  #   - All parameters now adjustable in function call
  # 2024-08-08 -
  #   - Added option for truncating each radial by distance.
  #     Probabilities beyond the truncation distance along each radial will be
  #     set to NA. This provides a means of excluding distant portions of a
  #     radial, that are on land or are in the acoustic shadow of islands.
  #   - Added option for "truncating" the SNR-detection function below a
  #     specified SNR. Below the truncation SNR probabilities will be set to 0.
  #     This option facilitates the use of SNR-detection functions that are
  #     mostly well behaved on some interval, but do not go all the way to zero
  #     over the range of SNR that can occur in the simulation. Such situations
  #     are most likely to arise with VGAM-based SNR detection functions, e.g.
  #     those generated by the VGAM package or callDensity::fitSNRvglm.

#*******************************************************************************
  #STEP 1(a) - Thin TL data

  # Number of measurements i.e. number of rows
  numTLrow<-dim(TLlookup)[1]

  # Use all TL columns from TLlookup file (1st column is range)
  no.profiles <- dim(TLlookup)[2]-1;

  # Check that transectTruncationDistances have same dimensions as TL
  numTruncDist <- length(truncationDistance)
  if  (is.null(numTruncDist) | numTruncDist == 1) { # Default or single value
    truncationDistance = rep(truncationDistance, no.profiles)
  }

  if (dim(truncationDistance)[2] != no.profiles){
    errorText <-"Error: Number of Transect Truncation Distances must be 1 or
                        the same as the number of TL transects. "
    stop(errorText)
  }

  # Convert to matrix:
  allTLlookup_h <- as.matrix(TLlookup) # to convert to matrix, with no header
  allTLlookup <- matrix(allTLlookup_h, ncol = ncol(allTLlookup_h),
                        dimnames=NULL)

  #Now thin data, so that there is only a TL value every 100m; Presently this
  #step just uses every nth measurement based on the spacing, but could instead
  #consider averaging or low-pass filtering and resampling to smooth the profile
  rangeStep.original <- (TLlookup$range_m[2]-TLlookup$range_m[1])
  spacing <- output.resolution.m/rangeStep.original

  #take 10000 measurements, and include 100m and 1000km
  subset<-seq(0,numTLrow,spacing)
  allTLlookupsubset<-allTLlookup[subset,]
  dim(allTLlookupsubset) #so now have approximately 10000 TL values per profile

  # number of measurements i.e. number of rows
  numTLrowsubset<-dim(allTLlookupsubset)[1]

  #delete the range column (first column) of allTLlookup - for matrix of TL only
  #will be used later in the loop
  allTLsubset<-allTLlookupsubset[,2:(no.profiles+1)]
  allTrunc<-matrix(data=0, nrow=dim( allTLsubset)[1], ncol=dim(allTLsubset)[2])
  for (i in 1:no.profiles){
    na.index <- allTLlookupsubset[,1]>truncationDistance[i]
    allTrunc[na.index,i]<-rep(NA)
  }


  #remove TLlookup to save memory
  rm(TLlookup)

  gc() #garbage collect function to free up memory after removing an object
  #*******************************************************************************

  # STEP 2
    # Reference: https://smolski.github.io/livroavancado/reglog.html#o-modelo
  # "Dispersion of" events "and" non-events "(binary response variable)":


  #*******************************************************************************
  # STEP 3
  #create a summary results matrix

  #Number of model parameters: 2 columns for SL mean/st.dev, 2 for NL$mean/st.dev,
  # plus 2 or more  columns for coef values
  no.core.params <- 4;

  no.coef <- length(coef(res.1))
  no.model.params <- no.core.params  + no.coef

  results1000sim<-matrix(NA,outerloop,no.model.params)

  # create a large results matrix for all pdets drawn from simulation
  # this is needed to calculate detection function plots
  # the dimensions of this will be 10000 columns by no.profiles*1000 rows
  # (there are 10000 range steps per profile, and there are no.profiles profiles
  # and 1000 iterations)
  resultsallpdets<-matrix(NA,numTLrowsubset,(no.profiles*outerloop))
  dim(resultsallpdets)

  #START OF THE OUTER LOOP
  for (i in 1:outerloop){
    #***************************************************************************
    #STEP 3(a) - bootstrap from the assumed SL distribuation and calculate new
    #mean and standard deviation. The number of samples bootstrapped reflects
    #the sample size used to create the original distribution (defined outside
    #the loop)
    SLdist<-stats::rnorm(SL$sampleSize, SL$mean, SL$sd)

    #calculate mean & stdev of set of resamples and save in main results matrix
    results1000sim[i,1]<-mean(SLdist)
    results1000sim[i,2]<-stats::sd(SLdist)
    #***************************************************************************
    #STEP 3(b) - bootstrap from the assumed noise distribuation and calculate
    #new mean and standard deviation. The number of samples bootstrapped
    #reflects the sample size used to create the original distribution (defined
    #outside the loop)

    # NL now calculated outside this function and specified as input parameter
    NLdist<-stats::rnorm(NL$sampleSize, NL$mean, NL$sd)

    #calculate mean & stdev of set of resamples and save in main results matrix
    results1000sim[i,3]<-mean(NLdist)
    results1000sim[i,4]<-stats::sd(NLdist)
    #***************************************************************************
    #STEP 3(c) - produce a random realisation of the detector characterisation
    #curve The random part involves taking a random set of the coefficients for
    #the model - these will be used later
    # vcov adopted in the place of Vp for GAM
    if (any(class(res.1)=='scam') | any(class(res.1)=='gam') ){
      br<-MASS::mvrnorm(1,as.vector(coef(res.1)),res.1$Vp)
    }else { # Class is a glm or vglm
      br<-MASS::mvrnorm(1,as.vector(coef(res.1)),vcov(res.1))
    }
    results1000sim[i,5:no.model.params]<-br # save in parameter matrix

    #********************************************************************************
    #STEP 4
    #For each combination of outer loop parameters:
    # For each stage of the TL lookup table (will be repeated 8 times, as 8
    # transects):
    #(a) draw 10000 source level values from a distribution generated in the
    #    outer loop and apply to all of the virtual calls along the transect
    #(b) draw 10000 noise level values from a distribution generated in the
    #    outer loop and apply to all of the range steps along the transect
    #(c) using the TL values, calculate the RL given SL for each virtual call
    #(d) calculate SNR of each virtual call, given the assigned noise levels
    #(e) predict the prob(detect) given the SNR for all SNR values (10 000 for
    #    each profile)

    #need some places to store the results
    allSL<-matrix(NA,numTLrowsubset,no.profiles)
    allNL<-matrix(NA,numTLrowsubset,no.profiles)
    allSNR<-matrix(NA,numTLrowsubset,no.profiles)
    allpdet<-matrix(NA,numTLrowsubset,no.profiles)

    #********************************************************************************
    # STEP 4(a)
    #draw 10000 values from the generated SL distribution and put in allSL column
    #repeat for all transects
    for (j in 1:no.profiles){
      allSL[,j]<-stats::rnorm(numTLrowsubset,results1000sim[i,1],results1000sim[i,2])
    }
    #********************************************************************************
    # STEP 4(b)
    #draw 10000 values from the generated NL distribution and put in allNL column
    #repeat for all transects
    for (j in 1:no.profiles){
      allNL[,j]<-stats::rnorm(numTLrowsubset,results1000sim[i,3],results1000sim[i,4])
    }
    #********************************************************************************
    # STEPS 4(c) and 4(d)
    #calculate SNR across the whole matrix (all 8 profiles)
    allSNR<- allSL-allTLsubset-allNL
    #********************************************************************************
    # STEPS 4(e): now for each profile, calculate predicted p(dets) for each range

    for (j in 1:no.profiles){
      # Create dataset using allSNR - this will be used to create predictions of p(dets) using the GAM
      newd<-data.frame(SNR=allSNR[,j])

      # Generate the lpmatrix for the values (if GAM) See 'mgcv' documentation
      # for more information about the predict function The lpmatrix contains
      # values which when multiplied by the coefficients, gives predicted values
      # (on the scale of the link function)
      #   Xp<-predict(res.1,newd,type="lpmatrix")
      #   #multiply the lpmatrix by the randomly chosen coefficient values. Use
      #   #inv.logit to put the predicted values on the scale of the response
      #   #function (i.e p(det) values)
      #   predmatrix<-boot::inv.logit(Xp %*% br)

      # Code above only seems to apply for standard MGCV GAMs Since we might
      # have a GAM, GLM, or SCAM we need to do something different. We will swap
      # the resampled coefficents into the res.1 model, then predict on the
      # scale of the response using the new dataset
      res.1.newcoeff<-res.1

      # VGLM/VGAMs handled a bit differently than GAMs, SCAMs, and GLMs
      # The differences here go beyond just syntax (tho this is different too).
      # VGAMs predict the response of all observers as a matrix with column
      # for each observer. So here we need to know which column to use.
      # WhichObserver will be included in the model under extra if it has been
      # fit via callDensity::fitSNRvglm. Otherwise, the default behaviour here
      # is to use the response of the last observer.
      if (any(class(res.1)=='vglm')){
        res.1.newcoeff@coefficients<-br
        predmatrix<-VGAM::predict(res.1.newcoeff,newdata=newd,type="response")

        index = ifelse(is.null(res.1.newcoeff@extra$whichObserver),
               dim(predmatrix)[2],
               which(colnames(predmatrix)==res.1.newcoeff@extra$whichObserver) )

        predmatrix<- predmatrix[,index]
      }else {
        res.1.newcoeff$coefficients<-br
        predmatrix<-predict(res.1.newcoeff,newd,type="response")
      }

      # set pdet<-0 when SNR is below the truncation threshold
      predmatrix[which(newd< snrTruncationThreshold)]<-0

      # set pdet<-NA when the distance along radial[,j] >= truncDistance[,j]
      predmatrix <- predmatrix + allTrunc[,j]

      # Check that: 0 <= predmatrix/probabilities <= 1 or else force them to be
      # This is needed because some models (looking at you VGAMs) produce
      # probabilities below zero or above 1 (e.g. -Inf and Inf)
      if (any( predmatrix < 0, na.rm = TRUE)){
        warning("Probabilities <0 encountered and set to zero")
        predmatrix[which(predmatrix < 0)]<-0
      }

      if (any( predmatrix > 1, na.rm = TRUE)){
        warning("Probabilities >1 encountered and set to 1")
        predmatrix[which(predmatrix > 1)]<-1
      }

      # Save the prob(det) in allpdet
      allpdet[,j]<-predmatrix

    }

    # Once pdet full with 5 columns, save allpdet in resultsallpdets
    resultsallpdets[,((no.profiles*i)-(no.profiles-1)):(no.profiles*i)]<-allpdet

    # print(i) # To keep track of the simulation

  }


  ##############################################################
  #SIMULATION ENDED - NOW PRODUCE THE RESULTS
  #Save originally selected parameters for the simulation
  utils::write.table(results1000sim,file=simResultsFile,row.names=F,col.names=F)
  #********************************************************************************

  # STEPS 5 - for all p(det), calculate a weighted average (weighting by distance)

  # Weight each prob(detect) by the range (i.e. multiply the p(det) by the range)
  # Then divide all weighted p(dets) for a transect by the sum of ALL ranges across
  # each transect.
  # want to save the 1000 weighted mean p(dets) for each transect
  # a bit easier to handle than resultsallpdets
  # alltransectspdetsweight1000<-matrix(NA,8,1000)
  alltransectspdetsweight1000<-matrix(NA,no.profiles,outerloop)
  #repeat all steps for each transect
  for (i in 1:no.profiles){
    # extract the data (10000x1000) for each transect
    transectallpdets<-resultsallpdets[,seq(i,dim(resultsallpdets)[2],no.profiles)]
    #check dimensions are ok
    #dim(transectallpdets)
    # Now want to weight the pdets using the ranges
    # define a matrix to save the results - will be numTLrowsubset x outerloop
    # for each transect
    transectallpdetsweight<-matrix(NA,numTLrowsubset,outerloop)
    #For each range step, multiply the 1000 x p(det) values by the range

    for (j in 1:numTLrowsubset){
      transectallpdetsweight[j,]<-allTLlookupsubset[j,1]*transectallpdets[j,]
    }
    #calculate average pdet for each transect:
    #So now need to divide each column (1000 iterations) by the sum of the ranges
    transectpdetsweight1000<-c()
    for (j in 1:outerloop){
      # Debugging code to check for NA
      # browser(expr = {any(is.na(transectallpdetsweight[,j]))})
      include <- !is.na(transectallpdetsweight[,j])
      transectpdetsweight1000[j]<-
        sum(transectallpdetsweight[include,j])/sum(allTLlookupsubset[include,1])
    }
    #This gives 1000 weighted average probability of detections for each iteration
    #Specific to one of the transects over 1000KM
    #Save in alltransectspdetsweight1000 matrix
    alltransectspdetsweight1000[i,]<-transectpdetsweight1000
  }
  #So now have a results matrix with no.profilesx1000 average weighted p(dets)
  #********************************************************************************
  # STEP 6 - Final p(det) calculations
  # STEP 6(a) Calculate a transect-specific mean p(det) for each transect
  transectavgweightpdets<-rowMeans(alltransectspdetsweight1000, na.rm = T)
  # STEP 6(b) Calculate an overall p(det) for the whole simulation

  # Shouldn't be any NA's here now, but we set na.rm=T anyway (maybe a whole transect is NA)
  # now that transects might be different lengths, we need to weight by length
  # when averaging across them.
  # overallpdet<-mean(transectavgweightpdets, na.rm=T)
  overallpdet<-weighted.mean(transectavgweightpdets, w=truncationDistance^2, na.rm=T)

  # STEP 6(c) Create a results matrix to hold p(det) values and st.error/st.dev values

  pdetsandvar<-matrix(NA,no.profiles+1,2)
  pdetsandvar[1:no.profiles,1]<-transectavgweightpdets
  pdetsandvar[no.profiles+1,1]<-overallpdet
  #********************************************************************************
  # STEP 7 - Variance calculations
  # STEP 7(a) Need to calculate the standard error of the transect-specific p(det) values
  st.errorPt<-sd(transectavgweightpdets,na.rm = T)/sqrt(no.profiles)
  # STEP 7(b) Need to calculate the standard deviation of the transect-specific p(det) values
  st.devPt<-c()
  for (i in 1:no.profiles){
    st.devPt[i]<-sd(alltransectspdetsweight1000[i,],na.rm = T)
  }
  # STEP 7(c) save the results in the results matrix and print out
  pdetsandvar[1:no.profiles,2]<-st.devPt
  pdetsandvar[no.profiles+1,2]<-st.errorPt
  utils::write.table(pdetsandvar,file=paFile, row.names=F,
                     col.names = c("Mean","SD"))


  ##############################################################

  #PLOTTING UNWEIGHTED P(DET) AGAINST RANGE TO DETERMINE APPROPRIATE VALUE FOR W
  #(These plots are loosely referred to as detection function plots)
  #for each profile, extract the correct columns of resultsallpdets

  # windows(record=TRUE)

  # par(mar=c(4,4,3,2))

  #         NW  N  NE  W   m  E   SW   S  SE
  #          2, 3, 6,  9,  8,  7,  4   1,  5
  #          1, 2, 3,  4,  5,  6,  7,  8   9
  if (season=='year' & FALSE){
    profile<-c(0,45,90,135,180,225,270,315)
    # windows(record=TRUE)
    # windows()
    # for(i in c(8, 1, 2,  7,  8,  3,  6,  5,  4)){

    graphics::layout(matrix(c(8, 1, 2, 7, 9, 3, 6, 5, 4), ncol=3, byrow=TRUE))
    for (i in seq(1,no.profiles,no.profiles/length(profile)) ){
      #extract the transect specific data (10000x1000) for each transect
      transectallpdets<-resultsallpdets[,seq(i,dim(resultsallpdets)[2],no.profiles)]
      #check dimensions are ok
      dim(transectallpdets)
      #save these results, so can recreate the work anytime
      # write.table(transectallpdets,file=paste("transectallpdets_",i,".txt"),row.names=F,col.names=F)
      #calculate the mean pdet (unweighted) for each range step
      transectavgdetfunc<-rowMeans(transectallpdets)
      #check that the dimensions are 10000 x 1
      length(transectavgdetfunc)
      #plot to take a look
      plot( (allTLlookupsubset[,1]/1000),transectavgdetfunc, cex=0.1, pch = 20, xlab =
              "Range(km)", ylab = "P(det)",
            main = paste(profile[((i-1)/3)+1],'degrees') )
    }

  }
  range_m <- allTLlookupsubset[,1]
  #Calculate an average plot, averaging across p(det) for each range step for all transects
  #Use the main results matrix, resultsallpdets
  combineddetfunc<-matrix(NA,no.profiles,dim(resultsallpdets)[1]) #10000
  for (i in 1:no.profiles){
    #extract the transect specific data (10000x1000) for each transect
    transectallpdets<-resultsallpdets[,seq(i,dim(resultsallpdets)[2],no.profiles)]
    #check dimensions are ok
    dim(transectallpdets)
    #calculate the mean pdet (unweighted) for each range step
    transectavgdetfunc<-rowMeans(transectallpdets, na.rm=T)
    #check that the dimensions are 10000 x 1
    length(transectavgdetfunc)
    combineddetfunc[i,]<-transectavgdetfunc
  }
  allDetFunctions <- allTLlookupsubset
  allDetFunctions[,2:(no.profiles+1)]<-t(combineddetfunc)
  allDetFunctions<- data.table::as.data.table(allDetFunctions)
  names(allDetFunctions)<-dimnames(allTLlookup_h)[[2]]
  system.time(data.table::fwrite(allDetFunctions, file=transectFile) )
  #average across transects to produce an overall det function
  averagecombineddetfunc<-colMeans(combineddetfunc, na.rm=T)

  #plot to take a look
  # plot((allTLlookupsubset[,1]/1000),averagecombineddetfunc, cex=0.1, pch = 20, xlab =
  #        "Range(km)", ylab = "P(det)", main = "All transects")


  ### If W ok, FINAL ###

  ##############################################################

  # Once have looked at the results, decide an appropriate truncation distance, w
  #Redo Steps 5-7 using the different truncation distance
  #NB: The truncation distance may well be different for each source level
  #In this analysis, a truncation distance of 500 km could be applied to the
  #simulation run with Sri Lankan call source levels.
  #Therefore Steps 5-7 were redone with a distance of 500km (below)
  #In the simulation using Antarctic call source levels, the truncation distance was left at 1000km.

}
