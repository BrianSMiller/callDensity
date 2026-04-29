#' Monte-Carlo simulation to predict average probability of detection in area
#'
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
#' @param detFun SNR-Detection function from fitSNRdetectionFunc. Must be an
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
#' @param parallel Logical. If TRUE, the Monte Carlo outer loop runs in
#'    parallel via future.apply::future_lapply, using whatever future plan
#'    is currently active (set via future::plan() before calling). If FALSE
#'    (default), runs serially - identical behaviour to the original
#'    implementation. Requires the future.apply package to be installed
#'    when TRUE; falls back to serial with a warning if not available.
#'
#' @export
pDetInArea <-
  function(detFun, SL, TLlookup,  NL, # Sonar equation inputs
           transectFile=NULL, simResultsFile=NULL, paFile=NULL,# output names
           output.resolution.m = 100, outerloop = 1000,
           truncationDistance=max(TL[,1]), snrTruncationThreshold=-Inf,
           parallel = FALSE) {
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

  if (parallel && !inherits(detFun, "vglm")) {
    message("Note: parallel=TRUE with non-VGLM models (e.g. scam/gam/glm) ",
            "may be slower than serial due to dispatch overhead. ",
            "Consider parallel=FALSE for lightweight models.")
  }

#*******************************************************************************
###################### STEP 1(a) - Thin TL data ################################

  # Number of measurements i.e. number of rows
  numTLrow<-dim(TLlookup)[1]

  # Use all TL columns from TLlookup file (1st column is range)
  no.profiles <- dim(TLlookup)[2]-1;

  # Check that transectTruncationDistances have same dimensions as TL
  numTruncDist <- length(truncationDistance)
  if  (is.null(numTruncDist) | numTruncDist == 1) { # Default or single value
    truncationDistance = rep(truncationDistance, no.profiles)
  }else if (dim(truncationDistance)[2] != no.profiles){
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
  output_range_m <- allTLlookupsubset[,1]
  allTLsubset<-allTLlookupsubset[,2:(no.profiles+1)]
  allTrunc<-matrix(data=0, nrow=dim( allTLsubset)[1], ncol=dim(allTLsubset)[2])
  for (i in 1:no.profiles){
    na.index <- output_range_m > truncationDistance[i]
    allTrunc[na.index,i]<-rep(NA)
  }

  #remove TLlookup to save memory
  rm(TLlookup)

  # gc() #garbage collect function to free up memory after removing an object
  #*****************************************************************************

  # STEP 2 #############################
    # Reference: https://smolski.github.io/livroavancado/reglog.html#o-modelo
  # "Dispersion of" events "and" non-events "(binary response variable)":


  #*****************************************************************************
  # STEP 3 #################################
  #create a summary results matricies

  #Number of model parameters: 2 columns for SL mean/st.dev, 2 for NL$mean/st.dev,
  # plus 2 or more  columns for coef values
  no.core.params <- 4;

  no.coef <- length(coef(detFun))
  no.model.params <- no.core.params  + no.coef

  # Matrix that holds bootstrapped parameters (inputs and not actually results)
  results1000sim<-matrix(NA,outerloop,no.model.params)

  # create a large results matrix for all pdets drawn from simulation
  # this is needed to calculate detection function plots
  # the dimensions of this will be 10000 columns by no.profiles*1000 rows
  # (there are 10000 range steps per profile, and there are no.profiles profiles
  # and 1000 iterations)
  resultsallpdets<-matrix(NA,numTLrowsubset,(no.profiles*outerloop))

  # Pre-allocate matrices outside the main loop
  allSL <- matrix(NA, numTLrowsubset, no.profiles)
  allNL <- matrix(NA, numTLrowsubset, no.profiles)
  allSNR <- matrix(NA, numTLrowsubset, no.profiles)
  allpdet <- matrix(NA, numTLrowsubset, no.profiles)

  ## START OF THE OUTER LOOP (PARALLELIZED) ####
  # The Monte Carlo outer loop is embarrassingly parallel: each iteration
  # draws its own SL/NL/coefficient samples and writes to a disjoint column
  # slice of `resultsallpdets`. We split the `outerloop` iterations into
  # chunks and dispatch them via future.apply::future_lapply. Each worker
  # returns its slice of `resultsallpdets` and the corresponding rows of
  # `results1000sim`, which we then stitch back into the pre-allocated
  # matrices so all downstream code (Steps 5-7) is unchanged.
  #
  # Parallel execution is opt-in via the `parallel` argument. When TRUE,
  # set the parallel plan *before* calling cde() / pDetInArea(), e.g.
  #   future::plan(future::multisession, workers = 8)
  # When FALSE (default), runs serially - identical to the original code.

  # Worker function: runs iterations `iters` and returns a list with:
  #   pdets       : numTLrowsubset x (no.profiles * length(iters)) matrix
  #   sim_params  : length(iters) x no.model.params matrix
  #   prob_warn   : logical, TRUE if any probability was clamped
  .pDetInAreaChunk <- function(iters,
                               detFun, SL, NL,
                               allTLsubset, allTrunc,
                               numTLrowsubset, no.profiles, no.model.params,
                               snrTruncationThreshold) {
    chunk_len   <- length(iters)
    chunk_pdets <- matrix(NA_real_, numTLrowsubset, no.profiles * chunk_len)
    chunk_sim   <- matrix(NA_real_, chunk_len, no.model.params)
    prob_warn   <- FALSE

    for (ci in seq_len(chunk_len)) {
      # Preserve original variable name `i` for readability vs. the serial code
      i <- ci

      # -----------------------------------------------------------------
      # STEP 3(a) - SL parametric bootstrap
      SLdist <- stats::rnorm(SL$sampleSize, SL$mean, SL$sd)
      chunk_sim[i, 1] <- mean(SLdist)
      chunk_sim[i, 2] <- stats::sd(SLdist)

      # STEP 3(b) - NL parametric bootstrap
      NLdist <- stats::rnorm(NL$sampleSize, NL$mean, NL$sd)
      chunk_sim[i, 3] <- mean(NLdist)
      chunk_sim[i, 4] <- stats::sd(NLdist)

      # STEP 3(c) - random realisation of the detector characterisation curve
      if (any(class(detFun) == 'scam') | any(class(detFun) == 'gam')) {
        br <- MASS::mvrnorm(1, as.vector(coef(detFun)), detFun$Vp)
      } else { # glm or vglm
        br <- MASS::mvrnorm(1, as.vector(coef(detFun)), vcov(detFun))
      }
      chunk_sim[i, 5:no.model.params] <- br

      # -----------------------------------------------------------------
      # STEP 4 - per-transect SL/NL/SNR draws and p(det) prediction
      allSL   <- matrix(NA, numTLrowsubset, no.profiles)
      allNL   <- matrix(NA, numTLrowsubset, no.profiles)
      allSNR  <- matrix(NA, numTLrowsubset, no.profiles)
      allpdet <- matrix(NA, numTLrowsubset, no.profiles)

      for (j in 1:no.profiles) {
        allSL[, j] <- stats::rnorm(numTLrowsubset, chunk_sim[i, 1], chunk_sim[i, 2])
      }
      for (j in 1:no.profiles) {
        allNL[, j] <- stats::rnorm(numTLrowsubset, chunk_sim[i, 3], chunk_sim[i, 4])
      }

      allSNR <- allSL - allTLsubset - allNL

      for (j in 1:no.profiles) {
        newd <- data.frame(SNR = allSNR[, j])

        detFun.newcoeff <- detFun

        if (inherits(detFun, 'vglm')) {
          detFun.newcoeff@coefficients <- br
          pred0 <- VGAM::predict(detFun.newcoeff, newdata = newd, type = "response",
                                 type.fitted = 'onempall0', na.action = na.pass)
          predmatrix <- VGAM::predict(detFun.newcoeff, newdata = newd, type = "response",
                                      na.action = na.pass)
          index <- ifelse(is.null(detFun.newcoeff@extra$whichObserver),
                          dim(predmatrix)[2],
                          which(colnames(predmatrix) == detFun.newcoeff@extra$whichObserver))
          predmatrix <- predmatrix[, index]
          predmatrix <- predmatrix * pred0
        } else {
          detFun.newcoeff$coefficients <- br
          predmatrix <- predict(detFun.newcoeff, newd, type = "response")
        }

        # SNR truncation and per-radial distance truncation
        predmatrix[which(newd$SNR < snrTruncationThreshold | is.na(allTrunc[, j]))] <- NA

        # Clamp probabilities to [0, 1]
        if (any(predmatrix < 0, na.rm = TRUE)) {
          prob_warn <- TRUE
          predmatrix[predmatrix < 0] <- 0
        }
        if (any(predmatrix > 1, na.rm = TRUE)) {
          prob_warn <- TRUE
          predmatrix[predmatrix > 1] <- 1
        }

        allpdet[, j] <- predmatrix
      }

      # Write this iteration's results into the chunk slice.
      col_start <- (no.profiles * i) - (no.profiles - 1)
      col_end   <- no.profiles * i
      chunk_pdets[, col_start:col_end] <- allpdet
    }

    list(pdets = chunk_pdets, sim_params = chunk_sim, prob_warn = prob_warn)
  }

  # Decide whether to run in parallel based on user request and availability.
  use_parallel <- isTRUE(parallel)
  if (use_parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    warning("parallel=TRUE requested but future.apply is not installed; ",
            "falling back to serial execution.")
    use_parallel <- FALSE
  }

  if (use_parallel) {
    # Strip the model object down to just what workers need (coef, vcov/Vp,
    # predict, class, extra). scam/gam models can be >1 GB because they store
    # the full training data, smooth matrices, etc. Workers only call coef(),
    # vcov(), and predict() with swapped coefficients, so we can safely remove
    # everything else. This typically shrinks a 1 GB scam to <1 MB.
    detFun_slim <- detFun
    if (inherits(detFun_slim, "vglm")) {
      # S4: clear large slots
      detFun_slim@model        <- data.frame()
      detFun_slim@x            <- matrix(0, 0, 0)
      detFun_slim@y            <- matrix(0, 0, 0)
      detFun_slim@residuals    <- matrix(0, 0, 0)
      detFun_slim@fitted.values <- matrix(0, 0, 0)
      detFun_slim@effects      <- numeric(0)
    } else {
      # S3 (scam, gam, glm): clear training data and large intermediates
      detFun_slim$data              <- NULL
      detFun_slim$model             <- NULL
      detFun_slim$y                 <- NULL
      detFun_slim$residuals         <- NULL
      detFun_slim$fitted.values     <- NULL
      detFun_slim$linear.predictors <- NULL
      detFun_slim$effects           <- NULL
      detFun_slim$weights           <- NULL
      detFun_slim$prior.weights     <- NULL
      detFun_slim$offset            <- NULL
      detFun_slim$hat               <- NULL
      detFun_slim$R                 <- NULL
      # scam-specific large matrices
      detFun_slim$X                 <- NULL
      detFun_slim$Vp.t              <- NULL
      detFun_slim$Ve.t              <- NULL
      detFun_slim$db.drho           <- NULL
      detFun_slim$Xt                <- NULL
    }

    # Split outerloop iterations into chunks (one chunk per worker).
    n_workers <- future::nbrOfWorkers()
    # Cap chunks at outerloop (don't make empty chunks if workers > outerloop)
    n_chunks  <- max(1L, min(as.integer(n_workers), outerloop))
    chunks    <- split(seq_len(outerloop),
                       cut(seq_len(outerloop), n_chunks, labels = FALSE))

    # Raise the globals size limit to accommodate large TL matrices.
    # Save and restore the old value so we don't affect the user's session.
    old_maxSize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
    on.exit(options(future.globals.maxSize = old_maxSize), add = TRUE)

    chunk_results <- future.apply::future_lapply(
      chunks,
      .pDetInAreaChunk,
      detFun                 = detFun_slim,
      SL                     = SL,
      NL                     = NL,
      allTLsubset            = allTLsubset,
      allTrunc               = allTrunc,
      numTLrowsubset         = numTLrowsubset,
      no.profiles            = no.profiles,
      no.model.params        = no.model.params,
      snrTruncationThreshold = snrTruncationThreshold,
      future.packages        = c("MASS", "VGAM"),
      future.seed            = TRUE
    )
  } else {
    # Serial: single chunk containing all iterations.
    chunks <- list(seq_len(outerloop))
    chunk_results <- list(
      .pDetInAreaChunk(
        iters                  = chunks[[1]],
        detFun                 = detFun,
        SL                     = SL,
        NL                     = NL,
        allTLsubset            = allTLsubset,
        allTrunc               = allTrunc,
        numTLrowsubset         = numTLrowsubset,
        no.profiles            = no.profiles,
        no.model.params        = no.model.params,
        snrTruncationThreshold = snrTruncationThreshold
      )
    )
  }

  # ---- Stitch chunk results back into the pre-allocated matrices --------
  # `chunks[[k]]` holds the original iteration indices for chunk k. We use
  # those indices to place each chunk's slice into the correct columns of
  # `resultsallpdets` and rows of `results1000sim`.
  prob_warnings <- FALSE
  for (k in seq_along(chunks)) {
    iters_k <- chunks[[k]]
    cr      <- chunk_results[[k]]

    # Rows of results1000sim correspond 1:1 with iters_k
    results1000sim[iters_k, ] <- cr$sim_params

    # Columns of resultsallpdets: iteration i uses columns
    # ((no.profiles*i)-(no.profiles-1)):(no.profiles*i)
    for (local_idx in seq_along(iters_k)) {
      i <- iters_k[local_idx]
      dst_cols <- ((no.profiles * i) - (no.profiles - 1)):(no.profiles * i)
      src_cols <- ((no.profiles * local_idx) - (no.profiles - 1)):(no.profiles * local_idx)
      resultsallpdets[, dst_cols] <- cr$pdets[, src_cols]
    }

    if (isTRUE(cr$prob_warn)) prob_warnings <- TRUE
  }

  if (prob_warnings) {
    warning("Probabilities outside [0,1] were encountered and clamped")
  }



  #*****************************************************************************
  # SIMULATION ENDED - NOW PRODUCE THE RESULTS #################################
  #Save originally selected parameters for the simulation
  if (!is.null(simResultsFile)){
    dir.create(dirname(results1000sim), showWarnings = FALSE)
    utils::write.table(results1000sim,file=simResultsFile,
                     row.names=F,col.names=F)
  }
  #*****************************************************************************

  # STEPS 5 - for all p(det), calculate a distance-weighted average ############

  # Weight each prob(detect) by the range (i.e. multiply the p(det) by the
  # range) Then divide all weighted p(dets) for a transect by the sum of ALL
  # ranges across each transect. want to save the 1000 weighted mean p(dets) for
  # each transect a bit easier to handle than resultsallpdets

  alltransectspdetsweight1000<-matrix(NA,no.profiles,outerloop)

  ## Repeat all steps for each transect (in loop) #############
  for (i in 1:no.profiles){
    # extract the data (10000x1000) for each transect
    transectallpdets<-resultsallpdets[,seq(i,dim(resultsallpdets)[2],no.profiles)]

    # Now want to weight the pdets using the ranges define a matrix to save the
    # results - will be numTLrowsubset x outerloop for each transect
    transectallpdetsweight<-matrix(NA,numTLrowsubset,outerloop)

    #For each range step, multiply the 1000 x p(det) values by the range
    for (j in 1:numTLrowsubset){
      transectallpdetsweight[j,]<-allTLlookupsubset[j,1]*transectallpdets[j,]
    }

    #calculate average pdet for each transect: So now need to divide each column
    #(1000 iterations) by the sum of the ranges
    transectpdetsweight1000<-c()
    for (j in 1:outerloop){
      # Debugging code to check for NA
      # browser(expr = {any(is.na(transectallpdetsweight[,j]))})
      include <- !is.na(transectallpdetsweight[,j])
      transectpdetsweight1000[j]<-
        sum(transectallpdetsweight[include,j])/sum(allTLlookupsubset[include,1])
    }
    #This gives 1000 weighted average probability of detections for each
    #iteration (Specific to one of the transects over the length w)

    #Save in alltransectspdetsweight1000 matrix
    alltransectspdetsweight1000[i,]<-transectpdetsweight1000
  }
  #So now have a results matrix with no.profilesx1000 average weighted p(dets)
  #*****************************************************************************
  # STEP 6 - Final p(det) calculations ####
  ## STEP 6(a) Calculate a transect-specific mean p(det) for each transect ####
  transectavgweightpdets<-rowMeans(alltransectspdetsweight1000, na.rm = T)

  ## STEP 6(b) Calculate an overall p(det) for the whole simulation ###########

  # Shouldn't be any NA's here now, but we set na.rm=T anyway (maybe a whole transect is NA)
  # now that transects might be different lengths, we need to weight by length
  # when averaging across them.
  # overallpdet<-mean(transectavgweightpdets, na.rm=T)
  overallpdet<-weighted.mean(transectavgweightpdets, w=truncationDistance^2,
                             na.rm=T)

  ## STEP 6(c) Create a results matrix to hold p(det) values and ##############
  # st.error/st.dev values

  pdetsandvar<-matrix(NA,no.profiles+1,2)
  pdetsandvar[1:no.profiles,1]<-transectavgweightpdets
  pdetsandvar[no.profiles+1,1]<-overallpdet
  #*****************************************************************************
  # STEP 7 - Variance calculations
  # STEP 7(a) Need to calculate the standard error of the transect-specific
  #    p(det) values
  st.errorPt<-sd(transectavgweightpdets,na.rm = T)/sqrt(no.profiles)

  # STEP 7(b) Need to calculate the standard deviation of the transect-specific
  #   p(det) values
  st.devPt<-c()
  for (i in 1:no.profiles){
    st.devPt[i]<-sd(alltransectspdetsweight1000[i,],na.rm = T)
  }
  # STEP 7(c) save the results in the results matrix and print out
  pdetsandvar[1:no.profiles,2]<-st.devPt
  pdetsandvar[no.profiles+1,2]<-st.errorPt
  if (!is.null(paFile)){
    dir.create(dirname(paFile), showWarnings = FALSE)
    utils::write.table(pdetsandvar,file=paFile, row.names=F,
                       col.names = c("Mean","SD"))
  }

  #****************************************************************************
  #Calculate an average across p(det) for each range step for all transects####
  #Use the main results matrix, resultsallpdets

  range_m <- allTLlookupsubset[,1]
  combineddetfunc<-matrix(NA,no.profiles,dim(resultsallpdets)[1]) #10000
  for (i in 1:no.profiles){
    #extract the transect specific data (10000x1000) for each transect
    transectallpdets<-resultsallpdets[,seq(i,dim(resultsallpdets)[2],no.profiles)]

    #calculate the mean pdet (unweighted) for each range step
    transectavgdetfunc<-rowMeans(transectallpdets, na.rm=T)
    combineddetfunc[i,]<-transectavgdetfunc
  }

  allDetFunctions <- allTLlookupsubset
  allDetFunctions[,2:(no.profiles+1)]<-t(combineddetfunc)
  allDetFunctions<- data.table::as.data.table(allDetFunctions)
  names(allDetFunctions)<-dimnames(allTLlookup_h)[[2]]

  if (!is.null(transectFile)){
    dir.create(dirname(transectFile), showWarnings = FALSE)
    data.table::fwrite(allDetFunctions, file=transectFile)
  }


  ## Average across transects to produce an overall det function ###############
  averagecombineddetfunc<-colMeans(combineddetfunc, na.rm=T)

  meanOfAllTransects = data.frame(
    range_m=output_range_m,
    pDet = averagecombineddetfunc
  )
  return(list(overall = overallpdet,
              perTransectMeanSD = pdetsandvar,
              meanOfAllTransects = meanOfAllTransects))

  # Once have looked at the results, decide an appropriate truncation distance,
  # w Redo Steps 5-7 using the different truncation distance
  ### If W ok, FINAL ###

}
