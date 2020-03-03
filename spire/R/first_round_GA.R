#' Title
#'
#' @param max  values from sigmoidal fit of infection data
#' @param slope  values from sigmoidal fit of infection data
#' @param midpoint  values from sigmoidal fit of infection data
#' @param lysis values from sigmoidal fit of infection data
#' @param MOI MOI that experiments were run at
#' @param pop_size the population size to be used by the GA
#' @param n_iterations number of generation to run the GA
#' @param n_replicates number of replicate parameter sets to generate
#' @param n_core number of cores to be used when evaluating final population
#' @param initpop if a previous population estimation is to be used to seed arnother round of the GA
#' @param precent_error an exceptable amount of error outside of the bounds of the supplied parameter distributions
#'
#' @param write a boolean variable that tells the simulation if you want to save the times when events in the simulation occur
#'
#' @return Population sorted by the best values
#' @export
#'

first_round_GA <-
  function(max, slope, midpoint, lysis,MOI,pop_size, n_iterations, n_replicates, min_time, max_time, initpop, n_cores,precent_error,write){
    #upper and lower bound for the GA
    #setting rules no less than -15, no more than 20
    #note that the GA will use all availbe cores, i just pass it the n_core param so that the call is the same as the abc
    low<-c(-5.9,-13,   -10.4, -6.8,-0.25,-10,      4.1,  11,  1,  0,  0, 0,0,0,0,0,0,0,0,0,0,0)
    high<-c(0.5, -5.1,   0.5,  0,   2.5, -0.8,   8,    17.5, 10,  0.95,.0001,1,1,1,1,1,1,1,1,1,1,1)

    min_time<-1
    max_time<-50

    min_shift<-0
    max_shift<-10

    min_var_time<-0
    max_var_time<-5

    min_var_shift<-0
    max_var_shift<-20


    low<-c(low,min_time,min_shift,min_var_time,min_var_shift)
    high<-c(high,max_time,max_shift,max_var_time,max_var_shift)

    savedParameters <- NULL
    savedParameters_var <- NULL
    savedJitteredParameters <- NULL
    savedParameters2 <- NULL
    savedVar <- NULL


    #parameters for sigmoidal
    est_max <- NULL
    est_slope1 <- NULL
    est_midpoint1 <- NULL

    #parameters for linear
    est_inter <- NULL
    est_slope <- NULL
    est_best <- NULL


    #function to do fitting
    CallSim <- function(x) {
      low<-c(-5.9,-13,   -10.4, -6.8,-0.25,-10,      4.1,  11,  1,  0,  0)
      high<-c(0.5, -5.1,   0.5,  0,   2.5, -0.8,   8,    17.5, 10,  0.95,.0001)

      Vect_max <- NULL
      Vect_slope <- NULL
      Vect_mid <- NULL
      Vect_best <- NULL
      Vect_ly <- NULL

      #calculate bounds for distributions
      max_ka <- max(max) + max(max) * precent_error
      min_ka <- min(max) - min(max) * precent_error

      max_b1 <- max(slope) + max(slope) * precent_error
      min_b1 <- min(slope) - min(slope) * precent_error


      max_m1 <- max(midpoint) + max(midpoint) * precent_error
      min_m1 <- min(midpoint) - min(midpoint) * precent_error


      #find out what the max lys time is that is under 24 hours.
      cp_lysis <- lysis
      cp_lysis[cp_lysis == 24] <- 0
      cut_off_lysis <- max(cp_lysis)

      #expriments were only run for 24 hours
      max_ly <- 24
      min_ly <- min(lysis) - min(lysis) * precent_error

      #store input from GA and add noise to each of these parameter sets
      param_list <- x
      for (j in 1:n_replicates) {
        #add noise to the parameters
        param_temp <- c(0, 0, 0,  0,     0,   0,  0,   0,  0,  0, 0)
        for (k in 1:11) {
          try_newparam <- rnorm(1,mean = param_list[k],sd = param_list[k + 11])
          while (try_newparam > high[k] || try_newparam < low[k]) {
            try_newparam <- rnorm(1,mean = param_list[k],sd = param_list[k + 11])
          }
          param_temp[k] <- try_newparam
        }

        #sample a lysis time
        ly_time<-rgamma(1,shape=param_list[23],rate=param_list[25])
        while ( ly_time <  min_ly-3){ #minus 3 because obs dont start till after 3 hours

          ly_time<-rgamma(1,shape=param_list[23],rate=param_list[25])
        }


        #experiments were only checked every 30 mins, so we want to run at 30 min intervales
        lysis_time <- round(ly_time / .5, digits = 0) * .5

        #sample a start time
        st_time<-rnorm(1,mean=param_list[24],sd=param_list[26])

        while(st_time < 0 || st_time > max_time){
          st_time<-rnorm(1,mean=param_list[24],sd=param_list[26])
        }


        #experiments were only checked every 30 mins, so we want to run at 30 min intervales
        start_time <- round(st_time / .5, digits = 0) * .5


        #if the start time and the lysis time are too long, adjust so we don't get values over 24
        if (start_time + lysis_time > cut_off_lysis) {
          lysis_time <- 24 - start_time
          ly_time <- 24 - st_time
        }

        #if the values are outside the bounds of our start and stop time don't even bother running the model and fitting it
        if( start_time+lysis_time >  max_ly || start_time+lysis_time < min_ly || lysis_time < start_time){
          return(Inf)
        }


        #create time sequence
        time = seq(start_time,start_time + lysis_time,0.5)

        #min lysis time *2 - number of time points until lysis, in real data they dont start obs until hours 3, so 3*2 points out
        if (length(time) < (min_ly * 2) - (3 * 2)) {
          return(Inf)
        }

        sim_min_time <- 0
        sim_max_time <- lysis_time

        #call the simulation
        ABC_data <- NULL
        ABC_data <-
          avepa(param_temp,MOI,sim_min_time, sim_max_time,write)

        #if the simulation returned dyanmics were something actually happened then do the model fitting
        #run the sim until infection happens
        if(length(unique(ABC_data))!=1){

          dataInput=data.frame(intensity=ABC_data,time=time)

          dataOutput = normalizeData(dataInput)
          dataInput2=dataOutput


          parameterVector_sig=multipleFitFunction(dataInput=dataInput2,model="sigmoidal",n_runs_min=20,n_runs_max=100)

          #make sure fit works
          if( parameterVector_sig$isThisaFit &  parameterVector_sig$maximum_Estimate < max_ka & parameterVector_sig$maximum_Estimate
              > min_ka & parameterVector_sig$slopeParam_Estimate < max_b1 & parameterVector_sig$slopeParam_Estimate > min_b1 &
              parameterVector_sig$midPoint_Estimate < max_m1 & parameterVector_sig$midPoint_Estimate > min_m1){


            Vect_max<-c(Vect_max,parameterVector_sig$maximum_Estimate)
            Vect_slope<-c(Vect_slope, parameterVector_sig$slopeParam_Estimate)
            Vect_mid<-c(Vect_mid, parameterVector_sig$midPoint_Estimate)
            #Vect_best<-c(Vect_best,best_model)
            Vect_ly<-c(Vect_ly,st_time+ly_time)
          }

        }


      }#end replicates

      #if there are saved values calculate the score
      if (length(Vect_max) != 0) {
        max_ks <- ks.test(Vect_max, max,exact = FALSE)

        slope_ks <- ks.test(Vect_slope, slope,exact = FALSE)

        midpoint_ks <- ks.test(Vect_mid, midpoint,exact = FALSE)

        ly_ks <- ks.test(Vect_ly, lysis,exact = FALSE)

        max_d <-
          sum(
            max_ks$statistic,slope_ks$statistic,midpoint_ks$statistic,ly_ks$statistic
          )

        return(max_d)
      }
      #if there are no saved values tell the GA it did poorly
      return(Inf)
    }#end function fitting


    parVar = c(
      "max","slope","midpoint","lysis","n_replicates","MOI","max_time","min_time","precent_error"
    )
    optimize_function_cmp <- cmpfun(CallSim)

    #ptm<-proc.time()
    round1<-DEoptim(optimize_function_cmp,lower=low,upper=high, control=list(NP=pop_size, initialpop=initpop, parallelType=1,packages=c(library(sicegar),library(spire),library(Rcpp),library(R.utils)),parVar=parVar,strategy=6, c=.15, itermax=n_iterations, storepopfrom=1, storepopfreq=1 ,trace=TRUE))

    return(round1)
  }


#helper function to sort the ending population and get the distributions of the max, slope, mid, and lysis point from it (I don't use this anymore)
score_pop <-
  function(pop,max, slope, midpoint, lysis,MOI,pop_size, n_iterations, n_replicates, min_time, max_time,n_cores,precent_error,write) {
    registerDoMC(n_cores)
    max_ly <- 24
    min_ly <- min(lysis) - min(lysis) * precent_error
    #find out what the max lys time is that is under 24 hours.
    cp_lysis <- lysis
    cp_lysis[cp_lysis == 24] <- 0
    cut_off_lysis <- max(cp_lysis)
    Vect_max <- NULL
    Vect_slope <- NULL
    Vect_mid <- NULL
    Vect_best <- NULL
    Vect_ly <- NULL
    max_ka <- max(max) + max(max) * precent_error
    min_ka <- min(max) - min(max) * precent_error

    max_b1 <- max(slope) + max(slope) * precent_error
    min_b1 <- min(slope) - min(slope) * precent_error


    max_m1 <- max(midpoint) + max(midpoint) * precent_error
    min_m1 <- min(midpoint) - min(midpoint) * precent_error

    low <-
      c(-2,-10,-10,-10,-10,-10,  2,   2,  1,  0, 0, 0,0,0,0,0,0,0,0,0,0,0)
    high <-
      c(6, 0,    1.5,   20,   20,    0,   20,   15,  10, 1, .0001,1,1,1,1,1,1,1,1,1,1,1)

    #get the last population for the GA
    last_pop <- pop$member$pop

    Pop <- NULL
    #for each member of the last population generate replicates
    for (i in 1:length(last_pop[,1])) {
      param_list <- last_pop[i,]
      max_score = 0
      Vect_max <- NULL
      Vect_slope <- NULL
      Vect_mid <- NULL
      Vect_ly <- NULL

      Res2 <- foreach(j = 1:n_replicates) %dopar% {
        score = 0
        param_temp <- c(0, 0, 0,  0,     0,   0,  0,   0,  0,  0, 0)
        for (k in 1:11) {
          try_newparam <- rnorm(1,mean = param_list[k],sd = param_list[k + 11])
          while (try_newparam > high[k] || try_newparam < low[k]) {
            try_newparam <- rnorm(1,mean = param_list[k],sd = param_list[k + 11])
          }
          param_temp[k] <- try_newparam
        }

        ly_time <- rgamma(1,shape = param_list[23],rate = param_list[25])
        while (ly_time <  min_ly - 3) {
          #minus 3 because obs dont start till after 3 hours
          ly_time <- rgamma(1,shape = param_list[23],rate = param_list[25])
        }

        lysis_time <- round(ly_time / .5, digits = 0) * .5

        st_time <- rnorm(1,mean = param_list[24],sd = param_list[26])

        while (st_time < 0 || st_time > max_time) {
          st_time <- rnorm(1,mean = param_list[24],sd = param_list[26])
        }

        start_time <- round(st_time / .5, digits = 0) * .5


        if (start_time + lysis_time > cut_off_lysis) {
          lysis_time <- 24 - start_time
          ly_time <- 24 - st_time
        }

        #stopping conditions if something is wrong, don't bother calling the model
        if (start_time + lysis_time >  max_ly ||
            start_time + lysis_time < min_ly || lysis_time < start_time) {
          score <- -Inf
        }

        if (score != -Inf) {
          time = seq(start_time,start_time + lysis_time,0.5)

          if (length(time) < (min_ly * 2) - (3 * 2)) {
            score = -Inf
          }
          #is things are ok, go ahead and run the model
          if (score != -Inf) {
            param_temp <- param_list[1:11]
            sim_min_time <- 0
            sim_max_time <- lysis_time


            ABC_data <- NULL
            ABC_data <-
              avepa(param_temp,MOI,sim_min_time, sim_max_time,write)

            if (length(unique(ABC_data)) != 1) {
              dataInput = data.frame(intensity = ABC_data,time = time)
              dataOutput = normalizeData(dataInput)
              dataInput2 = dataOutput

              parameterVector_sig = multipleFitFunction(
                dataInput = dataInput2,model = "sigmoidal",n_runs_min = 20,n_runs_max =
                  100
              )

              #check that fit is good and that estimates are within bounds
              if (parameterVector_sig$isThisaFit &
                  parameterVector_sig$maximum_Estimate < max_ka &
                  parameterVector_sig$maximum_Estimate
                  > min_ka &
                  parameterVector_sig$slopeParam_Estimate < max_b1 &
                  parameterVector_sig$slopeParam_Estimate > min_b1 &
                  parameterVector_sig$midPoint_Estimate < max_m1 &
                  parameterVector_sig$midPoint_Estimate > min_m1) {
                list(
                  parameterVector_sig$maximum_Estimate,parameterVector_sig$slopeParam_Estimate,parameterVector_sig$midPoint_Estimate, best_model,param_temp,st_time +
                    ly_time
                )

              }

            }

          }#end correct time

        }# end good score

      }# end par
      #save all distributions from par call
      maxa_dist <- do.call(cbind, lapply(Res2, "[[", 1))

      slope1_dist <- do.call(cbind, lapply(Res2, "[[", 2))

      midpoint1_dist <- do.call(cbind, lapply(Res2, "[[", 3))

      best_dist <- do.call(cbind, lapply(Res2, "[[", 4))

      jitter_parameters <- do.call(rbind,lapply(Res2, "[[",5))

      ly_dist <- do.call(cbind,lapply(Res2, "[[",6))


      #store the final scores and each of the distributions
      max_score <- NA
      if (length(maxa_dist) != 0) {
        max_ks <- ks.test(maxa_dist , max,exact = FALSE)

        slope_ks <- ks.test(slope1_dist , slope,exact = FALSE)

        midpoint_ks <- ks.test(midpoint1_dist, midpoint,exact = FALSE)

        ly_ks <- ks.test(ly_dist, lysis,exact = FALSE)

        max_d <-
          sum(
            max_ks$statistic,slope_ks$statistic,midpoint_ks$statistic,ly_ks$statistic
          )

        max_score <- max_d

        temp <- NULL
        temp$param_list <- as.vector(param_list)
        temp$score <- as.vector(max_score)
        temp$max <- as.vector(maxa_dist)
        temp$slope <- as.vector(slope1_dist)
        temp$mid <- as.vector(midpoint1_dist)
        temp$ly <- as.vector(ly_dist)
        ttt <- as.data.frame(t(temp))
        print(ttt)
        Pop <- rbind(Pop,ttt)
      }

    }

    #sort the final population
    df <- as.data.frame(lapply(Pop$score, unlist))
    Pop <- cbind(Pop,t(df))
    colnames(Pop)[7] <- "sort_score"
    Pop <- arrange(Pop, sort_score)
    return(Pop)


  }#end scoring function
