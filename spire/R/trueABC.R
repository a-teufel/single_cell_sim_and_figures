#' Title
#'
#' @param first_round_ABC_parameters parameter set to start from this can be from the GA or a previous round of ABC
#' @param max  values from sigmoidal fit of infection data
#' @param slope  values from sigmoidal fit of infection data
#' @param midpoint  values from sigmoidal fit of infection data
#' @param lysis values from sigmoidal fit of infection data
#' @param MOI MOI that experiments were run at
#' @param pop_size the population size to be used by the GA
#' @param n_saves number of save before the round of ABC is over
#' @param n_replicates number of replicate parameter sets to generate
#' @param n_core number of cores to be used when evaluating final population
#' @param initpop if a previous population estimation is to be used to seed arnother round of the GA
#' @param precent_error an exceptable amount of error outside of the bounds of the supplied parameter distributions
#'
#' @return Function returns parameter values accepted from ABC simulation and the max,slope, and midpoint values
#' @export


trueABC <-
  function(first_round_ABC_parameters,max, slope, midpoint, lysis,  MOI, n_saves, n_replicates, n_cores, accept_crit, precent_error,write) {
    #register how many cores are to be used when calling replicates
    registerDoMC(n_cores)


    #upper an lower bounds for ABC
    low<-c(-5.9,-13,-10.4,-10,-0.25,-10,4.1,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
    high<-c(0.5,0,.5,0,2.5,0,10,17.5,10,0.95,1,1,1,1,1,1,1,1,1,1,1,1)


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


    max_ka<-max(max)+max(max)*precent_error
    min_ka<-min(max)-min(max)*precent_error

    max_b1<-max(slope)+max(slope)*precent_error
    min_b1<-min(slope)-min(slope)*precent_error

    max_m1<-max(midpoint)+max(midpoint)*precent_error
    min_m1<-min(midpoint)-min(midpoint)*precent_error

    #find out what the max lys time is that is under 24 hours.
    cp_lysis <- lysis
    cp_lysis[cp_lysis == 24] <- 0
    cut_off_lysis <- max(cp_lysis)

    max_ly <- 24
    min_ly <- min(lysis) - min(lysis) * precent_error

    savedParameters <- first_round_ABC_parameters
    savedParameters_var <- NULL
    savedJitteredParameters <- NULL
    savedParameters2 <- NULL
    savedVar <- NULL

    #interal variance to use
    sigma_prime <- sig(first_round_ABC_parameters)

    #parameters for sigmoidal
    est_maxa <- NULL
    est_slope1 <- NULL
    est_midpoint1 <- NULL
    est_lys <- NULL

    #parameters for linear
    est_max_d<-NULL
    est_best<-NULL

    n_saved <- 0
    i <- 0

    while (n_saved < n_saves) {
      print("cut off is ")
      print(accept_crit)
      i <- i + 1
      print("round")
      print(i)
      print("numbers saved")
      print(n_saved)
      inputs <- nrow(savedParameters)

      #randomly pick a new parameter set from the population
      rand <- runif(26, 1, inputs)

      #create a parameter vector to send to the abc
      param_list_var <-
        c(
          savedParameters[rand[1],1],savedParameters[rand[2],2],savedParameters[rand[3],3],savedParameters[rand[4],4],savedParameters[rand[5],5],savedParameters[rand[6],6],savedParameters[rand[7],7],savedParameters[rand[8],8],
          savedParameters[rand[9],9],savedParameters[rand[10],10],savedParameters[rand[11],11],savedParameters[rand[12],12],savedParameters[rand[13],13],savedParameters[rand[14],14],savedParameters[rand[15],15],savedParameters[rand[16],16],savedParameters[rand[17],17],savedParameters[rand[18],18],
          savedParameters[rand[19],19],savedParameters[rand[20],20],savedParameters[rand[21],21],savedParameters[rand[22],22],savedParameters[rand[23],23],savedParameters[rand[24],24],savedParameters[rand[25],25],savedParameters[rand[26],26]
        )

      #add some variance to the parameters
      param_temp_var <- param_list_var

      for (k in 1:length(param_list_var)) {
        try_newparam <- rnorm(1,mean = param_list_var[k],sd = sigma_prime[k])

        while (try_newparam > high[k] || try_newparam < low[k]) {
          try_newparam <- rnorm(1,mean = param_list_var[k],sd = sigma_prime[k])
        }

        param_temp_var[k] <- try_newparam
      }

      #for the new parameter set generate replicates
      Res2 <- foreach(j = 1:n_replicates) %dopar% {

        param_temp <- c(0, 0, 0,  0,     0,   0,  0,   0,  0,  0, 0)
        for (k in 1:11) {
          try_newparam <-
            rnorm(1,mean = param_temp_var[k],sd = param_temp_var[k + 11])

          while (try_newparam > high[k] || try_newparam < low[k]) {
            try_newparam <-
              rnorm(1,mean = param_temp_var[k],sd = param_temp_var[k + 11])
          }
          param_temp[k] <- try_newparam
        }


        ly_time <-
          rgamma(1,shape = param_temp_var[23],rate = param_temp_var[25])
        while (ly_time <  min_ly - 3) {
          #minus 3 because obs dont start till after 3 hours
          ly_time <-
            rgamma(1,shape = param_temp_var[23],rate = param_temp_var[25])
        }


        lysis_time <- round(ly_time / .5, digits = 0) * .5

        st_time <-
          rnorm(1,mean = param_temp_var[24],sd = param_temp_var[26])

        while (st_time < 0 || st_time > max_time) {
          st_time <-
            rnorm(1,mean = param_temp_var[24],sd = param_temp_var[26])
        }

        if (st_time >= lysis_time) {
          st_time <- 0
        }

        start_time <- round(st_time / .5, digits = 0) * .5

        if (start_time + lysis_time > cut_off_lysis) {
          start_time <- 1
          st_time <- 1
          lysis_time <-24 - start_time
          ly_time<-24-st_time
        }

        time = seq(start_time,start_time + lysis_time,0.5)


        #if a valid time length is sampled then run the simulation
        if (length(time) > (min_ly * 2) - (3 * 2) &&
            start_time + lysis_time <=  max_ly &&
            start_time + lysis_time > min_ly) {
          sim_min_time <- 0
          sim_max_time <- lysis_time


          ABC_data<-NULL

          ABC_data<-avepa(param_temp,MOI,sim_min_time, sim_max_time,write)


          best_model<-NULL
          parameterVector_sig<-NULL
          parameterVector_line<-NULL
          param_temp<-c(param_temp,ly_time,start_time)

          #if simulation results are ok, then do model fitting
          if (length(unique(ABC_data)) != 1) {
            dataInput = data.frame(intensity = ABC_data,time = time)

            dataOutput = normalizeData(dataInput)
            dataInput2 = dataOutput

            parameterVector_sig = multipleFitFunction(
              dataInput = dataInput2,model = "sigmoidal",n_runs_min = 20,n_runs_max =
                100
            )
            print("model fit")

            #if data is resonable save the values
            if (parameterVector_sig$isThisaFit &
                parameterVector_sig$maximum_Estimate < max_ka &
                parameterVector_sig$maximum_Estimate
                > min_ka &
                parameterVector_sig$slopeParam_Estimate < max_b1 &
                parameterVector_sig$slopeParam_Estimate > min_b1 &
                parameterVector_sig$midPoint_Estimate < max_m1 &
                parameterVector_sig$midPoint_Estimate > min_m1) {
              list(parameterVector_sig$maximum_Estimate,parameterVector_sig$slopeParam_Estimate,parameterVector_sig$midPoint_Estimate, best_model,param_temp,st_time+ly_time)


            }
          }
        }
      }#end par
      print("ended")
      #get the results into a usuable from the parallel call
      maxa_dist <- do.call(cbind, lapply(Res2, "[[", 1))

      slope1_dist <- do.call(cbind, lapply(Res2, "[[", 2))

      midpoint1_dist <- do.call(cbind, lapply(Res2, "[[", 3))

      best_dist <- do.call(cbind, lapply(Res2, "[[", 4))

      jitter_parameters <- do.call(rbind,lapply(Res2, "[[",5))

      ly_dist <- do.call(cbind,lapply(Res2, "[[",6))

      #if the run was ok then score the set
      if (length(maxa_dist) != 0) {

        max_ks <- ks.test(maxa_dist, max,exact = FALSE)

        slope_ks <- ks.test(slope1_dist, slope,exact = FALSE)

        midpoint_ks <-
          ks.test(midpoint1_dist, midpoint,exact = FALSE)

        ly_ks <- ks.test(ly_dist, lysis,exact = FALSE)

        max_d <-
          sum(
            max_ks$statistic,slope_ks$statistic,midpoint_ks$statistic,ly_ks$statistic
          )

        #if the score is good enough then save that parameter set
        if (max_d <= accept_crit) {
          parm <- cbind(param_temp_var)
          savedParameters2 <- rbind(savedParameters2,parm)

          jittered_parameters <- cbind(jitter_parameters)
          savedJitteredParameters <-
            rbind(savedJitteredParameters,jitter_parameters)
          est_maxa <- cbind(est_maxa,maxa_dist)
          est_slope1 <- cbind(est_slope1,slope1_dist)
          est_midpoint1 <- cbind(est_midpoint1,midpoint1_dist)
          est_lys <- cbind(est_lys,ly_dist)
          est_max_d <- cbind(est_max_d,max_d)
          n_saved <- n_saved + 1
        }

      }#end first round of ABC
    }#end while

    trueabc <- NULL
    trueabc$est_ABC_parameters <- savedParameters2
    trueabc$est_maxa <- est_maxa
    trueabc$est_slope1 <- est_slope1
    trueabc$est_midpoint1 <- est_midpoint1
    trueabc$est_lys <- est_lys
    trueabc$est_max_d <- est_max_d
    trueabc$est_jitter_ABC_parameters <- savedJitteredParameters


    return(trueabc)
  }#function

