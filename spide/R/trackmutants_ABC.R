
#' Title
#'
#' @param second_round_ABC_parameters parameter set that generations are to be estimated from
#' @param MOI moi of experiments
#' @param n_replicates number of replicates of a given parameter value

#' @return generations
#' @export
trackmutants_ABC <-function(second_round_ABC_parameters, MOI, n_replicates){

  G_dist<-NULL
  G<-NULL
  low <-
    c(-2,-10,-10,-10,-10,-10,  2,   2,  1,  0, 0, 0,0,0,0,0,0,0,0,0,0,0)
  high <-
    c(6, 0,    1.5,   20,   20,    0,   20,   15,  10, 1, .0001,1,1,1,1,1,1,1,1,1,1,1)



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



   savedParameters<-second_round_ABC_parameters #just being lazy cuz I didn't want to change all the value names in my code fix later

 dist<-NULL
 ANS<-NULL
  #for each parameter set in the population simulate
  for(i in 1:nrow(savedParameters)){

    param_list<-as.numeric(savedParameters[i,])


    Res2 <- foreach(h = 1:n_replicates) %dopar% {
    #for each parameter set pick a variant
    param_temp<-c(0, 0, 0,  0,     0,   0,  0,   0,  0,  0)
	  for(k in 1:11){
	      try_newparam<-rnorm(1,mean=param_list[k],sd=param_list[k+11])
		  param_temp[k]<-try_newparam
	    }


            ly_time<-rgamma(1,shape=param_list[23],rate=param_list[25])
	    while ( ly_time < 0){
	      ly_time<-rgamma(1,shape=param_list[23],rate=param_list[25])
	      }


	   lysis_time<-round( ly_time/.5, digits = 0)*.5


           sim_min_time<-0
	   sim_max_time<-lysis_time

	   ABC_gens<-0

	  #if the sim time isn't less then an hour
	  if(sim_max_time > 2){
	    ABC_gens<-avetrack(param_temp,MOI, 0, sim_max_time)
	      }
	   list(ABC_gens)
	  }#end loop over all parameter sets

   G<-do.call(cbind, lapply(Res2, "[[", 1))

  dist<-cbind(dist,G)
  }
  return(dist)
  } #end function
