#include <Rcpp.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace Rcpp;

// Above this value for the master rate, time steps are estimated as their expected value; below this threshold
// they are drawn from an exponential distribution
//#define threshold (1000) what they had set.

#define threshold (1000)
//#define threshold (0)
// Fitted parameters for absorption time distribution
#define SHAPE (0.678)
#define RATE (0.02)


// [[Rcpp::export]]
SEXP pa(NumericVector p,double MOI,double min_time, double max_time, int rep, int write)
{


  //set up R seed so each call does not produce the same results. Remove if you want each call to be the same
  Environment g = Environment::global_env();
  Environment::Binding RandomSeed = g[".Random.seed"];

  clock_t start, end;
  start = clock();
  //FREE PARAMETERS
  //limiting constant for translation
  double c_trans;
  //constant for rate of compartmentaliztion
  double c_compart;
  //constant for 3CD-linear RNA reaction
  double c_circ;
  //constants for replication
  //Replication OF positives FROM negatives
  double c_repP;
  //Replication OF negatives FROM positives
  double c_repN;
  //encapsidation rate
  double c_pack;
  //number of possible compartments
  double MEMBRANE;
  //Limit constant for replication
  unsigned long long repLimit;
  //molecules of 3A consumed per compartment
  int consumption;
  //fraction of unencapsulated +strands that stay in their compartment
  double SAH;
  //scaler to convert gfp to RNA
  double scale;

  //FIXED PARAMETERS
  //polyproteins per capsid
  int nC = 60;
  //max time (minutes)
  int maxT = (max_time * 60)+5;
  double masterRate;
  //double rates [6];
  NumericVector rates(5);

  //set free variables from input
  c_trans = exp(p[0]);
  c_compart = exp(p[1]);
  c_circ = exp(p[2]);
  c_repP = exp(p[3]);
  c_repN = exp(p[4]);
  c_pack = exp(p[5]);
  MEMBRANE = exp(p[6]);
  repLimit = (unsigned long long)(exp(p[7]));
  consumption = round(exp(p[8]));
  SAH = p[9];
  scale=p[10];


  NumericVector meanReplicated(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector meanGFP(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector sampleTimes(((max_time*2)-(min_time*2))+1);
  NumericVector Timeout(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  for(int i = 0; i <= (((max_time)*2)-(min_time*2)); i++)
  {

    meanReplicated[i]=0;
    meanGFP[i]=0;
    Timeout[i]=0;
    sampleTimes[i]=60*((i+(min_time*2))*.5);
    // printf("same: %f",sampleTimes[i]);
  }

  NumericVector TimeoutRates(5);
  NumericVector SumRates(5);
  for(int i=0; i < 6; i++){
    TimeoutRates[i]=0;
    SumRates[i] =0;
  }

  //Concentration of protein used to make compartments
  long long protein = 0;
  int comparts = round(MEMBRANE);
  int packaged = 0;
  int CAP = 0;
  int RNA_P = 0;
  int RNA_N = 0;
  double t = 0;
  int minuteCounter = 0;
  int unattached = 0;
  int coms = 0;
  int tComs = 0;
  unsigned long long cds = 0;
  double tau;
  unsigned long long replicated = 0;
  unsigned long long GFP = 0;

  int buffer = 10000;

  // Keeps track of how many compartments have $index copies of the replication-competence protein
  //int* cs = (int*)(malloc(buffer * sizeof(int)));
  //int* cs=new int[buffer];
  NumericVector cs(buffer);
  for(int i = 0; i < buffer; i++) cs[i] = 0;

  // Used to efficiently scale replication rate to logistic limit
  double pacer = 1;
  double decrement = 1.0f / repLimit;

  // Assigns each infecting virus to a timepoint and sorts those into a schedule
  int popCounter = 0;
  int incoming;
  do
  {
    //incoming = gsl_ran_poisson(rgB, MOI);
    incoming=R::rpois(MOI);
    //printf("incoming: %f", incoming);

  }while(incoming == 0);

  //double* popTimes = malloc(sizeof(double) * incoming);
  //double* popTimes=new double[incoming];
  NumericVector popTimes(incoming);
  for(int i = 0; i < incoming; i++)
  {
    //popTimes[i] = gsl_ran_gamma(rgB, SHAPE, 1.0f / RATE);
    popTimes[i] = R::rgamma(SHAPE, 1.0f / RATE);
  }
  if(incoming > 1)
  {
    int changed = 0;
    do
    {
      changed = 0;
      for(int i = 1; i < incoming; i++)
      {
        if(popTimes[i] < popTimes[i-1])
        {
          changed++;
          double temp = popTimes[i];
          popTimes[i] = popTimes[i-1];
          popTimes[i-1] = temp;
        }
      }
    }while(changed > 0);

  }

  //Advance to the first infection time and, if less than maxT, begin the simulation
  t = popTimes[0];
  popCounter++;
  unattached = 1;



  //set free variables from input
  c_trans = exp(p[0]);
  c_compart = exp(p[1]);
  c_circ = exp(p[2]);
  c_repP = exp(p[3]);
  c_repN = exp(p[4]);
  c_pack = exp(p[5]);
  MEMBRANE = exp(p[6]);
  repLimit = (unsigned long long)(exp(p[7]));
  consumption = round(exp(p[8]));
  SAH = p[9];
  scale=p[10];
  //scale=1;


  std::ofstream file;
  std::stringstream name;
  if(write==1){
    name << "/stor/work/Wilke/ateufel/Polio_track_final_redux/events" << "_"<< c_trans << "_" << c_compart << "_" << c_circ << "_" << c_repP << "_" << c_repN << "_" << c_pack << "_" << MEMBRANE << "_" << repLimit << "_" << consumption << "_" << SAH << "_" << scale << "rep_" << rep << ".txt";
    //name << "events.txt";
    std::string s = name.str();
    file.open(s.c_str());
    file << "t,protein, comparts,packaged,CAP,RNA_P,RNA_N,replicated,GFP" << std::endl;
  }
  int Timeoutcnt =0;
  while(t < maxT)
  {


    //added to spot lots of system calls
    Timeoutcnt++;

    //if( Timeoutcnt % 100000 == 0){
    end = clock();
    //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
    if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
      // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
      List ret;
      ret["GFP"] = Timeout;
      ret["rates"] = TimeoutRates;
      return(ret);
      //return(Timeout);
      //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

    }
    //}
    assert(RNA_P >= 0);
    assert(unattached >= 0);
    assert(RNA_N >= 0);


    //taxonomy of events
    // 0 -- attached + strand begins translation
    // 1 -- floating + strand attached to membrane
    // 2 -- + strand circularized -- by compartment!
    // 3 -- circularized strand begins replicating -- global
    // 4 -- neg. strand begins replicating -- global
    // 5 -- floating + strand begins translation
    // printf("t: %f",t);
    int event = 0;

    rates[0] = c_trans  * tComs; //had pacer here
    if(protein > 0) rates[1] = protein * c_compart * unattached * (double)(comparts - coms) / comparts;
    else rates[1] = 0;
    rates[2] = c_circ * cds;
    rates[3] = c_repP * RNA_P * pacer;
    rates[4] = c_repN * RNA_N * pacer;
    rates[5] = c_trans* unattached; //had pacer here
    assert(rates[1] >= 0);
    masterRate = rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5];


    // std::cout << "rates: " << rates[0] << " " << rates[1] << " " << rates[2] << " " << rates[3] << " " << rates[4] << " " << rates[5] << std::endl;
    if(masterRate != 0)
    {

      //keep track of sum of rates to make sure every reaction was possible durring the course of the sim.
      for(int i=0; i < 6; i++){
        SumRates[i]=SumRates[i]+rates[i];
      }

      if(masterRate > threshold) tau = 1 / masterRate;
      //else tau = gsl_ran_exponential(rgB, 1 / masterRate);
      else tau=R::rexp(1/masterRate);
      assert(masterRate > 0);
      //figure out event
      //double r = gsl_rng_uniform(rgB) * masterRate;
      double r = R::runif(0,1)*masterRate;
      while(r >= rates[event])
      {
        r -= rates[event];
        event++;
        assert(event < 6);


        end = clock();
        //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
        if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
          // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
          List ret;
          ret["GFP"] = Timeout;
          ret["rates"] = TimeoutRates;
          return(ret);
        }
      }
    }
    else
    {
      t = maxT;
      event = -1;
      tau = 0;
    }
    end = clock();
    //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
    if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
      // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
      List ret;
      ret["GFP"] = Timeout;
      ret["rates"] = TimeoutRates;
      return(ret);
      //return(Timeout);
      //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

    }


    if(popCounter < incoming && (t + tau) > popTimes[popCounter])
    {
      unattached++;
      t = popTimes[popCounter];
      popCounter++;

      while(t > sampleTimes[minuteCounter] && minuteCounter <= ((max_time)*2 -(min_time)*2))
      {

        end = clock();
        //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
        if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
          // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
          List ret;
          ret["GFP"] = Timeout;
          ret["rates"] = TimeoutRates;
          return(ret);
          //return(Timeout);
          //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

        }

        //value that you used to use was .0000025
        meanReplicated[minuteCounter] += replicated *scale;
        meanGFP[minuteCounter] += GFP *scale;

        file << minuteCounter << "," <<  protein << "," << comparts << "," << packaged << "," << CAP << "," << RNA_P << ","<< RNA_N << "," << replicated << "," << GFP << std::endl;
        minuteCounter++;

      }

    }
    else
    {


      while((t+tau) > sampleTimes[minuteCounter] && minuteCounter <= ((max_time)*2 -(min_time)*2))
      {
        end = clock();
        //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
        if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
          // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
          List ret;
          ret["GFP"] = Timeout;
          ret["rates"] = TimeoutRates;
          return(ret);
          //return(Timeout);
          //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

        }

        meanReplicated[minuteCounter] += replicated *scale;
        meanGFP[minuteCounter] += GFP *scale;
        file << minuteCounter << "," <<  protein << "," << comparts << "," << packaged << "," << CAP << "," << RNA_P << ","<< RNA_N << "," << replicated << "," << GFP << std::endl;
        minuteCounter++;

      }
      end = clock();
      // std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
      if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
        // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
        List ret;
        ret["GFP"] = Timeout;
        ret["rates"] = TimeoutRates;
        return(ret);
        //return(Timeout);
        //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

      }


      //deal with event
      switch(event)
      {
      case 0:
      {
        // Attached + strand begins translation.
        // Chooses a translating compartment uniformly.
        //int target = (int)gsl_rng_uniform_int(rgB, tComs);
        int target=R::runif(0,tComs-1); //should convert to int rather then the doulbe it produces
        int cursor = 0;
        while(target >= cs[cursor])
        {
          target -= cs[cursor];
          cursor++;

          end = clock();
          //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
          if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
            // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
            List ret;
            ret["GFP"] = Timeout;
            ret["rates"] = TimeoutRates;
            return(ret);
            //return(Timeout);
            //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

          }

        }
        cs[cursor]--;
        assert( (cursor + 1) < buffer);
        cs[cursor+1]++;
        cds++;
        CAP++;
        protein++;
        GFP++;



      }
        break;
      case 1:
      {
        // Floating + strand attached to membrane
        unattached--;
        cs[0]++;
        coms++;
        tComs++;
        protein -= consumption;

      }
        break;
      case 2:
      {
        // + strand circularized -- by compartment!
        // Choose a compartment with probability proportional to cs[$index]
        //int r = (int)gsl_rng_uniform_int(rgB, cds);
        int r = R::runif(0,cds-1);
        int cursor = 1;
        while(r >= (cs[cursor] * cursor))
        {
          r -= cs[cursor] * cursor;
          cursor++;
          end = clock();
          //std::cout << "time:" <<((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
          if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
            // std::cout << "time out " << Timeoutcnt << "pac:" << pacer <<  std::endl;
            List ret;
            ret["GFP"] = Timeout;
            ret["rates"] = TimeoutRates;
            return(ret);
            //return(Timeout);
            //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

          }

        }
        cds -= cursor;
        cs[cursor]--;
        RNA_P++;
        //GFP++;
        tComs--;

      }
        break;
      case 3:
        // Circularized strand begins replicating -- global
        replicated++;
        pacer -= decrement;
        if(pacer < 0) pacer = 0;
        RNA_N++;
        //GFP++;
        break;
      case 4:
        // Neg. strand begins replicating -- global
        // Encapsidate or not
        //if(gsl_rng_uniform(rgB) < (1 - exp(-c_pack * CAP)))
        if(R::runif(0,1) < (1 - exp(-c_pack * CAP)))
        {
          packaged++;
          CAP -= nC;
        }
        else
        {
          //if(gsl_rng_uniform(rgB) < SAH) RNA_P++;
          if(R::runif(0,1) < SAH) RNA_P++;
          else unattached++;
        }
        replicated++;
        pacer -= decrement;
        if(pacer < 0) pacer = 0;
        break;
      case 5:
        // Floating + strand begins translation
        CAP++;
        protein++;
        GFP++;
        break;
      }
      t += tau;
    }

  }
  file.close();


  //std::cout << " cnt: " << Timeoutcnt << std::endl;
  List ret;
  ret["GFP"] = meanGFP;
  ret["rates"] = SumRates;
  //std::cout << "pacer " << pacer << "cut: " << (.5) <<  std::endl;
  //replication limit has to be reached, we are timing things up to lysis, so it should max out
  //if(pacer > .5){
  //std::cout << " max not reached " << pacer << std::endl;
  //  ret["GFP"] = Timeout;
  // ret["rates"] = TimeoutRates;
  //}
  //std::cout << " max reached " << pacer << std::endl;
  return( ret );

}


// [[Rcpp::export]]
NumericVector avepa(NumericVector p,double MOI,double min_time, double max_time, int write){

  //int reps=25;
  int reps=1;

  NumericVector meanGFP(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector NoRates(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6

  for(int i = 0; i <= (((max_time)*2)-(min_time*2)); i++){
    meanGFP[i]=0;
    NoRates[i]=0;
  }

  NumericVector rates(5);
  for(int i = 0; i < 6; i++){
    rates[i]=0;
  }


  int reps_suc=0;
  //std::cout <<"do a rep" << std::endl;
  //do replicates of gillepsy sim
  for(int i =0; i < reps; i++){

    NumericVector temp(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6
    NumericVector temp_rates(6);

    for(int j = 0; j <= (((max_time)*2)-(min_time*2)); j++){
      temp[j]=0;
    }

    for(int j = 0; j < 6; j++){
      temp_rates[j]=0;
    }

    //printf("replicates");
    List t =  pa(p, MOI, min_time,  max_time, i, write);
    temp = t["GFP"];
    temp_rates = t["rates"];
    //std::cout << t << std::endl;

    int time_out=0;
    for(int j=0; j < temp.length(); j++){
      //std::cout << temp[j] << ","; //<< std::endl;
      if(temp[j]==0){
        time_out++;
      }
    }
    //std::cout<<std::endl;
    //std::cout << "time out?" << time_out << " length: " << temp.length() << std::endl;
    //check to see if everything was zero, if so it timed out and we shouldnt increase the count of what we are diving the replicats of
    if(time_out<temp.length()){
      reps_suc++;
      //std::cout << "good call" << std::endl;
    }
    if(time_out==temp.length()){
      // std::cout << "bad call" << std::endl;
      return(temp);
    }


    //keep track of sum
    for(int j=0; j < meanGFP.length(); j++){

      //std::cout << "temp[j] " << temp[j] << std::endl;
      meanGFP[j]=meanGFP[j]+temp[j];
    }

    //for(int j=0; j < 6; j++){
    //std::cout << "rates[j] " << rates[j] << std::endl;
    //rates[j]=rates[j]+ temp_rates[j];
    //}

  }

  //check if nothing happend
  int no_growth=0;
  for(int i = 0; i <= (((max_time)*2)-(min_time*2)); i++){
    if( meanGFP[i] !=0){
      no_growth++;
    }
  }

  //no growth happened, all the GFP was zero
  if(no_growth == 0){
    return(meanGFP);

  }


  //std::cout << "aver regsults" << std::endl;
  //ave results
  for(int j=0; j < meanGFP.length(); j++){
    // std::cout << meanGFP[j]/((double)reps_suc)<<",";
    meanGFP[j]=meanGFP[j]/((double)reps_suc);

  }
  //std::cout << std::endl;
  // printf("rates");
  //check if rates are non zero
  //for(int j=0; j<6; j++){
  // std::cout << "rates:" << rates[j] << std::endl;
  //if(rates[j] <=0){
  //return(NoRates);
  //}
  //}

  //if(rates[4] < rates[5]){
  //return(NoRates);
  //}

  return(meanGFP);


}



// [[Rcpp::export]]
NumericVector sig(NumericMatrix p0){
  NumericVector sigma(26);
  int inputs=p0.nrow();
  for(int i = 0; i < 26; i++)
  {
    double wMean = 0;
    double wSD = 0;
    for(int j = 0; j < inputs; j++)
    {
      wMean += p0(j,i);
    }
    wMean /= inputs;
    for(int j = 0; j < inputs; j++)
    {
      wSD += pow(p0(j,i) - wMean, 2);
    }
    wSD /= inputs;
    wSD = sqrt(wSD);
    sigma[i] = 1. * wSD;
  }
  return(sigma);
}






// [[Rcpp::export]]
SEXP track_mutants(NumericVector p, double MOI, double min_time, double max_time)
{
  //set up R seed so each call does not produce the same results. Remove if you want each call to be the same
  Environment g = Environment::global_env();
  Environment::Binding RandomSeed = g[".Random.seed"];

  clock_t start, end;
  start = clock();
  //FREE PARAMETERS
  //limiting constant for translation
  double c_trans;
  //constant for rate of compartmentaliztion
  double c_compart;
  //constant for 3CD-linear RNA reaction
  double c_circ;
  //constants for replication
  //Replication OF positives FROM negatives
  double c_repP;
  //Replication OF negatives FROM positives
  double c_repN;
  //encapsidation rate
  double c_pack;
  //number of possible compartments
  double MEMBRANE;
  //Limit constant for replication
  unsigned long long repLimit;
  //molecules of 3A consumed per compartment
  int consumption;
  //fraction of unencapsulated +strands that stay in their compartment
  double SAH;

  double scale;

  //FIXED PARAMETERS

  //FIXED PARAMETERS
  //polyproteins per capsid
  int nC = 60;
  //max time (minutes)
  int maxT = (max_time * 60)+5;
  double masterRate;
  NumericVector rates(5);
  NumericVector TimeoutRates(5);
  NumericVector SumRates(5);

  for(int i =0; i < 6; i++){
    SumRates[i] = 0;
    TimeoutRates[i]=0;
  }

  //set free variables from input
  c_trans = exp(p[0]);
  c_compart = exp(p[1]);
  c_circ = exp(p[2]);
  c_repP = exp(p[3]);//
  c_repN = exp(p[4]);
  c_pack = exp(p[5]);
  MEMBRANE = exp(p[6]);
  repLimit = (unsigned long long)(exp(p[7]));
  consumption = round(exp(p[8]));
  SAH = p[9];
  scale=p[10];


  NumericVector meanReplicated(((max_time*2)-(min_time*2))+1);//15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector sampleTimes(((max_time*2)-(min_time*2))+1);
  NumericVector meanPluses(((max_time*2)-(min_time*2))+1);
  NumericVector meanMinuses(((max_time*2)-(min_time*2))+1);
  NumericVector meanCap(((max_time*2)-(min_time*2))+1);
  NumericVector meanVirions(((max_time*2)-(min_time*2))+1);
  NumericVector Timeout(((max_time*2)-(min_time*2))+1);
  for(int i = 0; i <= (((max_time)*2)-(min_time*2)); i++)
  {

    meanReplicated[i]=0;
    meanPluses[i] = 0;
    meanMinuses[i] = 0;
    meanCap[i] = 0;
    meanVirions[i] = 0;
    sampleTimes[i]=60*((i+(min_time*2))*.5);
    Timeout[i] =0;
    // printf("same: %f",sampleTimes[i]);
  }

  // Sum # of positive strands and packaged genomes and their generation numbers
  double meanPositive = 0;
  double meanPackaged = 0;
  double totalPositive = 0;
  double totalPackaged = 0;
  double firstPackaged = 0;

  NumericVector totalP(((max_time*2)-(min_time*2))+1);
  NumericVector transP(((max_time*2)-(min_time*2))+1);
  NumericVector meanRate(((max_time*2)-(min_time*2))+1);
  for(int i = 0; i <= (((max_time)*2)-(min_time*2)); i++)
  {
    totalP[i] = 0;
    transP[i] = 0;
    meanRate[i] = 0;
  }


  //Concentration of protein used to make compartments
  long long protein = 0;
  int comparts = round(MEMBRANE);
  int packaged = 0;
  int CAP = 0;
  int RNA_P = 0;
  int RNA_N = 0;
  double t = 0;
  int minuteCounter = 0;
  int unattached = 0;
  int coms = 0;
  int tComs = 0;
  unsigned long long cds = 0;
  double tau;
  unsigned long long replicated = 0;
  int packed = 0;

  //double* popTimes = malloc(sizeof(double) * incoming);
  //double* popTimes=new double[incoming];
  NumericVector cs(comparts);
  NumericVector cg(comparts);

  //int* cs = (int*)(malloc(comparts * sizeof(int)));
  //int* cg = (int*)(malloc(comparts * sizeof(int)));
  for(int i = 0; i < comparts; i++)
  {
    cs[i] = 0;
    cg[i] = 0;
  }


  //unsigned long long una_cat [1000];
  //unsigned long long pos_cat [1000];
  //unsigned long long neg_cat [1000];
  NumericVector una_cat(1000); //they have 3 data points, we have 24. 3*8=24
  NumericVector pos_cat(1000);
  NumericVector neg_cat(1000);
  for(int i = 0; i < 1000; i++)
  {
    una_cat[i] = 0;
    pos_cat[i] = 0;
    neg_cat[i] = 0;
  }
  //std::cout << "made it past un int" << std::endl;

  //What are cg and cn?
  //cg stores the generation number of the +sense strand that founds a compartment
  // the purpose of storing cg for each compartment is to deal with the lag between formation of a compartment
  // and its transition to replication competency. this needs to be retained.
  // However, I don't currently see the need to implement cn as a vector of length comparts --surely we can just
  // retain a record based on class? And do the same thing for cp (positives), and use those for replication decisions.
  // This has the added bonus of dealing with SAH positives.

  double pacer = 1;
  double decrement = 1.0f / repLimit;

  int popCounter = 0;
  int incoming;
  do
  {
    //incoming = gsl_ran_poisson(rgB, MOIs[code]);
    incoming=R::rpois(MOI);
    //printf("in loop");
    //std::cout << "in come" << incoming << std::endl;
  }while(incoming == 0);
  //double* popTimes = malloc(sizeof(double) * incoming);
  //double* popTimes=new double[incoming];
  //double* popTimes = malloc(sizeof(double) * incoming);
  //double* popTimes=new double[incoming];
  NumericVector popTimes(incoming);
  for(int i = 0; i < incoming; i++)
  {
    //popTimes[i] = gsl_ran_gamma(rgB, SHAPE, 1.0f / RATE);
    popTimes[i] = R::rgamma(SHAPE, 1.0f / RATE);
    //printf("in second loop");
  }
  if(incoming > 1)
  {
    int changed = 0;
    do
    {
      changed = 0;
      for(int i = 1; i < incoming; i++)
      {
        if(popTimes[i] < popTimes[i-1])
        {
          changed++;
          double temp = popTimes[i];
          popTimes[i] = popTimes[i-1];
          popTimes[i-1] = temp;
        }
      }
    }while(changed > 0);
    //printf("in third loop");
  }
  t = popTimes[0];
  popCounter++;
  unattached = 1;
  una_cat[0]++;
  //printf("u 0 0\n");
  int Timeoutcnt =0;
  while(t < maxT)
  {


    //added to spot lots of system calls
    Timeoutcnt++;

    //if( Timeoutcnt % 100000 == 0){
    end = clock();
    std::cout << "time is : " << ((float)(end - start) / (float)(CLOCKS_PER_SEC)) << std::endl;
    if( ((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5){ //check if sim is running for too long, if so stop
      List ret;
      std::cout << "time out" << std::endl;
      ret["Gen"] = 0.0;
      ret["rates"] = TimeoutRates;
      return(ret);
      //return(Timeout);
      //return List::create(Named("GFP") = NumericVector::create(Timeout),Named("rates") = NumericVector::create(rates));

    }

    //printf("in sim");
    //std::cout << "t: " << t << "t max:"<< maxT << std::endl;
    assert(RNA_P >= 0);
    assert(unattached >= 0);
    assert(RNA_N >= 0);
    //taxonomy of events
    // 0 -- attached + strand begins translation
    // 1 -- floating + strand attached to membrane
    // 2 -- + strand circularized -- by compartment!
    // 3 -- circularized strand begins replicating -- global
    // 4 -- neg. strand begins replicating -- global
    // 5 -- floating + strand begins translation

    int event = 0;

    rates[0] = c_trans * tComs;
    if(protein > 0) rates[1] = protein * c_compart * unattached * (double)(comparts - coms) / comparts;
    else rates[1] = 0;
    rates[2] = c_circ * cds;
    rates[3] = c_repP * RNA_P * pacer;
    rates[4] = c_repN * RNA_N * pacer;
    rates[5] = c_trans * unattached;
    assert(rates[1] >= 0);


    masterRate = rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5];
    //std::cout << "rates: " << rates[0] << " " << rates[1] << " " << rates[2] << " " << rates[3] << " " << rates[4] << " " << rates[5] << std::endl;
    if(masterRate != 0)
    {

      for(int i=0; i < 6; i++){
        SumRates[i]=SumRates[i]+rates[i];
      }

      //std::cout << "masterRate:" << masterRate << " threshold:" << threshold << std::endl;
      if(masterRate > threshold) tau = 1 / masterRate;
      // else tau = gsl_ran_exponential(rgB, 1 / masterRate);
      else tau=R::rexp(1/masterRate);
      assert(masterRate > 0);
      //figure out event
      // double r = gsl_rng_uniform(rgB) * masterRate;
      double r = R::runif(0,1)*masterRate;
      //std::cout << "rand: " << r << std::endl;
      while(r >= rates[event])
      {
        r -= rates[event];
        event++;
        assert(event < 6);
      }
    }
    else
    {
      t = maxT;
      event = -1;
      tau = 0;
    }

    if(popCounter < incoming && (t + tau) > popTimes[popCounter])
    {
      //std::cout << "update 1 " << std::endl;
      unattached++;
      una_cat[0]++;
      t = popTimes[popCounter];
      //printf("u %f 0\n", t);
      popCounter++;
      while(t > sampleTimes[minuteCounter] && minuteCounter <= ((max_time)*2 -(min_time)*2))
      {
        meanPluses[minuteCounter] += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        meanMinuses[minuteCounter] += RNA_N;
        meanCap[minuteCounter] += CAP;
        totalP[minuteCounter] += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        transP[minuteCounter] += unattached + tComs;
        meanVirions[minuteCounter] += packaged;
        meanRate[minuteCounter] += rates[0] + rates[5];
        //printf("# %f %f %f %f %f %f\n", rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]);
        minuteCounter++;
      }
    }
    else
    {
      //std::cout << "update 2" << std::endl;
      while((t+tau) > sampleTimes[minuteCounter] && minuteCounter <= ((max_time)*2 -(min_time)*2))
      {
        meanPluses[minuteCounter] += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        meanMinuses[minuteCounter] += RNA_N;
        meanCap[minuteCounter] += CAP;
        totalP[minuteCounter] += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        transP[minuteCounter] += unattached + tComs;
        meanVirions[minuteCounter] += packaged;
        meanRate[minuteCounter] += rates[0] + rates[5];
        //printf("# %f %f %f %f %f %f\n", rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]);
        minuteCounter++;
        //std::cout << "min count" << minuteCounter << std::endl;
      }
      //std::cout << "switch:" << event << std::endl;
      //deal with event
      switch(event)
      {
        //std::cout << "in switch:" << event << std::endl;
      case 0:
      {
        // std::cout << "case 0 " << std::endl;
        int target;
        do
        {

          target=R::runif(0,coms); //should convert to int rather then the doulbe it produces
          // target=R::runif(0,tComs); //should convert to int rather then the doulbe it produces


        }while(cs[target] == -1);
        cs[target]++;
        cds++;
        CAP++;
        protein++;
      }
        break;


      case 1:
      {
        //std::cout << "case 1" << std::endl;
        unsigned long long  r=R::runif(0,unattached-1); //should convert to int rather then the doulbe it produces
        int cursor = 0;
        //std::cout << "r: " << r << std::endl;
        while(r >= una_cat[cursor])
        {
          //std::cout << "r: " << r << " una-cat" << una_cat[cursor] <<std::endl;
          r -= una_cat[cursor];
          cursor++;
          assert(cursor < 1000);
        }
        cg[coms] = cursor;
        una_cat[cursor]--;
        unattached--;
        assert(cs[coms] == 0);
        coms++;
        tComs++;
        protein -= consumption;
      }
        break;
      case 2:
      {
        //std::cout << "case 2 " << std::endl;
        //int r = (int)gsl_rng_uniform_int(rgB, cds);
        int r = R::runif(0,cds-1);
        int cursor = 0;
        while(cs[cursor] == -1 || r >= cs[cursor])
        {
          if(cs[cursor] != -1) r -= cs[cursor];
          cursor++;
          assert(cursor < coms);
        }
        assert(cs[cursor] > 0);
        cds -= cs[cursor];
        cs[cursor] = -1;
        RNA_P++;
        pos_cat[cg[cursor]]++;
        //printf("p %f %d\n", t+tau, cg[cursor]);
        tComs--;
      }
        break;
      case 3:
      {
        //std::cout << "case 3" << std::endl;
        replicated++;
        pacer -= decrement;
        if(pacer < 0) pacer = 0;
        RNA_N++;
        //unsigned long long r = gsl_rng_uniform_int(rgB, RNA_P);
        unsigned long long  r=R::runif(0,RNA_P-1);
        int cursor = 0;
        while(r >= pos_cat[cursor])
        {
          //std::cout << "1r:" << r << std::endl;
          //std::cout << "1cursor" << cursor << std::endl;
          r -= pos_cat[cursor];
          cursor++;
          assert(cursor < 1000);
        }
        neg_cat[cursor + 1]++;
        //printf("n %f %d\n", t+tau, cursor + 1);
      }
        break;
      case 4:
      {
        //std::cout << "case 4 " << std::endl;
        //unsigned long long r = gsl_rng_uniform_int(rgB, RNA_N);
        unsigned long long  r=R::runif(0,RNA_N-1);
        int cursor = 0;
        while(r >= neg_cat[cursor])
        {
          r -= neg_cat[cursor];
          cursor++;
          assert(cursor < 1000);
        }
        //encapsidate or not
        //if(gsl_rng_uniform(rgB) < (1 - exp(-c_pack * CAP)))
        if(R::runif(0,1) < (1 - exp(-c_pack * CAP)))
        {
          if(packed == 0)
          {
            firstPackaged += t + tau;
            packed = 1;
          }
          totalPackaged++;
          //printf("%d %d \n", pCode, cursor + 1);
          meanPackaged += cursor + 1;
          packaged++;
          CAP -= nC;
          //printf("v %f %d\n", t+tau, cursor + 1);
        }
        else
        {
          //if(gsl_rng_uniform(rgB) < SAH)
          if(R::runif(0,1) < SAH)
          {
            RNA_P++;
            pos_cat[cursor + 1]++;
            //printf("p %f %d\n", t+tau, cursor + 1);
          }
          else
          {
            unattached++;
            una_cat[cursor + 1]++;
            //printf("u %f %d\n", t+tau, cursor + 1);
          }
        }
        replicated++;
        pacer -= decrement;
        if(pacer < 0) pacer = 0;
      }
        break;
      case 5:
        //std::cout <<"case 5" << std::endl;
        CAP++;
        protein++;
        break;
      }
      t += tau;
      //std::cout << "tau is:"<< tau << std::endl;
    }

  }
  //std::cout << "made it past sim" << std::endl;
  for(int i = 0; i < 1000; i++)
  {
    //std::cout << "made it poscat"  <<pos_cat[i] << std::endl;
    //std::cout << "made it una"  <<una_cat[i] << std::endl;
    totalPositive += pos_cat[i] + una_cat[i];
    meanPositive += (pos_cat[i] + una_cat[i]) * (i + 1);
  }
  totalPositive += totalPackaged;
  meanPositive += meanPackaged;
  //free(cs);
  //free(cg);
  //free(popTimes);
  //std::cout << "rates: " << rates[0] << " " << rates[1] << " " << rates[2] << " " << rates[3] << " " << rates[4] << " " << rates[5] << std::endl;
  meanPackaged /= totalPackaged;
  meanPositive /= totalPositive;
  // if(MODE == 1 || MODE == 3 || MODE == 4)
  printf(" %f %f %f %f %f \n",meanPackaged, meanPositive, totalPackaged, totalPositive, firstPackaged);

  if(NumericVector::is_na(meanPackaged)){

    meanPackaged = 0;
  }
  std::cout << " mean packaged: " << meanPackaged << std::endl;
  List ret;
  ret["Gen"] = meanPackaged;
  ret["rates"] = SumRates;



  return( ret );

}


// [[Rcpp::export]]
double avetrack(NumericVector p,double MOI,double min_time, double max_time){

  int reps=1;
  double track=0;
  double non_zero=0;

  NumericVector rates(5);
  for(int i = 0; i < 6; i++){
    rates[i]=0;
  }


  NumericVector temp_rates(5);

  //do replicates of gillepsy sim
  printf("ref");
  for(int i =0; i < reps; i++){

    //refill up temp rates
    for(int j = 0; j < 6; j++){
      temp_rates[j]=0;
    }

    std::cout << "calling function" << std::endl;
    List t =  track_mutants(p, MOI, min_time,  max_time);
    std::cout << "function called" << std::endl;
    double temp = t["Gen"];
    temp_rates = t["rates"];
    track=temp;

  }


  return(track);

}
