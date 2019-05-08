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

#define threshold (1000)
// Fitted parameters for absorption time distribution
#define SHAPE (0.678)
#define RATE (0.02)


// [[Rcpp::export]]
SEXP pa(NumericVector p, double MOI, double min_time, double max_time, int rep, int write)
{

  // set up R seed so each call does not produce the same results. Remove if you want each call to
  // be the same
  Environment g = Environment::global_env();
  Environment::Binding RandomSeed = g[".Random.seed"];

  clock_t start, end;
  start = clock();
  // FREE PARAMETERS
  // limiting constant for translation
  double c_trans;
  // constant for rate of compartmentaliztion
  double c_compart;
  // constant for 3CD-linear RNA reaction
  double c_circ;
  // constants for replication
  // Replication OF positives FROM negatives
  double c_repP;
  // Replication OF negatives FROM positives
  double c_repN;
  // encapsidation rate
  double c_pack;
  // number of possible compartments
  double MEMBRANE;
  // Limit constant for replication
  unsigned long long repLimit;
  // molecules of 3A consumed per compartment
  int consumption;
  // fraction of unencapsulated +strands that stay in their compartment
  double SAH;
  // scaler to convert gfp to RNA
  double scale;

  // FIXED PARAMETERS
  // polyproteins per capsid
  int nC = 60;
  // max time (minutes)
  int maxT = (max_time * 60) + 5;
  double masterRate;
  // double rates [6];
  NumericVector rates(5);

  // set free variables from input
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
  scale = p[10];


  NumericVector meanReplicated(((max_time * 2) - (min_time * 2))
                                 + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector meanGFP(((max_time * 2) - (min_time * 2))
                          + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector sampleTimes(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector Timeout(((max_time * 2) - (min_time * 2))
                          + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  for (int i = 0; i <= (((max_time)*2) - (min_time * 2)); i++)
  {

    meanReplicated[i] = 0;
    meanGFP[i] = 0;
    Timeout[i] = 0;
    sampleTimes[i] = 60 * ((i + (min_time * 2)) * .5);
  }

  NumericVector TimeoutRates(5);
  NumericVector SumRates(5);
  for (int i = 0; i < 6; i++)
  {
    TimeoutRates[i] = 0;
    SumRates[i] = 0;
  }

  // Concentration of protein used to make compartments
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

  NumericVector cs(buffer);
  for (int i = 0; i < buffer; i++)
    cs[i] = 0;

  // Used to efficiently scale replication rate to logistic limit
  double pacer = 1;
  double decrement = 1.0f / repLimit;

  // Assigns each infecting virus to a timepoint and sorts those into a schedule
  int popCounter = 0;
  int incoming;
  do
  {

    incoming = R::rpois(MOI);


  } while (incoming == 0);


  NumericVector popTimes(incoming);
  for (int i = 0; i < incoming; i++)
  {
    popTimes[i] = R::rgamma(SHAPE, 1.0f / RATE);
  }
  if (incoming > 1)
  {
    int changed = 0;
    do
    {
      changed = 0;
      for (int i = 1; i < incoming; i++)
      {
        if (popTimes[i] < popTimes[i - 1])
        {
          changed++;
          double temp = popTimes[i];
          popTimes[i] = popTimes[i - 1];
          popTimes[i - 1] = temp;
        }
      }
    } while (changed > 0);
  }

  // Advance to the first infection time and, if less than maxT, begin the simulation
  t = popTimes[0];
  popCounter++;
  unattached = 1;


  // set free variables from input
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
  scale = p[10];

  // for tracking event timing
  std::ofstream file;
  std::stringstream name;
  if (write == 1)
  {
    name << "/stor/work/Wilke/ateufel/Polio_track_final/events"
         << "_" << c_trans << "_" << c_compart << "_" << c_circ << "_" << c_repP << "_"
         << c_repN << "_" << c_pack << "_" << MEMBRANE << "_" << repLimit << "_" << consumption
         << "_" << SAH << "_" << scale << "rep_" << rep << ".txt";
    std::string s = name.str();
    file.open(s.c_str());
    file << "t,protein, comparts,packaged,CAP,RNA_P,RNA_N,replicated,GFP" << std::endl;
  }
  int Timeoutcnt = 0;
  while (t < maxT)
  {
    // added to spot lots of system calls
    Timeoutcnt++;

    if (Timeoutcnt % 100000 == 0)
    {
      end = clock();

      if (((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5)
      { // check if sim is running for too long, if so stop
        List ret;
        ret["GFP"] = Timeout;
        ret["rates"] = TimeoutRates;
        return (ret);
      }
    }
    assert(RNA_P >= 0);
    assert(unattached >= 0);
    assert(RNA_N >= 0);

    // taxonomy of events
    // 0 -- attached + strand begins translation
    // 1 -- floating + strand attached to membrane
    // 2 -- + strand circularized -- by compartment!
    // 3 -- circularized strand begins replicating -- global
    // 4 -- neg. strand begins replicating -- global
    // 5 -- floating + strand begins translation
    // printf("t: %f",t);
    int event = 0;

    rates[0] = c_trans * tComs;
    if (protein > 0)
      rates[1] = protein * c_compart * unattached * (double)(comparts - coms) / comparts;
    else
      rates[1] = 0;
    rates[2] = c_circ * cds;
    rates[3] = c_repP * RNA_P * pacer;
    rates[4] = c_repN * RNA_N * pacer;
    rates[5] = c_trans * unattached;
    assert(rates[1] >= 0);
    masterRate = rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5];

    if (masterRate != 0)
    {

      // keep track of sum of rates to make sure every reaction was possible durring the
      // course of the sim.
      for (int i = 0; i < 6; i++)
      {
        SumRates[i] = SumRates[i] + rates[i];
      }

      if (masterRate > threshold)
        tau = 1 / masterRate;

      else
        tau = R::rexp(1 / masterRate);
      assert(masterRate > 0);

      double r = R::runif(0, 1) * masterRate;
      while (r >= rates[event])
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


    if (popCounter < incoming && (t + tau) > popTimes[popCounter])
    {
      unattached++;
      t = popTimes[popCounter];
      popCounter++;

      while (t > sampleTimes[minuteCounter] && minuteCounter <= ((max_time)*2 - (min_time)*2))
      {

        meanReplicated[minuteCounter] += replicated * scale;
        meanGFP[minuteCounter] += GFP * scale;

        file << minuteCounter << "," << protein << "," << comparts << "," << packaged << ","
             << CAP << "," << RNA_P << "," << RNA_N << "," << replicated << "," << GFP
             << std::endl;
        minuteCounter++;
      }
    }
    else
    {


      while ((t + tau) > sampleTimes[minuteCounter]
               && minuteCounter <= ((max_time)*2 - (min_time)*2))
      {
        meanReplicated[minuteCounter] += replicated * scale;
        meanGFP[minuteCounter] += GFP * scale;
        file << minuteCounter << "," << protein << "," << comparts << "," << packaged << ","
             << CAP << "," << RNA_P << "," << RNA_N << "," << replicated << "," << GFP
             << std::endl;
        minuteCounter++;
      }

      // deal with event
      switch (event)
      {
      case 0:
      {
        // Attached + strand begins translation.
        // Chooses a translating compartment uniformly.
        // int target = (int)gsl_rng_uniform_int(rgB, tComs);
        int target = R::runif(
          0, tComs - 1); // should convert to int rather then the doulbe it produces
        int cursor = 0;
        while (target >= cs[cursor])
        {
          target -= cs[cursor];
          cursor++;
        }
        cs[cursor]--;
        assert((cursor + 1) < buffer);
        cs[cursor + 1]++;
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
        int r = R::runif(0, cds - 1);
        int cursor = 1;
        while (r >= (cs[cursor] * cursor))
        {
          r -= cs[cursor] * cursor;
          cursor++;
        }
        cds -= cursor;
        cs[cursor]--;
        RNA_P++;
        tComs--;
      }
        break;
      case 3:
        // Circularized strand begins replicating -- global
        replicated++;
        pacer -= decrement;
        if (pacer < 0)
          pacer = 0;
        RNA_N++;
        break;
      case 4:
        // Neg. strand begins replicating -- global
        // Encapsidate or not
        if (R::runif(0, 1) < (1 - exp(-c_pack * CAP)))
        {
          packaged++;
          CAP -= nC;
        }
        else
        {
          if (R::runif(0, 1) < SAH)
            RNA_P++;
          else
            unattached++;
        }
        replicated++;
        pacer -= decrement;
        if (pacer < 0)
          pacer = 0;
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


  List ret;
  ret["GFP"] = meanGFP;
  ret["rates"] = SumRates;

  // replication limit has to be reached, we are timing things up to lysis, so it should max out
  if (pacer > .5)
  {
    ret["GFP"] = Timeout;
    ret["rates"] = TimeoutRates;
  }
  return (ret);
}


// [[Rcpp::export]]
NumericVector avepa(NumericVector p, double MOI, double min_time, double max_time, int write)
{

  int reps = 25;

  NumericVector meanGFP(((max_time * 2) - (min_time * 2))
                          + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector NoRates(((max_time * 2) - (min_time * 2))
                          + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6

  for (int i = 0; i <= (((max_time)*2) - (min_time * 2)); i++)
  {
    meanGFP[i] = 0;
    NoRates[i] = 0;
  }

  NumericVector rates(5);
  for (int i = 0; i < 6; i++)
  {
    rates[i] = 0;
  }


  int reps_suc = 0;
  // do replicates of gillepsy sim
  for (int i = 0; i < reps; i++)
  {

    NumericVector temp(((max_time * 2) - (min_time * 2))
                         + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6
    NumericVector temp_rates(6);

    for (int j = 0; j <= (((max_time)*2) - (min_time * 2)); j++)
    {
      temp[j] = 0;
    }

    for (int j = 0; j < 6; j++)
    {
      temp_rates[j] = 0;
    }

    List t = pa(p, MOI, min_time, max_time, i, write);
    temp = t["GFP"];
    temp_rates = t["rates"];

    int time_out = 0;
    for (int j = 0; j < temp.length(); j++)
    {
      if (temp[j] == 0)
      {
        time_out++;
      }
    }

    // check to see if everything was zero, if so it timed out and we shouldnt increase the
    // count of what we are diving the replicats of
    if (time_out < temp.length())
    {
      reps_suc++;
    }
    if (time_out == temp.length())
    {
      return (temp);
    }
    // keep track of sum
    for (int j = 0; j < meanGFP.length(); j++)
    {
      meanGFP[j] = meanGFP[j] + temp[j];
    }
  }
  // check if nothing happend
  int no_growth = 0;
  for (int i = 0; i <= (((max_time)*2) - (min_time * 2)); i++)
  {
    if (meanGFP[i] != 0)
    {
      no_growth++;
    }
  }

  // no growth happened, all the GFP was zero
  if (no_growth == 0)
  {
    return (meanGFP);
  }

  // ave results
  for (int j = 0; j < meanGFP.length(); j++)
  {
    meanGFP[j] = meanGFP[j] / ((double)reps_suc);
  }


  return (meanGFP);
}

// [[Rcpp::export]]
NumericVector sig(NumericMatrix p0)
{
  NumericVector sigma(26);
  int inputs = p0.nrow();
  for (int i = 0; i < 26; i++)
  {
    double wMean = 0;
    double wSD = 0;
    for (int j = 0; j < inputs; j++)
    {
      wMean += p0(j, i);
    }
    wMean /= inputs;
    for (int j = 0; j < inputs; j++)
    {
      wSD += pow(p0(j, i) - wMean, 2);
    }
    wSD /= inputs;
    wSD = sqrt(wSD);
    sigma[i] = 1. * wSD;
  }
  return (sigma);
}


// [[Rcpp::export]]
SEXP track_mutants(NumericVector p, double MOI, double min_time, double max_time)
{
  // set up R seed so each call does not produce the same results. Remove if you want each call to
  // be the same
  Environment g = Environment::global_env();
  Environment::Binding RandomSeed = g[".Random.seed"];

  clock_t start, end;
  start = clock();
  // FREE PARAMETERS
  // limiting constant for translation
  double c_trans;
  // constant for rate of compartmentaliztion
  double c_compart;
  // constant for 3CD-linear RNA reaction
  double c_circ;
  // constants for replication
  // Replication OF positives FROM negatives
  double c_repP;
  // Replication OF negatives FROM positives
  double c_repN;
  // encapsidation rate
  double c_pack;
  // number of possible compartments
  double MEMBRANE;
  // Limit constant for replication
  unsigned long long repLimit;
  // molecules of 3A consumed per compartment
  int consumption;
  // fraction of unencapsulated +strands that stay in their compartment
  double SAH;

  double scale;

  // FIXED PARAMETERS

  // FIXED PARAMETERS
  // polyproteins per capsid
  int nC = 60;
  // max time (minutes)
  int maxT = (max_time * 60) + 5;
  double masterRate;
  NumericVector rates(5);
  NumericVector TimeoutRates(5);
  NumericVector SumRates(5);

  for (int i = 0; i < 6; i++)
  {
    SumRates[i] = 0;
    TimeoutRates[i] = 0;
  }

  // set free variables from input
  c_trans = exp(p[0]);
  c_compart = exp(p[1]);
  c_circ = exp(p[2]);
  c_repP = exp(p[3]); //
  c_repN = exp(p[4]);
  c_pack = exp(p[5]);
  MEMBRANE = exp(p[6]);
  repLimit = (unsigned long long)(exp(p[7]));
  consumption = round(exp(p[8]));
  SAH = p[9];
  scale = p[10];


  NumericVector meanReplicated(((max_time * 2) - (min_time * 2))
                                 + 1); // 15 hours, by the half hour starting at hour 3 thats why (24*2)-6
  NumericVector sampleTimes(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector meanPluses(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector meanMinuses(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector meanCap(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector meanVirions(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector Timeout(((max_time * 2) - (min_time * 2)) + 1);
  for (int i = 0; i <= (((max_time)*2) - (min_time * 2)); i++)
  {

    meanReplicated[i] = 0;
    meanPluses[i] = 0;
    meanMinuses[i] = 0;
    meanCap[i] = 0;
    meanVirions[i] = 0;
    sampleTimes[i] = 60 * ((i + (min_time * 2)) * .5);
    Timeout[i] = 0;
  }

  // Sum # of positive strands and packaged genomes and their generation numbers
  double meanPositive = 0;
  double meanPackaged = 0;
  double totalPositive = 0;
  double totalPackaged = 0;
  double firstPackaged = 0;

  NumericVector totalP(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector transP(((max_time * 2) - (min_time * 2)) + 1);
  NumericVector meanRate(((max_time * 2) - (min_time * 2)) + 1);
  for (int i = 0; i <= (((max_time)*2) - (min_time * 2)); i++)
  {
    totalP[i] = 0;
    transP[i] = 0;
    meanRate[i] = 0;
  }


  // Concentration of protein used to make compartments
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


  NumericVector cs(comparts);
  NumericVector cg(comparts);

  for (int i = 0; i < comparts; i++)
  {
    cs[i] = 0;
    cg[i] = 0;
  }


  NumericVector una_cat(1000);
  NumericVector pos_cat(1000);
  NumericVector neg_cat(1000);
  for (int i = 0; i < 1000; i++)
  {
    una_cat[i] = 0;
    pos_cat[i] = 0;
    neg_cat[i] = 0;
  }


  double pacer = 1;
  double decrement = 1.0f / repLimit;

  int popCounter = 0;
  int incoming;
  do
  {

    incoming = R::rpois(MOI);

  } while (incoming == 0);

  NumericVector popTimes(incoming);
  for (int i = 0; i < incoming; i++)
  {

    popTimes[i] = R::rgamma(SHAPE, 1.0f / RATE);
  }
  if (incoming > 1)
  {
    int changed = 0;
    do
    {
      changed = 0;
      for (int i = 1; i < incoming; i++)
      {
        if (popTimes[i] < popTimes[i - 1])
        {
          changed++;
          double temp = popTimes[i];
          popTimes[i] = popTimes[i - 1];
          popTimes[i - 1] = temp;
        }
      }
    } while (changed > 0);
  }
  t = popTimes[0];
  popCounter++;
  unattached = 1;
  una_cat[0]++;
  int Timeoutcnt = 0;
  while (t < maxT)
  {

    Timeoutcnt++;

    if (Timeoutcnt % 100000 == 0)
    {
      end = clock();

      if (((float)(end - start) / (float)(CLOCKS_PER_SEC)) > 5)
      { // check if sim is running for too long, if so stop
        List ret;

        ret["Gen"] = 0.0;
        ret["rates"] = TimeoutRates;
        return (ret);
      }
    }

    assert(RNA_P >= 0);
    assert(unattached >= 0);
    assert(RNA_N >= 0);
    // taxonomy of events
    // 0 -- attached + strand begins translation
    // 1 -- floating + strand attached to membrane
    // 2 -- + strand circularized -- by compartment!
    // 3 -- circularized strand begins replicating -- global
    // 4 -- neg. strand begins replicating -- global
    // 5 -- floating + strand begins translation

    int event = 0;

    rates[0] = c_trans * tComs;
    if (protein > 0)
      rates[1] = protein * c_compart * unattached * (double)(comparts - coms) / comparts;
    else
      rates[1] = 0;
    rates[2] = c_circ * cds;
    rates[3] = c_repP * RNA_P * pacer;
    rates[4] = c_repN * RNA_N * pacer;
    rates[5] = c_trans * unattached;
    assert(rates[1] >= 0);


    masterRate = rates[0] + rates[1] + rates[2] + rates[3] + rates[4] + rates[5];
    if (masterRate != 0)
    {

      for (int i = 0; i < 6; i++)
      {
        SumRates[i] = SumRates[i] + rates[i];
      }

      if (masterRate > threshold)
        tau = 1 / masterRate;

      else
        tau = R::rexp(1 / masterRate);
      assert(masterRate > 0);

      double r = R::runif(0, 1) * masterRate;

      while (r >= rates[event])
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

    if (popCounter < incoming && (t + tau) > popTimes[popCounter])
    {

      unattached++;
      una_cat[0]++;
      t = popTimes[popCounter];

      popCounter++;
      while (t > sampleTimes[minuteCounter] && minuteCounter <= ((max_time)*2 - (min_time)*2))
      {
        meanPluses[minuteCounter]
        += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        meanMinuses[minuteCounter] += RNA_N;
        meanCap[minuteCounter] += CAP;
        totalP[minuteCounter]
        += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        transP[minuteCounter] += unattached + tComs;
        meanVirions[minuteCounter] += packaged;
        meanRate[minuteCounter] += rates[0] + rates[5];

        minuteCounter++;
      }
    }
    else
    {

      while ((t + tau) > sampleTimes[minuteCounter]
               && minuteCounter <= ((max_time)*2 - (min_time)*2))
      {
        meanPluses[minuteCounter]
        += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        meanMinuses[minuteCounter] += RNA_N;
        meanCap[minuteCounter] += CAP;
        totalP[minuteCounter]
        += unattached + packaged + RNA_P + tComs + incoming - popCounter;
        transP[minuteCounter] += unattached + tComs;
        meanVirions[minuteCounter] += packaged;
        meanRate[minuteCounter] += rates[0] + rates[5];

        minuteCounter++;
      }

      switch (event)
      {

      case 0:
      {

        int target;
        do
        {

          target = R::runif(
            0, coms); // should convert to int rather then the doulbe it produces


        } while (cs[target] == -1);
        cs[target]++;
        cds++;
        CAP++;
        protein++;
      }
        break;


      case 1:
      {

        unsigned long long r = R::runif(0,
                                        unattached - 1); // should convert to int rather then the doulbe it produces
        int cursor = 0;

        while (r >= una_cat[cursor])
        {

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

        int r = R::runif(0, cds - 1);
        int cursor = 0;
        while (cs[cursor] == -1 || r >= cs[cursor])
        {
          if (cs[cursor] != -1)
            r -= cs[cursor];
          cursor++;
          assert(cursor < coms);
        }
        assert(cs[cursor] > 0);
        cds -= cs[cursor];
        cs[cursor] = -1;
        RNA_P++;
        pos_cat[cg[cursor]]++;

        tComs--;
      }
        break;
      case 3:
      {

        replicated++;
        pacer -= decrement;
        if (pacer < 0)
          pacer = 0;
        RNA_N++;

        unsigned long long r = R::runif(0, RNA_P - 1);
        int cursor = 0;
        while (r >= pos_cat[cursor])
        {

          r -= pos_cat[cursor];
          cursor++;
          assert(cursor < 1000);
        }
        neg_cat[cursor + 1]++;
      }
        break;
      case 4:
      {

        unsigned long long r = R::runif(0, RNA_N - 1);
        int cursor = 0;
        while (r >= neg_cat[cursor])
        {
          r -= neg_cat[cursor];
          cursor++;
          assert(cursor < 1000);
        }
        // encapsidate or not

        if (R::runif(0, 1) < (1 - exp(-c_pack * CAP)))
        {
          if (packed == 0)
          {
            firstPackaged += t + tau;
            packed = 1;
          }
          totalPackaged++;

          meanPackaged += cursor + 1;
          packaged++;
          CAP -= nC;
        }
        else
        {

          if (R::runif(0, 1) < SAH)
          {
            RNA_P++;
            pos_cat[cursor + 1]++;
          }
          else
          {
            unattached++;
            una_cat[cursor + 1]++;
          }
        }
        replicated++;
        pacer -= decrement;
        if (pacer < 0)
          pacer = 0;
      }
        break;
      case 5:

        CAP++;
        protein++;
        break;
      }
      t += tau;
    }
  }

  for (int i = 0; i < 1000; i++)
  {

    totalPositive += pos_cat[i] + una_cat[i];
    meanPositive += (pos_cat[i] + una_cat[i]) * (i + 1);
  }
  totalPositive += totalPackaged;
  meanPositive += meanPackaged;
  meanPackaged /= totalPackaged;
  meanPositive /= totalPositive;

  List ret;
  ret["Gen"] = meanPackaged;
  ret["rates"] = SumRates;

  if (pacer > .5)
  {

    ret["Gen"] = 0.0;
    ret["rates"] = TimeoutRates;
  }

  return (ret);
}


// [[Rcpp::export]]
double avetrack(NumericVector p, double MOI, double min_time, double max_time)
{

  int reps = 25;
  double track = 0;
  double non_zero = 0;

  NumericVector rates(5);
  for (int i = 0; i < 6; i++)
  {
    rates[i] = 0;
  }


  NumericVector temp_rates(5);

  // do replicates of gillepsy sim

  for (int i = 0; i < reps; i++)
  {

    // refill up temp rates
    for (int j = 0; j < 6; j++)
    {
      temp_rates[j] = 0;
    }


    List t = track_mutants(p, MOI, min_time, max_time);

    double temp = t["Gen"];
    temp_rates = t["rates"];


    if (temp > 2.0)
    {
      track = track + temp;
      non_zero++;
    }


    for (int k = 0; k < 6; k++)
    {

      rates[k] = rates[k] + temp_rates[k];
    }
  }

  if (non_zero == 0)
  {
    return (0.0);
  }

  return (track / non_zero);
}
