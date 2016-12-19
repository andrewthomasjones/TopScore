//#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <math.h>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <valarray>
#include <regex>
#include <fstream>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
typedef std::vector<double> stdvec;
typedef std::pair<int,double> mypair;
//calculate alpha
// [[Rcpp::export]]
double alphaCalc(double error){
    return(log((1-error)/error));
}

// [[Rcpp::export]]
vector<int> sort_indexes(vector<double> v) {

  // initialize original index locations
  vector<int> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}




// [[Rcpp::export]]
mat getData(NumericMatrix obsvData, NumericVector  labels,  int n_pair, NumericVector  weights, bool disjoint_mode){
    // creates Rcpp matrix from SEXP
    int p = obsvData.nrow(), n = obsvData.ncol();
    mat matData(obsvData.begin(), p, n, false);
    
    //creates vector in same way
    int vec_len = labels.size();
    int w_len = labels.size();
    vec lab(labels.begin(), vec_len , false);
    vec weigh(weights.begin(), w_len , false);
    int nClass = lab.max();
    
    if (vec_len != n | w_len != n) {           
        stop("Label vector must be same length as number of obsverations (columns)");
    }
    
    //split by class
    field<mat> dataByClass(nClass);
    vec obsv_cl = zeros<vec>(nClass);
    
    for(int i=0; i < nClass; i++){
      uvec ab = find(lab==(i+1));    
      dataByClass(i) = matData.cols(ab(0),ab(ab.n_rows-1));
      obsv_cl(i) = dataByClass(i).n_cols;
    }
        
    //initialise storage
    field<mat> probByClass(nClass);
    for(int i=0; i < nClass; i++){
       probByClass(i) = zeros<mat>(p,p);
    }
    
    //obsvation in class
      
     //initialise storage tie
    field<mat> tieBreak(nClass);
    for(int i=0; i < nClass; i++){
       tieBreak(i) = zeros<mat>(p,p);
    }
    
     //initialise storage rank
    field<umat> rankByClass(nClass);
    for(int i=0; i < nClass; i++){
       rankByClass(i) = zeros<umat>(p,p);
    }
    
    //get ranks
    for(int i=0; i < nClass; i++){
      for(int m=0; m< obsv_cl(i); m++){
          vec temp = dataByClass(i).col(m);
          stdvec z = conv_to< stdvec >::from(temp);
          vector<int> sort = sort_indexes(z);
          rankByClass(i).col(m) = arma::conv_to<arma::uvec>::from(sort);      
      }
    }

    //uvec test = (dataByClass(0).col(5) > dataByClass(0)(10,5));
    //test.print();
    
    
    for(int k=0; k< nClass; k++){
      
      //each predictor
      for(int i=0; i< p; i++){
                
            for(int m=0; m< obsv_cl(k); m++){
            //fill out a col at a time
           
            //for tie breaker - possibly omit for speed
            mat temp2 = arma::conv_to<arma::mat>::from(rankByClass(k).col(m));
            double temp3 = (double)(rankByClass(k)(i,m));
            mat temp4 =  (temp2-temp3); 
            tieBreak(k).col(i) = tieBreak(k).col(i) +  temp4 ;
            //for prob
            umat temp = (dataByClass(k).col(m) > dataByClass(k)(i,m));
            probByClass(k).col(i) =  probByClass(k).col(i) + weigh(m)*arma::conv_to<arma::mat>::from(temp)  ;
          
          }
          
          
      }
      
    }
    //calculate score and tie breaker
    mat delta = zeros<mat>(p,p);
    mat gamma = zeros<mat>(p,p);
    delta = ((arma::conv_to<arma::mat>::from(probByClass(0))/obsv_cl(0)) - (arma::conv_to<arma::mat>::from(probByClass(1))/obsv_cl(1)));
    gamma = ((arma::conv_to<arma::mat>::from(tieBreak(0))/obsv_cl(0)) - (arma::conv_to<arma::mat>::from(tieBreak(1))/obsv_cl(1)));
    
    //cut to triangle to avoid doubles/confusion
    mat L = trimatl(abs(delta));
    uword pair1, pair2;
    L.max(pair1, pair2);
    
    //get top pairs
    vec Lvec = vectorise(L);
    stdvec z = conv_to< stdvec >::from(Lvec);
    vector<int> indices (n_pair);
    std::fill (indices.begin(),indices.end(),-1);
    //sort and get indexes, sort_index in arma produces weird results!
    vector<int> test = sort_indexes(z);

    
    // prep output
    mat output(n_pair,4, fill::zeros);
    div_t divresult, divresult2;
    
    //list of already used pairs
    list<int> pairList;
        
    int i=0; int j=0;
    
    while(i<n_pair){
      list<double> scoreList,scoreListCopy;
      div_t div1,div2;

      //set score level to start
      div1 = div(test[j],p);
      double scoreLevel = abs(delta(div1.quot,div1.rem));
      
      //list pairs in tied groups
      int x=0;
      while((j+x) < test.size()){
        div2 = div(test[j+x],p);
        double newScore = abs(delta(div2.quot,div2.rem));
        if(abs(scoreLevel - newScore)<0.001){
          double tie = (gamma(div2.quot,div2.rem));
          scoreList.push_back(tie);
          x++;
        }else{
          break;
        }
      }
          
    //make copy of list so when delete bits can still find original
    int counter = scoreList.size();
    scoreListCopy = scoreList;
    while(!scoreList.empty()){
        //find max and then find again in original
        //this is not efficient here
        list<double>::iterator result;
        result = max_element(scoreList.begin(), scoreList.end());
        int choose = distance(scoreList.begin(), result);
        list<double>::iterator result2; 
        result2 = find(scoreListCopy.begin(),scoreListCopy.end(),*result);
        int choose2 = distance(scoreListCopy.begin(), result2); 
        //get pair 
        divresult2 = div(test[j+choose2],p);
        
        
        bool A = false;
        bool B = false;
        
        if(disjoint_mode == true){
          //check is a new disjoint pair
           A = (find(pairList.begin(), pairList.end(), divresult2.quot) != pairList.end());
           B = (find(pairList.begin(), pairList.end(), divresult2.rem) != pairList.end());
        }
        //std::cout<< "A: " << A <<   " B: " << B <<  std::endl;
        if(!A & !B){
          if(i==n_pair){
            //enough pairs already
            break;
          }else{
            output(i,0) = divresult2.quot;
            output(i,1) = divresult2.rem;
            output(i,2) = delta(output(i,0),output(i,1));
            output(i,3) = gamma(output(i,0),output(i,1));
            i++;
            pairList.push_back(divresult2.quot);
            pairList.push_back(divresult2.rem);
          }
       }
        
        scoreList.erase(result);
        
        
      }
    j+=counter;
    } 
    
    
    
    return(output);
}

// [[Rcpp::export]]
vec quickTSP(NumericVector  labels, NumericMatrix obsvData, NumericVector  weights){
    // creates Rcpp matrix from SEXP
    int p = obsvData.nrow(), n = obsvData.ncol();
    mat matData(obsvData.begin(), p, n, false);
    
    //creates vector in same way
    int vec_len = labels.size();
    int w_len = labels.size();
    vec lab(labels.begin(), vec_len , false);
    vec weigh(weights.begin(), w_len , false);
    int nClass = lab.max();
    
    if (vec_len != n | w_len != n) {           
        stop("Label vector must be same length as number of obsverations (columns)");
    }
    
    //split by class
    field<mat> dataByClass(nClass);
    vec obsv_cl = zeros<vec>(nClass);
    for(int i=0; i < nClass; i++){
      obsv_cl(i) = sum(lab==(i+1));
    }
    
    for(int i=0; i < nClass; i++){
      dataByClass(i) = zeros<mat>(p,obsv_cl(i));
    }  
    
    int count0 =0;int count1 =0;
    for(int i=0; i < n; i++){
      
      if(lab(i)==1){
        dataByClass(0).col(count0) = matData.col(i);
        count0++;
        
      }else{
        dataByClass(1).col(count1) =matData.col(i);
        count1++;
      }
    
    
    }
    
    
    
    

    //initialise storage
    field<mat> probByClass(nClass);
    for(int i=0; i < nClass; i++){
       probByClass(i) = zeros<mat>(p,p);
    }
    
        
    
    for(int k=0; k< nClass; k++){
      
      //each predictor
      for(int i=0; i< p; i++){
                
            for(int m=0; m< obsv_cl(k); m++){
            //fill out a col at a time
            //for prob
            umat temp = (dataByClass(k).col(m) > dataByClass(k)(i,m));
            probByClass(k).col(i) =  probByClass(k).col(i) + weigh(m)*arma::conv_to<arma::mat>::from(temp)  ;
          
          }
          
          
      }
      
    }
    //calculate score and tie breaker
    mat delta = zeros<mat>(p,p);
    delta = ((arma::conv_to<arma::mat>::from(probByClass(0))/obsv_cl(0)) - (arma::conv_to<arma::mat>::from(probByClass(1))/obsv_cl(1)));
    
    //cut to triangle to avoid doubles/confusion
    mat L = trimatl(abs(delta));
    uword pair1, pair2;
    L.max(pair1, pair2);
    
    //get top pairs
    vec Lvec = vectorise(L);
    stdvec z = conv_to< stdvec >::from(Lvec);
    //sort and get indexes, sort_index in arma produces weird results!
    vector<int> test = sort_indexes(z);

    
    // prep output
    vec output(3, fill::zeros);
    div_t divresult, divresult2;
    
    //list of already used pairs
    list<int> pairList;
        
    div_t div1;

    div1 = div(test[0],p);
     
    output(0) = div1.quot;
    output(1) = div1.rem;
    output(2) = delta(div1.quot,div1.rem);
 
    return(output);
}



