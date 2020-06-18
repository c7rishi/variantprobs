#include <Rcpp.h>


// modified from R package edgeR to include standard errors of the estimates

using namespace Rcpp;

double sq (double x) {
  return x * x;
}

// [[Rcpp::export]]
Rcpp::List good_turing_multinom(
    NumericVector Obs,
    NumericVector Freq,
    double confid_factor) {

  const int nrows=Obs.size();
  if (nrows!=Freq.size()) {
    stop("lengths of obs and freq vectors must match");
  }

  // Computing constant values.
  double bigN=0;
  double XYs=0, meanX=0, meanY=0, Xsquares=0;
  std::vector<double> log_obs(nrows);
  const int last=nrows-1;

  for (int i=0; i<nrows; ++i) {
    const double& o=Obs[i];
    const double& f=Freq[i];
    bigN+=o*f;

    // Computing log data.
    const double& x=(i==0 ? 0 : Obs[i-1]);
    const double& logO=(log_obs[i]=std::log(double(o)));
    const double logZ=std::log(2*f/double(i==last ? 2*(o-x) : Obs[i+1]-x));
    meanX+=logO;
    meanY+=logZ;
    XYs+=logO*logZ;
    // X2Ys+=logO*logO*logZ;
    Xsquares+=logO*logO;
  }

  meanX/=nrows;
  meanY/=nrows;
  Xsquares-=meanX*meanX*nrows;
  XYs-=meanX*meanY*nrows;
  const double slope=XYs/Xsquares;

  double Xcent2overYs = 0.0;
  for (int i = 0; i < nrows; ++i) {
    Xcent2overYs = sq(log_obs[i]-meanX)/Freq[i];
  }
  const double sd_slope = std::sqrt(Xcent2overYs)/Xsquares;

  // Rcout << "slope= " << slope
  //       << ", slope sdnum= " << sd_slope
  //       <<  std::endl;
  const bool slope_ok = (slope + 1 <= 0);

  // const double intrcpt=meanY-slope*meanX;

  // Computing other various bits and pieces.
  const double& PZero = ((nrows==0 || Obs[0]!=1) ? 0 : Freq[0]/double(bigN));

  // Collecting results.
  double bigNprime=0;
  bool indiffValsSeen=false;
  Rcpp::NumericVector outp(nrows), outpse(nrows);


  for (int i=0; i<nrows; ++i) {

    const double& o=Obs[i];
    const double& f=Freq[i];
    double& op=outp[i];
    double& opse=outpse[i];

    const double next_obs=o+1;
    const double y = next_obs*std::exp(slope*(std::log(double(next_obs))-log_obs[i])); // don't need intercept, cancels out.
    const double yse = o*pow((1+1/o), slope+1)*std::log(1+1/o)*sd_slope;
    if (i==last || Obs[i+1]!=next_obs) {
      indiffValsSeen=true;
    }

    if (!indiffValsSeen) {
      const double& next_n=Freq[i+1];
      const double x = next_obs*next_n/double(f);
      const double xse = x * std::sqrt(1.0/next_n + 1.0/double(f));
      if (std::abs(x - y) <= confid_factor * xse) { // Simplified expression.
        indiffValsSeen=true;
      } else {
        op=x;
        opse=xse;
      }
    }
    if (indiffValsSeen) {
      op=y;
      opse=yse;
    }

    // Rcout << "r=" << o << ", Nr=" << f
    //       << ", indiffValsSeen=" << indiffValsSeen << ", GTprob=" << op << ", SE=" << opse << std::endl;
    bigNprime+=op*f;
  }

  // Running through them to compute the remaining bit.
  const double factor=(1.0-PZero)/bigNprime;
  for (auto& op : outp) {
    op*=factor;
  }
  for (auto& opse : outpse) {
    opse*=factor;
  }

  return Rcpp::List::create(
    Rcpp::NumericVector::create(PZero),
    outp,
    outpse,
    slope_ok
  );
}
