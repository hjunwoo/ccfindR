// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
//#include "asa103.hpp"
#include <RcppEigen.h>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <cstdlib>
#include <iostream>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List vbnmf_update(const Eigen::MatrixXd &X, const Rcpp::List &wh, 
      const Rcpp::List &hyper, const Rcpp::NumericVector &fudge){

    double fud = fudge[0];
    int n = X.rows();
    int m = X.cols();
    Eigen::MatrixXd lw = wh["lw"];
    Eigen::MatrixXd lh = wh["lh"];
    Eigen::MatrixXd ew = wh["ew"];
    Eigen::MatrixXd eh = wh["eh"];
    int r = lw.cols();
    double aw = hyper["aw"];
    double ah = hyper["ah"];
    double bw = hyper["bw"];
    double bh = hyper["bh"];

    Eigen::MatrixXd wth = lw * lh;
    Eigen::MatrixXd xwh = X.array()/wth.array();
    Eigen::MatrixXd sw  = lw.array() * (xwh * lh.transpose()).array();
    Eigen::MatrixXd sh  = lh.array() * (lw.transpose() * xwh).array();

//  std::cout << lw(0,0) << " " << lw(0,1) << "\n";
//  std::cout << lw(1,0) << " " << lw(1,1) << "\n";
//  std::cout << lw(2,0) << " " << lw(2,1) << "\n";

//  std::cout << lh(0,0) << " " << lh(0,1) <<  lh(0,2) << "\n";
//  std::cout << lh(1,0) << " " << lh(1,1) <<  lh(1,2) << "\n";

//  std::cout << aw << " " << bw << "\n";
//  std::cout << eh(0,0) << " " << eh(0,1) <<  eh(0,2) << "\n";
//  std::cout << eh(1,0) << " " << eh(1,1) <<  eh(1,2) << "\n";

    Eigen::MatrixXd alw(n,r);
    alw = Eigen::MatrixXd::Constant(n,r,aw) + sw;
    Eigen::MatrixXd bew(n,r);
    bew = Eigen::MatrixXd::Constant(n,r,aw/bw);
    for(int i=0; i<n; i++)
      bew.row(i) = bew.row(i) + eh.rowwise().sum().transpose();
    ew = alw.array() / bew.array();

    Eigen::MatrixXd alh(r,m);
    alh = Eigen::MatrixXd::Constant(r,m,ah) + sh;
    Eigen::MatrixXd beh(r,m);
    beh = Eigen::MatrixXd::Constant(r,m,ah/bh);
    for(int j=0; j<m; j++)
      beh.col(j) = beh.col(j) + ew.colwise().sum().transpose();
    eh = alh.array() / beh.array();

//  std::cout << bew(0,0) << " " << bew(0,1) << "\n";
//  std::cout << bew(1,0) << " " << bew(1,1) << "\n";
//  std::cout << bew(2,0) << " " << bew(2,1) << "\n";

    for(int i=0;i<n;i++) for(int k=0;k<r;k++){
      double tmp = exp(gsl_sf_psi(alw(i,k)))/bew(i,k);
      lw(i,k) = (tmp > fud ? tmp : fud);
    }
    for(int k=0;k<r;k++) for(int j=0;j<m;j++){
      double tmp = exp(gsl_sf_psi(alh(k,j)))/beh(k,j);
      lh(k,j) = (tmp > fud ? tmp : fud);
    }
//  std::cout << lw(0,0) << " " << lw(0,1) << "\n";
//  std::cout << lw(1,0) << " " << lw(1,1) << "\n";
//  std::cout << lw(2,0) << " " << lw(2,1) << "\n";

//  std::cout << lh(0,0) << " " << lh(0,1) << " " << lh(0,2) << "\n";
//  std::cout << lh(1,0) << " " << lh(1,1) << " " << lh(1,2) << "\n";

    wth = lw * lh;
//  std::cout << wth(0,0) << " " << wth(0,1) << " " << wth(0,2) << "\n";
//  std::cout << wth(1,0) << " " << wth(1,1) << " " << wth(1,2) << "\n";
//  std::cout << wth(2,0) << " " << wth(2,1) << " " << wth(2,2) << "\n";

    Eigen::MatrixXd A = lw.array()*(lw.array().log());
    A = A * lh;
//  std::cout << A(0,0) << " " << A(0,1) << " " << A(0,2) << "\n";
//  std::cout << A(1,0) << " " << A(1,1) << " " << A(1,2) << "\n";
//  std::cout << A(2,0) << " " << A(2,1) << " " << A(2,2) << "\n";
    Eigen::MatrixXd B = lh.array()*(lh.array().log());
    B = lw * B;
//  std::cout << B(0,0) << " " << B(0,1) << " " << B(0,2) << "\n";
//  std::cout << B(1,0) << " " << B(1,1) << " " << B(1,2) << "\n";
//  std::cout << B(2,0) << " " << B(2,1) << " " << B(2,2) << "\n";
    Eigen::MatrixXd lwth = wth.array().log();
    Eigen::MatrixXd U1 = A + B;
    U1 = U1.array() / wth.array();
    U1 = U1 - lwth;
    U1 = X.array() * U1.array();
    U1 = -ew * eh - U1;
//  std::cout << U1(0,0) << " " << U1(0,1) << " " << U1(0,2) << "\n";
//  std::cout << U1(1,0) << " " << U1(1,1) << " " << U1(1,2) << "\n";
//  std::cout << U1(2,0) << " " << U1(2,1) << " " << U1(2,2) << "\n";
    double U=0;
    for(int i=0;i<n;i++) for(int j=0;j<m;j++)
      U += U1(i,j) - gsl_sf_lngamma(X(i,j)+1.0);
    double lga = -gsl_sf_lngamma(aw) + aw*log(aw/bw);

//  std::cout << gsl_sf_lngamma(alw(0,0)) << " " << gsl_sf_lngamma(alw(0,1)) << "\n";
//  std::cout << gsl_sf_lngamma(alw(1,0)) << " " << gsl_sf_lngamma(alw(1,1)) << "\n";
//  std::cout << gsl_sf_lngamma(alw(2,0)) << " " << gsl_sf_lngamma(alw(2,1)) << "\n";

//  std::cout << alw(0,0) << " " << alw(0,1) << "\n";
//  std::cout << alw(1,0) << " " << alw(1,1) << "\n";
//  std::cout << alw(2,0) << " " << alw(2,1) << "\n";

//  std::cout << log(bew(0,0)) << " " << log(bew(0,1)) << "\n";
//  std::cout << log(bew(1,0)) << " " << log(bew(1,1)) << "\n";
//  std::cout << log(bew(2,0)) << " " << log(bew(2,1)) << "\n";

//  std::cout << alh(0,0) << " " << alh(0,1) << " " << alh(0,2) << "\n";
//  std::cout << alh(1,0) << " " << alh(1,1) << " " << alh(0,2) << "\n";
//  double dU1=0;
//  double dU2=0;
    for(int i=0;i<n;i++) for(int k=0;k<r;k++){
      U += -(aw/bw)*ew(i,k) + lga + alw(i,k)*(1.0-log(bew(i,k)))+gsl_sf_lngamma(alw(i,k));
//    dU1 += -(aw/bw)*ew(i,k) + lga;
//    dU2 += alw(i,k)*(1.0-log(bew(i,k)))+gsl_sf_lngamma(alw(i,k));
    }
//  std::cout << dU1 << "\n";
//  std::cout << dU2 << "\n";
    lga = -gsl_sf_lngamma(ah) + ah*log(ah/bh);
    for(int k=0;k<r;k++) for(int j=0;j<m;j++)
      U += -(ah/bh)*eh(k,j) + lga + alh(k,j)*(1.0-log(beh(k,j)))+gsl_sf_lngamma(alh(k,j));
    U/=n*m;

    Rcpp::List z = Rcpp::List::create(Rcpp::Named("w")=ew,
                           Rcpp::Named("h")=eh,
                           Rcpp::Named("lw")=lw,
                           Rcpp::Named("lh")=lh,
                           Rcpp::Named("ew")=ew,
                           Rcpp::Named("eh")=eh,
                           Rcpp::Named("lkh")=U);
    return z;
}
