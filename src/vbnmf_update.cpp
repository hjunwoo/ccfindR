// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

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

    Eigen::MatrixXd alw(n,r);
    alw = Eigen::MatrixXd::Constant(n,r,aw) + sw;
    Eigen::MatrixXd bew(n,r);
    bew = Eigen::MatrixXd::Constant(n,r,aw/bw);
    for(int i=0; i<n; i++)
      bew.row(i) = bew.row(i) + eh.rowwise().sum().transpose();
    ew = alw.array() / bew.array();
    
    Eigen::MatrixXd dw = alw.array() / bew.array() / bew.array();
      
    Eigen::MatrixXd alh(r,m);
    alh = Eigen::MatrixXd::Constant(r,m,ah) + sh;
    Eigen::MatrixXd beh(r,m);
    beh = Eigen::MatrixXd::Constant(r,m,ah/bh);
    for(int j=0; j<m; j++)
      beh.col(j) = beh.col(j) + ew.colwise().sum().transpose();
    eh = alh.array() / beh.array();
    
    Eigen::MatrixXd dh = alh.array() / beh.array() / beh.array();

    for(int i=0;i<n;i++) for(int k=0;k<r;k++){
      double tmp = exp(gsl_sf_psi(alw(i,k)))/bew(i,k);
      lw(i,k) = (tmp > fud ? tmp : fud);
    }
    for(int k=0;k<r;k++) for(int j=0;j<m;j++){
      double tmp = exp(gsl_sf_psi(alh(k,j)))/beh(k,j);
      lh(k,j) = (tmp > fud ? tmp : fud);
    }

    wth = lw * lh;

    Eigen::MatrixXd A = lw.array()*(lw.array().log());
    A = A * lh;
    Eigen::MatrixXd B = lh.array()*(lh.array().log());
    B = lw * B;
    Eigen::MatrixXd lwth = wth.array().log();
    Eigen::MatrixXd U1 = A + B;
    U1 = U1.array() / wth.array();
    U1 = U1 - lwth;
    U1 = X.array() * U1.array();
    U1 = -ew * eh - U1;
    double U=0;
    for(int i=0;i<n;i++) for(int j=0;j<m;j++)
      U += U1(i,j) - gsl_sf_lngamma(X(i,j)+1.0);
    double lga = -gsl_sf_lngamma(aw) + aw*log(aw/bw);

    for(int i=0;i<n;i++) for(int k=0;k<r;k++){
      U += -(aw/bw)*ew(i,k) + lga + alw(i,k)*(1.0-log(bew(i,k)))+gsl_sf_lngamma(alw(i,k));
    }
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
                           Rcpp::Named("lkh")=U,
                           Rcpp::Named("dw")=dw,
                           Rcpp::Named("dh")=dh);
    return z;
}
