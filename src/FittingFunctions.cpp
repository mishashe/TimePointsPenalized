// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

arma::vec Thresholding(arma::vec p, double Threshold) {
  for (int s=0; s<p.size(); s++) {
    if (p(s) > 1.0-Threshold) p(s)=1.0-Threshold; 
    else if (p(s)<Threshold) p(s)=Threshold;
  }
  return p;
}

arma::vec sgn(arma::vec val) {
  arma::vec out(val.size());
  for (int i=0; i<val.size(); i++)
  {
    if (val(i)>0) {
      out(i)=1;
    }
    else if (val(i)<0) {
      out(i)=-1;
    } else {
      out(i)=0;
    }
  }
  return out;
}

// Calculates intercept in iterative way 
double CalculateDeltaIntercept(arma::vec y, arma::vec p, arma::vec w){
  int ns = y.size();
  double DeltaBeta0=1.0;
  double DeltaBeta0Prev=50.0;
  while(fabs(DeltaBeta0-DeltaBeta0Prev)/(1e-5+fabs(DeltaBeta0)+fabs(DeltaBeta0Prev))>1.0e-5) {
    double Sum1=0;
    double Sum2=0;
    DeltaBeta0Prev=DeltaBeta0;
    for (int s=0; s<ns; s++) {
      Sum1+=w(s)*y(s);
      Sum2+=w(s)/(1.0/p(s)-1.0+DeltaBeta0);
    }
    DeltaBeta0 = Sum1/Sum2;
  }
  return log(DeltaBeta0);
}

// Calculates first (b) and second (a) derivatives of the function to minimize 
void GetHessian(arma::mat x, arma::vec beta, arma::vec p, arma::vec y, double lam2, arma::vec w,
                arma::vec& b, arma::mat& a) {
  int ng = beta.size();
  int ns = y.size();
  a.fill(0.0);
  b.fill(0.0);
  
  // Calc with lambda1=0
  for (int g1=0; g1<ng; g1++)
  {
    if (g1>0) {b(g1) += 2.0*lam2*(beta(g1)-beta(g1-1)); a(g1,g1) += 2.0*lam2;}
    if (g1<(ng-1)) {b(g1) += 2.0*lam2*(beta(g1)-beta(g1+1)); a(g1,g1) += 2.0*lam2;}
    for (int s=0; s<ns; s++)
    {
      b(g1) += -x(s,g1)*(y(s)*(1.0-p(s))-(1.0-y(s))*p(s))*w(s);
      a(g1,g1) += x(s,g1)*x(s,g1)*p(s)*(1.0-p(s))*w(s);    
    }
    
    for (int g2=g1+1; g2<ng; g2++)
    {
      if (g2==g1-1) a(g1,g2) += -2.0*lam2;
      if (g2==g1+1) a(g1,g2) += -2.0*lam2;
      // for (int s=0; s<ns; s++)
      // {
      //   a(g1,g2) += x(s,g1)*x(s,g2)*p(s)*(1.0-p(s))*w(s);  
      // }
    }
    for (int g2=0; g2<g1; g2++) a(g1,g2) = a(g2,g1);
  }
  return;
}


void UpdateIntercept(arma::mat x0, arma::vec y, arma::vec tV, double lam1, double lam2, 
                     arma::vec beta, arma::vec& Intercept, arma::vec w, arma::vec IndFor0,
                     arma::vec IndTFor0, arma::vec& M, double& LLmin){
  arma::vec InterceptPrev = Intercept;
  arma::vec Mprev = M;
  double LLprev= LLmin;
  
  int m = beta.size();
  int nt = tV.size();
  int ns = y.size();
  int m0 = m/nt;
  arma::vec p(ns); 
  double LL;
  
  p = 1.0/(1.0+exp(-M));
  p = Thresholding(p, 1.0e-2);
  double F=0;
  for (int s=0;s<ns;s++) {
    F += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  }
  // Update group of intercepts
  for (int it=0;it<nt;it++){
    std::vector<double> yT,pT,wT;
    
    for (int s=0;s<ns;s++) {
      if (IndTFor0(s)==it) {
        yT.push_back(y(IndFor0(s)));
        pT.push_back(p(IndFor0(s)));
        wT.push_back(w(IndFor0(s)));
      }
    }
    double dIntercept = CalculateDeltaIntercept(yT, pT, wT);
    Intercept(it) = Intercept(it) + dIntercept;
    for (int s=0;s<ns;s++) {
      if (IndTFor0(s)==it) {
        M(s) = M(s) + dIntercept;
      }
    }
  }
  p = 1.0/(1.0+exp(-M));
  p = Thresholding(p, 1.0e-2);
  for (int s=0;s<ns;s++) {
    LLmin += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  }
  LLmin += - F;
  if (LLmin>LLprev) {
    Intercept = InterceptPrev;
    M = Mprev;
    LLmin = LLprev;
  }
  return;
}



void SingleGeneRound(arma::mat x0, arma::vec y, arma::vec tV, double lam1, double lam2,
                          arma::vec& beta, arma::vec& Intercept, arma::vec w, arma::vec IndFor0,
                          arma::vec IndTFor0, arma::vec& M, double& LLmin){
  int m = beta.size();
  int nt = tV.size();
  int ns = y.size();
  int m0 = m/nt;
  arma::vec b(nt);
  arma::mat a(nt,nt);
  arma::vec p(ns); 
  arma::uvec dGg(nt);
  arma::mat x0G(ns,nt);
  double LL;
  UpdateIntercept(x0, y, tV, lam1, lam2, beta, Intercept, w,  IndFor0, IndTFor0, M, LLmin);
  p = 1.0/(1.0+exp(-M));
  p = Thresholding(p, 1.0e-2);
  for (int g=0;g<m0;g++) {
    for (int it=0;it<nt;it++) {dGg(it) = it*m0 + g;}
    x0G.fill(0);
    for (int s=0;s<ns;s++) {x0G(s,IndTFor0(s)) = x0(IndFor0(s),g);}
    GetHessian(x0G, beta(dGg), p, y, lam2,w,b,a);

    for (int it=0;it<nt;it++) {
      double betaOld = beta(g+m0*it);
      double betaTry = (b(it) - a(it,it) * betaOld)/a(it,it);
      int sigma0 = (betaTry > 0 ? 1 : (betaTry < 0 ? -1 : 0 ) );
      double betaNew = (b(it) - a(it,it) * betaOld + sigma0*lam1)/a(it,it);
      int sigma = (betaNew > 0 ? 1 : (betaNew < 0 ? -1 : 0 ) );
      if (sigma0!=sigma) {
        betaNew=0;
      }
      if (betaNew!=betaOld) {
        arma::vec Mnew = M;
        for (int s=0;s<ns;s++) Mnew(s) += x0G(s,it) * (betaNew-betaOld);
        arma::vec pnew = 1.0/(1.0+exp(-Mnew));
        pnew = Thresholding(pnew, 1.0e-2);
        beta(g+it*m0) = betaNew;
        LL=0;
        for (int s=0;s<ns;s++) {
          LL += -y(s)*log(pnew(s))*w(s) - (1-y(s))*log(1.0-pnew(s))*w(s);
        }
        LL += lam1*(fabs(betaNew)-fabs(betaOld));
        if (it>0) {
          LL+= lam2*(betaNew-beta(g+(it-1)*m0))*(betaNew-beta(g+(it-1)*m0))
               -lam2*(betaOld-beta(g+(it-1)*m0))*(betaOld-beta(g+(it-1)*m0));
        }
        if (it<(nt-1)) {
          LL+= lam2*(betaNew-beta(g+(it+1)*m0))*(betaNew-beta(g+(it+1)*m0))
               -lam2*(betaOld-beta(g+(it+1)*m0))*(betaOld-beta(g+(it+1)*m0));
        }
        if (LL<LLmin) {
          LLmin = LL;
          p = pnew;
          M = Mnew;
        }
        else {
          beta(g+it*m0) = betaOld;
        }
      }
    }
  }
  return;
}

// Make one step for the t-group of a gene (soft threshold)
arma::vec glmnetSimple(arma::mat X, arma::vec Y, double lam1){
  int ng = Y.size();
  arma::vec beta0(ng);
  beta0 = -solve( X, Y);
  arma::vec sign0(ng);
  sign0 = sgn(beta0);
  arma::vec beta(ng);
  beta = -solve( X, Y + sign0*lam1);
  arma::vec sign(ng);
  sign = sgn(beta);
  for (int g=0;g<ng;g++) {
    if (sign0(g)!=sign(g)) {
      return(beta*0);
    }
  }
  return beta;
}


void GroupRound(arma::mat x0, arma::vec y, arma::vec tV, double lam1, double lam2, 
                     arma::vec& beta, arma::vec& Intercept, arma::vec w, arma::vec IndFor0,
                     arma::vec IndTFor0, arma::vec& M, double& LLmin){
  int m = beta.size();
  int nt = tV.size();
  int ns = y.size();
  int m0 = m/nt;
  arma::uvec dGg(nt);
  arma::mat x0G(ns,nt);
  arma::vec b(nt);
  arma::mat a(nt,nt);
  arma::vec p(ns); 
  arma::vec betaOld(nt);
  arma::vec betaNew(nt);
  double LL;
  UpdateIntercept(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0, IndTFor0, M, LLmin);
  p = 1.0/(1.0+exp(-M));
  p = Thresholding(p, 1.0e-2);
  for (int g=0;g<m0;g++) {
    for (int it=0;it<nt;it++) {
      dGg(it) = it*m0 + g;
    }
    betaOld = beta(dGg);
    x0G.fill(0);
    
    for (int s=0;s<ns;s++) {
      x0G(s,IndTFor0(s)) = x0(IndFor0(s),g);
    }
    GetHessian(x0G, betaOld, p, y, lam2,w,b,a);
    //Rcout<<a(0,0)<<" "<<a(0,1)<<" "<<a(1,0)<<" "<<a(1,1)<<" "<<a(1,2)<<" "<<a(2,1)<<" "<<std::endl;
    //Rcout<<b(0)<<" "<<b(1)<<" "<<b(2)<<" "<<b(3)<<" "<<std::endl;
    
    betaNew = glmnetSimple(a,b - a * betaOld,lam1);
    //Rcout<<betaNew(1)<<std::endl;
    // return;
    if (any(betaNew!=betaOld)) {
      arma::vec Mnew = M + x0G * (betaNew-betaOld);
      arma::vec pnew = 1.0/(1.0+exp(-Mnew));
      pnew = Thresholding(pnew, 1.0e-2);
      beta(dGg) = betaNew;
      LL=0;
      for (int s=0;s<ns;s++) {
        LL += -y(s)*log(pnew(s))*w(s) - (1-y(s))*log(1.0-pnew(s))*w(s);
      }
      for (int it=0;it<nt;it++) {
        LL += lam1*fabs(betaNew(it))-lam1*fabs(betaOld(it));
      }
      for (int it=1;it<nt;it++) {
        LL += lam2*(betaNew(it)-betaNew(it-1))*(betaNew(it)-betaNew(it-1))-lam2*(betaOld(it)
              -betaOld(it-1))*(betaOld(it)-betaOld(it-1));
      }
      for (int it=0;it<(nt-1);it++) {
        LL += lam2*(betaNew(it)-betaNew(it+1))*(betaNew(it)-betaNew(it+1))-lam2*(betaOld(it)
              -betaOld(it+1))*(betaOld(it)-betaOld(it+1));
      }
      if (LL<LLmin) {
        LLmin = LL;
        p = pnew;
        M = Mnew;
      }
      else {
        beta(dGg) = betaOld;
      }
    }
  }
  return;
}

// [[Rcpp::export]]
List Fit(arma::mat x0, arma::vec y, arma::vec tV, double lam1, double lam2,
                   arma::vec beta, arma::vec Intercept, arma::vec w, arma::vec IndFor0,
                   arma::vec IndTFor0){
  IndFor0 = IndFor0-1; 
  IndTFor0 = IndTFor0-1;
  int nt = tV.size();
  w=w/accu(w);
  int m = beta.size();
  int ns = y.size();
  int m0 = m/nt;
  arma::vec M(ns);
  for (int s=0;s<ns;s++) {
    for (int it=0;it<nt;it++) {
      if (IndTFor0(s)==it) {
        M(s) = M(s) + Intercept(it);
      }
    }
    for (int g=0;g<m0;g++) {
        M(s) = M(s) + x0(IndFor0(s),g) * beta(g+IndTFor0(s)*m0);
    }
  }
  arma::vec p(ns);
  p = 1.0/(1.0+exp(-M));
  p = Thresholding(p, 1.0e-2);

  double LL=0.0;   

  for (int s=0;s<ns;s++) {
    LL += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  }
  for (int g=1;g<m0;g++) {
    for (int it=0;it<nt;it++) {
      LL += lam1*fabs(beta(g+it*m0));
    }
  }
  for (int g=1;g<m0;g++) {
    for (int it=0;it<nt-1;it++) {
      LL+= lam2*(beta(g+(it+1)*m0)-beta(g+it*m0))*(beta(g+(it+1)*m0)-beta(g+it*m0));
    }
  }
  double LLprev; LLprev = -50.1654465;
  arma::vec betaPrev = -(beta+0.001);
  //Rcout<<abs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)<<std::endl;
  bool dontstop = TRUE;
  do {
    Rcout<<" LL s "<<LL<<" "<<LLprev<<" "<<fabs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)<<" "<<any(sgn(beta) != sgn(betaPrev))<<std::endl;
    Rcout<<" LL s 2 "<<LL<<" "<<LLprev<<" "<<fabs(LL-LLprev) <<" "<<sqrt(LLprev*LLprev+LL*LL)<<" "<<any(sgn(beta) != sgn(betaPrev))<<std::endl;
    Rcout<<(fabs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)>1.0e-5 | any(sgn(beta) != sgn(betaPrev)))<<std::endl;
    LLprev = LL+0.00000001;
    betaPrev = beta;
    //SingleGeneRound(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0,IndTFor0, M, LL);
    //GroupRound(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0,IndTFor0, M, LL);
    Rcout<<" LL "<<LL<<" "<<LLprev<<" "<<fabs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)<<" "<<any(sgn(beta) != sgn(betaPrev))<<std::endl;
    Rcout<<" LL 2 "<<LL<<" "<<LLprev<<" "<<fabs(LL-LLprev) <<" "<<sqrt(LLprev*LLprev+LL*LL)<<" "<<any(sgn(beta) != sgn(betaPrev))<<std::endl;
    Rcout<<(fabs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)>1.0e-5 | any(sgn(beta) != sgn(betaPrev)))<<std::endl;
    dontstop = (fabs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)>1.0e-5 | any(sgn(beta) != sgn(betaPrev)));
    Rcout<<dontstop<<std::endl;
  }
  while (dontstop)
  return(List::create(Named("beta") = beta, Named("Intercept") = Intercept, Named("nG") = accu(beta!=0)));
}





