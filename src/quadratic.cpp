#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <iostream>
using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

NumericVector concatenate(const NumericVector &x, const NumericVector &y) {
  NumericVector result = x;
  int n = y.size();
  for (int i = 0; i < n; i++) {
    result.push_back(y[i]);
  }
  return result;
}
  
IntegerVector which(const LogicalVector &a) {
  IntegerVector result;
  int n = a.size();
  for (int i = 0; i < n; i++) {
    if (a[i]) result.push_back(i);
  }
  return result;
}

IntegerVector order(const NumericVector &x) {
  const int n = x.size();
  std::vector<std::pair<double, int> > pairs;

  for(size_t i = 0; i < n; i++) {
    pairs.push_back(std::pair<double, int>(x[i], i));
  }
  
  std::sort(pairs.begin(), pairs.end());
  
  Rcpp::IntegerVector result;
  for(size_t i = 0; i < n; i++) {
    result.push_back(pairs[i].second);
  }
  return result;
}

// [[Rcpp::export]]
List quadratic_1D_mu(NumericVector mu, NumericVector lambda = NumericVector::create(), NumericVector weights = 1) {
  int m = mu.length();
//  cout << m << endl;
  
  if (lambda.length() == 0) {
    lambda.push_back(R_PosInf);
  }
  
  if (lambda.length() == 1) {
    lambda = rep(lambda[0], m);
  }
  
  if (weights.length() == 1) {
    weights = rep(weights[0], m);
  }
  
  double w_sum = 0, w_mu_sum = 0, w_mu2_sum = 0, w_lambda_sum = 0;
  vector<pair<double, int> > end_pts;
  end_pts.reserve(m * 2);
  bool included[m];
  for (int i = 0; i < m; i++) {
    if (weights[i] == 0) continue;
    if (lambda[i] == R_PosInf) { //untruncated
      w_sum += weights[i];
      w_mu_sum += weights[i] * mu[i];
      w_mu2_sum += weights[i] * mu[i] * mu[i];
      included[i] = true;
    } else {
      included[i]  = false;
      w_lambda_sum += weights[i] * lambda[i];
      if (lambda[i] > 0) { //truncated
        end_pts.push_back(make_pair(mu[i] - sqrt(lambda[i]), i));
        end_pts.push_back(make_pair(mu[i] + sqrt(lambda[i]), i));
      }
    } 
  }
  sort(end_pts.begin(), end_pts.end());
  double par = (w_sum > 0) ? w_mu_sum / w_sum : 0;
  double value = w_mu2_sum - w_sum * par * par + w_lambda_sum;
  for (int i = 0; i < (int)end_pts.size(); i++) {
    int p = end_pts[i].second;
    included[p] = !included[p];
    double sign = included[p] * 2 - 1;
    w_sum += sign * weights[p];
    w_mu_sum += sign * weights[p] * mu[p];
    w_mu2_sum += sign * weights[p] * mu[p] * mu[p];
    w_lambda_sum -= sign * weights[p] * lambda[p];
    if (w_sum > 0) {
      double x = w_mu_sum / w_sum;
      double f = w_mu2_sum - w_sum * x * x + w_lambda_sum;
      if (f < value) {
        value = f;
        par = x;
      }
    }
  }
  return List::create(Named("par") = par, Named("value") = value);
}

// [[Rcpp::export]]
List quadratic_1D_coef(NumericMatrix coef, NumericVector lambda = NumericVector::create()) {
  const NumericVector& a = coef.column(0);
  NumericVector mu = -coef.column(1) / (2*a);
  NumericVector offset = coef.column(2) - mu * mu * a;

  if (lambda.length() == 0) {
    lambda.push_back(R_PosInf);
  }
  
  if (lambda.length() == 1) {
    lambda = rep(lambda[0], mu.size());
  }
  
  lambda = (lambda - offset) / a;

  List result = quadratic_1D_mu(mu, lambda, a);
  
  return List::create(Named("par") = result["par"], Named("value") = (double)result["value"] + sum(offset));
}

// [[Rcpp::export]]
List quadratic_cpp(SEXP A, NumericVector b = NumericVector::create(), NumericVector x0 = NumericVector::create(), NumericVector lambda = NumericVector::create(), int max_it = 1000, double tol = 1e-3, bool trace = false) {
    bool sparse = !Rf_isMatrix(A);
    int n, p;
    if (sparse) {
        S4 mat(A);
        IntegerVector dim = mat.slot("Dim");
        n = dim[0];
        p = dim[1];
    } else {
        NumericMatrix mat(A);
        n = mat.nrow();
        p = mat.ncol();
    }
    
    if (b.length() == 0) {
        b.push_back(0);
    }
    if (b.length() == 1) {
        b = rep(b[0], n);
    }
    if (x0.length() == 0) {
        x0.push_back(0);
    }
    if (x0.length() == 1) {
        x0 = rep(x0[0], p);
    }
    if (lambda.length() == 0) {
        lambda.push_back(R_PosInf);
    }
    if (lambda.length() == 1) {
        lambda = rep(lambda[0], n);
    }
    
    NumericVector x = clone(x0);
    NumericVector Axb = sparse ? wrap(as<MappedSparseMatrix<double> >(A) * as<Map<VectorXd> >(x) + as<Map<VectorXd> >(b)) : wrap(as<Map<MatrixXd> >(A) * as<Map<VectorXd> >(x) + as<Map<VectorXd> >(b));
    double obj = R_PosInf;
    
    for (int it = 1; it <= max_it; it++) {
        NumericVector x_old = Rcpp::clone(x);
        double obj_old = obj;
        
        for (int i = 0; i < p; i++) {
            if (sparse) {
                S4 mat(A);
                IntegerVector index = mat.slot("i");
                IntegerVector pointer = mat.slot("p");
                NumericVector coef = mat.slot("x");
                int start = pointer[i];
                int end = pointer[i+1];
                vector<double> mu, temp_Axb, lam, weights;
                mu.resize(end - start);
                temp_Axb.resize(end - start);
                lam.resize(end - start);
                weights.resize(end - start);
                for (int j = start; j < end; j++) {
                    int row = index[j];
                    temp_Axb[j - start] = Axb[row] - coef[j] * x[i];
                    mu[j - start] = -temp_Axb[j - start] / coef[j];
                    lam[j - start] = lambda[row] / coef[j] / coef[j];
                    weights[j - start] = coef[j] * coef[j];
                }
                List sol = quadratic_1D_mu(wrap(mu), wrap(lam), wrap(weights));
                x[i] = sol[0];
                for (int j = start; j < end; j++) {
                    int row = index[j];
                    Axb[row] = temp_Axb[j - start] + coef[j] * x[i];
                }
            } else {
                NumericMatrix mat(A);
                NumericVector coef = mat.column(i);
                NumericVector temp_Axb = Axb - coef * x[i];
                List sol = quadratic_1D_mu(-temp_Axb / coef, lambda / coef / coef, coef * coef);
                x[i] = sol[0];
                Axb = temp_Axb + coef * x[i];
            }
        }
        
        double dxnorm = Rcpp::max(Rcpp::abs(x - x_old));
        
        if (trace) {
            obj = Rcpp::sum(Rcpp::pmin(Axb * Axb, lambda));
            if (it == 1) cout << "iteration\tdelta_x\tobj_value\tdelta_obj_value\n";
            cout << it << "\t" << dxnorm << "\t" << obj << "\t" << obj - obj_old << endl;
            if (obj > obj_old) cout << "WARNING: objective values increased!\n";
        }
        
        if (dxnorm < tol) break;
    }
    
    obj = Rcpp::sum(Rcpp::pmin(Axb * Axb, lambda));
    return List::create(Named("par") = x, Named("value") = obj);
}
