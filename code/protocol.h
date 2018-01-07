#ifndef __PROTOCOL_H_
#define __PROTOCOL_H_

#include "gwasiter.h"
#include "mpc.h"
#include "util.h"
#include <vector>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <chrono>

using namespace NTL;
using namespace std;

using msec = chrono::milliseconds;
using get_time = chrono::steady_clock;

#define ABS(a) (((a)<0)?-(a):(a))

auto clock_start = get_time::now();

void tic() {
  clock_start = get_time::now();
}

int toc() {
  auto clock_end = get_time::now();
  int duration = chrono::duration_cast<msec>(clock_end - clock_start).count();
  cout << "Elapsed time is " << duration / 1000.0 << " secs" << endl;
  return duration;
}

bool DecrComp(const pair<int, double> &a, const pair<int, double> &b) {
  return a.second > b.second;
}

string cache(int pid, string desc) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX << "_" << desc << ".bin";
  return oss.str();
}

string cache(int pid, int index) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX << "_" << index << ".bin";
  return oss.str();
}

string outname(string desc) {
  ostringstream oss;
  oss << Param::OUTPUT_FILE_PREFIX << "_" << desc << ".txt";
  return oss.str();
}

bool logireg_protocol(MPCEnv& mpc, int pid) {
  SetNumThreads(Param::NUM_THREADS);
  cout << AvailableThreads() << " threads created" << endl;

  int ntop = 100;

  int n0 = Param::NUM_INDS;
  int m0 = Param::NUM_SNPS;
  int k = Param::NUM_DIM_TO_REMOVE;

  cout << "n0: " << n0 << ", " << "m0: " << m0 << endl;

  // Shared variables
  string s;
  fstream fs;
  ofstream ofs;
  ifstream ifs;
  streampos strpos;
  int ind;
  Vec<ZZ_p> tmp_vec;
  Mat<ZZ_p> tmp_mat;

  //mpc.ProfilerPushState("main");

  ZZ_p fp_one = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);

  Vec<ZZ_p> pheno;
  Init(pheno, n0);

  Mat<ZZ_p> cov;
  Init(cov, n0, Param::NUM_COVS);

  if (!exists(cache(pid, "input_geno")) || !exists(cache(pid, "input_pheno_cov"))) {
    cout << "Initial data sharing results not found:" << endl;
    cout << "\t" << cache(pid, "input_geno") << endl;
    cout << "\t" << cache(pid, "input_pheno_cov") << endl;
    return false;
  }

  cout << "Initial data sharing results found" << endl;

  ifs.open(cache(pid, "input_pheno_cov").c_str(), ios::binary);
  mpc.ReadFromFile(pheno, ifs, n0);
  mpc.ReadFromFile(cov, ifs, n0, Param::NUM_COVS);
  ifs.close();

  cout << "Phenotypes and covariates loaded" << endl;

  Vec<ZZ_p> gkeep1;
  Init(gkeep1, m0);
  
  cout << "Using locus missing rate filter from a previous run" << endl;
  
  if (pid == 2) {
    ifs.open(outname("gkeep1").c_str());
    for (int i = 0; i < m0; i++) {
      ifs >> gkeep1[i];
    }
    ifs.close();

    mpc.SendVec(gkeep1, 0);
    mpc.SendVec(gkeep1, 1);
  } else {
    mpc.ReceiveVec(gkeep1, 2, m0);
  }

  uint m1 = conv<uint>(Sum(gkeep1));
  cout << "n0: " << n0 << ", " << "m1: " << m1 << endl;

  Vec<ZZ_p> ikeep;
  Init(ikeep, n0);

  cout << "Using individual missing rate/het rate filters from a previous run" << endl;
  
  if (pid == 2) {
    ifs.open(outname("ikeep").c_str());
    for (int i = 0; i < n0; i++) {
      ifs >> ikeep[i];
    }
    ifs.close();

    mpc.SendVec(ikeep, 0);
    mpc.SendVec(ikeep, 1);
  } else {
    mpc.ReceiveVec(ikeep, 2, n0);
  }

  uint n1 = conv<uint>(Sum(ikeep));
  cout << "n1: " << n1 << ", " << "m1: " << m1 << endl;

  cout << "Filtering phenotypes and covariates" << endl;
  mpc.Filter(pheno, ikeep, n1);
  mpc.FilterRows(cov, ikeep, n1);

  Vec<ZZ_p> gkeep2;
  Init(gkeep2, m1);

  cout << "Using MAF/HWE filters from a previous run" << endl;
  
  if (pid == 2) {
    ifs.open(outname("gkeep2").c_str());
    for (int i = 0; i < m1; i++) {
      ifs >> gkeep2[i];
    }
    ifs.close();

    mpc.SendVec(gkeep2, 0);
    mpc.SendVec(gkeep2, 1);
  } else {
    mpc.ReceiveVec(gkeep2, 2, m1);
  }

  uint m2 = conv<uint>(Sum(gkeep2));
  cout << "n1: " << n1 << ", " << "m2: " << m2 << endl;

  cout << "Using CA statistics from a previous run" << endl;
  
  Vec<ZZ_p> gkeep3;
  Init(gkeep3, m2);

  if (pid == 2) {
    vector<pair<int, double> > cavec(m2);

    ifs.open(outname("assoc").c_str());
    double val;
    for (int i = 0; i < m2; i++) {
      ifs >> val;
      cavec[i] = make_pair(i, val * val);
    }
    ifs.close();

    sort(cavec.begin(), cavec.end(), DecrComp);

    cout << "Selected top " << ntop << " candidates" << endl;
    cout << "Top 5 CA stats: " << cavec[0].second;
    for (int i = 1; i < 5; i++) {
      cout << ", " << cavec[i].second;
    }
    cout << endl;

    for (int i = 0; i < ntop; i++) {
      gkeep3[cavec[i].first] = 1;
    }
    mpc.SendVec(gkeep3, 0);
    mpc.SendVec(gkeep3, 1);
  } else {
    mpc.ReceiveVec(gkeep3, 2, m2);
  }

  Mat<ZZ_p> V;
  Init(V, k, n1);

  cout << "Using eigenvectors from a previous run" << endl;
  ifs.open(cache(pid, "eigen").c_str(), ios::binary);
  mpc.ReadFromFile(V, ifs, k, n1);
  ifs.close();

  // Concatenate covariate matrix and jointly orthogonalize
  mpc.Transpose(cov);
  V.SetDims(k + Param::NUM_COVS, n1);
  if (pid > 0) {
    for (int i = 0; i < Param::NUM_COVS; i++) {
      V[k + i] = cov[i] * fp_one;
    }
  }
  cov.kill();

  mpc.OrthonormalBasis(V, V);

  Vec<ZZ_p> V_mean;
  Init(V_mean, V.NumRows());
  ZZ_p fp_denom = DoubleToFP(1 / ((double) V.NumCols()), Param::NBIT_K, Param::NBIT_F);
  for (int i = 0; i < V_mean.length(); i++) {
    V_mean[i] = Sum(V[i]) * fp_denom;
  }
  mpc.Trunc(V_mean);
  for (int i = 0; i < V_mean.length(); i++) {
    AddScalar(V[i], -V_mean[i]);
  }

  Vec<ZZ_p> V_var;
  mpc.InnerProd(V_var, V);
  mpc.Trunc(V_var);
  V_var *= fp_denom;
  mpc.Trunc(V_var);

  Vec<ZZ_p> V_stdinv, dummy_vec;
  mpc.FPSqrt(dummy_vec, V_stdinv, V_var);
  
  for (int i = 0; i < V_mean.length(); i++) {
    mpc.MultMat(V[i], V[i], V_stdinv[i]);
  }
  mpc.Trunc(V);

  Vec<bool> gkeep;
  gkeep.SetLength(m0);
  for (int j = 0; j < m0; j++) {
    gkeep[j] = gkeep1[j] == 1;
  }

  ind = 0;
  for (int j = 0; j < m0; j++) {
    if (gkeep[j]) {
      gkeep[j] = gkeep2[ind] == 1;
      ind++;
    }
  }

  ind = 0;
  for (int j = 0; j < m0; j++) {
    if (gkeep[j]) {
      gkeep[j] = gkeep3[ind] == 1;
      ind++;
    }
  }

  Mat<ZZ_p> X, X_mask;
  if (exists(cache(pid, "logi_input"))) {
    cout << "logi_input cache found" << endl;
    ifs.open(cache(pid, "logi_input").c_str(), ios::in | ios::binary);
    mpc.BeaverReadFromFile(X, X_mask, ifs, ntop, n1);
    ifs.close();
  } else {
    X.SetDims(n1, ntop);
    X_mask.SetDims(n1, ntop);

    ifs.open(cache(pid, "input_geno").c_str(), ios::binary);
    if (pid > 0) {
      mpc.ImportSeed(10, ifs);
    } else {
      for (int p = 1; p <= 2; p++) {
        mpc.ImportSeed(10 + p, ifs);
      }
    }

    cout << "Collecting genotypes for top candidates" << endl;
    ind = -1;
    int batch_size = n1 / 10;
    tic();
    for (int cur = 0; cur < n1; cur++) {
      ind++;

      if ((cur + 1) % batch_size == 0 || cur == n1 - 1) {
        cout << cur + 1 << "/" << n1 << ", "; toc();
        tic();
      }

      Mat<ZZ_p> g0, g0_mask;
      Vec<ZZ_p> miss0, miss0_mask;

      while (ikeep[ind] != 1) {
        if (pid > 0) {
          mpc.SkipData(ifs, 3, m0); // g
          mpc.SkipData(ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          }
        }
        ind++;
      }

      if (pid > 0) {
        mpc.ReadFromFile(g0, ifs, 3, m0); // g
        mpc.ReadFromFile(miss0, ifs, m0); // miss

        mpc.SwitchSeed(10);
        mpc.RandMat(g0_mask, 3, m0);
        mpc.RandVec(miss0_mask, m0);
        mpc.RestoreSeed();
      } else {
        Init(g0, 3, m0);
        Init(g0_mask, 3, m0);
        Init(miss0, m0);
        Init(miss0_mask, m0);

        for (int p = 1; p <= 2; p++) {
          mpc.SwitchSeed(10 + p);
          mpc.RandMat(tmp_mat, 3, m0);
          mpc.RandVec(tmp_vec, m0);
          mpc.RestoreSeed();

          g0_mask += tmp_mat;
          miss0_mask += tmp_vec;
        }
      }
      
      Mat<ZZ_p> g, g_mask;
      Vec<ZZ_p> miss, miss_mask;
      g.SetDims(3, ntop);
      miss.SetLength(ntop);
      g_mask.SetDims(3, ntop);
      miss_mask.SetLength(ntop);
      int ind2 = 0;
      for (int j = 0; j < m0; j++) {
        if (gkeep[j]) {
          for (int k = 0; k < 3; k++) {
            g[k][ind2] = g0[k][j];
            g_mask[k][ind2] = g0_mask[k][j];
          }
          miss[ind2] = miss0[j];
          miss_mask[ind2] = miss0_mask[j];
          ind2++;
        }
      }

      X[cur] = g[1] + 2 * g[2];
      X_mask[cur] = g_mask[1] + 2 * g_mask[2];
    }

    mpc.Transpose(X); // ntop-by-n1
    transpose(X_mask, X_mask);

    fs.open(cache(pid, "logi_input").c_str(), ios::out | ios::binary);
    mpc.BeaverWriteToFile(X, X_mask, fs);
    fs.close();
  }

  // Shuffle
  Mat<ZZ_p> xrt, xmt;
  transpose(xrt, X);
  transpose(xmt, X_mask);

  transpose(V, V);

  Vec<long> indices;
  indices.SetLength(n1);
  for (int i = 0; i < n1; i++) {
    indices[i] = i + 1;
  }

  mpc.SwitchSeed(-1);
  for (int i = 0; i < n1-1; i++) {
    long chosen = RandomBnd(n1-i);
    Vec<ZZ_p> tmp;
    if (chosen > 0) {
      tmp = xrt[i];
      xrt[i] = xrt[i + chosen];
      xrt[i + chosen] = tmp;

      tmp = xmt[i];
      xmt[i] = xmt[i + chosen];
      xmt[i + chosen] = tmp;

      tmp = V[i];
      V[i] = V[i + chosen];
      V[i + chosen] = tmp;

      ZZ_p tmpy;
      tmpy = pheno[i];
      pheno[i] = pheno[i + chosen];
      pheno[i + chosen] = tmpy;

      long t = indices[i];
      indices[i] = indices[i + chosen];
      indices[i + chosen] = t;
    }
  }
  mpc.RestoreSeed();

  transpose(X, xrt);
  transpose(X_mask, xmt);
  transpose(V, V);
  xrt.kill();
  xmt.kill();

  Mat<ZZ_p> V_mask;
  mpc.BeaverPartition(V_mask, V);

  Vec<ZZ_p> pheno_mask;
  mpc.BeaverPartition(pheno_mask, pheno);

  Vec<ZZ_p> b0;
  Mat<ZZ_p> bv;
  Vec<ZZ_p> bx;
  mpc.ParallelLogisticRegression(b0, bv, bx, X, X_mask, V, V_mask, pheno, pheno_mask, 500);

  fs.open(cache(pid, "logireg_final_coeff").c_str(), ios::out | ios::binary);
  if (pid > 0) {
    mpc.WriteToFile(b0, fs);
    mpc.WriteToFile(bv, fs);
    mpc.WriteToFile(bx, fs);
  }
  fs.close();

  mpc.RevealSym(bx);
  if (pid == 2) {
    Vec<double> bx_double;
    FPToDouble(bx_double, bx, Param::NBIT_K, Param::NBIT_F);
    ofs.open(outname("logi_coeff").c_str(), ios::out);
    for (int i = 0; i < bx_double.length(); i++) {
      ofs << bx_double[i] << endl;
    }
    ofs.close();
    cout << "Result written to " << outname("logi_coeff") << endl;
  }

  return true;
}

bool data_sharing_protocol(MPCEnv& mpc, int pid) {
  int n = Param::NUM_INDS;

  fstream fs;

  Vec<ZZ_p> pheno;
  Init(pheno, n);

  Mat<ZZ_p> cov;
  Init(cov, n, Param::NUM_COVS);

  fs.open(cache(pid, "input_geno").c_str(), ios::out | ios::binary);
  if (pid > 0) {
    mpc.ExportSeed(fs, 0);
  } else {
    for (int p = 1; p <= 2; p++) {
      mpc.ExportSeed(fs, p);
    }
  }

  GwasIterator git(mpc, pid);

  git.Init(true, true);

  long bsize = n / 10;

  cout << "Begin processing:" << endl;

  tic();
  for (int i = 0; i < n; i++) {
    Mat<ZZ_p> g;
    Vec<ZZ_p> miss, p;

    git.GetNextGMP(g, miss, p);

    if (pid > 0) {
      pheno[i] = p[0];
      for (int j = 0; j < Param::NUM_COVS; j++) {
        cov[i][j] = p[1 + j];
      }
    }

    // In practice this would be done in one batch
    Mat<ZZ_p> g_mask;
    Vec<ZZ_p> miss_mask;
    mpc.BeaverPartition(g_mask, g);
    mpc.BeaverPartition(miss_mask, miss);

    if (pid > 0) {
      // Note: g_mask and miss_mask can be recovered from PRG and
      // need not be written
      mpc.WriteToFile(g, fs);
      mpc.WriteToFile(miss, fs);
    }

    if ((i + 1) % bsize == 0 || i == n - 1) {
      cout << "\t" << i+1 << " / " << n << ", "; toc(); tic();
    }
  }

  git.Terminate();

  fs.close();

  cout << "Finished writing Beaver partitioned genotype data" << endl;

  if (Param::DEBUG) {
    cout << "pheno" << endl;
    mpc.Print(pheno, 5);
    cout << "cov" << endl;
    mpc.Print(cov[0], 5);
  }

  fs.open(cache(pid, "input_pheno_cov").c_str(), ios::out | ios::binary);
  mpc.WriteToFile(pheno, fs);
  mpc.WriteToFile(cov, fs);
  fs.close();

  cout << "Finished writing phenotype and covariate data" << endl;

  return true;
}

bool gwas_protocol(MPCEnv& mpc, int pid) {
  SetNumThreads(Param::NUM_THREADS);
  cout << AvailableThreads() << " threads created" << endl;

  int n0 = Param::NUM_INDS;
  int m0 = Param::NUM_SNPS;
  int k = Param::NUM_DIM_TO_REMOVE;
  int kp = k + Param::NUM_OVERSAMPLE;

  cout << "n0: " << n0 << ", " << "m0: " << m0 << endl;

  // Shared variables
  string s;
  fstream fs;
  ofstream ofs;
  ifstream ifs;
  streampos strpos;
  int ind;
  Vec<ZZ_p> tmp_vec;
  Mat<ZZ_p> tmp_mat;

  mpc.ProfilerPushState("main");

  // Read in SNP list
  Vec<ZZ> snp_pos;
  Init(snp_pos, m0);

  ifs.open(Param::SNP_POS_FILE.c_str());
  if (!ifs.is_open()) {
    cout << "Could not open SNP_POS_FILE: " << Param::SNP_POS_FILE << endl;
    return false;
  }

  for (int i = 0; i < m0; i++) {
    long chrom, pos;
    ifs >> chrom >> pos;
    snp_pos[i] = ZZ(chrom) * 1e9 + ZZ(pos);
  }
  ifs.close();

  /* Useful constants */
  ZZ_p two(2), twoinv;
  inv(twoinv, two);

  ZZ_p fp_one = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);

  Vec<ZZ_p> pheno;
  Init(pheno, n0);

  Mat<ZZ_p> cov;
  Init(cov, n0, Param::NUM_COVS);

  if (!exists(cache(pid, "input_geno")) || !exists(cache(pid, "input_pheno_cov"))) {
    cout << "Initial data sharing results not found:" << endl;
    cout << "\t" << cache(pid, "input_geno") << endl;
    cout << "\t" << cache(pid, "input_pheno_cov") << endl;
    return false;
  }

  cout << "Initial data sharing results found" << endl;

  ifs.open(cache(pid, "input_pheno_cov").c_str(), ios::binary);
  mpc.ReadFromFile(pheno, ifs, n0);
  mpc.ReadFromFile(cov, ifs, n0, Param::NUM_COVS);
  ifs.close();

  cout << "Phenotypes and covariates loaded" << endl;

  if (Param::DEBUG) {
    cout << "pheno" << endl;
    mpc.Print(pheno, 5);
    cout << "cov" << endl;
    mpc.Print(cov[0], 5);
  }

  mpc.ProfilerPushState("qc");

  mpc.ProfilerPushState("snp_miss");

  Vec<ZZ_p> gkeep1;
  Init(gkeep1, m0);

  if (Param::SKIP_QC) {
    for (int i = 0; i < m0; i++) {
      gkeep1[i] = 1;
    }
    cout << "Locus missing rate filter skipped" << endl;
  } else {

    bool history;
    if (pid == 2) {
      history = exists(outname("gkeep1"));
      mpc.SendBool(history, 0);
      mpc.SendBool(history, 1);
    } else {
      // ask P2 if gkeep1 has been computed before
      history = mpc.ReceiveBool(2);
    }

    if (history) {
      cout << "Using locus missing rate filter from a previous run" << endl;
   
      if (pid == 2) {
        ifs.open(outname("gkeep1").c_str());
        for (int i = 0; i < m0; i++) {
          ifs >> gkeep1[i];
        }
        ifs.close();

        mpc.SendVec(gkeep1, 0);
        mpc.SendVec(gkeep1, 1);
      } else {
        mpc.ReceiveVec(gkeep1, 2, m0);
      }
    } else {
      Vec<ZZ_p> gmiss;
      Init(gmiss, m0);
      
      if (exists(cache(pid, "gmiss"))) {

        cout << "Locus missing rate cache found" << endl;

        ifs.open(cache(pid, "gmiss").c_str(), ios::binary);
        mpc.ReadFromFile(gmiss, ifs, m0);
        ifs.close();

      } else {

        cout << "Taking a pass to calculate locus missing rates:" << endl;

        if (pid > 0) {
          ifs.open(cache(pid, "input_geno").c_str(), ios::binary);

          mpc.ImportSeed(10, ifs);

          long bsize = n0 / 10;

          tic();
          for (int i = 0; i < n0; i++) {
            Vec<ZZ_p> miss, miss_mask;

            // Load stored Beaver partition
            mpc.SwitchSeed(10);
            mpc.RandMat(tmp_mat, 3, m0); // g_mask
            mpc.RandVec(miss_mask, m0);
            mpc.RestoreSeed();

            if (pid == 2) {
              mpc.SkipData(ifs, 3, m0);
              mpc.ReadFromFile(miss, ifs, m0);
            }

            // Recover secret shares from Beaver partition
            if (pid == 1) {
              miss = miss_mask;
            } else {
              miss += miss_mask;
            }

            // Add to running sum
            gmiss += miss;

            if ((i + 1) % bsize == 0 || i == n0 - 1) {
              cout << "\t" << i+1 << " / " << n0 << ", "; toc(); tic();
            }
          }

          ifs.close();
        }

        fs.open(cache(pid, "gmiss").c_str(), ios::out | ios::binary);
        mpc.WriteToFile(gmiss, fs);
        fs.close();

        cout << "Wrote results to cache" << endl;

      }

      if (Param::DEBUG) {
        cout << "gmiss" << endl;
        mpc.Print(gmiss, 5);
      }
      cout << "Locus missing rate filter ... "; tic();

      ZZ_p gmiss_ub = ZZ_p((long) (n0 * Param::GMISS_UB));
      mpc.LessThanPublic(gkeep1, gmiss, gmiss_ub);
      cout << "done. "; toc();
      
      mpc.RevealSym(gkeep1);

      if (pid == 2) {
        mpc.SendVec(gkeep1, 0);
      } else if (pid == 0) {
        mpc.ReceiveVec(gkeep1, 2, m0);
      }

      if (pid == 2) {
        ofs.open(outname("gkeep1").c_str());
        for (int i = 0; i < gkeep1.length(); i++) {
          ofs << gkeep1[i] << endl;
        }
        ofs.close();
      }

    }
  }

  uint m1 = conv<uint>(Sum(gkeep1));
  cout << "n0: " << n0 << ", " << "m1: " << m1 << endl;

  cout << "Filtering SNP position vector" << endl;
  FilterVec(snp_pos, gkeep1);

  mpc.ProfilerPopState(true); // snp_miss

  mpc.ProfilerPushState("ind_miss/het");

  Vec<ZZ_p> ikeep;
  Init(ikeep, n0);

  if (Param::SKIP_QC) {
    for (int i = 0; i < n0; i++) {
      ikeep[i] = 1;
    }
    cout << "Individual missing rate/het rate filters skipped" << endl;
  } else {

    bool history;
    if (pid == 2) {
      history = exists(outname("ikeep"));
      mpc.SendBool(history, 0);
      mpc.SendBool(history, 1);
    } else {
      // ask P2 if ikeep has been computed before
      history = mpc.ReceiveBool(2);
    }

    if (history) {
      cout << "Using individual missing rate/het rate filters from a previous run" << endl;
   
      if (pid == 2) {
        ifs.open(outname("ikeep").c_str());
        for (int i = 0; i < n0; i++) {
          ifs >> ikeep[i];
        }
        ifs.close();

        mpc.SendVec(ikeep, 0);
        mpc.SendVec(ikeep, 1);
      } else {
        mpc.ReceiveVec(ikeep, 2, n0);
      }
    } else {

      Vec<ZZ_p> imiss, ihet;
      Init(imiss, n0);
      Init(ihet, n0);

      if (exists(cache(pid, "imiss_ihet"))) {
        cout << "Individual missing rate and het rate cache found" << endl;

        ifs.open(cache(pid, "imiss_ihet").c_str(), ios::binary);
        mpc.ReadFromFile(imiss, ifs, n0);
        mpc.ReadFromFile(ihet, ifs, n0);
        ifs.close();

      } else {

        cout << "Taking a pass to calculate individual missing rates and het rates:" << endl;

        mpc.ProfilerPushState("data_scan");

        if (pid > 0) {
          ifs.open(cache(pid, "input_geno").c_str(), ios::binary);

          mpc.ImportSeed(10, ifs);

          long bsize = n0 / 10;

          tic();
          for (int i = 0; i < n0; i++) {
            Mat<ZZ_p> g, g_mask;
            Vec<ZZ_p> miss, miss_mask;

            // Load stored Beaver partition
            mpc.SwitchSeed(10);
            mpc.RandMat(g_mask, 3, m0);
            mpc.RandVec(miss_mask, m0);
            mpc.RestoreSeed();

            if (pid == 2) {
              mpc.ReadFromFile(g, ifs, 3, m0);
              mpc.ReadFromFile(miss, ifs, m0);
            }

            // Recover secret shares from Beaver partition
            if (pid == 1) {
              g = g_mask;
              miss = miss_mask;
            } else {
              g += g_mask;
              miss += miss_mask;
            }

            // Add to running sum
            for (int j = 0; j < m0; j++) {
              if (gkeep1[j] == 1) {
                imiss[i] += miss[j];
                ihet[i] += g[1][j];
              }
            }

            if ((i + 1) % bsize == 0 || i == n0 - 1) {
              cout << "\t" << i+1 << " / " << n0 << ", "; toc(); tic();
            }
          }
          ifs.close();
        }

        fs.open(cache(pid, "imiss_ihet").c_str(), ios::out | ios::binary);
        mpc.WriteToFile(imiss, fs);
        mpc.WriteToFile(ihet, fs);
        fs.close();

        cout << "Wrote results to cache" << endl;

        mpc.ProfilerPopState(false); // data_scan
      }

      mpc.ProfilerPushState("miss_filt");

      // Individual missingness filter
      cout << "Individual missing rate filter ... "; tic();
      ZZ_p imiss_ub = ZZ_p((long) (m1 * Param::IMISS_UB));
      mpc.LessThanPublic(ikeep, imiss, imiss_ub);
      cout << "done. "; toc();

      mpc.ProfilerPopState(true); // miss_filt
      mpc.ProfilerPushState("het_filt");

      // Individual heterozygosity filter
      cout << "Individual heterozygosity rate filter ... "; tic();
      ZZ_p ihet_ub_frac = DoubleToFP(Param::HET_UB, Param::NBIT_K, Param::NBIT_F);
      ZZ_p ihet_lb_frac = DoubleToFP(Param::HET_LB, Param::NBIT_K, Param::NBIT_F);

      // Number of observed SNPs per individual
      Vec<ZZ_p> m1_obs;
      Init(m1_obs, n0);
      if (pid > 0) {
        for (int i = 0; i < n0; i++) {
          m1_obs[i] = -imiss[i];
          if (pid == 1) {
            m1_obs[i] += m1;
          }
        }
      }

      Vec<ZZ_p> ihet_ub, ihet_lb;
      Init(ihet_ub, n0);
      Init(ihet_lb, n0);

      if (pid > 0) {
        for (int i = 0; i < n0; i++) {
          ihet_ub[i] = m1_obs[i] * ihet_ub_frac;
          ihet_lb[i] = m1_obs[i] * ihet_lb_frac;
          ihet[i] *= fp_one;
        }
      }

      Vec<ZZ_p> het_filt;
      mpc.LessThan(het_filt, ihet, ihet_ub);
      mpc.NotLessThan(tmp_vec, ihet, ihet_lb);
      mpc.MultElem(het_filt, het_filt, tmp_vec);

      mpc.MultElem(ikeep, ikeep, het_filt);
      het_filt.kill();
      cout << "done. "; toc();

      mpc.ProfilerPopState(true); // het_filt

      // Reveal samples to be filtered
      mpc.RevealSym(ikeep);

      if (pid == 2) {
        mpc.SendVec(ikeep, 0);
      } else if (pid == 0) {
        mpc.ReceiveVec(ikeep, 2, n0);
      }

      if (pid == 2) {
        ofs.open(outname("ikeep"));
        for (int i = 0; i < ikeep.length(); i++) {
          ofs << ikeep[i] << endl;
        }
        ofs.close();
      }
    }
  }

  mpc.ProfilerPopState(true); // ind_miss/het

  uint n1 = conv<uint>(Sum(ikeep));

  cout << "n1: " << n1 << ", " << "m1: " << m1 << endl;

  cout << "Filtering phenotypes and covariates" << endl;
  mpc.Filter(pheno, ikeep, n1);
  mpc.FilterRows(cov, ikeep, n1);

  Vec<ZZ_p> ctrl;
  mpc.FlipBit(ctrl, pheno);

  Vec<ZZ_p> ctrl_mask;
  mpc.BeaverPartition(ctrl_mask, ctrl);

  Vec<ZZ_p> dosage_sum;
  Vec<ZZ_p> gmiss, gmiss_ctrl, dosage_sum_ctrl;
  Mat<ZZ_p> g_count_ctrl;
  ZZ_p n1_ctrl(0);

  Init(gmiss, m1);
  Init(gmiss_ctrl, m1);
  Init(dosage_sum, m1);
  Init(dosage_sum_ctrl, m1);
  Init(g_count_ctrl, 3, m1);

  mpc.ProfilerPushState("data_scan");

  if (exists(cache(pid, "geno_stats"))) {
    cout << "Genotype statistics cache found" << endl;

    ifs.open(cache(pid, "geno_stats").c_str(), ios::binary);
    mpc.ReadFromFile(gmiss, ifs, m1);
    mpc.ReadFromFile(gmiss_ctrl, ifs, m1);
    mpc.ReadFromFile(dosage_sum, ifs, m1);
    mpc.ReadFromFile(dosage_sum_ctrl, ifs, m1);
    mpc.ReadFromFile(g_count_ctrl, ifs, 3, m1);
    mpc.ReadFromFile(n1_ctrl, ifs);
    ifs.close();
  } else {
    cout << "Taking a pass to calculate genotype statistics:" << endl;

    ifs.open(cache(pid, "input_geno").c_str(), ios::binary);
    if (pid > 0) {
      mpc.ImportSeed(10, ifs);
    } else {
      for (int p = 1; p <= 2; p++) {
        mpc.ImportSeed(10 + p, ifs);
      }
    }

    long report_bsize = n1 / 10;

    long bsize = Param::PITER_BATCH_SIZE;

    // Containers for batching the computation
    Vec< Mat<ZZ_p> > g, g_mask;
    Mat<ZZ_p> dosage, dosage_mask;
    Mat<ZZ_p> miss, miss_mask;
    Vec<ZZ_p> ctrl_vec, ctrl_mask_vec;
    g.SetLength(3);
    g_mask.SetLength(3);
    dosage.SetDims(bsize, m1);
    dosage_mask.SetDims(bsize, m1);
    miss.SetDims(bsize, m1);
    miss_mask.SetDims(bsize, m1);
    for (int k = 0; k < 3; k++) {
      g[k].SetDims(bsize, m1);
      g_mask[k].SetDims(bsize, m1);
    }
    ctrl_vec.SetLength(bsize);
    ctrl_mask_vec.SetLength(bsize);

    ind = -1;
    tic();
    for (int i = 0; i < n1; i++) {
      ind++;

      mpc.ProfilerPushState("file_io/rng");

      Mat<ZZ_p> g0, g0_mask;
      Vec<ZZ_p> miss0, miss0_mask;

      while (ikeep[ind] != 1) {
        if (pid > 0) {
          mpc.SkipData(ifs, 3, m0); // g
          mpc.SkipData(ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          }
        }
        ind++;
      }

      if (pid > 0) {
        mpc.ReadFromFile(g0, ifs, 3, m0); // g
        mpc.ReadFromFile(miss0, ifs, m0); // miss

        mpc.SwitchSeed(10);
        mpc.RandMat(g0_mask, 3, m0);
        mpc.RandVec(miss0_mask, m0);
        mpc.RestoreSeed();
      } else {
        Init(g0, 3, m0);
        Init(g0_mask, 3, m0);
        Init(miss0, m0);
        Init(miss0_mask, m0);

        for (int p = 1; p <= 2; p++) {
          mpc.SwitchSeed(10 + p);
          mpc.RandMat(tmp_mat, 3, m0);
          mpc.RandVec(tmp_vec, m0);
          mpc.RestoreSeed();

          g0_mask += tmp_mat;
          miss0_mask += tmp_vec;
        }
      }
      
      mpc.ProfilerPopState(false); // file_io/rng

      // Filter out loci that failed missing rate filter
      int ind2 = 0;
      for (int j = 0; j < m0; j++) {
        if (gkeep1[j] == 1) {
          for (int k = 0; k < 3; k++) {
            g[k][i % bsize][ind2] = g0[k][j];
            g_mask[k][i % bsize][ind2] = g0_mask[k][j];
          }
          miss[i % bsize][ind2] = miss0[j];
          miss_mask[i % bsize][ind2] = miss0_mask[j];
          ind2++;
        }
      }

      dosage[i % bsize] = g[1][i % bsize] + 2 * g[2][i % bsize];
      dosage_mask[i % bsize] = g_mask[1][i % bsize] + 2 * g_mask[2][i % bsize];

      ctrl_vec[i % bsize] = ctrl[i];
      ctrl_mask_vec[i % bsize] = ctrl_mask[i];

      // Update running sums
      if (pid > 0) {
        n1_ctrl += ctrl_mask[i];
        gmiss += miss_mask[i % bsize];
        dosage_sum += dosage_mask[i % bsize];

        if (pid == 1) {
          n1_ctrl += ctrl[i];
          gmiss += miss[i % bsize];
          dosage_sum += dosage[i % bsize];
        }
      }

      if (i % bsize == bsize - 1 || i == n1 - 1) {
        if (i % bsize < bsize - 1) {
          int new_bsize = (i % bsize) + 1;
          for (int k = 0; k < 3; k++) {
            g[k].SetDims(new_bsize, m1);
            g_mask[k].SetDims(new_bsize, m1);
          }
          dosage.SetDims(new_bsize, m1);
          dosage_mask.SetDims(new_bsize, m1);
          miss.SetDims(new_bsize, m1);
          miss_mask.SetDims(new_bsize, m1);
          ctrl_vec.SetLength(new_bsize);
          ctrl_mask_vec.SetLength(new_bsize);
        }

        mpc.BeaverMult(gmiss_ctrl, ctrl_vec, ctrl_mask_vec, miss, miss_mask);
        mpc.BeaverMult(dosage_sum_ctrl, ctrl_vec, ctrl_mask_vec, dosage, dosage_mask);
        for (int k = 0; k < 3; k++) {
          mpc.BeaverMult(g_count_ctrl[k], ctrl_vec, ctrl_mask_vec, g[k], g_mask[k]);
        }
      }

      if ((i + 1) % report_bsize == 0 || i == n1 - 1) {
        cout << "\t" << i+1 << " / " << n1 << ", "; toc(); tic();
      }
    }

    ifs.close();

    mpc.BeaverReconstruct(gmiss_ctrl);
    mpc.BeaverReconstruct(dosage_sum_ctrl);
    mpc.BeaverReconstruct(g_count_ctrl);

    // Write to cache
    fs.open(cache(pid, "geno_stats").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(gmiss, fs);
    mpc.WriteToFile(gmiss_ctrl, fs);
    mpc.WriteToFile(dosage_sum, fs);
    mpc.WriteToFile(dosage_sum_ctrl, fs);
    mpc.WriteToFile(g_count_ctrl, fs);
    mpc.WriteToFile(n1_ctrl, fs);
    fs.close();

    cout << "Wrote results to cache" << endl;
  }

  mpc.ProfilerPopState(true); // data_scan

  mpc.ProfilerPushState("maf/hwe");

  if (Param::DEBUG) {
    cout << "gmiss" << endl;
    mpc.Print(gmiss, 5);
    cout << "gmiss_ctrl" << endl;
    mpc.Print(gmiss_ctrl, 5);
    cout << "dosage_sum" << endl;
    mpc.Print(dosage_sum, 5);
    cout << "dosage_sum_ctrl" << endl;
    mpc.Print(dosage_sum_ctrl, 5);
    cout << "g_count_ctrl" << endl;
    for (int i = 0; i < 3; i++) {
      mpc.Print(g_count_ctrl[i], 5);
    }
  }

  mpc.ProfilerPushState("maf");

  // SNP MAF filter
  cout << "Locus minor allele frequency (MAF) filter ... " << endl;
  ZZ_p maf_lb = DoubleToFP(Param::MAF_LB, Param::NBIT_K, Param::NBIT_F);
  ZZ_p maf_ub = DoubleToFP(Param::MAF_UB, Param::NBIT_K, Param::NBIT_F);

  Vec<ZZ_p> dosage_tot, dosage_tot_ctrl;
  if (pid > 0) {
    dosage_tot = -gmiss;
    dosage_tot_ctrl = -gmiss_ctrl;
    mpc.AddPublic(dosage_tot, ZZ_p(n1));
    mpc.Add(dosage_tot_ctrl, n1_ctrl);
    dosage_tot *= 2;
    dosage_tot_ctrl *= 2;
  } else {
    dosage_tot.SetLength(m1);
    dosage_tot_ctrl.SetLength(m1);
  }

  cout << "Calculating MAFs ... " << endl; tic();
  Vec<ZZ_p> maf, maf_ctrl;
  if (exists(cache(pid, "maf"))) {
    cout << "maf cache found" << endl;
    ifs.open(cache(pid, "maf").c_str(), ios::binary);
    mpc.ReadFromFile(maf, ifs, dosage_tot.length());
    mpc.ReadFromFile(maf_ctrl, ifs, dosage_tot_ctrl.length());
    ifs.close();
  } else {
    mpc.ProfilerPushState("div");
    mpc.FPDiv(maf, dosage_sum, dosage_tot); 
    mpc.FPDiv(maf_ctrl, dosage_sum_ctrl, dosage_tot_ctrl); 
    mpc.ProfilerPopState(false); // div

    fs.open(cache(pid, "maf").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(maf, fs);
    mpc.WriteToFile(maf_ctrl, fs);
    fs.close();
  }
  cout << "done. "; toc();

  Vec<ZZ_p> Maf, Maf_ctrl; // MAJOR allele freq
  if (pid > 0) {
    Maf = -maf;
    Maf_ctrl = -maf_ctrl;
    mpc.AddPublic(Maf, fp_one);
    mpc.AddPublic(Maf_ctrl, fp_one);
  } else {
    Maf.SetLength(m1);
    Maf_ctrl.SetLength(m1);
  }

  // Variance based on Bernoulli distribution over each allele
  Vec<ZZ_p> g_var_bern;
  mpc.MultElem(g_var_bern, maf, Maf);
  mpc.Trunc(g_var_bern);

  mpc.ProfilerPopState(true); // maf

  if (Param::DEBUG) {
    cout << "maf" << endl;
    mpc.PrintFP(maf, 5);
    cout << "maf_ctrl" << endl;
    mpc.PrintFP(maf_ctrl, 5);
  }

  Vec<ZZ_p> gkeep2;
  Init(gkeep2, m1);

  if (Param::SKIP_QC) {
    for (int i = 0; i < m1; i++) {
      gkeep2[i] = 1;
    }
    cout << "SNP MAF/HWE filters skipped" << endl;
  } else {
    bool history;
    if (pid == 2) {
      history = exists(outname("gkeep2"));
      mpc.SendBool(history, 0);
      mpc.SendBool(history, 1);
    } else {
      // ask P2 if gkeep2 has been computed before
      history = mpc.ReceiveBool(2);
    }

    if (history) {
      cout << "Using MAF/HWE filters from a previous run" << endl;
   
      if (pid == 2) {
        ifs.open(outname("gkeep2").c_str());
        for (int i = 0; i < m1; i++) {
          ifs >> gkeep2[i];
        }
        ifs.close();

        mpc.SendVec(gkeep2, 0);
        mpc.SendVec(gkeep2, 1);
      } else {
        mpc.ReceiveVec(gkeep2, 2, m1);
      }
    } else {
      
      mpc.ProfilerPushState("maf_filt");

      mpc.LessThanPublic(gkeep2, maf, maf_ub);
      mpc.NotLessThanPublic(tmp_vec, maf, maf_lb);
      mpc.MultElem(gkeep2, gkeep2, tmp_vec);

      mpc.ProfilerPopState(true); // maf_filt

      mpc.ProfilerPushState("hwe_filt");

      cout << "Locus Hardy-Weinberg equilibrium (HWE) filter ... " << endl; tic();
      ZZ_p hwe_ub = DoubleToFP(Param::HWE_UB, Param::NBIT_K, Param::NBIT_F); // p < 1e-7
      
      // Calculate expected genotype distribution in control group
      Mat<ZZ_p> g_exp_ctrl;
      Init(g_exp_ctrl, 3, m1);

      mpc.MultElem(g_exp_ctrl[0], Maf_ctrl, Maf_ctrl); // AA
      mpc.MultElem(g_exp_ctrl[1], Maf_ctrl, maf_ctrl); // Aa
      if (pid > 0) {
        g_exp_ctrl[1] *= 2;
      }
      mpc.MultElem(g_exp_ctrl[2], maf_ctrl, maf_ctrl); // aa

      for (int i = 0; i < 3; i++) {
        mpc.MultElem(g_exp_ctrl[i], g_exp_ctrl[i], dosage_tot_ctrl);
      }
      g_exp_ctrl *= twoinv; // dosage_tot_ctrl is twice the # individuals we actually want

      mpc.Trunc(g_exp_ctrl);

      cout << "\tCalculated expected genotype counts, "; toc(); tic();

      Vec<ZZ_p> hwe_chisq; 
      Init(hwe_chisq, m1);

      if (exists(cache(pid, "hwe"))) {
        cout << "HWE cache found" << endl;
        ifs.open(cache(pid, "hwe").c_str(), ios::binary);
        mpc.ReadFromFile(hwe_chisq, ifs, m1);
        ifs.close();
      } else {
        for (int i = 0; i < 3; i++) {
          Vec<ZZ_p> diff;
          if (pid > 0) {
            diff = fp_one * g_count_ctrl[i] - g_exp_ctrl[i];
          } else {
            diff.SetLength(m1);
          }

          mpc.MultElem(diff, diff, diff); // square
          mpc.Trunc(diff);

          mpc.ProfilerPushState("div");
          mpc.FPDiv(tmp_vec, diff, g_exp_ctrl[i]);
          mpc.ProfilerPopState(false); // div
          hwe_chisq += tmp_vec;

          cout << "\tChi-square test (" << i+1 << "/3), "; toc(); tic();
        }

        fs.open(cache(pid, "hwe").c_str(), ios::out | ios::binary);
        mpc.WriteToFile(hwe_chisq, fs);
        fs.close();
      }

      if (Param::DEBUG) {
        cout << "hwe" << endl;
        mpc.PrintFP(hwe_chisq, 5);
      }
      
      Vec<ZZ_p> hwe_filt;
      mpc.LessThanPublic(hwe_filt, hwe_chisq, hwe_ub);
      mpc.MultElem(gkeep2, gkeep2, hwe_filt);
      hwe_filt.kill();

      // Reveal which SNPs to discard 
      mpc.RevealSym(gkeep2);
        
      if (pid == 2) {
        mpc.SendVec(gkeep2, 0);
      } else if (pid == 0) {
        mpc.ReceiveVec(gkeep2, 2, m1);
      }

      if (pid == 2) {
        ofs.open(outname("gkeep2").c_str());
        for (int i = 0; i < gkeep2.length(); i++) {
          ofs << gkeep2[i] << endl;
        }
        ofs.close();
      }

      mpc.ProfilerPopState(true); // hwe_filt
    }
  }

  uint m2 = conv<uint>(Sum(gkeep2));
  cout << "n1: " << n1 << ", " << "m2: " << m2 << endl;

  cout << "Filtering genotype statistics" << endl;
  mpc.Filter(g_var_bern, gkeep2, m2);
  mpc.Filter(maf, gkeep2, m2);
  FilterVec(snp_pos, gkeep2);

  gmiss.kill();
  gmiss_ctrl.kill();
  dosage_sum.kill();
  dosage_sum_ctrl.kill();
  g_count_ctrl.kill();

  mpc.ProfilerPopState(false); // maf/hwe
  mpc.ProfilerPopState(true); // qc
  mpc.ProfilerPushState("std_param");

  Vec<ZZ_p> g_std_bern_inv;
  if (exists(cache(pid, "stdinv_bern"))) {
    cout << "Genotype standard deviation cache found" << endl;

    ifs.open(cache(pid, "stdinv_bern").c_str(), ios::binary);
    mpc.ReadFromFile(g_std_bern_inv, ifs, g_var_bern.length());
    ifs.close();

  } else {
    cout << "Calculating genotype standard deviations (inverse)" << endl;

    mpc.ProfilerPushState("sqrt");
    mpc.FPSqrt(tmp_vec, g_std_bern_inv, g_var_bern);
    mpc.ProfilerPopState(false); // sqrt

    fs.open(cache(pid, "stdinv_bern").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(g_std_bern_inv, fs);
    fs.close();
  }

  if (Param::DEBUG) {
    cout << "g_std_bern_inv" << endl;
    mpc.PrintFP(g_std_bern_inv, 5);
  }

  Vec<ZZ_p> g_mean;
  if (pid > 0) {
    g_mean = 2 * maf;
  } else {
    g_mean.SetLength(m2);
  }

  mpc.ProfilerPopState(true); // std_param

  cout << "Starting population stratification analysis" << endl;

  mpc.ProfilerPushState("pop_strat");
  mpc.ProfilerPushState("select_snp");

  Vec<int8_t> selected; // 1 selected, 0 unselected, -1 TBD
  Vec<bool> to_process;
  selected.SetLength(m2);
  to_process.SetLength(m2);

  for (int i = 0; i < m2; i++) {
    selected[i] = -1;
  }

  ZZ dist_thres(Param::LD_DIST_THRES);
 
  ZZ prev(-1);
  for (int i = 0; i < m2; i++) {
    selected[i] = 0;
    if (prev < 0 || snp_pos[i] - prev > dist_thres) {
      selected[i] = 1;
      prev = snp_pos[i];
    }
  }

  // At this point "selected" contains the SNP filter for PCA, shared across all parties
  uint32_t m3 = 0;
  for (int i = 0; i < selected.length(); i++) {
    if (selected[i] == 1) {
      m3++;
    }
  }

  cout << "SNP selection complete: " << m3 << " / " << m2 << " selected" << endl;
  mpc.ProfilerPopState(false); // select_snp
  mpc.ProfilerPushState("reduce_file");

  // Cache the reduced G for PCA
  if (exists(cache(pid, "pca_input"))) {
    cout << "pca_input cache found" << endl;
  } else {
    Vec<bool> gkeep3;
    gkeep3.SetLength(m0);
    for (int j = 0; j < m0; j++) {
      gkeep3[j] = gkeep1[j] == 1;
    }

    ind = 0;
    for (int j = 0; j < m0; j++) {
      if (gkeep3[j]) {
        gkeep3[j] = gkeep2[ind] == 1;
        ind++;
      }
    }

    ind = 0;
    for (int j = 0; j < m0; j++) {
      if (gkeep3[j]) {
        gkeep3[j] = selected[ind] == 1;
        ind++;
      }
    }

    ifs.open(cache(pid, "input_geno").c_str(), ios::binary);
    if (pid > 0) {
      mpc.ImportSeed(10, ifs);
    } else {
      for (int p = 1; p <= 2; p++) {
        mpc.ImportSeed(10 + p, ifs);
      }
    }

    long bsize = n1 / 10;

    cout << "Caching input data for PCA:" << endl;

    fs.open(cache(pid, "pca_input").c_str(), ios::out | ios::binary);

    ind = -1;
    tic();
    for (int i = 0; i < n1; i++) {
      ind++;

      mpc.ProfilerPushState("file_io/rng");

      Mat<ZZ_p> g0, g0_mask;
      Vec<ZZ_p> miss0, miss0_mask;

      while (ikeep[ind] != 1) {
        if (pid > 0) {
          mpc.SkipData(ifs, 3, m0); // g
          mpc.SkipData(ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          }
        }
        ind++;
      }

      if (pid > 0) {
        mpc.ReadFromFile(g0, ifs, 3, m0); // g
        mpc.ReadFromFile(miss0, ifs, m0); // miss

        mpc.SwitchSeed(10);
        mpc.RandMat(g0_mask, 3, m0);
        mpc.RandVec(miss0_mask, m0);
        mpc.RestoreSeed();
      } else {
        Init(g0, 3, m0);
        Init(g0_mask, 3, m0);
        Init(miss0, m0);
        Init(miss0_mask, m0);

        for (int p = 1; p <= 2; p++) {
          mpc.SwitchSeed(10 + p);
          mpc.RandMat(tmp_mat, 3, m0);
          mpc.RandVec(tmp_vec, m0);
          mpc.RestoreSeed();

          g0_mask += tmp_mat;
          miss0_mask += tmp_vec;
        }
      }
      
      mpc.ProfilerPopState(false); // file_io/rng

      // Filter out loci that failed missing rate filter
      Mat<ZZ_p> g, g_mask;
      Vec<ZZ_p> miss, miss_mask;
      g.SetDims(3, m3);
      g_mask.SetDims(3, m3);
      miss.SetLength(m3);
      miss_mask.SetLength(m3);
      int ind2 = 0;
      for (int j = 0; j < m0; j++) {
        if (gkeep3[j]) {
          for (int k = 0; k < 3; k++) {
            g[k][ind2] = g0[k][j];
            g_mask[k][ind2] = g0_mask[k][j];
          }
          miss[ind2] = miss0[j];
          miss_mask[ind2] = miss0_mask[j];
          ind2++;
        }
      }

      Vec<ZZ_p> dosage, dosage_mask;
      dosage = g[1] + 2 * g[2];
      dosage_mask = g_mask[1] + 2 * g_mask[2];

      mpc.BeaverWriteToFile(dosage, dosage_mask, fs);
      mpc.BeaverWriteToFile(miss, miss_mask, fs);

      if ((i + 1) % bsize == 0 || i == n1 - 1) {
        cout << "\t" << i+1 << " / " << n1 << ", "; toc(); tic();
      }
    }

    ifs.close();
    fs.close();
  }

  mpc.ProfilerPopState(false); // reduce_file

  Vec<ZZ_p> g_mean_pca = g_mean;
  mpc.Filter(g_mean_pca, selected, m3);

  Vec<ZZ_p> g_stdinv_pca = g_std_bern_inv;
  mpc.Filter(g_stdinv_pca, selected, m3);

  Vec<ZZ_p> g_mean_pca_mask, g_stdinv_pca_mask;
  mpc.BeaverPartition(g_mean_pca_mask, g_mean_pca);
  mpc.BeaverPartition(g_stdinv_pca_mask, g_stdinv_pca);

  /* Pass 2: Random sketch */
  Mat<ZZ_p> Y_cur;
  Init(Y_cur, kp, m3);

  if (exists(cache(pid, "sketch"))) {

    cout << "sketch cache found" << endl;
    ifs.open(cache(pid, "sketch").c_str(), ios::in | ios::binary);
    ifs >> kp;
    mpc.ReadFromFile(Y_cur, ifs, kp, m3);
    ifs.close();

  } else {

    mpc.ProfilerPushState("rand_proj");
    
    Mat<ZZ_p> Y_cur_adj;
    Init(Y_cur_adj, kp, m3);

    Vec<int> bucket_count;
    bucket_count.SetLength(kp);
    for (int i = 0; i < kp; i++) {
      bucket_count[i] = 0;
    }

    ifs.open(cache(pid, "pca_input").c_str(), ios::in | ios::binary);
    for (int cur = 0; cur < n1; cur++) {
      // Count sketch (use global PRG)
      mpc.SwitchSeed(-1);
      long bucket_index = RandomBnd(kp);
      long rand_sign = RandomBnd(2) * 2 - 1;
      mpc.RestoreSeed();

      Vec<ZZ_p> g, g_mask, miss, miss_mask;
      mpc.BeaverReadFromFile(g, g_mask, ifs, m3);
      mpc.BeaverReadFromFile(miss, miss_mask, ifs, m3);

      // Flip miss bits so it points to places where g_mean should be subtracted
      mpc.BeaverFlipBit(miss, miss_mask);

      // Update running sum
      if (pid > 0) {
        Y_cur[bucket_index] += rand_sign * g_mask;
        if (pid == 1) {
          Y_cur[bucket_index] += rand_sign * g;
        }
      }

      // Update adjustment factor
      miss *= rand_sign;
      miss_mask *= rand_sign;
      mpc.BeaverMultElem(Y_cur_adj[bucket_index], miss, miss_mask, g_mean_pca, g_mean_pca_mask);

      bucket_count[bucket_index]++;
    }
    ifs.close();

    // Subtract the adjustment factor
    mpc.BeaverReconstruct(Y_cur_adj);
    if (pid > 0) {
      Y_cur = fp_one * Y_cur - Y_cur_adj;
    }
    Y_cur_adj.kill();

    if (Param::DEBUG) {
      cout << "Y_cur" << endl;
      mpc.PrintFP(Y_cur[0], 5);
      cout << "g_mean_pca" << endl;
      mpc.PrintBeaverFP(g_mean_pca, g_mean_pca_mask, 10);
      cout << "g_stdinv_pca" << endl;
      mpc.PrintBeaverFP(g_stdinv_pca, g_stdinv_pca_mask, 10);
    }

    // Get rid of empty buckets and normalize nonempty ones
    int empty_slot = 0;
    for (int i = 0; i < kp; i++) {
      if (bucket_count[i] > 0) {
        ZZ_p fp_count_inv = DoubleToFP(1 / ((double) bucket_count[i]), Param::NBIT_K, Param::NBIT_F);
        Y_cur[empty_slot] = Y_cur[i] * fp_count_inv;
        empty_slot++;
      }
    }
    kp = empty_slot;
    Y_cur.SetDims(kp, m3);
    mpc.Trunc(Y_cur);

    mpc.ProfilerPopState(true); // rand_proj

    fs.open(cache(pid, "sketch").c_str(), ios::out | ios::binary);
    fs << kp;
    if (pid > 0) {
      mpc.WriteToFile(Y_cur, fs);
    }
    fs.close();
  }

  mpc.ProfilerPushState("power_iter");

  Mat<ZZ_p> Y_cur_mask;
  mpc.BeaverPartition(Y_cur_mask, Y_cur);

  cout << "Initial sketch obtained, starting power iteration (num iter = " << Param::NUM_POWER_ITER << ")" << endl;
  tic();

  Mat<ZZ_p> gQ;

  if (exists(cache(pid, "piter"))) {

    cout << "piter cache found" << endl;
    ifs.open(cache(pid, "piter").c_str(), ios::in | ios::binary);
    mpc.ReadFromFile(gQ, ifs, n1, kp);
    ifs.close();
    
  } else {

    // Divide by standard deviation
    Mat<ZZ_p> Y;
    Init(Y, kp, m3);

    for (int i = 0; i < kp; i++) {
      mpc.BeaverMultElem(Y[i], Y_cur[i], Y_cur_mask[i], g_stdinv_pca, g_stdinv_pca_mask);
    }
    Y_cur.kill();
    Y_cur_mask.kill();

    mpc.BeaverReconstruct(Y);
    mpc.Trunc(Y);

    /* Calculate orthonormal bases of Y */
    Mat<ZZ_p> Q;
    mpc.ProfilerPushState("qr_m");
    mpc.OrthonormalBasis(Q, Y);
    mpc.ProfilerPopState(false); // qr_m
    Y.kill();

    Mat<ZZ_p> gQ_adj;
    Mat<ZZ_p> Q_mask;
    Mat<ZZ_p> Q_scaled, Q_scaled_mask;
    Mat<ZZ_p> Q_scaled_gmean, Q_scaled_gmean_mask;

    /* Power iteration */
    for (int pit = 0; pit <= Param::NUM_POWER_ITER; pit++) {
      /* This section is ran before each iteration AND once after all iterations */
      mpc.BeaverPartition(Q_mask, Q);

      // Normalize Q by standard deviations
      Init(Q_scaled, kp, m3);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(Q_scaled[i], Q[i], Q_mask[i], g_stdinv_pca, g_stdinv_pca_mask);
      }
      mpc.BeaverReconstruct(Q_scaled);
      mpc.Trunc(Q_scaled);

      mpc.BeaverPartition(Q_scaled_mask, Q_scaled);

      // Pre-multiply with g_mean to simplify calculation of centering matrix
      Init(Q_scaled_gmean, kp, m3);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(Q_scaled_gmean[i], Q_scaled[i], Q_scaled_mask[i],
                           g_mean_pca, g_mean_pca_mask);
      }
      mpc.BeaverReconstruct(Q_scaled_gmean);
      mpc.Trunc(Q_scaled_gmean);

      mpc.Transpose(Q_scaled); // m3-by-kp
      transpose(Q_scaled_mask, Q_scaled_mask); // m3-by-kp, unlike mpc.Transpose, P0 also transposes
      mpc.Transpose(Q_scaled_gmean); // m3-by-kp
      mpc.BeaverPartition(Q_scaled_gmean_mask, Q_scaled_gmean);

      Mat<ZZ_p> g, g_mask, miss, miss_mask;

      long bsize = Param::PITER_BATCH_SIZE;

      Init(g, bsize, m3);
      Init(g_mask, bsize, m3);
      Init(miss, bsize, m3);
      Init(miss_mask, bsize, m3);
      
      /* Pass 1 */
      Init(gQ, n1, kp);
      Init(gQ_adj, n1, kp);

      mpc.ProfilerPushState("data_scan0");

      mpc.ProfilerPushState("file_io");
      ifs.open(cache(pid, "pca_input").c_str(), ios::in | ios::binary);
      for (int cur = 0; cur < n1; cur++) {
        mpc.BeaverReadFromFile(g[cur % bsize], g_mask[cur % bsize], ifs, m3);
        mpc.BeaverReadFromFile(miss[cur % bsize], miss_mask[cur % bsize], ifs, m3);
        mpc.BeaverFlipBit(miss[cur % bsize], miss_mask[cur % bsize]);

        if (cur % bsize == bsize - 1) {
          mpc.ProfilerPopState(false); // file_io

          Init(tmp_mat, bsize, kp);
          mpc.BeaverMult(tmp_mat, g, g_mask, Q_scaled, Q_scaled_mask);
          for (int i = 0; i < bsize; i++) {
            gQ[cur-(bsize-1)+i] = tmp_mat[i];
          }

          Init(tmp_mat, bsize, kp);
          mpc.BeaverMult(tmp_mat, miss, miss_mask, Q_scaled_gmean, Q_scaled_gmean_mask);
          for (int i = 0; i < bsize; i++) {
            gQ_adj[cur-(bsize-1)+i] = tmp_mat[i];
          }

          mpc.ProfilerPushState("file_io");
        }
      }
      ifs.close();
      mpc.ProfilerPopState(false); // file_io

      long remainder = n1 % bsize;
      if (remainder > 0) {
        g.SetDims(remainder, m3);
        g_mask.SetDims(remainder, m3);
        miss.SetDims(remainder, m3);
        miss_mask.SetDims(remainder, m3);

        Init(tmp_mat, remainder, kp);
        mpc.BeaverMult(tmp_mat, g, g_mask, Q_scaled, Q_scaled_mask);
        for (int i = 0; i < remainder; i++) {
          gQ[n1-remainder+i] = tmp_mat[i];
        }

        Init(tmp_mat, remainder, kp);
        mpc.BeaverMult(tmp_mat, miss, miss_mask, Q_scaled_gmean, Q_scaled_gmean_mask);
        for (int i = 0; i < remainder; i++) {
          gQ_adj[n1-remainder+i] = tmp_mat[i];
        }

      }

      mpc.BeaverReconstruct(gQ);
      mpc.BeaverReconstruct(gQ_adj);
      if (pid > 0) {
        gQ -= gQ_adj;
      }

      mpc.ProfilerPopState(false); // data_scan1

      if (pit == Param::NUM_POWER_ITER) { // Quit if all iterations are performed
        break;
      }

      mpc.Transpose(gQ); // kp-by-n1
      mpc.ProfilerPushState("qr_n");
      mpc.OrthonormalBasis(Q, gQ);
      mpc.ProfilerPopState(false); // qr_n
      mpc.Transpose(Q); // n1-by-kp

      mpc.BeaverPartition(Q_mask, Q);

      Init(gQ, kp, m3);
      Init(gQ_adj, kp, m3);

      Init(g, bsize, m3);
      Init(g_mask, bsize, m3);
      Init(miss, bsize, m3);
      Init(miss_mask, bsize, m3);

      Mat<ZZ_p> Qsub, Qsub_mask;
      Init(Qsub, bsize, kp);
      Init(Qsub_mask, bsize, kp);

      mpc.ProfilerPushState("data_scan2");

      // Pass 2
      mpc.ProfilerPushState("file_io");
      ifs.open(cache(pid, "pca_input").c_str(), ios::in | ios::binary);
      for (int cur = 0; cur < n1; cur++) {
        mpc.BeaverReadFromFile(g[cur % bsize], g_mask[cur % bsize], ifs, m3);
        mpc.BeaverReadFromFile(miss[cur % bsize], miss_mask[cur % bsize], ifs, m3);
        mpc.BeaverFlipBit(miss[cur % bsize], miss_mask[cur % bsize]);

        Qsub[cur % bsize] = Q[cur];
        Qsub_mask[cur % bsize] = Q_mask[cur];

        if (cur % bsize == bsize - 1) {
          mpc.ProfilerPopState(false); // file_io

          mpc.Transpose(Qsub);
          transpose(Qsub_mask, Qsub_mask);

          mpc.BeaverMult(gQ, Qsub, Qsub_mask, g, g_mask);
          mpc.BeaverMult(gQ_adj, Qsub, Qsub_mask, miss, miss_mask);

          Qsub.SetDims(bsize, kp);
          Qsub_mask.SetDims(bsize, kp);

          mpc.ProfilerPushState("file_io");
        }
      }
      ifs.close();
      mpc.ProfilerPopState(false); // file_io

      remainder = n1 % bsize;
      if (remainder > 0) {
        g.SetDims(remainder, m3);
        g_mask.SetDims(remainder, m3);
        miss.SetDims(remainder, m3);
        miss_mask.SetDims(remainder, m3);
        Qsub.SetDims(remainder, kp);
        Qsub_mask.SetDims(remainder, kp);
        
        mpc.Transpose(Qsub);
        transpose(Qsub_mask, Qsub_mask);

        mpc.BeaverMult(gQ, Qsub, Qsub_mask, g, g_mask);
        mpc.BeaverMult(gQ_adj, Qsub, Qsub_mask, miss, miss_mask);
      }

      Qsub.kill();
      Qsub_mask.kill();
      g.kill();
      g_mask.kill();
      miss.kill();
      miss_mask.kill();

      mpc.BeaverReconstruct(gQ);
      mpc.BeaverReconstruct(gQ_adj);

      mpc.ProfilerPopState(false); // data_scan2

      Mat<ZZ_p> gQ_adj_mask;
      mpc.BeaverPartition(gQ_adj_mask, gQ_adj);

      Mat<ZZ_p> gQ_adj_gmean;
      Init(gQ_adj_gmean, kp, m3);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(gQ_adj_gmean[i], gQ_adj[i], gQ_adj_mask[i],
                           g_mean_pca, g_mean_pca_mask);
      }
      mpc.BeaverReconstruct(gQ_adj_gmean);
      mpc.Trunc(gQ_adj_gmean);

      if (pid > 0) {
        gQ -= gQ_adj_gmean;
      }
      gQ_adj_gmean.kill();

      Mat<ZZ_p> gQ_mask;
      mpc.BeaverPartition(gQ_mask, gQ);

      Mat<ZZ_p> gQ_scaled;
      gQ_scaled.SetDims(kp, m3);
      clear(gQ_scaled);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(gQ_scaled[i], gQ[i], gQ_mask[i], g_stdinv_pca, g_stdinv_pca_mask);
      }
      mpc.BeaverReconstruct(gQ_scaled);
      mpc.Trunc(gQ_scaled);

      mpc.ProfilerPushState("qr_m");
      mpc.OrthonormalBasis(Q, gQ_scaled);
      mpc.ProfilerPopState(false); // qr_m

      cout << "Iter " << pit + 1 << " complete, "; toc();
      tic();
    }

    fs.open(cache(pid, "piter").c_str(), ios::out | ios::binary);
    if (pid > 0) {
      mpc.WriteToFile(gQ, fs);
    }
    fs.close();

  }

  mpc.ProfilerPopState(true); // power_iter
  cout << "Power iteration complete" << endl;

  Mat<ZZ_p> Z = gQ;
  gQ.kill();

  cout << "Data projected to subspace" << endl;
  if (Param::DEBUG) {
    cout << "Z" << endl;
    mpc.PrintFP(Z[0], 5);
  }

  Mat<ZZ_p> V;
  Init(V, k, n1);

  /* Eigendecomposition */
  if (exists(cache(pid, "eigen"))) {

    cout << "eigen cache found" << endl;
    ifs.open(cache(pid, "eigen").c_str(), ios::binary);
    mpc.ReadFromFile(V, ifs, k, n1);
    ifs.close();

  } else {

    ZZ_p fp_m2_inv = DoubleToFP(1 / ((double) m2), Param::NBIT_K, Param::NBIT_F);
    Z *= fp_m2_inv;
    mpc.Trunc(Z);

    mpc.Transpose(Z); // kp-by-n1

    Mat<ZZ_p> Z_mask;
    mpc.BeaverPartition(Z_mask, Z);

    /* Form covariance matrix */
    Mat<ZZ_p> Z_gram;
    Init(Z_gram, kp, kp);
    for (int i = 0; i < kp; i++) {
      mpc.BeaverMult(Z_gram[i], Z, Z_mask, Z[i], Z_mask[i]);
    }
    mpc.BeaverReconstruct(Z_gram);
    mpc.Trunc(Z_gram);

    cout << "Constructed reduced eigenvalue problem" << endl;

    if (Param::DEBUG) {
      cout << "Z_gram" << endl;
      mpc.PrintFP(Z_gram[0], 5);
    }

    mpc.ProfilerPushState("eigen_solve");

    Mat<ZZ_p> U;
    Vec<ZZ_p> L;
    mpc.EigenDecomp(U, L, Z_gram);
    Z_gram.kill();

    // Select top eigenvectors and eigenvalues
    U.SetDims(k, kp);
    L.SetLength(k);

    cout << "Selected K eigenvectors" << endl;
    mpc.ProfilerPopState(false); // eigen_solve

    if (Param::DEBUG) {
      mpc.PrintFP(U[0], 5);
    }

    // Recover singular vectors
    Mat<ZZ_p> U_mask;
    mpc.BeaverPartition(U_mask, U);

    mpc.BeaverMultMat(V, U, U_mask, Z, Z_mask);
    U.kill();
    U_mask.kill();
    Z_mask.kill();
    mpc.BeaverReconstruct(V);
    mpc.Trunc(V);

    fs.open(cache(pid, "eigen").c_str(), ios::out | ios::binary);
    if (pid > 0) {
      mpc.WriteToFile(V, fs);
    }
    fs.close();

  }

  Z.kill();

  mpc.ProfilerPopState(true); // pop_strat

  mpc.ProfilerPushState("assoc_test");
  mpc.ProfilerPushState("covar");

  // Concatenate covariate matrix and jointly orthogonalize
  mpc.Transpose(cov);
  V.SetDims(k + Param::NUM_COVS, n1);
  if (pid > 0) {
    for (int i = 0; i < Param::NUM_COVS; i++) {
      V[k + i] = cov[i] * fp_one;
    }
  }
  cov.kill();
  mpc.OrthonormalBasis(V, V);

  Mat<ZZ_p> V_mask;
  mpc.BeaverPartition(V_mask, V);

  cout << "Bases for top singular vectors and covariates calculated" << endl;
  mpc.ProfilerPopState(false); // covar

  if (Param::DEBUG) {
    mpc.PrintBeaverFP(V[0], V_mask[0], 5);
  }

  /* Pass 4: Calculate GWAS statistics */

  Vec<ZZ_p> pheno_mask;
  mpc.BeaverPartition(pheno_mask, pheno);

  Vec<ZZ_p> Vp;
  Init(Vp, k + Param::NUM_COVS);
  mpc.BeaverMult(Vp, V, V_mask, pheno, pheno_mask);
  mpc.BeaverReconstruct(Vp);

  Vec<ZZ_p> Vp_mask;
  mpc.BeaverPartition(Vp_mask, Vp);
  
  Vec<ZZ_p> VVp;
  Init(VVp, n1);
  mpc.BeaverMult(VVp, Vp, Vp_mask, V, V_mask);
  mpc.BeaverReconstruct(VVp);
  mpc.Trunc(VVp);

  Vec<ZZ_p> VVp_mask;
  mpc.BeaverPartition(VVp_mask, VVp);

  Vec<ZZ_p> p_hat, p_hat_mask;
  p_hat = fp_one * pheno - VVp;
  p_hat_mask = fp_one * pheno_mask - VVp_mask;

  Vp.kill();
  Vp_mask.kill();
  VVp.kill();
  VVp_mask.kill();

  cout << "Phenotypes corrected" << endl;

  Vec<ZZ_p> V_sum, V_sum_mask;
  Init(V_sum, k + Param::NUM_COVS);
  Init(V_sum_mask, k + Param::NUM_COVS);
  for (int i = 0; i < k + Param::NUM_COVS; i++) {
    for (int j = 0; j < n1; j++) {
      V_sum[i] += V[i][j];
      V_sum_mask[i] += V_mask[i][j];
    }
  }

  Vec<ZZ_p> u;
  Init(u, n1);
  mpc.BeaverMult(u, V_sum, V_sum_mask, V, V_mask);
  mpc.BeaverReconstruct(u);
  mpc.Trunc(u);
  if (pid > 0) {
    u *= -1;
    mpc.AddPublic(u, fp_one);
  }

  Vec<ZZ_p> u_mask;
  mpc.BeaverPartition(u_mask, u);

  if (Param::DEBUG) {
    cout << "u" << endl;
    mpc.PrintBeaverFP(u, u_mask, 10);
  }

  cout << "Allocating sx, sxx, sxp, B ... ";

  Vec<ZZ_p> sx, sxx, sxp;
  Mat<ZZ_p> B;
  Init(sx, m2);
  Init(sxx, m2);
  Init(sxp, m2);
  Init(B, k + Param::NUM_COVS, m2);

  cout << "done.";

  mpc.ProfilerPushState("data_scan");

  if (exists(cache(pid, "gwas_stats"))) {
    cout << "GWAS statistics cache found" << endl;
    ifs.open(cache(pid, "gwas_stats").c_str(), ios::binary);
    mpc.ReadFromFile(sx, ifs, m2);
    mpc.ReadFromFile(sxx, ifs, m2);
    mpc.ReadFromFile(sxp, ifs, m2);
    mpc.ReadFromFile(B, ifs, k + Param::NUM_COVS, m2);
    ifs.close();

  } else {

    ifs.open(cache(pid, "input_geno").c_str(), ios::binary);
    if (pid > 0) {
      mpc.ImportSeed(10, ifs);
    } else {
      for (int p = 1; p <= 2; p++) {
        mpc.ImportSeed(10 + p, ifs);
      }
    }

    long bsize = Param::PITER_BATCH_SIZE;

    cout << "Allocating batch variables ... ";

    Mat<ZZ_p> dosage, dosage_mask;
    Init(dosage, bsize, m2);
    Init(dosage_mask, bsize, m2);

    Vec<ZZ_p> u_vec, u_mask_vec, p_hat_vec, p_hat_mask_vec;
    Init(u_vec, bsize);
    Init(u_mask_vec, bsize);
    Init(p_hat_vec, bsize);
    Init(p_hat_mask_vec, bsize);

    mpc.Transpose(V); // n1-by-(k + NUM_COVS)
    transpose(V_mask, V_mask);

    Mat<ZZ_p> V_sub, V_mask_sub;
    Init(V_sub, bsize, k + Param::NUM_COVS);
    Init(V_mask_sub, bsize, k + Param::NUM_COVS);

    cout << "done." << endl;

    Vec<bool> gkeep3;
    gkeep3.SetLength(m0);
    for (int j = 0; j < m0; j++) {
      gkeep3[j] = gkeep1[j] == 1;
    }

    ind = 0;
    for (int j = 0; j < m0; j++) {
      if (gkeep3[j]) {
        gkeep3[j] = gkeep2[ind] == 1;
        ind++;
      }
    }

    ind = -1;
    tic();
    mpc.ProfilerPushState("file_io/rng");
    cout << "GWAS pass:" << endl;
    for (int cur = 0; cur < n1; cur++) {
      ind++;

      Mat<ZZ_p> g0, g0_mask;
      Vec<ZZ_p> miss0, miss0_mask;

      while (ikeep[ind] != 1) {
        if (pid > 0) {
          mpc.SkipData(ifs, 3, m0); // g
          mpc.SkipData(ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          }
        }
        ind++;
      }

      if (pid > 0) {
        mpc.ReadFromFile(g0, ifs, 3, m0); // g
        mpc.ReadFromFile(miss0, ifs, m0); // miss

        mpc.SwitchSeed(10);
        mpc.RandMat(g0_mask, 3, m0);
        mpc.RandVec(miss0_mask, m0);
        mpc.RestoreSeed();
      } else {
        Init(g0, 3, m0);
        Init(g0_mask, 3, m0);
        Init(miss0, m0);
        Init(miss0_mask, m0);

        for (int p = 1; p <= 2; p++) {
          mpc.SwitchSeed(10 + p);
          mpc.RandMat(tmp_mat, 3, m0);
          mpc.RandVec(tmp_vec, m0);
          mpc.RestoreSeed();

          g0_mask += tmp_mat;
          miss0_mask += tmp_vec;
        }
      }
      
      Mat<ZZ_p> g, g_mask;
      Vec<ZZ_p> miss, miss_mask;
      g.SetDims(3, m2);
      miss.SetLength(m2);
      g_mask.SetDims(3, m2);
      miss_mask.SetLength(m2);
      int ind2 = 0;
      for (int j = 0; j < m0; j++) {
        if (gkeep3[j]) {
          for (int k = 0; k < 3; k++) {
            g[k][ind2] = g0[k][j];
            g_mask[k][ind2] = g0_mask[k][j];
          }
          miss[ind2] = miss0[j];
          miss_mask[ind2] = miss0_mask[j];
          ind2++;
        }
      }

      dosage[cur % bsize] = g[1] + 2 * g[2];
      dosage_mask[cur % bsize] = g_mask[1] + 2 * g_mask[2];

      u_vec[cur % bsize] = u[cur];
      u_mask_vec[cur % bsize] = u_mask[cur];
      p_hat_vec[cur % bsize] = p_hat[cur];
      p_hat_mask_vec[cur % bsize] = p_hat_mask[cur];

      V_sub[cur % bsize] = V[cur];
      V_mask_sub[cur % bsize] = V_mask[cur];

      if (cur % bsize == bsize - 1) {
        mpc.ProfilerPopState(false); // file_io/rng

        mpc.BeaverMult(sx, u_vec, u_mask_vec, dosage, dosage_mask);
        mpc.BeaverMult(sxp, p_hat_vec, p_hat_mask_vec, dosage, dosage_mask);

        Mat<ZZ_p> sxx_tmp;
        Init(sxx_tmp, bsize, m2);
        mpc.BeaverMultElem(sxx_tmp, dosage, dosage_mask, dosage, dosage_mask);
        for (int b = 0; b < bsize; b++) {
          sxx += sxx_tmp[b];
        }
        sxx_tmp.kill();

        mpc.Transpose(V_sub); // (k + NUM_COVS)-by-bsize
        transpose(V_mask_sub, V_mask_sub);

        mpc.BeaverMult(B, V_sub, V_mask_sub, dosage, dosage_mask);

        cout << "\t" << cur+1 << " / " << n1 << ", "; toc(); tic();

        Init(dosage, bsize, m2);
        Init(dosage_mask, bsize, m2);
        Init(V_sub, bsize, k + Param::NUM_COVS);
        Init(V_mask_sub, bsize, k + Param::NUM_COVS);

        mpc.ProfilerPushState("file_io/rng");
      }
    }
    ifs.close();
    mpc.ProfilerPopState(false); // file_io/rng

    long remainder = n1 % bsize;
    if (remainder > 0) {
      dosage.SetDims(remainder, m2);
      dosage_mask.SetDims(remainder, m2);
      u_vec.SetLength(remainder);
      u_mask_vec.SetLength(remainder);
      p_hat_vec.SetLength(remainder);
      p_hat_mask_vec.SetLength(remainder);
      V_sub.SetDims(remainder, k + Param::NUM_COVS);
      V_mask_sub.SetDims(remainder, k + Param::NUM_COVS);

      mpc.BeaverMult(sx, u_vec, u_mask_vec, dosage, dosage_mask);
      mpc.BeaverMult(sxp, p_hat_vec, p_hat_mask_vec, dosage, dosage_mask);

      Mat<ZZ_p> sxx_tmp;
      Init(sxx_tmp, remainder, m2);
      mpc.BeaverMultElem(sxx_tmp, dosage, dosage_mask, dosage, dosage_mask);
      for (int b = 0; b < remainder; b++) {
        sxx += sxx_tmp[b];
      }
      sxx_tmp.kill();

      mpc.Transpose(V_sub); // (k + NUM_COVS)-by-remainder
      transpose(V_mask_sub, V_mask_sub);

      mpc.BeaverMult(B, V_sub, V_mask_sub, dosage, dosage_mask);

      cout << "\t" << n1 << " / " << n1 << ", "; toc(); tic();
    }

    mpc.BeaverReconstruct(sx);
    mpc.BeaverReconstruct(sxp);
    mpc.BeaverReconstruct(sxx);
    mpc.BeaverReconstruct(B);
    sxx *= fp_one;

    fs.open(cache(pid, "gwas_stats").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(sx, fs);
    mpc.WriteToFile(sxx, fs);
    mpc.WriteToFile(sxp, fs);
    mpc.WriteToFile(B, fs);
    fs.close();

    cout << "Wrote results to cache" << endl;
  }

  mpc.ProfilerPopState(true); // data_scan

  if (Param::DEBUG) {
    cout << "sx" << endl;
    mpc.PrintFP(sx, 3);
    cout << "sxp" << endl;
    mpc.PrintFP(sxp, 3);
    cout << "sxx" << endl;
    mpc.PrintFP(sxx, 3);
    cout << "B" << endl;
    mpc.PrintFP(B, 3, 3);
  }

  mpc.Transpose(B); // m2-by-(k + Param::NUM_COVS)

  Vec<ZZ_p> BB;
  mpc.InnerProd(BB, B); // m2
  mpc.Trunc(BB);
  if (pid > 0) {
    sxx -= BB;
  }

  ZZ_p sp(0);
  if (pid > 0) {
    for (int i = 0; i < n1; i++) {
      sp += p_hat_mask[i];
      if (pid == 1) {
        sp += p_hat[i];
      }
    }
  }

  ZZ_p spp(0);
  mpc.BeaverInnerProd(spp, p_hat, p_hat_mask);
  mpc.BeaverReconstruct(spp);

  ZZ_p fp_n1_inv = DoubleToFP(1 / ((double) n1), Param::NBIT_K, Param::NBIT_F);
  sx *= fp_n1_inv;
  sp *= fp_n1_inv;

  mpc.Trunc(sx);
  mpc.Trunc(sp);
  mpc.Trunc(spp);

  Vec<ZZ_p> sx_mask;
  mpc.BeaverPartition(sx_mask, sx);

  ZZ_p sp_mask;
  mpc.BeaverPartition(sp_mask, sp);

  Vec<ZZ_p> spsx, sx2;
  ZZ_p sp2(0);
  Init(spsx, m2);
  Init(sx2, m2);

  mpc.BeaverMult(spsx, sx, sx_mask, sp, sp_mask);
  mpc.BeaverMult(sp2, sp, sp_mask, sp, sp_mask);
  mpc.BeaverMultElem(sx2, sx, sx_mask, sx, sx_mask);

  mpc.BeaverReconstruct(spsx);
  mpc.BeaverReconstruct(sp2);
  mpc.BeaverReconstruct(sx2);

  spsx *= n1;
  sp2 *= n1;
  sx2 *= n1;

  mpc.Trunc(spsx);
  mpc.Trunc(sp2);
  mpc.Trunc(sx2);

  Vec<ZZ_p> numer, denom;
  Init(numer, m2);
  Init(denom, m2 + 1);
  if (pid > 0) {
    numer = sxp - spsx;
    for (int i = 0; i < m2; i++) {
      denom[i] = sxx[i] - sx2[i];
    }
    denom[m2] = spp - sp2;
  }

  Vec<ZZ_p> denom1_sqrt_inv;
  if (exists(cache(pid, "denom_inv"))) {
    cout << "denom_inv cache found" << endl;
    ifs.open(cache(pid, "denom_inv").c_str(), ios::binary);
    mpc.ReadFromFile(denom1_sqrt_inv, ifs, denom.length());
    ifs.close();
  } else {
    mpc.ProfilerPushState("sqrt");
    mpc.FPSqrt(tmp_vec, denom1_sqrt_inv, denom);
    mpc.ProfilerPopState(false); // sqrt

    fs.open(cache(pid, "denom_inv").c_str(), ios::out | ios::binary);
    if (pid > 0) {
      mpc.WriteToFile(denom1_sqrt_inv, fs);
    }
    fs.close();
  }

  denom.kill();
  tmp_vec.kill();

  ZZ_p denom2_sqrt_inv = denom1_sqrt_inv[m2]; // p term
  denom1_sqrt_inv.SetLength(m2); // truncate

  Vec<ZZ_p> z;
  mpc.MultElem(z, numer, denom1_sqrt_inv);
  mpc.Trunc(z);

  mpc.MultMat(z, z, denom2_sqrt_inv);
  mpc.Trunc(z);

  mpc.ProfilerPopState(false); // assoc_test

  cout << "Association statistics calculated" << endl;
  mpc.RevealSym(z);
  if (pid == 2) {
    Vec<double> z_double;
    FPToDouble(z_double, z, Param::NBIT_K, Param::NBIT_F);
    ofs.open(outname("assoc").c_str(), ios::out);
    for (int i = 0; i < z_double.length(); i++) {
      ofs << z_double[i] << endl;
    }
    ofs.close();
    cout << "Result written to " << outname("assoc") << endl;
  }

  mpc.ProfilerPopState(true); // main

  return true;
}

#endif
