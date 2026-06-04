/* deltar.h: header file for deltar.cc */

#ifndef __deltarSMEFT__
#define __deltarSMEFT__

#include "deltar.h"

// define indices for the Wilson coeffs

#define CphiD     30
#define CphiWB    31    
#define Cphil1    32 
#define Cphil3    33   
#define Cphie     34  
#define Cphiq1    35  
#define Cphiq3    36  
#define Cphiu     37
#define Cphid     38
#define Cll       39
#define Clq1      40
#define Clq3      41
#define Cee       42
#define Ceu       43
#define Ced       44  
#define Cle       45
#define Clu       46
#define Cld       47
#define Cqe       48

namespace griffin {

// Delta r predicted in the SMEFT (at LO) and the SM (at NNLO+)
class dr_SMEFTLO : public dr_SMNNLO {
  public:
  using dr_SMNNLO::dr_SMNNLO;
  double resSMEFT(void) const;
  Cplx result(void) const{
    return(dr_SMNNLO::result()+resSMEFT());
  }
};

/* input class that computes MW from Gmu */
class invalGmuSMEFT : public invalGmu {
  protected:
    vector<vector<vector<vector<vector<double>>>>> dataSMEFTG;
    void compute(void);
  public:
    using invalGmu::invalGmu;
    using invalGmu::set;
    using invalGmu::get;
    invalGmuSMEFT(const int sizealloc) : dataSMEFTG(sizealloc, vector<vector<vector<vector<double>>>>(3,vector<vector<vector<double>>>(3,vector<vector<double>>(3,vector<double>(3,NAN))))) {};
    invalGmuSMEFT(void) : invalGmuSMEFT(SIZE1) {};
    invalGmuSMEFT(const inval& copyfrom) : invalGmu(copyfrom), dataSMEFTG(SIZE1, vector<vector<vector<vector<double>>>>(3,vector<vector<vector<double>>>(3,vector<vector<double>>(3,vector<double>(3,NAN))))) {};

    //set Wilson coefficient values
    void set(const int idx, const double val) override
    {
      if(idx < data.size() & idx > 29)
        dataSMEFTG[idx][0][0][0][0] = val;
      else if(idx < data.size()){
        data[idx] = val;
      } else
      {
        cerr << "Input value index " << idx << " outside of range" << endl;
        exit(1);
      }
      compute();
    }

    void set(const int idx, const double val, const int gen1, const int gen2)
    {
      if(idx < data.size() & idx > 29)
        dataSMEFTG[idx][gen1-1][gen2-1][0][0] = val;
      else if(idx < data.size()){
        data[idx] = val;
      } else
      {
        cerr << "Input value index " << idx << " outside of range" << endl;
        exit(1);
      }
        compute();
    }
    
    void set(const int idx, const double val, const int gen1, const int gen2, const int gen3, const int gen4)
    {
      if(idx < data.size() & idx > 29)
      {
        dataSMEFTG[idx][gen1-1][gen2-1][gen3-1][gen4-1] = val;
      }
      else if(idx < data.size())
      {
        data[idx] = val;
      } else
      {
        cerr << "Input value index " << idx << " outside of range" << endl;
        exit(1);
      }
      compute();
    }

    //get Wilson coefficient values
    double get(const int idx) const override
    {
      if(idx < data.size() & idx > 29)
      {
        if(isfinite(dataSMEFTG[idx][0][0][0][0]))
          return(dataSMEFTG[idx][0][0][0][0]);
      } else
      {
        if(isfinite(data[idx]))
          return(data[idx]); 
      }
      cerr << "Invalid or undefined input value for index " << idx << endl;
      exit(1);
    }

    double get(const int idx, const int gen1, const int gen2) const override
    {
      if(idx < data.size() & idx > 29)
      {
        if(isfinite(dataSMEFTG[idx][gen1-1][gen2-1][0][0]))
          return(dataSMEFTG[idx][gen1-1][gen2-1][0][0]);
      } else
      {
        if(isfinite(data[idx]))
          return(data[idx]); 
      }
      cerr << "Invalid or undefined input value for index " << idx << endl;
      exit(1);
    }
    
    double get(const int idx, const int gen1, const int gen2, const int gen3, const int gen4) const override
    {
      if(idx < data.size() & idx > 29)
      {
        if(isfinite(dataSMEFTG[idx][gen1-1][gen2-1][gen3-1][gen4-1]))
          return(dataSMEFTG[idx][gen1-1][gen2-1][gen3-1][gen4-1]);
      } else
      {
        if(isfinite(data[idx]))
          return(data[idx]); 
      }
      cerr << "Invalid or undefined input value for index " << idx << endl;
      exit(1);
    }
};

} // namespace

#endif // __deltarSMEFT__
