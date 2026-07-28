/* SMEFTval.h: define input class that includes values for Wilsonian SMEFT coeffs */

#ifndef __SMEFTval__
#define __SMEFTval__

#include "SMval.h"

namespace griffin {

class SMEFTval : public SMval {
protected:
  vector<vector<vector<vector<vector<double>>>>> dataSMEFT;
public:
  using SMval::SMval;
  using SMval::set;
  using SMval::get;
  SMEFTval(const int sizealloc) : dataSMEFT(sizealloc, vector<vector<vector<vector<double>>>>(3,vector<vector<vector<double>>>(3,vector<vector<double>>(3,vector<double>(3,NAN))))) {};
  SMEFTval(void) : SMEFTval(SIZE1) {};
  SMEFTval(const inval& copyfrom) : SMval(copyfrom), dataSMEFT(SIZE1, vector<vector<vector<vector<double>>>>(3,vector<vector<vector<double>>>(3,vector<vector<double>>(3,vector<double>(3,NAN))))) {};

  // set Wilson coefficient values
  void set(const int idx, const double val) override
  {
    if(idx < data.size() & idx > 29)
      dataSMEFT[idx][0][0][0][0] = val;
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
      dataSMEFT[idx][gen1-1][gen2-1][0][0] = val;
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
      dataSMEFT[idx][gen1-1][gen2-1][gen3-1][gen4-1] = val;
    else if(idx < data.size()){
      data[idx] = val;
    } else
    {
      cerr << "Input value index " << idx << " outside of range" << endl;
      exit(1);
    }
    compute();
  }

  // get Wilson coefficient values
  double get(const int idx) const override
  {
    if(idx < data.size() & idx > 29)
    {
      if(isfinite(dataSMEFT[idx][0][0][0][0]))
        return(dataSMEFT[idx][0][0][0][0]);
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
      if(isfinite(dataSMEFT[idx][gen1-1][gen2-1][0][0]))
        return(dataSMEFT[idx][gen1-1][gen2-1][0][0]);
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
      if(isfinite(dataSMEFT[idx][gen1-1][gen2-1][gen3-1][gen4-1]))
        return(dataSMEFT[idx][gen1-1][gen2-1][gen3-1][gen4-1]);
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

#endif // __SMEFTval__
