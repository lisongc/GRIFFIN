/* SMvalGMwMz.h: header file for SMvalGMwMz.cc */

#include "SMval.h"

namespace griffin {

/* input class for Gmu-MW-MZ scheme, which computes alpha from these inputs */
class SMvalGMwMz : public SMval {
protected:
  void compute(void);
public:
  using SMval::SMval;
  SMvalGMwMz(void) : SMval() {};
  SMvalGMwMz(const inval& copyfrom) : SMval(copyfrom) { compute(); };
};

} // namespace
