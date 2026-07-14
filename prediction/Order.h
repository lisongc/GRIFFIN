#ifndef GRIFFIN_PREDICTION_ORDER_H
#define GRIFFIN_PREDICTION_ORDER_H

#include <stdexcept>

namespace griffin {
namespace prediction {

// NNLOPlus denotes the content of GRIFFIN's existing highest-order SM
// classes: NNLO together with selected known higher-order terms.
enum class LoopOrder {
  LO = 0,
  NLO = 1,
  NNLOPlus = 2
};

inline const char* loopOrderName(LoopOrder order)
{
  switch(order)
  {
    case LoopOrder::LO: return "LO";
    case LoopOrder::NLO: return "NLO";
    case LoopOrder::NNLOPlus: return "NNLOPlus";
  }
  return "unknown order";
}

inline int orderIndex(LoopOrder order)
{
  return static_cast<int>(order);
}

enum class ContributionSector {
  ReferenceSM,
  ModelCorrection
};

struct OrderTag {
  ContributionSector sector;
  LoopOrder loopOrder;
  int qcdOrder;
  // 0: no inverse-scale expansion, 1: 1/Lambda^2, 2: 1/Lambda^4, ...
  int inverseScalePower;

  static OrderTag referenceSM(LoopOrder order, int qcd = 0)
  {
    OrderTag tag = {ContributionSector::ReferenceSM,order,qcd,0};
    return tag;
  }

  static OrderTag model(LoopOrder order, int inversePower, int qcd = 0)
  {
    OrderTag tag = {
      ContributionSector::ModelCorrection,order,qcd,inversePower
    };
    return tag;
  }
};

inline bool operator==(const OrderTag& left, const OrderTag& right)
{
  return left.sector == right.sector
      && left.loopOrder == right.loopOrder
      && left.qcdOrder == right.qcdOrder
      && left.inverseScalePower == right.inverseScalePower;
}

class Accuracy {
private:
  LoopOrder smOrder_;
  LoopOrder modelOrder_;
  bool includeModel_;
  int maxInverseScalePower_;

public:
  Accuracy(LoopOrder smOrder, LoopOrder modelOrder,
           bool includeModel, int maxInverseScalePower)
    : smOrder_(smOrder), modelOrder_(modelOrder),
      includeModel_(includeModel),
      maxInverseScalePower_(maxInverseScalePower)
  {
    if(maxInverseScalePower < 0)
      throw std::invalid_argument("inverse-scale accuracy cannot be negative");
  }

  static Accuracy smOnly(LoopOrder smOrder)
  {
    return Accuracy(smOrder,LoopOrder::LO,false,0);
  }

  static Accuracy smPlusModel(LoopOrder smOrder, LoopOrder modelOrder,
                              int maxInverseScalePower = 1)
  {
    return Accuracy(smOrder,modelOrder,true,maxInverseScalePower);
  }

  static Accuracy SMLO()
  {
    return smOnly(LoopOrder::LO);
  }

  static Accuracy SMNLO()
  {
    return smOnly(LoopOrder::NLO);
  }

  static Accuracy SMNNLOPlus()
  {
    return smOnly(LoopOrder::NNLOPlus);
  }

  static Accuracy SMLO_plus_ModelLO()
  {
    return smPlusModel(LoopOrder::LO,LoopOrder::LO,1);
  }

  static Accuracy SMNLO_plus_ModelLO()
  {
    return smPlusModel(LoopOrder::NLO,LoopOrder::LO,1);
  }

  static Accuracy SMNNLOPlus_plus_ModelLO()
  {
    return smPlusModel(LoopOrder::NNLOPlus,LoopOrder::LO,1);
  }

  LoopOrder smOrder() const { return smOrder_; }
  LoopOrder modelOrder() const { return modelOrder_; }
  bool includesModel() const { return includeModel_; }
  int maxInverseScalePower() const { return maxInverseScalePower_; }

  bool accepts(const OrderTag& tag) const
  {
    if(tag.sector == ContributionSector::ReferenceSM)
      return orderIndex(tag.loopOrder) <= orderIndex(smOrder_);
    return includeModel_
        && orderIndex(tag.loopOrder) <= orderIndex(modelOrder_)
        && tag.inverseScalePower <= maxInverseScalePower_;
  }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_ORDER_H
