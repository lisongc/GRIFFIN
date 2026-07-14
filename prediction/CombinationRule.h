#ifndef GRIFFIN_PREDICTION_COMBINATION_RULE_H
#define GRIFFIN_PREDICTION_COMBINATION_RULE_H

namespace griffin {
namespace prediction {

enum class CombinationKind {
  LinearEFT,
  AdditiveCorrection,
  FullModelDifference,
  ExactModel
};

class CombinationRule {
private:
  CombinationKind kind_;

  explicit CombinationRule(CombinationKind kind) : kind_(kind) {}

public:
  static CombinationRule linearEFT()
  {
    return CombinationRule(CombinationKind::LinearEFT);
  }

  static CombinationRule additiveCorrection()
  {
    return CombinationRule(CombinationKind::AdditiveCorrection);
  }

  static CombinationRule fullModelDifference()
  {
    return CombinationRule(CombinationKind::FullModelDifference);
  }

  static CombinationRule exactModel()
  {
    return CombinationRule(CombinationKind::ExactModel);
  }

  CombinationKind kind() const { return kind_; }

  const char* name() const
  {
    switch(kind_)
    {
      case CombinationKind::LinearEFT: return "LinearEFT";
      case CombinationKind::AdditiveCorrection: return "AdditiveCorrection";
      case CombinationKind::FullModelDifference:
        return "FullModelDifference";
      case CombinationKind::ExactModel: return "ExactModel";
    }
    return "unknown combination rule";
  }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_COMBINATION_RULE_H
