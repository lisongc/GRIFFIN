#ifndef GRIFFIN_PREDICTION_PERTURBATIVE_SERIES_H
#define GRIFFIN_PREDICTION_PERTURBATIVE_SERIES_H

#include <string>
#include <vector>

#include "prediction/Order.h"

namespace griffin {
namespace prediction {

template<typename Value>
struct PerturbativeTerm {
  OrderTag order;
  Value value;
  std::string label;
};

template<typename Value>
class PerturbativeSeries {
private:
  std::vector<PerturbativeTerm<Value> > terms_;

public:
  void add(const OrderTag& order, const Value& value,
           const std::string& label)
  {
    PerturbativeTerm<Value> term = {order,value,label};
    terms_.push_back(term);
  }

  const std::vector<PerturbativeTerm<Value> >& terms() const
  {
    return terms_;
  }

  Value selectedSum(const Accuracy& accuracy) const
  {
    Value result = Value();
    for(typename std::vector<PerturbativeTerm<Value> >::const_iterator it =
          terms_.begin(); it != terms_.end(); ++it)
      if(accuracy.accepts(it->order))
        result += it->value;
    return result;
  }

  Value exact(const OrderTag& order) const
  {
    Value result = Value();
    for(typename std::vector<PerturbativeTerm<Value> >::const_iterator it =
          terms_.begin(); it != terms_.end(); ++it)
      if(it->order == order)
        result += it->value;
    return result;
  }

  Value sumThrough(LoopOrder order) const
  {
    Value result = Value();
    for(typename std::vector<PerturbativeTerm<Value> >::const_iterator it =
          terms_.begin(); it != terms_.end(); ++it)
      if(orderIndex(it->order.loopOrder) <= orderIndex(order))
        result += it->value;
    return result;
  }

  bool hasLoopOrder(LoopOrder order) const
  {
    for(typename std::vector<PerturbativeTerm<Value> >::const_iterator it =
          terms_.begin(); it != terms_.end(); ++it)
      if(it->order.loopOrder == order)
        return true;
    return false;
  }

  bool empty() const { return terms_.empty(); }
};

} // namespace prediction
} // namespace griffin

#endif // GRIFFIN_PREDICTION_PERTURBATIVE_SERIES_H
