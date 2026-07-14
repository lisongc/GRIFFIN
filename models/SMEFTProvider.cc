#include "models/SMEFTProvider.h"

#include <stdexcept>

#include "ff0.h"
#include "smeft/EWPOZ_smeft.h"
#include "smeft/deltar_smeft.h"
#include "smeft/xsc_smeft.h"

namespace griffin {
namespace prediction {

PerturbativeSeries<Cplx> SMEFTProvider::pseudoObservableSeries(
    PseudoObservable observable, Fermion fermion) const
{
  const inval& input = context_.input();
  const EWInputScheme scheme = context_.inputScheme();
  Cplx correction = 0.0;

  switch(observable)
  {
    case PseudoObservable::SW:
      correction = smeft::SW_SMEFTLO(
          fermion,input,*eftInput_,scheme).resD6Tree();
      break;
    case PseudoObservable::FA:
      correction = smeft::FA_SMEFTLO(
          fermion,input,*eftInput_,scheme).resD6Tree();
      break;
    case PseudoObservable::FV:
      correction = smeft::FV_SMEFTLO(
          fermion,input,*eftInput_,scheme).resD6Tree();
      break;
    case PseudoObservable::DeltaR:
      correction = smeft::dr_SMEFTLO(
          input,*eftInput_,scheme).resD6Tree();
      break;
    case PseudoObservable::WMass:
      if(scheme != EWInputScheme::GmuAlphaMZ)
        correction = 0.0;
      else
        correction = smeft::derivedMWShiftD6(input,*eftInput_).deltaMW;
      break;
    default:
      throw std::invalid_argument("unsupported SMEFT pseudo-observable");
  }

  PerturbativeSeries<Cplx> series;
  series.add(OrderTag::model(LoopOrder::LO,1),correction,
             "linear dimension-six tree correction");
  return series;
}

PerturbativeSeries<Cplx> SMEFTProvider::amplitudeSeries(
    const ScatteringPoint& point, const CurrentChannel& channel) const
{
  const inval& input = context_.input();
  FA_SMLO faInitial(point.initial,input);
  FA_SMLO faFinal(point.final,input);
  SW_SMLO swInitial(point.initial,input);
  SW_SMLO swFinal(point.final,input);

  smeft::mat_SMEFTLO amplitude(
      point.initial,point.final,channel.initial,channel.final,
      faInitial,faFinal,swInitial,swFinal,
      point.s,point.cosTheta,input,*eftInput_,context_.inputScheme());

  PerturbativeSeries<Cplx> series;
  series.add(OrderTag::model(LoopOrder::LO,1),
             amplitude.resultD6Tree(),
             "linear dimension-six tree amplitude");
  return series;
}

} // namespace prediction
} // namespace griffin
