#include "models/SMProvider.h"

#include <stdexcept>

#include "EWPOZ.h"
#include "EWPOZ2.h"
#include "deltar.h"
#include "xscnnlo.h"

namespace griffin {
namespace prediction {
namespace {

Cplx smPseudoObservableAt(PseudoObservable observable, Fermion fermion,
                          LoopOrder order, const inval& input)
{
  switch(observable)
  {
    case PseudoObservable::SW:
      if(order == LoopOrder::LO) return SW_SMLO(fermion,input).result();
      if(order == LoopOrder::NLO) return SW_SMNLO(fermion,input).result();
      return SW_SMNNLO(fermion,input).result();

    case PseudoObservable::FA:
      if(order == LoopOrder::LO) return FA_SMLO(fermion,input).result();
      if(order == LoopOrder::NLO) return FA_SMNLO(fermion,input).result();
      return FA_SMNNLO(fermion,input).result();

    case PseudoObservable::FV:
      if(order == LoopOrder::LO) return FV_SMLO(fermion,input).result();
      if(order == LoopOrder::NLO) return FV_SMNLO(fermion,input).result();
      return FV_SMNNLO(fermion,input).result();

    case PseudoObservable::DeltaR:
      if(order == LoopOrder::LO) return 0.0;
      if(order == LoopOrder::NLO) return dr_SMNLO(input).result();
      return dr_SMNNLO(input).result();

    case PseudoObservable::WMass:
      // The context mass is an input or a separately derived reference value.
      // Its derivation accuracy is carried by the input object, not recomputed
      // by changing the requested observable accuracy.
      return input.get(InputPar::WMassComplex);
  }
  throw std::invalid_argument("unsupported SM pseudo-observable");
}

Cplx smAmplitudeAt(const ScatteringPoint& point,
                   const CurrentChannel& channel,
                   LoopOrder order, const inval& input)
{
  if(order == LoopOrder::LO)
  {
    FA_SMLO faInitial(point.initial,input);
    FA_SMLO faFinal(point.final,input);
    SW_SMLO swInitial(point.initial,input);
    SW_SMLO swFinal(point.final,input);
    matel amplitude(point.initial,point.final,channel.initial,channel.final,
                    faInitial,faFinal,swInitial,swFinal,
                    point.s,point.cosTheta,input);
    return amplitude.result();
  }

  if(order == LoopOrder::NLO)
  {
    FA_SMNLO faInitial(point.initial,input);
    FA_SMNLO faFinal(point.final,input);
    SW_SMNLO swInitial(point.initial,input);
    SW_SMNLO swFinal(point.final,input);
    matel amplitude(point.initial,point.final,channel.initial,channel.final,
                    faInitial,faFinal,swInitial,swFinal,
                    point.s,point.cosTheta,input);
    return amplitude.result();
  }

  FA_SMNNLO faInitial(point.initial,input);
  FA_SMNNLO faFinal(point.final,input);
  SW_SMNNLO swInitial(point.initial,input);
  SW_SMNNLO swFinal(point.final,input);
  mat_SMNNLO amplitude(point.initial,point.final,
                       channel.initial,channel.final,
                       faInitial,faFinal,swInitial,swFinal,
                       point.s,point.cosTheta,input);
  return amplitude.result();
}

} // namespace

PerturbativeSeries<Cplx> SMProvider::pseudoObservableSeries(
    PseudoObservable observable, Fermion fermion) const
{
  PerturbativeSeries<Cplx> series;
  const inval& input = context_.input();

  if(observable == PseudoObservable::WMass)
  {
    series.add(OrderTag::referenceSM(LoopOrder::LO),
               smPseudoObservableAt(observable,fermion,LoopOrder::LO,input),
               "context reference MW");
    return series;
  }

  const Cplx lo = smPseudoObservableAt(
      observable,fermion,LoopOrder::LO,input);
  const Cplx nlo = smPseudoObservableAt(
      observable,fermion,LoopOrder::NLO,input);
  const Cplx nnlo = smPseudoObservableAt(
      observable,fermion,LoopOrder::NNLOPlus,input);

  series.add(OrderTag::referenceSM(LoopOrder::LO),lo,
             "SM LO");
  series.add(OrderTag::referenceSM(LoopOrder::NLO),nlo-lo,
             "SM NLO increment");
  series.add(OrderTag::referenceSM(LoopOrder::NNLOPlus),nnlo-nlo,
             "SM NNLOPlus increment");
  return series;
}

PerturbativeSeries<Cplx> SMProvider::amplitudeSeries(
    const ScatteringPoint& point, const CurrentChannel& channel) const
{
  const inval& input = context_.input();
  const Cplx lo = smAmplitudeAt(point,channel,LoopOrder::LO,input);
  const Cplx nlo = smAmplitudeAt(point,channel,LoopOrder::NLO,input);
  const Cplx nnlo = smAmplitudeAt(
      point,channel,LoopOrder::NNLOPlus,input);

  PerturbativeSeries<Cplx> series;
  series.add(OrderTag::referenceSM(LoopOrder::LO),lo,
             "SM LO amplitude");
  series.add(OrderTag::referenceSM(LoopOrder::NLO),nlo-lo,
             "SM NLO amplitude increment");
  series.add(OrderTag::referenceSM(LoopOrder::NNLOPlus),nnlo-nlo,
             "SM NNLOPlus amplitude increment");
  return series;
}

} // namespace prediction
} // namespace griffin
