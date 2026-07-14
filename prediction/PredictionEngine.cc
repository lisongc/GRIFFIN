#include "prediction/PredictionEngine.h"

#include <stdexcept>
#include <string>

#include "xscmassless.h"

namespace griffin {
namespace prediction {
namespace {

template<typename Value>
void requireOrder(const PerturbativeSeries<Value>& series,
                  LoopOrder order, const char* providerName,
                  const char* quantity)
{
  if(!series.hasLoopOrder(order))
    throw std::domain_error(
        std::string(providerName)+" does not provide "+quantity+" at "+
        loopOrderName(order));
}

} // namespace

void PredictionEngine::validateContexts() const
{
  if(reference_->outputKind() != ProviderOutputKind::ReferenceTheory)
    throw std::invalid_argument(
        "PredictionEngine reference provider is not a reference theory");

  if(model_ && !reference_->context().compatibleWith(model_->context()))
    throw std::invalid_argument(
        "reference and model providers use different calculation contexts");
}

PredictionResult<Cplx> PredictionEngine::pseudoObservable(
    PseudoObservable observable, Fermion fermion,
    const Accuracy& accuracy) const
{
  if(observable == PseudoObservable::SW && isNeutrino(fermion))
    throw std::domain_error("SW is undefined for a neutral fermion");

  const PerturbativeSeries<Cplx> referenceSeries =
      reference_->pseudoObservableSeries(observable,fermion);
  if(observable != PseudoObservable::WMass)
    requireOrder(referenceSeries,accuracy.smOrder(),reference_->name(),
                 pseudoObservableName(observable));
  const Cplx referenceValue = referenceSeries.selectedSum(accuracy);

  PerturbativeSeries<Cplx> modelSeries;
  Cplx modelCorrection = 0.0;
  if(accuracy.includesModel())
  {
    if(!model_)
      throw std::invalid_argument(
          "model accuracy requested without a model provider");
    modelSeries = model_->pseudoObservableSeries(observable,fermion);
    requireOrder(modelSeries,accuracy.modelOrder(),model_->name(),
                 pseudoObservableName(observable));

    switch(combination_.kind())
    {
      case CombinationKind::LinearEFT:
        if(model_->outputKind() != ProviderOutputKind::AdditiveCorrection)
          throw std::invalid_argument(
              "LinearEFT requires an additive correction provider");
        if(accuracy.maxInverseScalePower() != 1)
          throw std::invalid_argument(
              "LinearEFT requires truncation at exactly 1/Lambda^2");
        modelCorrection = modelSeries.selectedSum(accuracy);
        break;

      case CombinationKind::AdditiveCorrection:
        if(model_->outputKind() != ProviderOutputKind::AdditiveCorrection)
          throw std::invalid_argument(
              "AdditiveCorrection requires an additive provider");
        modelCorrection = modelSeries.selectedSum(accuracy);
        break;

      case CombinationKind::FullModelDifference:
        if(model_->outputKind() != ProviderOutputKind::FullModel)
          throw std::invalid_argument(
              "FullModelDifference requires a full-model provider");
        modelCorrection = modelSeries.sumThrough(accuracy.modelOrder())
            - referenceSeries.sumThrough(accuracy.modelOrder());
        break;

      case CombinationKind::ExactModel:
        if(model_->outputKind() != ProviderOutputKind::FullModel)
          throw std::invalid_argument(
              "ExactModel requires a full-model provider");
        modelCorrection = modelSeries.sumThrough(accuracy.modelOrder())
            - referenceValue;
        break;
    }
  }

  return PredictionResult<Cplx>(
      referenceValue,modelCorrection,referenceSeries,modelSeries);
}

PredictionResult<double> PredictionEngine::differentialCrossSection(
    const ScatteringPoint& point, const Accuracy& accuracy) const
{
  const MasslessFermionCrossSectionBuilder builder(
      point.initial,point.final,point.s,point.cosTheta,point.unitConversion);
  if(accuracy.includesModel()
      && combination_.kind() != CombinationKind::LinearEFT)
    throw std::domain_error(
        "the current cross-section builder supports only LinearEFT model corrections");
  if(accuracy.includesModel() && accuracy.modelOrder() != LoopOrder::LO)
    throw std::domain_error(
        "the current cross-section builder has no model NLO amplitudes or real radiation");

  PerturbativeSeries<Cplx> referenceAmplitudeSeries[4];
  PerturbativeSeries<Cplx> modelAmplitudeSeries[4];
  MasslessFermionCrossSectionBuilder::Amplitudes smLO = {};
  MasslessFermionCrossSectionBuilder::Amplitudes smNLO = {};
  MasslessFermionCrossSectionBuilder::Amplitudes smNNLO = {};
  MasslessFermionCrossSectionBuilder::Amplitudes modelLO = {};

  for(std::size_t channel = 0; channel < builder.channelCount(); ++channel)
  {
    const CurrentChannel currentChannel(
        builder.initialCurrent(channel),builder.finalCurrent(channel));
    referenceAmplitudeSeries[channel] =
        reference_->amplitudeSeries(point,currentChannel);
    requireOrder(referenceAmplitudeSeries[channel],accuracy.smOrder(),
                 reference_->name(),"current amplitude");

    smLO[channel] = referenceAmplitudeSeries[channel].sumThrough(
        LoopOrder::LO);
    smNLO[channel] = referenceAmplitudeSeries[channel].sumThrough(
        LoopOrder::NLO);
    smNNLO[channel] = referenceAmplitudeSeries[channel].sumThrough(
        LoopOrder::NNLOPlus);
    modelLO[channel] = 0.0;

    if(accuracy.includesModel())
    {
      if(!model_)
        throw std::invalid_argument(
            "model accuracy requested without a model provider");
      if(model_->outputKind() != ProviderOutputKind::AdditiveCorrection)
        throw std::invalid_argument(
            "LinearEFT cross sections require an additive model provider");
      modelAmplitudeSeries[channel] =
          model_->amplitudeSeries(point,currentChannel);
      requireOrder(modelAmplitudeSeries[channel],LoopOrder::LO,
                   model_->name(),"current amplitude");
      modelLO[channel] = modelAmplitudeSeries[channel].exact(
          OrderTag::model(LoopOrder::LO,1));
    }
  }

  const double sigmaLO = builder.squaredPrediction(smLO);
  const double sigmaNLO = builder.squaredPrediction(smNLO);
  const double sigmaNNLO = builder.squaredPrediction(smNNLO);

  PerturbativeSeries<double> referenceSeries;
  referenceSeries.add(OrderTag::referenceSM(LoopOrder::LO),sigmaLO,
                      "SM LO cumulative cross section");
  referenceSeries.add(OrderTag::referenceSM(LoopOrder::NLO),
                      sigmaNLO-sigmaLO,
                      "SM NLO cumulative-prediction increment");
  referenceSeries.add(OrderTag::referenceSM(LoopOrder::NNLOPlus),
                      sigmaNNLO-sigmaNLO,
                      "SM NNLOPlus cumulative-prediction increment");
  const double referenceValue = referenceSeries.selectedSum(accuracy);

  PerturbativeSeries<double> modelSeries;
  double modelCorrection = 0.0;
  if(accuracy.includesModel())
  {
    modelCorrection = builder.linearInterference(smLO,modelLO);
    modelSeries.add(OrderTag::model(LoopOrder::LO,1),modelCorrection,
                    "tree-SM / tree-D6 interference");
  }

  return PredictionResult<double>(
      referenceValue,modelCorrection,referenceSeries,modelSeries);
}

} // namespace prediction
} // namespace griffin
