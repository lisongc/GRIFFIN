#ifndef GRIFFIN_PREDICTION_EW_INPUT_SCHEME_H
#define GRIFFIN_PREDICTION_EW_INPUT_SCHEME_H

namespace griffin {

// Electroweak input schemes belong to the complete prediction, not to one
// particular BSM model.  The SMEFT namespace re-exports this type for source
// compatibility with the original SMEFT API.
enum class EWInputScheme {
  GmuAlphaMZ,
  AlphaMWMZ
};

inline const char* ewInputSchemeName(EWInputScheme scheme)
{
  switch(scheme)
  {
    case EWInputScheme::GmuAlphaMZ:
      return "{Gmu, alpha, MZ}";
    case EWInputScheme::AlphaMWMZ:
      return "{alpha, MW, MZ}";
  }
  return "unknown electroweak input scheme";
}

} // namespace griffin

#endif // GRIFFIN_PREDICTION_EW_INPUT_SCHEME_H
