#ifndef GRIFFIN_PREDICTION_EW_INPUT_SCHEME_H
#define GRIFFIN_PREDICTION_EW_INPUT_SCHEME_H

namespace griffin {

// Shared electroweak input-scheme identifiers. Use griffin::EWInputScheme;
// there is no separate SMEFT-namespace alias.
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
