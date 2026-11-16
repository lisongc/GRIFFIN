/*
  File: tools.cc
  Purpose: Implement decay-width utilities and the angular conversion of
           supplied amplitudes into a differential cross section.
  Original authors: Lisong Chen and Ayres Freitas
  Integration contributors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
  The original author/revision record is preserved below.
*/

/*-----------------------------------------------------------------------------
tools.cc
Lisong Chen (lic114@pitt.edu), Ayres Freitas (afreitas@pitt.edu)
last revision: 19 Sep 2022
-------------------------------------------------------------------------------
some useful functions for computing decay widths and cross-sections
-----------------------------------------------------------------------------*/

#include "griffin/sm/tools.h"
#include "griffin/sm/ff0.h"

#include <cmath>
#include <stdexcept>

namespace griffin {

double differentialCrossSectionFromAmplitudes(
    const std::array<Cplx,4>& amplitudes, Fermion initial, Fermion final,
    double invariantMassSquared, double cosTheta, double unitConversion)
{
  if(initial != Fermion::electron)
    throw std::invalid_argument(
        "massless annihilation cross section currently requires an electron initial state");
  double color;
  switch(final)
  {
    case Fermion::muon:
    case Fermion::tau:
      color = 1.0;
      break;
    case Fermion::d:
    case Fermion::u:
    case Fermion::s:
    case Fermion::c:
    case Fermion::b:
      color = 3.0;
      break;
    case Fermion::electron:
      throw std::domain_error(
          "Bhabha scattering requires dedicated s/t-channel and scalar/pseudoscalar structures");
    case Fermion::nuElectron:
    case Fermion::nuMuon:
    case Fermion::nuTau:
      throw std::domain_error(
          "neutrino final states require a dedicated process implementation");
    default:
      throw std::invalid_argument(
          "unsupported final state for the massless annihilation cross section");
  }
  if(!std::isfinite(invariantMassSquared) || !std::isfinite(cosTheta)
      || !std::isfinite(unitConversion))
    throw std::invalid_argument("cross-section kinematics must be finite");
  if(invariantMassSquared <= 0.0 || cosTheta < -1.0 || cosTheta > 1.0)
    throw std::domain_error(
        "cross-section kinematics are outside the physical domain");
  if(unitConversion <= 0.0)
    throw std::domain_error("cross-section unit conversion must be positive");

  const double even = 1.0+cosTheta*cosTheta;
  const double odd = 4.0*cosTheta;
  const Cplx squared = even*(
      amplitudes[0]*conj(amplitudes[0])
      + amplitudes[1]*conj(amplitudes[1])
      + amplitudes[2]*conj(amplitudes[2])
      + amplitudes[3]*conj(amplitudes[3]))
      + odd*(amplitudes[0]*conj(amplitudes[3])
             + amplitudes[2]*conj(amplitudes[1]));
  const double prefactor = color*invariantMassSquared/(32.0*Pi)*unitConversion;
  return prefactor*real(squared);
}

// Note: partzwidth will modify the FA and FV objects passed to it,
// by setting the fermion type and input
double partzwidth(FA_SMLO& fa, FV_SMLO& fv, const int type, 
                  const inval& input, int scheme)
{
  double alp = input.get(al)/Pi * (1+ input.get(Delal)), 
         alsp = input.get(als)/Pi;
  double x = sqr(input.get(MZ)/input.get(MT));
  double mf, fac, y, RV, RA;
  switch(type)
  {
    case ELE: mf = input.get(ME);  alsp = 0;  fac = 1;  break;
    case MUO: mf = input.get(MM);  alsp = 0;  fac = 1;  break;
    case TAU: mf = input.get(ML);  alsp = 0;  fac = 1;  break;
    case NUE:
    case NUM:
    case NUT: mf = 0;  alsp = 0;  fac = 1;  break;
    case UQU: mf = input.get(MU);  fac = 3;  break;
    case CQU: mf = input.get(MC);  fac = 3;  break;
    case DQU: mf = input.get(MD);  fac = 3;  break;
    case SQU: mf = input.get(MS);  fac = 3;  break;
    case BQU: mf = input.get(MB);  fac = 3;  break;
    default: return(0);
  }
  // QCD+QED corrections from Phys.Rep.277,189; PRL 101,012002; PRL 108,222003;
  //  PLB 287,209
  y = sqr(mf/input.get(MZ));
  RV = 1+ 0.75*sqr(Qf[type])*alp + alsp
       + (-1.1529539789360381 + x*(0.13037037037037036 
           - 0.02962962962962963*log(x)))*sqr(Qf[type])*sqr(alp)
       - 0.25*sqr(Qf[type])*alp*alsp
       + (1.409230409109778 + x*(0.06518518518518518
           - 0.014814814814814815*log(x)))*sqr(alsp)
       - 12.767064735769054*powint(alsp,3) - 80.0075*powint(alsp,4);
  RA = RV - 2*I3f[type]*(
        (-3.083333 + 0.086420*x + 0.013167*x*x + log(x))*sqr(alsp)
	+ (-15.98773 + 3.72222*log(x) + 1.91667*sqr(log(x)))*powint(alsp,3)
	+ (49.0309 - 17.6637*log(x) + 14.6597*sqr(log(x)) 
	   + 3.6736*powint(log(x),3))*powint(alsp,4));
  RV += 12*y*alsp - 6*sqr(y);
  RA += -6*y - 22*y*alsp + 6*sqr(y);
  fa.setftype(type);
  fa.setinput(input);
  fv.setftype(type);
  fv.setinput(input);
  if(scheme == RUNWIDTHSCHEME)
    fac *= sqrt(1+sqr(input.get(GamZ)/input.get(MZ)));
  return(input.get(MZ)*fac/(12*Pi) 
  	 * (RA*real(fa.result()) + RV*real(fv.result())));
}

// Note: zwidth will modify the FA and FV objects passed to it,
// by setting the fermion type and input
double zwidth(FA_SMLO& fa, FV_SMLO& fv, const inval& input, int scheme)
{
  int gen;
  double res=0;
  for(gen=0; gen<=2; gen++)
    res += partzwidth(fa, fv, ELE+2*gen, input, scheme) 
          + partzwidth(fa, fv, NUE+2*gen, input, scheme)
	  + partzwidth(fa, fv, UQU+2*gen, input, scheme) 
          + partzwidth(fa, fv, DQU+2*gen, input, scheme);
  return(res);
}

} // namespace
