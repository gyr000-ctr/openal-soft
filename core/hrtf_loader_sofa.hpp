#ifndef CORE_HRTF_LOADER_SOFA_HPP
#define CORE_HRTF_LOADER_SOFA_HPP

#include <iosfwd>
#include <memory>


struct HrtfStore;

/* Loads an AES69 SOFA (Spatially Oriented Format for Acoustics) file as an
 * HRTF data set, in the SimpleFreeFieldHRIR convention. Only available when
 * built with runtime SOFA support (HAVE_MYSOFA).
 *
 * Throws one of the SofaLoadError subclasses declared in hrtf_loader.hpp on
 * invalid, unsupported, incomplete, or oversized data.
 */
auto LoadSofaHrtf(std::istream &stream) -> std::unique_ptr<HrtfStore>;

#endif /* CORE_HRTF_LOADER_SOFA_HPP */
