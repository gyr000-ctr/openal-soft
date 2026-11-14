#ifndef CORE_HRTF_LOADER_SOFA_HPP
#define CORE_HRTF_LOADER_SOFA_HPP

#include <iosfwd>
#include <memory>

struct HrtfStore;

/* Loads a SOFA (Spatially Oriented Format for Acoustics) file as an HRTF
 * data set. Only available when built with MySOFA support. Throws
 * std::runtime_error on invalid or unsupported data.
 */
auto LoadSofaHrtf(std::istream &stream) -> std::unique_ptr<HrtfStore>;

#endif /* CORE_HRTF_LOADER_SOFA_HPP */
