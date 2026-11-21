#ifndef CORE_HRTF_LOADER_INTERNAL_HPP
#define CORE_HRTF_LOADER_INTERNAL_HPP

#include <cstdint>
#include <memory>
#include <span>

#include "mixer/hrtfdefs.h"

/* NOTE: HrtfStore must be defined before including this header, either by
 * importing the hrtfstore module or including core/hrtfstore.hpp.
 *
 * Helpers shared between the native .mhr loader (hrtf_loader.cpp) and the
 * runtime SOFA loader (hrtf_loader_sofa.cpp).
 */

struct HrtfStore;


namespace al::hrtf_loader {

/* Data set limits must be the same as or more flexible than those defined in
 * the makemhr utility.
 */
constexpr auto MinFdCount = 1u;
constexpr auto MaxFdCount = 16u;

constexpr auto MinFdDistance = 50u;
constexpr auto MaxFdDistance = 2500u;

constexpr auto MinEvCount = 5u;
constexpr auto MaxEvCount = 181u;

constexpr auto MinAzCount = 1u;
constexpr auto MaxAzCount = 255u;


auto CreateHrtfStore(std::uint32_t rate, std::uint8_t irSize,
    std::span<HrtfStore::Field const> fields, std::span<HrtfStore::Elevation const> elevs,
    std::span<HrirArray const> coeffs, std::span<u8x2 const> delays)
    -> std::unique_ptr<HrtfStore>;

void MirrorLeftHrirs(std::span<HrtfStore::Elevation const> elevs, std::span<HrirArray> coeffs,
    std::span<u8x2> delays);

} // namespace al::hrtf_loader

#endif /* CORE_HRTF_LOADER_INTERNAL_HPP */
