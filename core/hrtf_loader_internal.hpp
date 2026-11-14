#ifndef CORE_HRTF_LOADER_INTERNAL_HPP
#define CORE_HRTF_LOADER_INTERNAL_HPP

#include <cstdint>
#include <memory>
#include <span>

#include "mixer/hrtfdefs.h"

/* NOTE: HrtfStore must be defined before including this header, either by
 * importing the hrtfstore module or including core/hrtfstore.hpp.
 */

/* Shared between the .mhr loader (hrtf_loader.cpp) and the runtime SOFA
 * loader (hrtf_loader_sofa.cpp).
 */

auto CreateHrtfStore(std::uint32_t rate, std::uint8_t irSize,
    std::span<HrtfStore::Field const> fields, std::span<HrtfStore::Elevation const> elevs,
    std::span<HrirArray const> coeffs, std::span<u8x2 const> delays)
    -> std::unique_ptr<HrtfStore>;

void MirrorLeftHrirs(std::span<HrtfStore::Elevation const> elevs, std::span<HrirArray> coeffs,
    std::span<u8x2> delays);

#endif /* CORE_HRTF_LOADER_INTERNAL_HPP */
