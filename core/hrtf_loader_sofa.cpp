#include "config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <istream>
#include <iterator>
#include <limits>
#include <memory>
#include <ranges>
#include <span>
#include <string_view>
#include <vector>

#include "hrtf_loader_sofa.hpp"

#include "hrtf_loader.hpp"
#include "mysofa.h"
#include "opthelpers.h"
#include "zudl.hpp"

#if HAVE_CXXMODULES
import common.std;

import fmtlib;
import format;
import gsl;
import hrtfstore;
import logging;
#include "gsl_macros.hpp"
#else

#include <stdexcept>

#include "alformat.hpp"
#include "core/hrtfstore.hpp"
#include "gsl/gsl"
#include "logging.h"
#endif

#include "hrtf_loader_internal.hpp"


namespace {

using namespace std::string_view_literals;

using namespace al::hrtf_loader;

/* Bounds to keep memory use and load time in check when reading arbitrary
 * SOFA files at runtime. The file size cap bounds both the read buffer and
 * the parse time; the remaining caps bound the intermediate and final data
 * structures. They are intentionally more permissive than the .mhr limits.
 */
constexpr auto MaxSofaFileSize = 64_uz * 1024 * 1024;
constexpr auto MaxSofaMeasurements = 65536u;
constexpr auto MaxSofaIrPoints = 4096u;
constexpr auto MaxSofaIrSamples = 1_uz << 24; /* Total M*R*N samples. */
constexpr auto MaxIrTotal = 32768_uz; /* Total HRIRs across all fields. */

constexpr auto MinSofaSampleRate = 4000.0f;
constexpr auto MaxSofaSampleRate = 384000.0f;

/* Samples kept ahead of the IR peak when windowing over-long responses. */
constexpr auto IrPeakPrepend = 8_uz;

/* Minimum acceptable HRIR peak amplitude, to catch null/silent responses. */
constexpr auto MinIrPeak = 1.0e-6f;

/* Maximum acceptable HRIR peak amplitude. SOFA HRIRs are conventionally
 * normalized to unit range, so anything beyond this indicates garbage data.
 */
constexpr auto MaxIrPeak = 100.0f;


struct MySofaDeleter { void operator()(MYSOFA_HRTF *sofa) const { mysofa_free(sofa); } };
using MySofaHrtfPtr = std::unique_ptr<MYSOFA_HRTF, MySofaDeleter>;

[[nodiscard]]
auto SofaErrorStr(int const err) noexcept -> const char*
{
    switch(err)
    {
    case MYSOFA_OK: return "OK";
    case MYSOFA_INVALID_FORMAT: return "Invalid format";
    case MYSOFA_UNSUPPORTED_FORMAT: return "Unsupported format";
    case MYSOFA_INTERNAL_ERROR: return "Internal error";
    case MYSOFA_NO_MEMORY: return "Out of memory";
    case MYSOFA_READ_ERROR: return "Read error";
    }
    return "Unknown";
}

[[nodiscard]]
auto FindAttribute(MYSOFA_ATTRIBUTE const *attrs, std::string_view const name) noexcept
    -> gsl::czstring
{
    while(attrs)
    {
        if(attrs->name && name == attrs->name)
            return attrs->value;
        attrs = attrs->next;
    }
    return nullptr;
}


auto GetSampleRate(MYSOFA_HRTF const *const sofaHrtf) -> unsigned
{
    auto const *srate_attrs = sofaHrtf->DataSamplingRate.attributes;
    auto const *const srate_dim = FindAttribute(srate_attrs, "DIMENSION_LIST"sv);
    if(!srate_dim)
        throw SofaMissingDataError{"Missing sample rate dimensions"};
    if(srate_dim != "I"sv)
    {
        throw SofaDimensionError{al::format("Unsupported sample rate dimensions: {}",
            srate_dim)};
    }

    auto const *const srate_units = FindAttribute(srate_attrs, "Units"sv);
    if(!srate_units)
        throw SofaMissingDataError{"Missing sample rate unit type"};
    if(srate_units != "hertz"sv)
    {
        throw SofaDimensionError{al::format("Unsupported sample rate unit type: {}",
            srate_units)};
    }

    if(!sofaHrtf->DataSamplingRate.values || sofaHrtf->DataSamplingRate.elements < 1)
        throw SofaMissingDataError{"Missing sample rate value"};

    /* I dimensions guarantees 1 element, so just extract it. */
    auto const srate = sofaHrtf->DataSamplingRate.values[0];
    if(!std::isfinite(srate) || srate < MinSofaSampleRate || srate > MaxSofaSampleRate)
    {
        throw SofaDimensionError{al::format("Sample rate out of range: {:f} (expected {} to "
            "{}hz)", srate, MinSofaSampleRate, MaxSofaSampleRate)};
    }
    return gsl::narrow_cast<unsigned>(std::lround(srate));
}

enum class DelayType : std::uint8_t {
    None,
    I_R, /* [1][Channels] */
    M_R, /* [HRIRs][Channels] */
};

auto GetDelayType(MYSOFA_HRTF const *const sofaHrtf) -> DelayType
{
    auto const *delay_attrs = sofaHrtf->DataDelay.attributes;
    auto const *const delay_dim = FindAttribute(delay_attrs, "DIMENSION_LIST"sv);
    if(!delay_dim)
    {
        TRACE("Missing delay dimensions, assuming no delays");
        return DelayType::None;
    }
    if(delay_dim == "I,R"sv)
        return DelayType::I_R;
    if(delay_dim == "M,R"sv)
        return DelayType::M_R;
    throw SofaDimensionError{al::format("Unsupported delay dimensions: {}", delay_dim)};
}

void CheckIrData(MYSOFA_HRTF const *const sofaHrtf)
{
    auto const *ir_attrs = sofaHrtf->DataIR.attributes;
    auto const *const ir_dim = FindAttribute(ir_attrs, "DIMENSION_LIST"sv);
    if(!ir_dim)
        throw SofaMissingDataError{"Missing IR dimensions"};
    if(ir_dim != "M,R,N"sv)
        throw SofaDimensionError{al::format("Unsupported IR dimensions: {}", ir_dim)};
}


using double3 = std::array<double, 3>;

/* Produces a sorted array of unique elements from a particular axis of the
 * triplets array. The filters are used to focus on particular coordinates of
 * other axes as necessary. The epsilons are used to constrain the equality of
 * unique elements.
 */
auto GetUniquelySortedElems(std::span<double3 const> const aers, unsigned const axis,
    std::array<const double*, 3> const &filters, std::array<double, 3> const &epsilons)
    -> std::vector<double>
{
    auto elems = std::vector<double>{};
    for(const double3 &aer : aers)
    {
        const double elem{aer[axis]};

        auto j = 0u;
        for(j = 0;j < 3;j++)
        {
            if(filters[j] && std::abs(aer[j] - *filters[j]) > epsilons[j])
                break;
        }
        if(j < 3)
            continue;

        auto iter = elems.begin();
        for(;iter != elems.end();++iter)
        {
            const double delta{elem - *iter};
            if(delta > epsilons[axis]) continue;
            if(delta >= -epsilons[axis]) break;

            iter = elems.emplace(iter, elem);
            break;
        }
        if(iter == elems.end())
            elems.emplace_back(elem);
    }
    return elems;
}

/* Given a list of azimuths, this will produce the smallest step size that can
 * uniformly cover the list. Ideally this will be over half, but in degenerate
 * cases this can fall to a minimum of 5 (the lower limit).
 */
auto GetUniformAzimStep(double const epsilon, std::span<double const> const elems) -> double
{
    if(elems.size() < 5) return 0.0;

    /* Get the maximum count possible, given the first two elements. It would
     * be impossible to have more than this since the first element must be
     * included.
     */
    auto count = static_cast<unsigned>(std::ceil(360.0 / (elems[1]-elems[0])));
    count = std::min(count, 255u);

    for(;count >= 5;--count)
    {
        /* Given the stepping value for this number of elements, check each
         * multiple to ensure there's a matching element.
         */
        const double step{360.0 / count};
        bool good{true};
        size_t idx{1u};
        for(auto mult=1u;mult < count && good;++mult)
        {
            const double target{step*mult + elems[0]};
            while(idx < elems.size() && target-elems[idx] > epsilon)
                ++idx;
            good &= (idx < elems.size()) && !(std::abs(target-elems[idx++]) > epsilon);
        }
        if(good)
            return step;
    }
    return 0.0;
}

/* Given a list of elevations, this will produce the smallest step size that
 * can uniformly cover the list. Ideally this will be over half, but in
 * degenerate cases this can fall to a minimum of 5 (the lower limit).
 */
auto GetUniformElevStep(double const epsilon, std::vector<double> &elems) -> double
{
    if(elems.size() < 5) return 0.0;

    /* Reverse the elevations so it increments starting with -90 (flipped from
     * +90). This makes it easier to work out a proper stepping value.
     */
    std::ranges::reverse(elems);
    for(auto &v : elems) v *= -1.0;

    auto count = static_cast<unsigned>(std::ceil(180.0 / (elems[1]-elems[0])));
    count = std::min(count, 255u);

    double ret{0.0};
    for(;count >= 5;--count)
    {
        const double step{180.0 / count};
        bool good{true};
        size_t idx{1u};
        /* Elevations don't need to match all multiples if there's not enough
         * elements to check. Missing elevations can be synthesized.
         */
        for(auto mult=1u;mult <= count && idx < elems.size() && good;++mult)
        {
            const double target{step*mult + elems[0]};
            while(idx < elems.size() && target-elems[idx] > epsilon)
                ++idx;
            good &= !(idx < elems.size()) || !(std::abs(target-elems[idx++]) > epsilon);
        }
        if(good)
        {
            ret = step;
            break;
        }
    }
    /* Re-reverse the elevations to restore the correct order. */
    for(auto &v : elems) v *= -1.0;
    std::ranges::reverse(elems);

    return ret;
}


/* Per-field layout info. */
struct SofaField {
    double mDistance;
    unsigned mEvCount;
    unsigned mEvStart; /* First elevation row with measured data. */
    unsigned mEvEnd; /* One past the last elevation row with measured data. */
    std::vector<unsigned> mAzCounts;
};

/* Attempts to produce a compatible layout. Most data sets tend to be uniform
 * and have the same major axis as used by OpenAL Soft's HRTF model. Fields
 * are returned in ascending distance order.
 */
auto GetCompatibleLayout(std::span<double3 const> const aers) -> std::vector<SofaField>
{
    const auto radii = GetUniquelySortedElems(aers, 2, {}, {0.1, 0.1, 0.001});
    if(radii.empty())
        throw SofaMissingDataError{"No source positions"};
    if(radii.size() > MaxFdCount)
    {
        throw SofaDimensionError{al::format("Incompatible layout ({} radii, max {})",
            radii.size(), MaxFdCount)};
    }

    auto fds = std::vector<SofaField>{};
    fds.reserve(radii.size());

    for(const double dist : radii)
    {
        if(!std::isfinite(dist) || dist < 0.01 || dist > 100.0)
        {
            throw SofaDimensionError{al::format("Unsupported field distance: {:f} meters",
                dist)};
        }

        auto elevs = GetUniquelySortedElems(aers, 1, {nullptr, nullptr, &dist},
            {0.1, 0.1, 0.001});

        /* Remove elevations that don't have a valid set of azimuths. */
        std::erase_if(elevs, [&dist,&aers](const double ev) -> bool
        {
            const auto azims = GetUniquelySortedElems(aers, 0, {nullptr, &ev, &dist},
                {0.1, 0.1, 0.001});

            if(std::abs(ev) > 89.999)
                return azims.size() != 1;
            if(azims.empty() || !(std::abs(azims[0]) < 0.1))
                return true;
            return GetUniformAzimStep(0.1, azims) <= 0.0;
        });

        const double step{GetUniformElevStep(0.1, elevs)};
        if(step <= 0.0)
        {
            if(elevs.empty())
            {
                WARN("No usable elevations on field distance {:f}", dist);
            }
            else
                WARN("Non-uniform elevations on field distance {:f}", dist);
            continue;
        }

        const auto evCount = static_cast<unsigned>(std::round(180.0 / step)) + 1;
        if(evCount < MinEvCount || evCount > MaxEvCount)
        {
            WARN("Unsupported elevation count {} on field distance {:f}", evCount, dist);
            continue;
        }

        SofaField field{};
        field.mDistance = dist;
        field.mEvCount = evCount;
        field.mEvStart = evCount;
        field.mEvEnd = 0u;
        field.mAzCounts.resize(evCount, 0u);

        for(auto ei=0u;ei < evCount;++ei)
        {
            const double ev{-90.0 + ei*step};
            const auto azims = GetUniquelySortedElems(aers, 0, {nullptr, &ev, &dist},
                {0.1, 0.1, 0.001});
            if(azims.empty())
                continue;

            if(std::abs(ev) > 89.999)
            {
                /* Poles must have a singular azimuth. */
                if(azims.size() != 1)
                {
                    throw SofaDimensionError{al::format(
                        "Non-singular pole on field distance {:f}", dist)};
                }
                field.mAzCounts[ei] = 1;
            }
            else
            {
                const double azstep{GetUniformAzimStep(0.1, azims)};
                if(azstep <= 0.0)
                {
                    throw SofaDimensionError{al::format(
                        "Non-uniform azimuths on elevation {:f}, field distance {:f}", ev, dist)};
                }
                const auto azCount = static_cast<unsigned>(std::round(360.0 / azstep));
                if(azCount < MinAzCount || azCount > MaxAzCount)
                {
                    throw SofaDimensionError{al::format(
                        "Unsupported azimuth count {} on elevation {:f}, field distance {:f}",
                        azCount, ev, dist)};
                }
                field.mAzCounts[ei] = azCount;
            }

            field.mEvStart = std::min(field.mEvStart, ei);
            field.mEvEnd = ei + 1;
        }

        if(field.mEvStart >= field.mEvEnd)
        {
            WARN("No usable measurements on field distance {:f}", dist);
            continue;
        }

        /* Give the rows without measurements the azimuth count of the nearest
         * measured row, except the poles which get a singular azimuth. Their
         * responses are synthesized from the nearest measured row.
         */
        for(auto ei=0u;ei < evCount;++ei)
        {
            if(field.mAzCounts[ei] != 0)
                continue;
            if(ei == 0 || ei == evCount-1)
                field.mAzCounts[ei] = 1;
            else if(ei < field.mEvStart)
                field.mAzCounts[ei] = field.mAzCounts[field.mEvStart];
            else
                field.mAzCounts[ei] = field.mAzCounts[field.mEvEnd-1];
        }

        fds.emplace_back(std::move(field));
    }

    if(fds.empty())
        throw SofaDimensionError{"Incompatible layout (no usable fields)"};
    return fds;
}


/* Read the entire stream, bounding the total size to keep memory use and load
 * time in check.
 */
auto ReadSofaData(std::istream &stream) -> std::vector<char>
{
    auto data = std::vector<char>{};
    auto chunk = std::array<char, 16384>{};
    while(stream)
    {
        stream.read(chunk.data(), chunk.size());
        const auto got = gsl::narrow_cast<std::size_t>(stream.gcount());
        if(data.size()+got > MaxSofaFileSize)
        {
            throw SofaFileTooLargeError{al::format("SOFA data is too large (over {} bytes)",
                MaxSofaFileSize)};
        }
        data.insert(data.end(), chunk.begin(), chunk.begin()+static_cast<std::ptrdiff_t>(got));
    }
    if(data.empty())
        throw SofaMissingDataError{"No SOFA data"};
    return data;
}

} // namespace


auto LoadSofaHrtf(std::istream &stream) -> std::unique_ptr<HrtfStore>
{
    const auto data = ReadSofaData(stream);

    auto err = int{};
    auto sofaHrtf = MySofaHrtfPtr{mysofa_load_data(data.data(), data.size(), &err)};
    if(!sofaHrtf)
    {
        throw SofaLoadError{al::format("Could not parse SOFA data: {} ({})", SofaErrorStr(err),
            err)};
    }
    auto *sofa = sofaHrtf.get();

    /* NOTE: Some valid SOFA files fail this check, so only warn. */
    if(const auto chk = mysofa_check(sofa); chk != MYSOFA_OK)
        WARN("Supposedly malformed SOFA data: {} ({})", SofaErrorStr(chk), chk);

    /* Make sure the basic dimensions are sane before touching the data. */
    if(sofa->E != 1)
    {
        throw SofaDimensionError{al::format("{} emitters not supported (expected 1)",
            sofa->E)};
    }
    if(sofa->R < 1 || sofa->R > 2)
    {
        throw SofaDimensionError{al::format("{} receivers not supported (expected 1 or 2)",
            sofa->R)};
    }
    if(sofa->C != 3)
    {
        throw SofaDimensionError{al::format("{} coordinate dimensions not supported "
            "(expected 3)", sofa->C)};
    }
    if(sofa->I != 1)
    {
        throw SofaDimensionError{al::format("{} singleton dimensions not supported "
            "(expected 1)", sofa->I)};
    }
    if(sofa->M < 1 || sofa->M > MaxSofaMeasurements)
    {
        throw SofaDimensionError{al::format("{} measurements not supported (1 to {})",
            sofa->M, MaxSofaMeasurements)};
    }
    if(sofa->N < MinIrLength || sofa->N > MaxSofaIrPoints)
    {
        throw SofaDimensionError{al::format("Unsupported HRIR size, N={} ({} to {})", sofa->N,
            MinIrLength, MaxSofaIrPoints)};
    }
    if(std::uint64_t{sofa->M}*sofa->R*sofa->N > MaxSofaIrSamples)
    {
        throw SofaDimensionError{al::format("Too much IR data ({}x{}x{} samples)", sofa->M,
            sofa->R, sofa->N)};
    }
    if(!sofa->DataIR.values || !sofa->SourcePosition.values)
        throw SofaMissingDataError{"Missing IR or source position data"};

    /* Make sure the data arrays are actually as large as the dimensions
     * claim, so nothing is read out of bounds.
     */
    if(std::uint64_t{sofa->DataIR.elements} < std::uint64_t{sofa->M}*sofa->R*sofa->N)
        throw SofaMissingDataError{"IR data is smaller than M*R*N"};
    if(std::uint64_t{sofa->SourcePosition.elements} < std::uint64_t{sofa->M}*sofa->C)
        throw SofaMissingDataError{"Source position data is smaller than M*C"};

    const auto rate = GetSampleRate(sofa);
    if(rate > MaxHrtfSampleRate)
    {
        throw SofaDimensionError{al::format("Sample rate is too large (max: {}hz)",
            MaxHrtfSampleRate)};
    }

    const auto delayType = GetDelayType(sofa);
    if(delayType == DelayType::I_R)
    {
        if(!sofa->DataDelay.values
            || std::uint64_t{sofa->DataDelay.elements} < std::uint64_t{sofa->I}*sofa->R)
            throw SofaMissingDataError{"Missing or undersized delay data"};
    }
    else if(delayType == DelayType::M_R)
    {
        if(!sofa->DataDelay.values
            || std::uint64_t{sofa->DataDelay.elements} < std::uint64_t{sofa->M}*sofa->R)
            throw SofaMissingDataError{"Missing or undersized delay data"};
    }

    CheckIrData(sofa);

    /* Convert the source positions to spherical coordinates, validating them
     * in the process.
     */
    mysofa_tocartesian(sofa);

    auto aers = std::vector(std::size_t{sofa->M}, double3{});
    {
        const auto srcPosValues = std::span{sofa->SourcePosition.values,
            std::size_t{sofa->M}*3};
        for(const auto si : std::views::iota(0_uz, std::size_t{sofa->M}))
        {
            auto pos = std::array{srcPosValues[3*si], srcPosValues[3*si + 1],
                srcPosValues[3*si + 2]};
            if(!std::isfinite(pos[0]) || !std::isfinite(pos[1]) || !std::isfinite(pos[2]))
                throw SofaMissingDataError{al::format("Non-finite source position [{}]", si)};
            mysofa_c2s(pos.data());
            aers[si] = {pos[0], pos[1], pos[2]};
        }
    }

    const auto fds = GetCompatibleLayout(aers);

    /* Set up the HRTF grid. Fields are stored farthest-first, so iterate the
     * layout (which is nearest-first) in reverse.
     */
    const auto channels = std::size_t{sofa->R};
    auto fields = std::vector<HrtfStore::Field>{};
    auto elevs = std::vector<HrtfStore::Elevation>{};
    auto fieldbases = std::vector<std::size_t>{}; /* Per-field base elevation index. */
    auto irTotal = 0_uz;
    for(const auto &field : fds | std::views::reverse)
    {
        fields.emplace_back(HrtfStore::Field{.distance=gsl::narrow_cast<float>(field.mDistance),
            .evCount=gsl::narrow<std::uint8_t>(field.mEvCount)});
        fieldbases.emplace_back(elevs.size());
        for(const auto azCount : field.mAzCounts)
        {
            if(irTotal+azCount > MaxIrTotal)
                throw SofaDimensionError{"Too many HRIRs in data set"};
            elevs.emplace_back(HrtfStore::Elevation{.azCount=gsl::narrow<std::uint16_t>(azCount),
                .irOffset=gsl::narrow<std::uint16_t>(irTotal)});
            irTotal += azCount;
        }
    }
    /* fieldbases was built in farthest-first order matching fields/elevs; flip
     * it to layout (nearest-first) order for filling.
     */
    std::ranges::reverse(fieldbases);

    auto coeffs = std::vector(irTotal, HrirArray{});
    auto delays = std::vector(irTotal, u8x2{});
    auto filled = std::vector(irTotal, false);
    /* Per-cell delays in (fractional) samples, converted to the final
     * fixed-point format after any needed normalization.
     */
    auto sampleDelays = std::vector(irTotal, std::array{0.0, 0.0});

    const auto irSize = std::min(sofa->N, HrirLength);

    /* Load each measurement into its matching grid cell. */
    const auto irValues = std::span{sofa->DataIR.values, std::size_t{sofa->M}*sofa->R*sofa->N};
    auto delayValues = std::span<const float>{};
    if(delayType == DelayType::I_R)
        delayValues = std::span{sofa->DataDelay.values, std::size_t{sofa->I}*sofa->R};
    else if(delayType == DelayType::M_R)
        delayValues = std::span{sofa->DataDelay.values, std::size_t{sofa->M}*sofa->R};

    for(const auto si : std::views::iota(0_uz, std::size_t{sofa->M}))
    {
        auto aer = aers[si];
        if(std::abs(aer[1]) >= 89.999)
            aer[0] = 0.0;
        else
            aer[0] = std::fmod(360.0 - aer[0], 360.0);

        const auto fielditer = std::ranges::find_if(fds, [&aer](SofaField const &fld)
        { return std::abs(aer[2] - fld.mDistance) < 0.001; });
        if(fielditer == fds.end())
            continue;
        const auto &field = *fielditer;

        const auto evscale = 180.0 / static_cast<double>(field.mEvCount-1);
        const auto ef = (90.0 + aer[1]) / evscale;
        const auto ei = static_cast<long>(std::lround(ef));
        if(ei < 0 || ei >= static_cast<long>(field.mEvCount)
            || std::abs(ef - static_cast<double>(ei))*evscale >= 0.1)
            continue;

        const auto azcount = field.mAzCounts[ei];
        const auto azscale = 360.0 / static_cast<double>(azcount);
        const auto af = aer[0] / azscale;
        const auto ai = static_cast<long>(std::lround(af));
        if(std::abs(af - static_cast<double>(ai))*azscale >= 0.1)
            continue;

        const auto fi = gsl::narrow_cast<std::size_t>(std::distance(fds.begin(), fielditer));
        const auto ei_global = fieldbases[fi] + static_cast<unsigned>(ei);
        const auto idx = std::size_t{elevs[ei_global].irOffset}
            + static_cast<std::size_t>(ai % static_cast<long>(azcount));
        if(filled[idx])
        {
            throw SofaDimensionError{al::format(
                "Multiple measurements near [ a={:f}, e={:f}, r={:f} ]", aer[0], aer[1], aer[2])};
        }
        filled[idx] = true;

        for(const auto ti : std::views::iota(0_uz, channels))
        {
            const auto ir = irValues.subspan((si*sofa->R + ti)*sofa->N, sofa->N);

            /* Validate the response and locate its peak for windowing. */
            auto peak = 0.0f;
            auto peakpos = 0_uz;
            for(const auto n : std::views::iota(0_uz, ir.size()))
            {
                const auto val = ir[n];
                if(!std::isfinite(val))
                {
                    throw SofaMissingDataError{al::format(
                        "Non-finite IR sample in measurement {}", si)};
                }
                if(const auto mag = std::abs(val); mag > peak)
                {
                    peak = mag;
                    peakpos = n;
                }
            }
            if(!(peak > MinIrPeak))
                throw SofaMissingDataError{al::format("Silent IR in measurement {}", si)};
            if(!(peak < MaxIrPeak))
            {
                throw SofaDimensionError{al::format("Abnormal IR amplitude in measurement {}: "
                    "{}", si, peak)};
            }

            /* Window over-long responses around the response onset, keeping
             * the skipped samples as extra delay.
             */
            auto onset = 0_uz;
            if(sofa->N > HrirLength)
                onset = (peakpos > IrPeakPrepend) ? peakpos - IrPeakPrepend : 0_uz;

            for(const auto n : std::views::iota(0_uz, std::size_t{irSize}))
                coeffs[idx][n][ti] = (onset+n < ir.size()) ? ir[onset+n] : 0.0f;

            auto delay = static_cast<double>(onset);
            if(delayType == DelayType::I_R)
                delay += gsl::narrow_cast<double>(delayValues[ti]) * rate;
            else if(delayType == DelayType::M_R)
                delay += gsl::narrow_cast<double>(delayValues[si*sofa->R + ti]) * rate;
            if(!std::isfinite(delay))
                throw SofaMissingDataError{al::format("Non-finite delay in measurement {}", si)};
            sampleDelays[idx][ti] = delay;
        }
        if(channels == 1)
            sampleDelays[idx][1] = sampleDelays[idx][0];
    }

    sofaHrtf.reset();

    /* Make sure every measured elevation row is fully populated. */
    for(const auto fi : std::views::iota(0_uz, fds.size()))
    {
        const auto &field = fds[fi];
        for(auto ei=field.mEvStart;ei < field.mEvEnd;++ei)
        {
            const auto base = std::size_t{elevs[fieldbases[fi]+ei].irOffset};
            for(const auto ai : std::views::iota(0_uz, std::size_t{field.mAzCounts[ei]}))
            {
                if(!filled[base+ai])
                {
                    throw SofaMissingDataError{al::format("Missing measurement for source "
                        "[{}, {}, {}]", fi, ei, ai)};
                }
            }
        }
    }

    /* Synthesize the elevation rows without measurements by copying the
     * nearest measured row's responses (nearest azimuth wins).
     */
    for(const auto fi : std::views::iota(0_uz, fds.size()))
    {
        const auto &field = fds[fi];
        const auto fb = fieldbases[fi];
        for(const auto ei : std::views::iota(0u, field.mEvCount))
        {
            if(ei >= field.mEvStart && ei < field.mEvEnd)
                continue;
            const auto src_ei = std::clamp(ei, field.mEvStart, field.mEvEnd-1);
            const auto src_azcount = std::size_t{elevs[fb+src_ei].azCount};
            const auto src_base = std::size_t{elevs[fb+src_ei].irOffset};
            const auto dst_azcount = std::size_t{elevs[fb+ei].azCount};
            const auto dst_base = std::size_t{elevs[fb+ei].irOffset};
            for(const auto ai : std::views::iota(0_uz, dst_azcount))
            {
                const auto az = 360.0*gsl::narrow_cast<double>(ai)
                    / gsl::narrow_cast<double>(dst_azcount);
                const auto src_ai = std::min(src_azcount-1,
                    static_cast<std::size_t>(std::lround(
                        az*gsl::narrow_cast<double>(src_azcount)/360.0)) % src_azcount);
                coeffs[dst_base+ai] = coeffs[src_base+src_ai];
                sampleDelays[dst_base+ai] = sampleDelays[src_base+src_ai];
            }
        }
    }

    /* Normalize the responses so the largest coefficient is 1. */
    auto maxval = 0.0f;
    for(const auto &hrir : coeffs)
    {
        for(const auto &sample : hrir)
            maxval = std::max({maxval, std::abs(sample[0]), std::abs(sample[1])});
    }
    if(maxval > 0.0f)
    {
        const auto scale = 1.0f / maxval;
        for(auto &hrir : coeffs)
        {
            for(auto &sample : hrir)
            {
                sample[0] *= scale;
                sample[1] *= scale;
            }
        }
    }

    /* Validate the delays and shift them down by the smallest (a constant
     * offset for the whole set is inaudible; the mixer only uses the relative
     * interaural delay, and this is what makemhr does too). If the spread
     * still exceeds the max, scale them down to fit (essentially shrinking the
     * head radius; not ideal, but better than rejecting the data set).
     */
    auto delay_min = std::numeric_limits<double>::max();
    auto delay_max = 0.0;
    for(const auto &d : sampleDelays)
    {
        for(const auto v : d)
        {
            if(v < -0.001)
            {
                throw SofaDimensionError{al::format("Invalid negative delay: {:f}", v)};
            }
            delay_min = std::min(delay_min, v);
            delay_max = std::max(delay_max, v);
        }
    }
    if(delay_min > 0.0)
    {
        for(auto &d : sampleDelays)
        {
            std::ranges::transform(d, d.begin(), [delay_min](double const v)
            { return v - delay_min; });
        }
        delay_max -= delay_min;
    }
    if(delay_max > double{MaxHrirDelay})
    {
        WARN("Delay spread exceeds max ({:f} > {}), scaling down", delay_max, MaxHrirDelay);
        const auto delay_scale = double{MaxHrirDelay} / delay_max;
        for(auto &d : sampleDelays)
        {
            std::ranges::transform(d, d.begin(), [delay_scale](double const v)
            { return v * delay_scale; });
        }
    }
    for(const auto i : std::views::iota(0_uz, delays.size()))
    {
        for(const auto ti : std::views::iota(0_uz, 2_uz))
        {
            const auto d = std::round(sampleDelays[i][ti]*double{HrirDelayFracOne});
            delays[i][ti] = gsl::narrow<std::uint8_t>(std::clamp(d, 0.0,
                double{MaxHrirDelay*HrirDelayFracOne}));
        }
    }

    /* Mirror the left ear responses to the right ear for mono data sets. */
    if(channels == 1)
        MirrorLeftHrirs(elevs, coeffs, delays);

    TRACE("Loaded SOFA data set: {}hz, {}-sample filter, {} field(s), {} HRIRs", rate, irSize,
        fields.size(), irTotal);
    return CreateHrtfStore(rate, gsl::narrow<std::uint8_t>(irSize), fields, elevs, coeffs,
        delays);
}
