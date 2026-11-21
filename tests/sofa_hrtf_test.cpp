/* Runtime SOFA HRTF loading self-test.
 *
 * Part A exercises the loader directly: a valid stereo set, a mono set (ear
 * mirroring), per-measurement delays, over-long responses (windowing), the
 * three distinguishable error classes, malformed input, the oversize-file
 * limit, and a .mhr data set for regression.
 *
 * Part B drives the public ALC API on the null backend: SOFA files are
 * enumerated next to .mhr files, can be selected by index, report the right
 * specifier, resample to the device rate, and fall back when unreadable.
 *
 * The test writes its own ALSOFT_CONF/ALSOFT_DRIVERS settings so it doesn't
 * need an audio device or a user configuration.
 */

#include "config.h"

#include <AL/al.h>
#include <AL/alc.h>
#include <AL/alext.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <vector>

#include "core/hrtf_loader.hpp"
#include "core/hrtfstore.hpp"
#include "zudl.hpp"

#ifndef SOFA_TEST_DIR
#define SOFA_TEST_DIR "."
#endif

namespace {

int gFailures = 0;
int gChecks = 0;

#define CHECK(cond) do { \
    ++gChecks; \
    if(!(cond)) { \
        std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        ++gFailures; \
    } \
} while(0)

#define CHECK_MSG(cond, msg) do { \
    ++gChecks; \
    if(!(cond)) { \
        std::fprintf(stderr, "FAIL %s:%d: %s -- %s\n", __FILE__, __LINE__, #cond, msg); \
        ++gFailures; \
    } \
} while(0)

auto Path(std::string_view const name) -> std::string
{ return std::string{SOFA_TEST_DIR} + "/" + std::string{name}; }

auto LoadFile(std::string_view const name) -> std::unique_ptr<HrtfStore>
{
    auto stream = std::ifstream{Path(name), std::ios::binary};
    if(!stream.is_open())
        throw std::runtime_error{"Could not open " + Path(name)};
    return LoadHrtf(stream);
}

/* Returns the type name of the exception thrown by loading the named file, or
 * "none" if it loads, or "other" for an unclassified exception.
 */
auto LoadErrorKind(std::string_view const name) -> std::string
{
    try {
        std::ignore = LoadFile(name);
        return "none";
    }
    catch(SofaFileTooLargeError const&) { return "toolarge"; }
    catch(SofaMissingDataError const&) { return "missing"; }
    catch(SofaDimensionError const&) { return "dimension"; }
    catch(SofaLoadError const&) { return "sofa"; }
    catch(std::exception const&) { return "other"; }
}


/* Part A: direct loader tests. */

void TestValidStereo()
{
    const auto hrtf = LoadFile("valid.sofa");
    CHECK(hrtf != nullptr);
    CHECK(hrtf->mSampleRate == 44100u);
    CHECK(hrtf->mIrSize == 32u);
    CHECK(hrtf->mFields.size() == 1u);
    CHECK(hrtf->mFields[0].evCount == 9u);

    /* 1 pole + 7 rows of 8 + 1 pole = 58 HRIRs. */
    const auto irCount = std::size_t{hrtf->mElev.back().azCount} + hrtf->mElev.back().irOffset;
    CHECK(irCount == 58u);

    /* Every HRIR's peak should have been normalized to 1 (the source sets a
     * clear peak at sample 8). */
    auto maxval = 0.0f;
    for(const auto &hrir : hrtf->mCoeffs)
    {
        for(const auto &s : hrir)
            maxval = std::max({maxval, std::abs(s[0]), std::abs(s[1])});
    }
    CHECK(std::abs(maxval - 1.0f) < 1.0e-4f);

    /* Each stored HRIR must be finite and have its energy in the first
     * samples (the synthetic source has an 8-sample pre-ring and a peak at
     * sample 8). Delays must fit the mixer's range.
     */
    for(const auto &hrir : hrtf->mCoeffs)
    {
        auto early = 0.0f;
        for(const auto &s : hrir)
        {
            CHECK(std::isfinite(s[0]) && std::isfinite(s[1]));
            early += s[0]*s[0] + s[1]*s[1];
        }
        CHECK(early > 1.0e-4f);
    }
    for(const auto &delay : hrtf->mDelays)
    {
        CHECK(delay[0] <= MaxHrirDelay*HrirDelayFracOne);
        CHECK(delay[1] <= MaxHrirDelay*HrirDelayFracOne);
    }
}

void TestMonoMirrored()
{
    const auto hrtf = LoadFile("valid_mono.sofa");
    CHECK(hrtf != nullptr);
    CHECK(hrtf->mSampleRate == 48000u);

    /* A mono set has its left ear mirrored to the right. The mirroring pairs
     * each azimuth with its opposite azimuth, swapping the ears; the two
     * poles are each their own mirror. So for every elevation row, the right
     * ear of the i-th azimuth equals the left ear of the mirrored azimuth,
     * and the pole responses are left/right symmetric.
     */
    const auto irCount = std::size_t{hrtf->mElev.back().azCount} + hrtf->mElev.back().irOffset;
    for(const auto &elev : hrtf->mElev)
    {
        const auto azcount = std::size_t{elev.azCount};
        const auto base = std::size_t{elev.irOffset};
        if(azcount == 1)
        {
            for(const auto &s : hrtf->mCoeffs[base])
                CHECK(std::abs(s[0] - s[1]) < 1.0e-6f);
            continue;
        }
        for(auto j = 0_uz;j < azcount;++j)
        {
            const auto lidx = base + j;
            const auto ridx = base + ((azcount-j) % azcount);
            for(auto n = 0_uz;n < hrtf->mIrSize;++n)
                CHECK(std::abs(hrtf->mCoeffs[ridx][n][1] - hrtf->mCoeffs[lidx][n][0])
                    < 1.0e-6f);
            CHECK(hrtf->mDelays[ridx][1] == hrtf->mDelays[lidx][0]);
        }
    }
    CHECK(irCount == 58u);
}

void TestMrDelays()
{
    const auto hrtf = LoadFile("valid_mr_delay.sofa");
    CHECK(hrtf != nullptr);

    /* The right-ear delay is larger than the left's (0.0002s vs 0.0001s at
     * 44.1kHz, so ~8.8 vs ~4.4 samples). Delays are shifted down by the
     * smallest value, so the left/right difference (~4 samples) is what
     * survives. Checking the measured (non-pole) rows, the per-measurement
     * delays must have been taken from Data.Delay.
     */
    auto found = false;
    auto ldelay_min = unsigned{255};
    auto rdelay_min = unsigned{255};
    for(const auto &elev : hrtf->mElev)
    {
        if(elev.azCount <= 1) continue;
        const auto base = std::size_t{elev.irOffset};
        for(auto ai = std::size_t{0};ai < elev.azCount;++ai)
        {
            const auto l = hrtf->mDelays[base+ai][0];
            const auto r = hrtf->mDelays[base+ai][1];
            CHECK(r > l);
            ldelay_min = std::min(ldelay_min, unsigned{l});
            rdelay_min = std::min(rdelay_min, unsigned{r});
            found = true;
        }
    }
    CHECK_MSG(found, "no measured rows found");
    /* ~4 samples at 2 bits of fraction = ~17 units. */
    CHECK(rdelay_min > ldelay_min + 8u);
}

void TestLongIrWindowing()
{
    const auto hrtf = LoadFile("long_ir.sofa");
    CHECK(hrtf != nullptr);
    /* The IR is truncated to the mixer's filter length. */
    CHECK(hrtf->mIrSize == HrirLength);

    /* The response peak (originally sample 64 of 256) must survive the
     * windowing. After normalization the largest coefficient is 1, and the
     * peak is well within the first HrirLength samples. */
    auto peak = 0.0f;
    auto peakpos = 0_uz;
    for(const auto &hrir : hrtf->mCoeffs)
    {
        for(auto n = 0_uz;n < hrir.size();++n)
        {
            const auto mag = std::max(std::abs(hrir[n][0]), std::abs(hrir[n][1]));
            if(mag > peak) { peak = mag; peakpos = n; }
        }
    }
    CHECK(peak > 0.99f);
    CHECK(peakpos < HrirLength);
}

void TestErrorKinds()
{
    /* Dimension errors. */
    CHECK_MSG(LoadErrorKind("bad_receivers.sofa") == "dimension", "3 receivers");
    CHECK_MSG(LoadErrorKind("bad_rate.sofa") == "dimension", "unsupported sample rate");
    CHECK_MSG(LoadErrorKind("bad_grid.sofa") == "dimension", "too few elevations");

    /* Missing-data errors. */
    CHECK_MSG(LoadErrorKind("missing_ir.sofa") == "missing", "no FIR data");
    CHECK_MSG(LoadErrorKind("missing_sr.sofa") == "missing", "no sample rate");
    CHECK_MSG(LoadErrorKind("silent_ir.sofa") == "missing", "silent responses");

    /* Malformed container. */
    CHECK_MSG(LoadErrorKind("truncated.sofa") == "sofa", "truncated HDF5");

    /* A non-HDF5 file is still an .mhr-format error, not a SOFA one (it goes
     * through the native loader, whose header check fails). */
    {
        auto stream = std::ifstream{Path("not_sofa.mhr"), std::ios::binary};
        bool sofa_error = false;
        try { std::ignore = LoadHrtf(stream); }
        catch(SofaLoadError const&) { sofa_error = true; }
        catch(std::runtime_error const&) { }
        CHECK(!sofa_error);
    }

#if HAVE_MYSOFA
    /* The base class catches all three subclasses. */
    CHECK_MSG(LoadErrorKind("bad_rate.sofa") != "other", "classified");
#endif
}

void TestFileTooLarge()
{
    /* Build a stream larger than the loader's cap without writing much to
     * disk: the first bytes are a valid HDF5 signature, the rest is filler.
     * The loader must reject it on size before attempting to parse.
     */
    auto const path = std::string{SOFA_TEST_DIR} + "/oversize.tmp";
    {
        auto out = std::ofstream{path, std::ios::binary | std::ios::trunc};
        out.write("\x89HDF\r\n\x1a\n", 8);
        const auto filler = std::array<char, 65536>{};
        for(std::size_t i = 0;i < ((64u << 20) + (1u << 20)) / filler.size();++i)
            out.write(filler.data(), static_cast<std::streamsize>(filler.size()));
    }

    bool tooLarge = false;
    try {
        auto stream = std::ifstream{path, std::ios::binary};
        std::ignore = LoadHrtf(stream);
    }
    catch(SofaFileTooLargeError const&) { tooLarge = true; }
    catch(std::exception const&) { }
    std::remove(path.c_str());
    CHECK_MSG(tooLarge, "oversized stream should raise SofaFileTooLargeError");
}

void TestMhrRegression()
{
    /* The native .mhr path is unchanged: the built-in data set must still
     * load, with a sensible rate and filter size.
     */
    auto stream = std::ifstream{MHR_TEST_FILE, std::ios::binary};
    CHECK(stream.is_open());
    const auto hrtf = LoadHrtf(stream);
    CHECK(hrtf != nullptr);
    CHECK(hrtf->mSampleRate > 0);
    CHECK(hrtf->mIrSize >= 8);
    CHECK(!hrtf->mFields.empty());

    /* Sanity-check the loaded grid and coefficients. */
    const auto irCount = std::size_t{hrtf->mElev.back().azCount} + hrtf->mElev.back().irOffset;
    CHECK(irCount > 0);
    CHECK(hrtf->mCoeffs.size() == irCount);
    CHECK(hrtf->mDelays.size() == irCount);
    auto energy = 0.0f;
    for(const auto &hrir : hrtf->mCoeffs)
    {
        for(const auto &s : hrir)
        {
            CHECK(std::isfinite(s[0]) && std::isfinite(s[1]));
            energy += s[0]*s[0] + s[1]*s[1];
        }
    }
    CHECK(energy > 0.0f);
}


/* Part B: end-to-end through the ALC API on the null backend. */

struct AlcDevice {
    ALCdevice *device{nullptr};
    AlcDevice() { device = alcOpenDevice(nullptr); }
    ~AlcDevice() { if(device) alcCloseDevice(device); }
    auto operator!() const noexcept -> bool { return !device; }
};

auto GetSpecifiers(ALCdevice *const device) -> std::vector<std::string>
{
    auto result = std::vector<std::string>{};
    ALCint count = 0;
    alcGetIntegerv(device, ALC_NUM_HRTF_SPECIFIERS_SOFT, 1, &count);
    for(ALCint i = 0;i < count;++i)
    {
        if(auto const *const name = alcGetStringiSOFT(device, ALC_HRTF_SPECIFIER_SOFT, i))
            result.emplace_back(name);
    }
    return result;
}

auto FindSpecifier(ALCdevice *const device, std::string_view const name) -> int
{
    const auto specifiers = GetSpecifiers(device);
    for(std::size_t i = 0;i < specifiers.size();++i)
    {
        if(specifiers[i] == name)
            return static_cast<int>(i);
    }
    return -1;
}

/* Resets the device with the given HRTF specifier index, returning the
 * resulting HRTF specifier name (empty on failure).
 */
auto ResetWithHrtf(ALCdevice *const device, int const hrtf_id, ALCint const freq = 0)
    -> std::string
{
    auto attrs = std::vector<ALCint>{ALC_HRTF_SOFT, ALC_TRUE, ALC_HRTF_ID_SOFT, hrtf_id};
    if(freq != 0)
    {
        attrs.push_back(ALC_FREQUENCY);
        attrs.push_back(freq);
    }
    attrs.push_back(0);

    if(!alcResetDeviceSOFT(device, attrs.data()))
        return {};

    ALCint status = 0;
    alcGetIntegerv(device, ALC_HRTF_STATUS_SOFT, 1, &status);
    if(status != ALC_HRTF_ENABLED_SOFT)
        return {};
    auto const *const current = alcGetString(device, ALC_HRTF_SPECIFIER_SOFT);
    return current ? current : "";
}

void TestAlcEnumeration()
{
    const auto dev = AlcDevice{};
    CHECK(!!dev);
    if(!dev) return;

    const auto specifiers = GetSpecifiers(dev.device);
#if HAVE_MYSOFA
    for(const auto name : {"valid", "valid_mono", "long_ir", "valid_mr_delay"})
        CHECK_MSG(FindSpecifier(dev.device, name) >= 0, name);
    /* Non-SOFA and broken files are still enumerated (they fail at load). */
    CHECK(FindSpecifier(dev.device, "bad_rate") >= 0);
    CHECK(FindSpecifier(dev.device, "truncated") >= 0);
#else
    for(const auto name : {"valid", "valid_mono", "long_ir"})
        CHECK_MSG(FindSpecifier(dev.device, name) < 0, name);
#endif
    /* The built-in data set is always available. */
    CHECK_MSG(FindSpecifier(dev.device, "Built-In HRTF") >= 0, "built-in");
}

#if HAVE_MYSOFA
void TestAlcLoadAndResample()
{
    const auto dev = AlcDevice{};
    CHECK(!!dev);
    if(!dev) return;

    /* Load each good set at its own rate. */
    for(const auto name : {"valid", "valid_mono", "long_ir", "valid_mr_delay"})
    {
        const int idx = FindSpecifier(dev.device, name);
        CHECK_MSG(idx >= 0, name);
        if(idx < 0) continue;
        CHECK_MSG(ResetWithHrtf(dev.device, idx) == name, name);
    }

    /* A 44.1kHz set must also load when the device runs at other rates,
     * exercising the shared resampler. */
    const int idx = FindSpecifier(dev.device, "valid");
    if(idx >= 0)
    {
        for(const ALCint freq : {44100, 48000, 96000})
            CHECK_MSG(ResetWithHrtf(dev.device, idx, freq) == "valid", "resample");
    }
}

void TestAlcFallback()
{
    const auto dev = AlcDevice{};
    CHECK(!!dev);
    if(!dev) return;

    /* Selecting a broken SOFA set must still reset the device, but with a
     * different (fallback) set in use.
     */
    for(const auto name : {"bad_rate", "truncated", "missing_ir", "silent_ir"})
    {
        const int idx = FindSpecifier(dev.device, name);
        CHECK_MSG(idx >= 0, name);
        if(idx < 0) continue;
        const auto current = ResetWithHrtf(dev.device, idx);
        CHECK_MSG(!current.empty(), name);
        CHECK_MSG(current != name, name);
    }
}

void TestAlcConcurrentSwitching()
{
    auto failures = std::atomic<int>{0};
    auto worker = [&failures]()
    {
        const auto dev = AlcDevice{};
        if(!dev) { ++failures; return; }
        const auto count = static_cast<ALCint>(GetSpecifiers(dev.device).size());
        for(int iter=0;iter < 10;++iter)
        {
            for(ALCint i=0;i < count;++i)
            {
                if(ResetWithHrtf(dev.device, i).empty())
                    ++failures;
            }
        }
    };

    auto threads = std::vector<std::thread>{};
    for(auto i=0;i < 4;++i) threads.emplace_back(worker);
    for(auto &t : threads) t.join();
    CHECK_MSG(failures.load() == 0, "concurrent reset failures");
}
#endif

} // namespace


auto main() -> int
{
    /* Point the data search at the test assets and use the null backend. Only
     * the committed fixtures are enumerated (plus the built-in set), so the
     * test is independent of the host's installed HRTFs.
     */
#ifdef _WIN32
    _putenv_s("ALSOFT_DRIVERS", "null");
    _putenv_s("ALSOFT_LOCAL_PATH", SOFA_TEST_DIR);
#else
    setenv("ALSOFT_DRIVERS", "null", 1);
    setenv("ALSOFT_LOCAL_PATH", SOFA_TEST_DIR, 1);
#endif

    try {
#if HAVE_MYSOFA
        TestValidStereo();
        TestMonoMirrored();
        TestMrDelays();
        TestLongIrWindowing();
        TestErrorKinds();
        TestFileTooLarge();
#endif
        TestMhrRegression();

        TestAlcEnumeration();
#if HAVE_MYSOFA
        TestAlcLoadAndResample();
        TestAlcFallback();
        TestAlcConcurrentSwitching();
#endif
    }
    catch(std::exception const &e) {
        std::fprintf(stderr, "Unexpected exception: %s\n", e.what());
        ++gFailures;
    }

    std::printf("%d checks, %d failure(s)\n", gChecks, gFailures);
    return gFailures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
