/* Standalone acceptance test for the ITU-R BS.1770-4 loudness meter.
 *
 * Uses no audio device (null backend): it synthesizes a 1 kHz 0 dBFS sine
 * and feeds it directly to the meter. The key assertion is that a
 * full-scale 1 kHz sine measures -3.01 LUFS. Also covers stereo and 5.1
 * multichannel, and distinguishes silence, invalid sample rate, invalid
 * channel count, and insufficient data. Runs in well under 30 seconds.
 */

#include "core/bs1770.h"

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>


namespace {

int gFailures{0};

void check(bool const ok, const char *const what)
{
    if(ok)
        std::printf("ok   : %s\n", what);
    else
    {
        ++gFailures;
        std::fprintf(stderr, "FAIL : %s\n", what);
    }
}

auto makeSine(double const freq, double const rate, unsigned const channels,
    double const seconds) -> std::vector<float>
{
    const auto frames = static_cast<std::size_t>(rate*seconds);
    auto out = std::vector<float>(frames*channels);
    for(auto i = std::size_t{0}; i < frames; ++i)
    {
        const auto v = static_cast<float>(
            std::sin(2.0*std::numbers::pi*freq*static_cast<double>(i)/rate));
        for(auto c = 0u; c < channels; ++c)
            out[i*channels + c] = v;
    }
    return out;
}

auto statusName(BS1770Status const s) -> const char*
{
    switch(s)
    {
    case BS1770Status::Ok: return "Ok";
    case BS1770Status::Silence: return "Silence";
    case BS1770Status::InvalidSampleRate: return "InvalidSampleRate";
    case BS1770Status::InvalidChannels: return "InvalidChannels";
    case BS1770Status::InsufficientData: return "InsufficientData";
    }
    return "?";
}

} // namespace

int main()
{
    constexpr auto rate = 48000.0;

    /* Mono 1 kHz 0 dBFS sine, 10 seconds: must read -3.01 LUFS. */
    const auto mono = makeSine(1000.0, rate, 1, 10.0);
    const auto monoRes = BS1770Meter::measure(rate, 1, mono);
    std::printf("mono 1kHz 0dBFS: status=%s LUFS=%.3f LKFS=%.3f\n",
        statusName(monoRes.mStatus), monoRes.mLUFS, monoRes.mLKFS);
    check(monoRes.mStatus == BS1770Status::Ok, "mono: status Ok");
    check(std::fabs(monoRes.mLUFS - (-3.01)) <= 0.1, "mono: -3.01 LUFS (+/-0.1)");
    check(monoRes.mLKFS == monoRes.mLUFS, "mono: LKFS equals LUFS");

    /* Stereo, same tone on both channels: +3.01 LU over mono. */
    const auto stereo = makeSine(1000.0, rate, 2, 10.0);
    const auto stereoRes = BS1770Meter::measure(rate, 2, stereo);
    std::printf("stereo 1kHz 0dBFS: status=%s LUFS=%.3f\n",
        statusName(stereoRes.mStatus), stereoRes.mLUFS);
    check(stereoRes.mStatus == BS1770Status::Ok, "stereo: status Ok");
    check(std::fabs(stereoRes.mLUFS - (monoRes.mLUFS + 3.01)) <= 0.05,
        "stereo: mono + 3.01 LU");

    /* 5.1 (FL, FR, FC, LFE, BL, BR), tone on all channels: LFE excluded,
     * surrounds weighted 1.41 -> sum of weights 1+1+1+0+1.41+1.41 = 5.82.
     */
    const auto surround = makeSine(1000.0, rate, 6, 10.0);
    const auto surroundRes = BS1770Meter::measure(rate, 6, surround);
    const auto expect51 = monoRes.mLUFS + 10.0*std::log10(5.82);
    std::printf("5.1 1kHz 0dBFS: status=%s LUFS=%.3f (expected %.3f)\n",
        statusName(surroundRes.mStatus), surroundRes.mLUFS, expect51);
    check(surroundRes.mStatus == BS1770Status::Ok, "5.1: status Ok");
    check(std::fabs(surroundRes.mLUFS - expect51) <= 0.05,
        "5.1: surround weights and LFE exclusion");

    /* 44.1 kHz sample rate: coefficients are regenerated per rate. */
    const auto mono44 = makeSine(1000.0, 44100.0, 1, 10.0);
    const auto mono44Res = BS1770Meter::measure(44100.0, 1, mono44);
    std::printf("mono 44.1kHz: status=%s LUFS=%.3f\n",
        statusName(mono44Res.mStatus), mono44Res.mLUFS);
    check(mono44Res.mStatus == BS1770Status::Ok, "44.1kHz: status Ok");
    check(std::fabs(mono44Res.mLUFS - (-3.01)) <= 0.1, "44.1kHz: -3.01 LUFS (+/-0.1)");

    /* Streaming in odd-sized chunks must match the one-shot result. */
    {
        auto meter = BS1770Meter{};
        check(meter.init(rate, 1) == BS1770Status::Ok, "stream: init Ok");
        constexpr auto chunk = std::size_t{997};
        for(auto pos = std::size_t{0}; pos < mono.size(); pos += chunk)
        {
            const auto len = std::min(chunk, mono.size() - pos);
            meter.process(std::span<const float>{mono.data() + pos, len});
        }
        const auto streamRes = meter.integrated();
        std::printf("streamed mono: status=%s LUFS=%.3f\n",
            statusName(streamRes.mStatus), streamRes.mLUFS);
        check(streamRes.mStatus == BS1770Status::Ok, "stream: status Ok");
        check(std::fabs(streamRes.mLUFS - monoRes.mLUFS) <= 1e-6,
            "stream: chunked feeding matches one-shot");
    }

    /* Digital silence: all blocks fall below the -70 LUFS absolute gate. */
    const auto silence = std::vector<float>(static_cast<std::size_t>(rate*5.0), 0.0f);
    const auto silenceRes = BS1770Meter::measure(rate, 1, silence);
    std::printf("silence: status=%s\n", statusName(silenceRes.mStatus));
    check(silenceRes.mStatus == BS1770Status::Silence, "silence: distinct status");

    /* Invalid sample rates are reported distinctly. */
    check(BS1770Meter::measure(0.0, 1, mono).mStatus == BS1770Status::InvalidSampleRate,
        "sample rate 0: InvalidSampleRate");
    check(BS1770Meter::measure(-48000.0, 1, mono).mStatus == BS1770Status::InvalidSampleRate,
        "sample rate -48000: InvalidSampleRate");
    check(BS1770Meter::measure(1.0e9, 1, mono).mStatus == BS1770Status::InvalidSampleRate,
        "sample rate 1e9: InvalidSampleRate");

    /* Invalid channel count is reported distinctly. */
    check(BS1770Meter::measure(rate, 0, mono).mStatus == BS1770Status::InvalidChannels,
        "0 channels: InvalidChannels");

    /* Less than one 400 ms gating block: insufficient data. */
    const auto short_ = makeSine(1000.0, rate, 1, 0.3);
    const auto shortRes = BS1770Meter::measure(rate, 1, short_);
    std::printf("300ms input: status=%s\n", statusName(shortRes.mStatus));
    check(shortRes.mStatus == BS1770Status::InsufficientData,
        "300ms input: InsufficientData");

    if(gFailures != 0)
    {
        std::fprintf(stderr, "%d check(s) FAILED\n", gFailures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
