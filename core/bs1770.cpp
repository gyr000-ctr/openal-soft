#include "bs1770.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>


namespace {

/* Number of 100 ms sub-blocks per 400 ms gating block (75% overlap). */
constexpr std::size_t SubBlocksPerBlock{4};

} // namespace

auto BS1770Meter::init(double const sampleRate, unsigned const channels) -> BS1770Status
{
    if(channels == 0 || channels > MaxChannels)
        return BS1770Status::InvalidChannels;
    auto weights = std::vector<double>(channels);
    if(!defaultChannelWeights(channels, weights))
        return BS1770Status::InvalidChannels;
    return init(sampleRate, weights);
}

auto BS1770Meter::init(double const sampleRate, std::span<const double> const channelWeights)
    -> BS1770Status
{
    /* Start from a clean slate so a failed init leaves a usable object. */
    mChannels.clear();
    mBlocks.clear();
    mSubBlocks.clear();
    mSubBlockLen = 0;
    mSubBlockPos = 0;
    mSubBlockSum = 0.0;

    if(!std::isfinite(sampleRate) || sampleRate < MinSampleRate || sampleRate > MaxSampleRate)
        return BS1770Status::InvalidSampleRate;
    if(channelWeights.empty() || channelWeights.size() > MaxChannels)
        return BS1770Status::InvalidChannels;
    if(!std::ranges::all_of(channelWeights,
        [](double const w) { return std::isfinite(w) && w >= 0.0; }))
        return BS1770Status::InvalidChannels;

    /* K-weighting pre-filter, stage 1: high shelf. The constants are the
     * analog-prototype parameters of the BS.1770-4 filter; the bilinear
     * transform below reproduces the spec's 48 kHz reference coefficients
     * (b = [1.53512485958697, -2.69169618940638, 1.19839281085285],
     *  a = [1, -1.69065929318241, 0.73248077421585]) to within 1e-13.
     */
    auto shelf = BiquadState{};
    {
        constexpr auto f0 = 1681.974450955533;
        constexpr auto G = 3.999843853973347;
        constexpr auto Q = 0.7071752369554196;
        const auto K = std::tan(std::numbers::pi*f0 / sampleRate);
        const auto Vh = std::pow(10.0, G/20.0);
        const auto Vb = std::pow(Vh, 0.499666774155);
        const auto a0 = 1.0 + K/Q + K*K;
        shelf.mB0 = (Vh + Vb*K/Q + K*K) / a0;
        shelf.mB1 = 2.0*(K*K - Vh) / a0;
        shelf.mB2 = (Vh - Vb*K/Q + K*K) / a0;
        shelf.mA1 = 2.0*(K*K - 1.0) / a0;
        shelf.mA2 = (1.0 - K/Q + K*K) / a0;
    }

    /* K-weighting pre-filter, stage 2: high-pass. As in the spec's
     * reference coefficients, the numerator stays exactly [1, -2, 1]
     * (a = [1, -1.99004745483398, 0.99007225036621] at 48 kHz).
     */
    auto highpass = BiquadState{};
    {
        constexpr auto f0 = 38.13547087602444;
        constexpr auto Q = 0.5003270373238773;
        const auto K = std::tan(std::numbers::pi*f0 / sampleRate);
        const auto a0 = 1.0 + K/Q + K*K;
        highpass.mB0 = 1.0;
        highpass.mB1 = -2.0;
        highpass.mB2 = 1.0;
        highpass.mA1 = 2.0*(K*K - 1.0) / a0;
        highpass.mA2 = (1.0 - K/Q + K*K) / a0;
    }

    mChannels.resize(channelWeights.size());
    for(auto i = std::size_t{0}; i < channelWeights.size(); ++i)
        mChannels[i] = Channel{shelf, highpass, channelWeights[i]};

    mSubBlockLen = std::max<std::size_t>(1,
        static_cast<std::size_t>(sampleRate*(BlockDuration/SubBlocksPerBlock) + 0.5));
    return BS1770Status::Ok;
}

void BS1770Meter::reset() noexcept
{
    for(auto &ch : mChannels)
    {
        ch.mShelf.mZ1 = ch.mShelf.mZ2 = 0.0;
        ch.mHighPass.mZ1 = ch.mHighPass.mZ2 = 0.0;
    }
    mBlocks.clear();
    mSubBlocks.clear();
    mSubBlockPos = 0;
    mSubBlockSum = 0.0;
}

void BS1770Meter::process(std::span<const float> const samples)
{
    const auto chcount = mChannels.size();
    if(chcount == 0 || mSubBlockLen == 0)
        return;

    const auto frames = samples.size() / chcount;
    const auto *in = samples.data();
    for(auto f = std::size_t{0}; f < frames; ++f)
    {
        for(auto &ch : mChannels)
        {
            const auto y = ch.mHighPass.step(ch.mShelf.step(static_cast<double>(*in++)));
            mSubBlockSum += ch.mWeight * y*y;
        }
        if(++mSubBlockPos == mSubBlockLen)
        {
            mSubBlocks.push_back(mSubBlockSum);
            mSubBlockSum = 0.0;
            mSubBlockPos = 0;
            if(mSubBlocks.size() == SubBlocksPerBlock)
            {
                /* A 400 ms gating block is the last four 100 ms sub-blocks;
                 * its energy is the weighted channel sum of mean squares.
                 */
                const auto total = std::accumulate(mSubBlocks.begin(), mSubBlocks.end(), 0.0);
                mBlocks.push_back(total / static_cast<double>(mSubBlockLen*SubBlocksPerBlock));
                mSubBlocks.pop_front();
            }
        }
    }
}

auto BS1770Meter::integrated() const -> Result
{
    static constexpr auto NegInf = -std::numeric_limits<double>::infinity();

    if(mBlocks.empty())
        return {BS1770Status::InsufficientData, NegInf, NegInf};

    /* Absolute gate at -70 LUFS, expressed as a block energy. */
    const auto absGate = std::pow(10.0, (AbsoluteGate - LoudnessOffset)/10.0);

    auto sum = 0.0;
    auto count = std::size_t{0};
    for(const auto energy : mBlocks)
    {
        if(energy >= absGate)
        {
            sum += energy;
            ++count;
        }
    }
    if(count == 0)
        return {BS1770Status::Silence, NegInf, NegInf};

    /* Relative gate, 10 LU below the absolute-gated mean loudness. */
    const auto relGate = (sum/static_cast<double>(count)) * std::pow(10.0, RelativeGate/10.0);

    auto gatedSum = 0.0;
    auto gatedCount = std::size_t{0};
    for(const auto energy : mBlocks)
    {
        if(energy >= absGate && energy >= relGate)
        {
            gatedSum += energy;
            ++gatedCount;
        }
    }

    const auto loudness = LoudnessOffset
        + 10.0*std::log10(gatedSum/static_cast<double>(gatedCount));
    return {BS1770Status::Ok, loudness, loudness};
}

auto BS1770Meter::measure(double const sampleRate, unsigned const channels,
    std::span<const float> const samples) -> Result
{
    static constexpr auto NegInf = -std::numeric_limits<double>::infinity();

    auto meter = BS1770Meter{};
    if(const auto status = meter.init(sampleRate, channels); status != BS1770Status::Ok)
        return {status, NegInf, NegInf};
    meter.process(samples);
    return meter.integrated();
}

auto BS1770Meter::defaultChannelWeights(unsigned const channels, std::span<double> out) -> bool
{
    if(channels == 0 || channels > MaxChannels || out.size() < channels)
        return false;

    std::ranges::fill_n(out.begin(), channels, 1.0);
    switch(channels)
    {
    case 4: /* Quad: FL, FR, BL, BR. */
        out[2] = out[3] = 1.41;
        break;
    case 5: /* 5.0: FL, FR, FC, BL, BR. */
        out[3] = out[4] = 1.41;
        break;
    case 6: /* 5.1: FL, FR, FC, LFE, BL, BR. */
        out[3] = 0.0;
        out[4] = out[5] = 1.41;
        break;
    case 7: /* 6.1: FL, FR, FC, LFE, BC, SL, SR. */
        out[3] = 0.0;
        out[4] = out[5] = out[6] = 1.41;
        break;
    case 8: /* 7.1: FL, FR, FC, LFE, BL, BR, SL, SR. */
        out[3] = 0.0;
        out[4] = out[5] = out[6] = out[7] = 1.41;
        break;
    default: /* Mono, stereo, 3-channel, and unknown layouts: all 1.0. */
        break;
    }
    return true;
}
