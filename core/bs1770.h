#ifndef CORE_BS1770_H
#define CORE_BS1770_H

#include <cstddef>
#include <deque>
#include <limits>
#include <span>
#include <vector>

/**
 * Loudness metering according to ITU-R BS.1770-4.
 *
 * Implements the complete BS.1770-4 algorithm: K-weighting (a high-shelf
 * stage followed by a high-pass stage), 400 ms gating blocks with 75%
 * overlap, the -70 LUFS absolute gate, and the -10 LU relative gate. The
 * result is reported as gated integrated loudness in LUFS; LKFS (the ITU
 * name for the identical scale) is reported alongside it for interfaces
 * that expect that term.
 *
 * The meter is a standalone, offline component. It only observes copies of
 * sample data handed to it and is not part of, nor has any effect on, the
 * render or mixing paths.
 */

/** Result state of a loudness measurement, distinguishing failure causes. */
enum class BS1770Status {
    Ok,                /**< The measurement is valid. */
    Silence,           /**< All gating blocks are below the -70 LUFS absolute gate. */
    InvalidSampleRate, /**< The sample rate is outside the supported range. */
    InvalidChannels,   /**< The channel count or channel weights are not usable. */
    InsufficientData,  /**< Less than one 400 ms gating block has been processed. */
};

class BS1770Meter {
public:
    /** The result of an integrated loudness measurement. */
    struct Result {
        BS1770Status mStatus{BS1770Status::InsufficientData};
        /** Gated integrated loudness in LUFS, or -inf when not Ok. */
        double mLUFS{-std::numeric_limits<double>::infinity()};
        /** Same value as mLUFS; LKFS is the ITU name for the LUFS scale. */
        double mLKFS{-std::numeric_limits<double>::infinity()};
    };

private:
    struct BiquadState {
        /* Transfer function coefficients (a0 pre-applied). */
        double mB0{0.0}, mB1{0.0}, mB2{0.0}, mA1{0.0}, mA2{0.0};
        /* Direct form II transposed delay elements. */
        double mZ1{0.0}, mZ2{0.0};

        auto step(double const in) noexcept -> double
        {
            const auto out = in*mB0 + mZ1;
            mZ1 = in*mB1 - out*mA1 + mZ2;
            mZ2 = in*mB2 - out*mA2;
            return out;
        }
    };

    struct Channel {
        BiquadState mShelf, mHighPass;
        double mWeight{1.0};
    };

    std::vector<Channel> mChannels;
    /* Weighted mean-square energies of completed 400 ms gating blocks. */
    std::vector<double> mBlocks;
    /* Weighted square sums of the most recent 100 ms sub-blocks. */
    std::deque<double> mSubBlocks;
    /* Frames per 100 ms sub-block, and position within the current one. */
    std::size_t mSubBlockLen{0};
    std::size_t mSubBlockPos{0};
    /* Accumulated weighted square sum of the current sub-block. */
    double mSubBlockSum{0.0};

public:
    /* BS.1770-4 algorithm constants. */
    static constexpr double AbsoluteGate{-70.0}; /**< Absolute gate, LUFS. */
    static constexpr double RelativeGate{-10.0}; /**< Relative gate, LU below the abs-gated mean. */
    static constexpr double BlockDuration{0.4};  /**< Gating block length, seconds. */
    static constexpr double BlockOverlap{0.75};  /**< Overlap between gating blocks. */
    static constexpr double LoudnessOffset{-0.691};

    static constexpr double MinSampleRate{8000.0};
    static constexpr double MaxSampleRate{768000.0};
    static constexpr std::size_t MaxChannels{64};

    BS1770Meter() = default;

    /**
     * Initializes the meter for the given sample rate and channel count,
     * using the default BS.1770 channel weights (see defaultChannelWeights).
     * Returns BS1770Status::Ok on success, or the reason for failure.
     */
    auto init(double sampleRate, unsigned channels) -> BS1770Status;

    /**
     * Initializes the meter with explicit per-channel weights (e.g. 1.0 for
     * front channels, 1.41 for surrounds, 0.0 to exclude a channel such as
     * LFE). Returns BS1770Status::Ok on success, or the reason for failure.
     */
    auto init(double sampleRate, std::span<const double> channelWeights) -> BS1770Status;

    /** Clears all measurement history, keeping the current configuration. */
    void reset() noexcept;

    /**
     * Feeds interleaved float PCM frames to the meter. The frame count is
     * samples.size() / channel count; a trailing partial frame is ignored.
     * Does nothing if the meter is not initialized.
     */
    void process(std::span<const float> samples);

    /**
     * Computes the gated integrated loudness of everything fed so far.
     * Returns a Result with status Ok and the loudness in LUFS/LKFS, or a
     * status distinguishing Silence from InsufficientData.
     */
    [[nodiscard]]
    auto integrated() const -> Result;

    /**
     * One-shot convenience: initializes a meter with the default channel
     * weights, processes the given interleaved samples, and returns the
     * gated integrated loudness.
     */
    [[nodiscard]]
    static auto measure(double sampleRate, unsigned channels, std::span<const float> samples)
        -> Result;

    /**
     * Fills `out` with the default BS.1770 channel weights for the given
     * channel count, assuming OpenAL channel ordering (FL, FR, FC, LFE, BL,
     * BR, SL, SR): 1.0 for front channels, 1.41 for surrounds, and 0.0 for
     * LFE. Unknown layouts default to 1.0 for all channels. Returns false if
     * the channel count is not usable.
     */
    static auto defaultChannelWeights(unsigned channels, std::span<double> out) -> bool;
};

#endif /* CORE_BS1770_H */
