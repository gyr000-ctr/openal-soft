/* Runtime SOFA HRTF loading tests.
 *
 * Exercises loading SOFA data sets through the public ALC API: enumeration,
 * selection, resampling to the device rate, fallback on corrupt files, and
 * concurrent switching. The data sets live in tests/sofa (see
 * generate-sofa-tests.py in that directory).
 *
 * ALSOFT_LOCAL_PATH is pointed at the test data directory so only the test
 * data sets (plus the built-in HRTF) are enumerated.
 */

#include <gtest/gtest.h>

#include <AL/al.h>
#include <AL/alc.h>
#include <AL/alext.h>

#include <atomic>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

namespace {

#ifndef SOFA_TEST_DIR
#define SOFA_TEST_DIR "."
#endif

/* Data sets expected to load successfully. */
const std::vector<std::string> kGoodSets{
    "stereo-44100",
    "stereo-48000",
    "stereo-96000",
    "stereo-48000-sofa2",
    "mono-48000",
    "mr-delay",
};

/* Data sets expected to fail loading, falling back to another data set. */
const std::vector<std::string> kCorruptSets{
    "corrupt-longir",
    "corrupt-nonfinite",
    "corrupt-silent",
    "corrupt-huge",
    "corrupt-truncated",
    "corrupt-not-sofa",
};

class SofaTest : public ::testing::Test {
protected:
    ALCdevice *mDevice{nullptr};

    void SetUp() override
    {
        /* Point the data file search at the test assets and use the null
         * backend so the tests run anywhere.
         */
#ifdef _WIN32
        _putenv_s("ALSOFT_DRIVERS", "null");
        _putenv_s("ALSOFT_LOCAL_PATH", SOFA_TEST_DIR);
#else
        setenv("ALSOFT_DRIVERS", "null", 1);
        setenv("ALSOFT_LOCAL_PATH", SOFA_TEST_DIR, 1);
#endif
        mDevice = alcOpenDevice(nullptr);
        ASSERT_TRUE(mDevice != nullptr);
    }

    void TearDown() override
    {
        if(mDevice)
            alcCloseDevice(mDevice);
    }

    auto GetSpecifiers() -> std::vector<std::string>
    {
        auto result = std::vector<std::string>{};
        ALCint count = 0;
        alcGetIntegerv(mDevice, ALC_NUM_HRTF_SPECIFIERS_SOFT, 1, &count);
        for(ALCint i = 0;i < count;++i)
        {
            if(const ALCchar *name = alcGetStringiSOFT(mDevice, ALC_HRTF_SPECIFIER_SOFT, i))
                result.emplace_back(name);
        }
        return result;
    }

    auto FindSpecifier(const std::string &name) -> int
    {
        const auto specifiers = GetSpecifiers();
        for(size_t i = 0;i < specifiers.size();++i)
        {
            if(specifiers[i] == name)
                return static_cast<int>(i);
        }
        return -1;
    }

    /* Resets the device with the given HRTF specifier index, returning the
     * resulting current HRTF name (empty on failure).
     */
    auto ResetWithHrtf(int hrtf_id, ALCint freq = 0) -> std::string
    {
        ALCint attrs[] = {
            ALC_HRTF_SOFT, ALC_TRUE,
            ALC_HRTF_ID_SOFT, hrtf_id,
            ALC_FREQUENCY, freq,
            0, 0
        };
        if(!alcResetDeviceSOFT(mDevice, attrs))
            return {};

        ALCint status = 0;
        alcGetIntegerv(mDevice, ALC_HRTF_STATUS_SOFT, 1, &status);
        EXPECT_EQ(status, ALC_HRTF_ENABLED_SOFT);

        const ALCchar *current = alcGetString(mDevice, ALC_HRTF_SPECIFIER_SOFT);
        return current ? current : "";
    }
};

TEST_F(SofaTest, Enumeration)
{
    const auto specifiers = GetSpecifiers();
    ASSERT_FALSE(specifiers.empty());

#if HAVE_MYSOFA
    /* All test data sets (good and corrupt) should be enumerated. */
    for(const auto &name : kGoodSets)
        EXPECT_NE(FindSpecifier(name), -1) << "missing specifier: " << name;
    for(const auto &name : kCorruptSets)
        EXPECT_NE(FindSpecifier(name), -1) << "missing specifier: " << name;
#else
    /* Without MySOFA, SOFA files must not be enumerated. */
    for(const auto &name : kGoodSets)
        EXPECT_EQ(FindSpecifier(name), -1) << "unexpected specifier: " << name;
    for(const auto &name : kCorruptSets)
        EXPECT_EQ(FindSpecifier(name), -1) << "unexpected specifier: " << name;
#endif
    /* The built-in data set is always available as a fallback. */
    EXPECT_NE(FindSpecifier("Built-In HRTF"), -1);
}

#if HAVE_MYSOFA
TEST_F(SofaTest, LoadGoodSets)
{
    for(const auto &name : kGoodSets)
    {
        const int idx = FindSpecifier(name);
        ASSERT_NE(idx, -1) << "missing specifier: " << name;
        EXPECT_EQ(ResetWithHrtf(idx), name) << "failed to load: " << name;
    }
}

TEST_F(SofaTest, CorruptSetsFallBack)
{
    for(const auto &name : kCorruptSets)
    {
        const int idx = FindSpecifier(name);
        ASSERT_NE(idx, -1) << "missing specifier: " << name;
        /* The reset itself must succeed, but a different (fallback) data set
         * must be in use afterwards.
         */
        const auto current = ResetWithHrtf(idx);
        EXPECT_FALSE(current.empty()) << "reset failed for: " << name;
        EXPECT_NE(current, name) << "corrupt set was loaded: " << name;
    }
}

TEST_F(SofaTest, ResampleToDeviceRate)
{
    const int idx = FindSpecifier("stereo-48000");
    ASSERT_NE(idx, -1);

    /* The data set is at 48kHz; the device rates exercise the resampler in
     * both directions (the 96kHz set is covered by LoadGoodSets at the
     * default device rate).
     */
    for(const ALCint freq : {44100, 48000, 96000})
        EXPECT_EQ(ResetWithHrtf(idx, freq), "stereo-48000") << "device rate: " << freq;
}

TEST_F(SofaTest, RapidSwitching)
{
    const auto specifiers = GetSpecifiers();
    ASSERT_FALSE(specifiers.empty());

    for(int iter = 0;iter < 10;++iter)
    {
        for(size_t i = 0;i < specifiers.size();++i)
        {
            const auto current = ResetWithHrtf(static_cast<int>(i));
            EXPECT_FALSE(current.empty()) << "switch failed at iter=" << iter << " id=" << i;
        }
    }
}

TEST_F(SofaTest, ConcurrentSwitching)
{
    /* Concurrent initialization/switching on separate devices, to shake out
     * races in the shared HRTF cache.
     */
    std::atomic<int> failures{0};
    auto worker = [&failures]()
    {
        ALCdevice *device = alcOpenDevice(nullptr);
        if(!device)
        {
            ++failures;
            return;
        }
        ALCint count = 0;
        alcGetIntegerv(device, ALC_NUM_HRTF_SPECIFIERS_SOFT, 1, &count);
        for(int iter = 0;iter < 30;++iter)
        {
            for(ALCint i = 0;i < count;++i)
            {
                ALCint attrs[] = {ALC_HRTF_SOFT, ALC_TRUE, ALC_HRTF_ID_SOFT, i, 0};
                if(!alcResetDeviceSOFT(device, attrs))
                    ++failures;
            }
        }
        alcCloseDevice(device);
    };

    std::thread t1(worker), t2(worker), t3(worker), t4(worker);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    EXPECT_EQ(failures.load(), 0);
}
#endif /* HAVE_MYSOFA */

TEST_F(SofaTest, DefaultStillWorks)
{
    /* The built-in data set must remain selectable regardless of MySOFA. */
    const int idx = FindSpecifier("Built-In HRTF");
    ASSERT_NE(idx, -1);
    EXPECT_EQ(ResetWithHrtf(idx), "Built-In HRTF");
}

} // namespace
