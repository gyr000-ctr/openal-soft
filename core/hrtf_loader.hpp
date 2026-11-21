#ifndef CORE_HRTF_LOADER_HPP
#define CORE_HRTF_LOADER_HPP

#include <iosfwd>
#include <memory>
#include <stdexcept>


struct HrtfStore;

/* Runtime loading of HRTF data sets. Both the native .mhr formats and AES69
 * SOFA files are accepted (the latter when built with ALSOFT_RUNTIME_SOFA);
 * the stream format is detected from its leading bytes.
 *
 * Throws std::runtime_error (or a subclass) on invalid or unsupported data.
 */
auto LoadHrtf(std::istream &stream) -> std::unique_ptr<HrtfStore>;


/* Base type for SOFA loading failures. */
class SofaLoadError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/* The file is valid SOFA, but its dimensions or layout don't fit the HRTF
 * model (bad emitter/receiver/IR dimensions, non-uniform measurement grid,
 * out-of-range sample rate or distances, and so on).
 */
class SofaDimensionError final : public SofaLoadError {
public:
    using SofaLoadError::SofaLoadError;
};

/* The file is structurally valid, but required measurement data is absent or
 * unusable (no IR/source positions/sample rate, missing measurements for a
 * declared grid cell, non-finite or silent responses).
 */
class SofaMissingDataError final : public SofaLoadError {
public:
    using SofaLoadError::SofaLoadError;
};

/* The file exceeds the size limit imposed for runtime loading, checked before
 * the contents are parsed.
 */
class SofaFileTooLargeError final : public SofaLoadError {
public:
    using SofaLoadError::SofaLoadError;
};

#endif /* CORE_HRTF_LOADER_HPP */
