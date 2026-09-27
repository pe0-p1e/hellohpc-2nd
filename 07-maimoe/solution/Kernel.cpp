#include "maimoe/kernel_api.hpp"

#include "MaiMoeEngine.hpp"
#include "maimoe/render.hpp"
#include "maimoe/sha256.hpp"
#include "maimoe/simai.hpp"
#include "maimoe/state.hpp"
#include "maimoe/types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <string_view>
#include <type_traits>

namespace maimoe::kernel {
namespace {

constexpr std::array<std::uint8_t, 8> kFrameMagic = {
    'M', 'M', 'F', 'R', 'M', 0, 0, 3
};
constexpr std::array<std::uint8_t, 4> kResultMagic = {
    'M', 'M', 'R', '3'
};

template <typename UInt>
inline void store_le(std::uint8_t* dst, UInt value) noexcept {
    static_assert(std::is_unsigned_v<UInt>);
    for (std::size_t i = 0; i < sizeof(UInt); ++i) {
        dst[i] = static_cast<std::uint8_t>(value & static_cast<UInt>(0xffU));
        value >>= 8U;
    }
}

template <typename UInt>
inline void hash_le(Sha256& hasher, UInt value) {
    static_assert(std::is_unsigned_v<UInt>);
    std::array<std::uint8_t, sizeof(UInt)> bytes{};
    store_le(bytes.data(), value);
    hasher.update(bytes);
}

inline void hash_lp32(Sha256& hasher, std::span<const std::uint8_t> value) {
    hash_le(hasher, static_cast<std::uint32_t>(value.size()));
    hasher.update(value);
}

inline void hash_lp32(Sha256& hasher, std::string_view value) {
    hash_le(hasher, static_cast<std::uint32_t>(value.size()));
    hasher.update(value);
}

[[nodiscard]] Digest valid_result_digest(
    std::uint64_t chart_id,
    const Counts& counts,
    std::span<const std::uint8_t, kSerializedStateBytes> state,
    const Digest& frame_begin,
    const Digest& frame_end) {
    Sha256 hasher;
    hash_lp32(hasher, std::string_view("MaiMoeResult/v3"));
    hash_le(hasher, chart_id);
    hash_le(hasher, counts.total);
    hash_le(hasher, counts.tap);
    hash_le(hasher, counts.slide);
    hash_le(hasher, counts.hold);
    hash_le(hasher, counts.touch);
    hash_le(hasher, counts.break_count);
    hash_le(hasher, counts.warning);
    hash_lp32(hasher, std::span<const std::uint8_t>(state.data(), state.size()));
    hash_lp32(hasher, std::span<const std::uint8_t>(frame_begin.data(), frame_begin.size()));
    hash_lp32(hasher, std::span<const std::uint8_t>(frame_end.data(), frame_end.size()));
    return hasher.finish();
}

[[nodiscard]] std::array<std::uint8_t, 16> canonical_error_fields(
    const SemanticError& error) {
    if (!valid_semantic_error_code(error.code) || error.location.line == 0U ||
        error.location.column == 0U ||
        error.location.column > std::numeric_limits<std::uint16_t>::max()) {
        throw FormatError("semantic error cannot be encoded");
    }
    std::array<std::uint8_t, 16> fields{};
    store_le<std::uint16_t>(fields.data() + 0U,
                            static_cast<std::uint16_t>(error.code));
    store_le<std::uint16_t>(fields.data() + 2U,
                            static_cast<std::uint16_t>(error.location.column));
    store_le<std::uint32_t>(fields.data() + 4U, error.location.line);
    store_le<std::uint32_t>(fields.data() + 8U, error.detail);
    store_le<std::uint32_t>(fields.data() + 12U, 0U);
    return fields;
}

[[nodiscard]] Digest error_frame_seed(
    std::uint64_t chart_id,
    FrameRole role,
    std::span<const std::uint8_t, 16> fields) {
    Sha256 hasher;
    hash_lp32(hasher, std::string_view("MaiMoeErrorFrame/v3"));
    hash_le(hasher, chart_id);
    const std::array<std::uint8_t, 1> role_byte = {
        static_cast<std::uint8_t>(role)
    };
    hasher.update(role_byte);
    hash_lp32(hasher, std::span<const std::uint8_t>(fields.data(), fields.size()));
    return hasher.finish();
}

inline void fill_error_payload(
    std::span<std::uint8_t, kOutputFramePayloadBytes> payload,
    const Digest& seed) noexcept {
    static_assert(kOutputFramePayloadBytes % Digest{}.size() == 0U);
    for (std::size_t offset = 0; offset < payload.size(); offset += seed.size()) {
        std::memcpy(payload.data() + offset, seed.data(), seed.size());
    }
}

[[nodiscard]] Digest error_result_digest(
    std::uint64_t chart_id,
    std::span<const std::uint8_t, 16> fields,
    const Digest& frame_begin,
    const Digest& frame_end) {
    Sha256 hasher;
    hash_lp32(hasher, std::string_view("MaiMoeError/v3"));
    hash_le(hasher, chart_id);
    hash_lp32(hasher, std::span<const std::uint8_t>(fields.data(), fields.size()));
    hash_lp32(hasher, std::span<const std::uint8_t>(frame_begin.data(), frame_begin.size()));
    hash_lp32(hasher, std::span<const std::uint8_t>(frame_end.data(), frame_end.size()));
    return hasher.finish();
}

inline void write_frame_header(
    std::uint8_t* frame,
    std::uint64_t chart_id,
    FrameRole role,
    engine::ResultStatus status,
    const Digest& digest) noexcept {
    std::memcpy(frame, kFrameMagic.data(), kFrameMagic.size());
    store_le<std::uint16_t>(frame + 8U,
                            static_cast<std::uint16_t>(kOutputFrameHeaderBytes));
    store_le<std::uint16_t>(frame + 10U, kOutputFrameWidth);
    store_le<std::uint16_t>(frame + 12U, kOutputFrameHeight);
    frame[14U] = 2U;
    frame[15U] = static_cast<std::uint8_t>(role);
    store_le<std::uint64_t>(frame + 16U, chart_id);
    frame[24U] = static_cast<std::uint8_t>(status);
    frame[25U] = 0U;
    store_le<std::uint16_t>(frame + 26U, 0U);
    store_le<std::uint32_t>(frame + 28U, kOutputFramePayloadBytes);
    std::memcpy(frame + 32U, digest.data(), digest.size());
}

inline void write_result(
    std::uint8_t* result,
    std::uint64_t chart_id,
    engine::ResultStatus status,
    const Counts& counts,
    const std::uint8_t* error_fields,
    const Digest& frame_begin,
    const Digest& frame_end,
    const Digest& result_digest) noexcept {
    std::memcpy(result, kResultMagic.data(), kResultMagic.size());
    result[4U] = static_cast<std::uint8_t>(status);
    result[5U] = 0U;
    store_le<std::uint16_t>(result + 6U,
                            static_cast<std::uint16_t>(kOutputResultFileBytes));
    store_le<std::uint64_t>(result + 8U, chart_id);
    store_le<std::uint32_t>(result + 16U, counts.total);
    store_le<std::uint32_t>(result + 20U, counts.tap);
    store_le<std::uint32_t>(result + 24U, counts.slide);
    store_le<std::uint32_t>(result + 28U, counts.hold);
    store_le<std::uint32_t>(result + 32U, counts.touch);
    store_le<std::uint32_t>(result + 36U, counts.break_count);
    store_le<std::uint32_t>(result + 40U, counts.warning);
    if (error_fields == nullptr) {
        std::memset(result + 44U, 0, 16U);
    } else {
        std::memcpy(result + 44U, error_fields, 16U);
    }
    std::memcpy(result + 60U, frame_begin.data(), frame_begin.size());
    std::memcpy(result + 92U, frame_end.data(), frame_end.size());
    std::memcpy(result + 124U, result_digest.data(), result_digest.size());
    std::memset(result + 156U, 0, 4U);
}

}  // namespace

#if defined(__GNUC__)
__attribute__((hot))
#endif
void process_chart(
    std::uint64_t chart_id,
    std::string_view chart_text,
    std::span<std::uint8_t, kEncodedChartBytes> output) {
    const ParsedChart chart = parse_maidata(chart_text);

    auto begin_payload = std::span<std::uint8_t, kOutputFramePayloadBytes>(
        output.data() + kOutputFrameHeaderBytes, kOutputFramePayloadBytes);
    auto end_payload = std::span<std::uint8_t, kOutputFramePayloadBytes>(
        output.data() + kOutputFrameFileBytes + kOutputFrameHeaderBytes,
        kOutputFramePayloadBytes);

    Digest begin_digest{};
    Digest end_digest{};
    Digest result_digest{};
    Counts counts{};
    engine::ResultStatus status = engine::ResultStatus::Valid;
    std::array<std::uint8_t, 16> error_fields{};
    const std::uint8_t* result_error_fields = nullptr;

    if (const SemanticError* error = select_semantic_error(chart);
        error != nullptr) [[unlikely]] {
        status = engine::ResultStatus::Error;
        error_fields = canonical_error_fields(*error);
        result_error_fields = error_fields.data();

        const Digest begin_seed = error_frame_seed(
            chart_id, FrameRole::Begin,
            std::span<const std::uint8_t, 16>(error_fields));
        const Digest end_seed = error_frame_seed(
            chart_id, FrameRole::End,
            std::span<const std::uint8_t, 16>(error_fields));
        fill_error_payload(begin_payload, begin_seed);
        fill_error_payload(end_payload, end_seed);

        begin_digest = sha256(std::span<const std::uint8_t>(
            begin_payload.data(), begin_payload.size()));
        end_digest = sha256(std::span<const std::uint8_t>(
            end_payload.data(), end_payload.size()));
        result_digest = error_result_digest(
            chart_id,
            std::span<const std::uint8_t, 16>(error_fields),
            begin_digest, end_digest);
    } else [[likely]] {
        counts = chart.counts;
        const EndpointStates states = evolve_endpoint_states(chart_id, chart);

        render_output_rgba_into(
            chart_id, chart, states.frame_begin, kFirstSampleFrame, begin_payload);
        render_output_rgba_into(
            chart_id, chart, states.frame_end, kLastSampleFrame, end_payload);

        std::array<std::uint8_t, kSerializedStateBytes> serialized_state{};
        serialize_state_into(states.frame_end, serialized_state);

        begin_digest = sha256(std::span<const std::uint8_t>(
            begin_payload.data(), begin_payload.size()));
        end_digest = sha256(std::span<const std::uint8_t>(
            end_payload.data(), end_payload.size()));
        result_digest = valid_result_digest(
            chart_id, counts, serialized_state, begin_digest, end_digest);
    }

    write_frame_header(
        output.data(), chart_id, FrameRole::Begin, status, begin_digest);
    write_frame_header(
        output.data() + kOutputFrameFileBytes,
        chart_id, FrameRole::End, status, end_digest);
    write_result(
        output.data() + 2U * kOutputFrameFileBytes,
        chart_id, status, counts, result_error_fields,
        begin_digest, end_digest, result_digest);
}

}  // namespace maimoe::kernel
