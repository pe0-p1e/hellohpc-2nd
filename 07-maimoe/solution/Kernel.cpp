#include "maimoe/kernel_api.hpp"

#include "MaiMoeEngine.hpp"
#include "maimoe/render.hpp"
#include "maimoe/sha256.hpp"
#include "maimoe/simai.hpp"
#include "maimoe/state.hpp"
#include "maimoe/types.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__aarch64__)
#include <arm_neon.h>
#endif

namespace maimoe::kernel {
namespace {

constexpr std::array<std::uint8_t, 8> kFrameMagic = {
    'M', 'M', 'F', 'R', 'M', 0, 0, 3
};
constexpr std::array<std::uint8_t, 4> kResultMagic = {
    'M', 'M', 'R', '3'
};

#if defined(__aarch64__)

alignas(16) constexpr std::array<std::uint32_t, 64> kSha256RoundConstantsArm = {
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
    0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
    0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
    0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
    0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
    0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
    0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
    0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
    0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
    0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
    0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
    0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
    0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
    0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
    0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
    0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U,
};

__attribute__((target("+crypto"), noinline))
void sha256_arm_blocks(std::uint32_t state[8],
                       const std::uint8_t* data,
                       std::size_t length) noexcept {
    uint32x4_t state0 = vld1q_u32(state + 0);
    uint32x4_t state1 = vld1q_u32(state + 4);

    while (length >= 64U) {
        const uint32x4_t save0 = state0;
        const uint32x4_t save1 = state1;

        uint32x4_t msg0 = vreinterpretq_u32_u8(
            vrev32q_u8(vld1q_u8(data + 0U)));
        uint32x4_t msg1 = vreinterpretq_u32_u8(
            vrev32q_u8(vld1q_u8(data + 16U)));
        uint32x4_t msg2 = vreinterpretq_u32_u8(
            vrev32q_u8(vld1q_u8(data + 32U)));
        uint32x4_t msg3 = vreinterpretq_u32_u8(
            vrev32q_u8(vld1q_u8(data + 48U)));

        uint32x4_t tmp0 = vaddq_u32(
            msg0, vld1q_u32(kSha256RoundConstantsArm.data() + 0));
        uint32x4_t tmp1;
        uint32x4_t tmp2;

        msg0 = vsha256su0q_u32(msg0, msg1);
        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg1, vld1q_u32(kSha256RoundConstantsArm.data() + 4));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);
        msg0 = vsha256su1q_u32(msg0, msg2, msg3);

        msg1 = vsha256su0q_u32(msg1, msg2);
        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg2, vld1q_u32(kSha256RoundConstantsArm.data() + 8));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);
        msg1 = vsha256su1q_u32(msg1, msg3, msg0);

        msg2 = vsha256su0q_u32(msg2, msg3);
        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg3, vld1q_u32(kSha256RoundConstantsArm.data() + 12));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);
        msg2 = vsha256su1q_u32(msg2, msg0, msg1);

        msg3 = vsha256su0q_u32(msg3, msg0);
        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg0, vld1q_u32(kSha256RoundConstantsArm.data() + 16));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);
        msg3 = vsha256su1q_u32(msg3, msg1, msg2);

        msg0 = vsha256su0q_u32(msg0, msg1);
        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg1, vld1q_u32(kSha256RoundConstantsArm.data() + 20));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);
        msg0 = vsha256su1q_u32(msg0, msg2, msg3);

        msg1 = vsha256su0q_u32(msg1, msg2);
        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg2, vld1q_u32(kSha256RoundConstantsArm.data() + 24));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);
        msg1 = vsha256su1q_u32(msg1, msg3, msg0);

        msg2 = vsha256su0q_u32(msg2, msg3);
        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg3, vld1q_u32(kSha256RoundConstantsArm.data() + 28));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);
        msg2 = vsha256su1q_u32(msg2, msg0, msg1);

        msg3 = vsha256su0q_u32(msg3, msg0);
        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg0, vld1q_u32(kSha256RoundConstantsArm.data() + 32));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);
        msg3 = vsha256su1q_u32(msg3, msg1, msg2);

        msg0 = vsha256su0q_u32(msg0, msg1);
        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg1, vld1q_u32(kSha256RoundConstantsArm.data() + 36));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);
        msg0 = vsha256su1q_u32(msg0, msg2, msg3);

        msg1 = vsha256su0q_u32(msg1, msg2);
        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg2, vld1q_u32(kSha256RoundConstantsArm.data() + 40));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);
        msg1 = vsha256su1q_u32(msg1, msg3, msg0);

        msg2 = vsha256su0q_u32(msg2, msg3);
        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg3, vld1q_u32(kSha256RoundConstantsArm.data() + 44));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);
        msg2 = vsha256su1q_u32(msg2, msg0, msg1);

        msg3 = vsha256su0q_u32(msg3, msg0);
        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg0, vld1q_u32(kSha256RoundConstantsArm.data() + 48));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);
        msg3 = vsha256su1q_u32(msg3, msg1, msg2);

        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg1, vld1q_u32(kSha256RoundConstantsArm.data() + 52));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);

        tmp2 = state0;
        tmp0 = vaddq_u32(
            msg2, vld1q_u32(kSha256RoundConstantsArm.data() + 56));
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);

        tmp2 = state0;
        tmp1 = vaddq_u32(
            msg3, vld1q_u32(kSha256RoundConstantsArm.data() + 60));
        state0 = vsha256hq_u32(state0, state1, tmp0);
        state1 = vsha256h2q_u32(state1, tmp2, tmp0);

        tmp2 = state0;
        state0 = vsha256hq_u32(state0, state1, tmp1);
        state1 = vsha256h2q_u32(state1, tmp2, tmp1);

        state0 = vaddq_u32(state0, save0);
        state1 = vaddq_u32(state1, save1);

        data += 64U;
        length -= 64U;
    }

    vst1q_u32(state + 0, state0);
    vst1q_u32(state + 4, state1);
}

[[nodiscard]] Digest sha256_arm(std::span<const std::uint8_t> bytes) noexcept {
    alignas(16) std::uint32_t state[8] = {
        0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
        0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U,
    };

    const std::size_t full_bytes = bytes.size() & ~std::size_t{63U};
    if (full_bytes != 0U) {
        sha256_arm_blocks(state, bytes.data(), full_bytes);
    }

    alignas(16) std::array<std::uint8_t, 128> tail{};
    const std::size_t remaining = bytes.size() - full_bytes;
    if (remaining != 0U) {
        std::memcpy(tail.data(), bytes.data() + full_bytes, remaining);
    }
    tail[remaining] = 0x80U;

    const std::size_t padded_bytes = remaining < 56U ? 64U : 128U;
    const std::uint64_t bit_length =
        static_cast<std::uint64_t>(bytes.size()) * 8ULL;
    for (std::size_t i = 0U; i < 8U; ++i) {
        tail[padded_bytes - 1U - i] =
            static_cast<std::uint8_t>(bit_length >> (8U * i));
    }
    sha256_arm_blocks(state, tail.data(), padded_bytes);

    Digest digest{};
    for (std::size_t i = 0U; i < 8U; ++i) {
        digest[4U * i + 0U] =
            static_cast<std::uint8_t>(state[i] >> 24U);
        digest[4U * i + 1U] =
            static_cast<std::uint8_t>(state[i] >> 16U);
        digest[4U * i + 2U] =
            static_cast<std::uint8_t>(state[i] >> 8U);
        digest[4U * i + 3U] =
            static_cast<std::uint8_t>(state[i]);
    }
    return digest;
}

[[nodiscard]] inline Digest fast_frame_sha256(
    std::span<const std::uint8_t> bytes) {
    static int arm_sha_mode = 0;
    if (arm_sha_mode > 0) [[likely]] {
        return sha256_arm(bytes);
    }
    if (arm_sha_mode < 0) [[unlikely]] {
        return sha256(bytes);
    }

    const Digest hardware = sha256_arm(bytes);
    const Digest software = sha256(bytes);
    if (hardware == software) {
        arm_sha_mode = 1;
        return hardware;
    }
    arm_sha_mode = -1;
    return software;
}

#else

[[nodiscard]] inline Digest fast_frame_sha256(
    std::span<const std::uint8_t> bytes) {
    return sha256(bytes);
}

#endif

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

struct SourceState {
    bool have_tick = false;
    std::uint32_t tick = 0;
    std::uint32_t active_end = 0;
};

struct CompactEvent {
    EventType type = EventType::Tap;
    std::uint32_t tick = 0;
    std::uint32_t end_tick = 0;
    std::uint8_t lane = 0;
    std::uint8_t target_lane = 0;
    std::uint8_t touch_sensor = 0;
    std::uint8_t strength = 0;
    std::uint8_t path_id = 0;
    SourceLocation location{};
    std::uint64_t word = 0;
};

struct FastMetadata {
    std::uint32_t whole_bpm_milli = 120000U;
};

struct FastChart {
    FastMetadata metadata;
    std::vector<CompactEvent> events;
    Counts counts;
    std::uint32_t warning_count = 0;
    std::uint32_t max_tick = 0;
};

struct FastParser {
    FastChart chart;
    std::uint64_t tick = 0;
    std::uint32_t division = 4;
    std::uint32_t bpm_milli = 0;
    std::uint32_t parsed_events = 0;
    std::uint32_t warning_count = 0;
    bool initial_bpm_control = true;
    bool initial_division_control = true;
    bool conventional_division_pending = false;
    bool ended = false;
    std::array<SourceState, 8> button_sources{};
    std::array<SourceState, kTouchSensorCount> touch_sources{};
};

template <typename UInt>
[[nodiscard]] bool parse_uint_fast(std::string_view text, UInt& value) noexcept {
    static_assert(std::is_unsigned_v<UInt>);
    if (text.empty() || (text.size() > 1U && text.front() == '0')) {
        return false;
    }
    value = 0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    return error == std::errc{} && end == text.data() + text.size();
}

[[nodiscard]] bool parse_bpm_fast(std::string_view text,
                                  std::uint32_t& result) noexcept {
    const std::size_t dot = text.find('.');
    if (dot != std::string_view::npos &&
        text.find('.', dot + 1U) != std::string_view::npos) {
        return false;
    }
    const std::string_view whole_text = text.substr(0U, dot);
    const std::string_view fraction =
        dot == std::string_view::npos ? std::string_view{} : text.substr(dot + 1U);
    if (whole_text.empty() ||
        (whole_text.size() > 1U && whole_text.front() == '0') ||
        (dot != std::string_view::npos &&
         (fraction.empty() || fraction.size() > 3U))) {
        return false;
    }

    std::uint32_t whole = 0;
    if (!parse_uint_fast(whole_text, whole)) {
        return false;
    }

    std::uint32_t fractional = 0;
    if (!fraction.empty()) {
        if (fraction.back() == '0') {
            return false;
        }
        for (const char digit : fraction) {
            if (digit < '0' || digit > '9') {
                return false;
            }
            fractional =
                fractional * 10U + static_cast<std::uint32_t>(digit - '0');
        }
        if (fraction.size() == 1U) {
            fractional *= 100U;
        } else if (fraction.size() == 2U) {
            fractional *= 10U;
        }
    }
    if (whole >
        (std::numeric_limits<std::uint32_t>::max() - fractional) / 1000U) {
        return false;
    }
    result = whole * 1000U + fractional;
    return true;
}

[[nodiscard]] bool metadata_value_fast(std::string_view line,
                                       std::string_view prefix) noexcept {
    if (!line.starts_with(prefix)) {
        return false;
    }
    const std::string_view value = line.substr(prefix.size());
    return !value.empty() && value.front() != ' ' && value.back() != ' ';
}

[[nodiscard]] bool touch_sensor_fast(std::string_view text,
                                     std::uint8_t& sensor) noexcept {
    if (text == "C") {
        sensor = 16U;
        return true;
    }
    if (text.size() != 2U || text[1] < '1' || text[1] > '8') {
        return false;
    }
    const std::uint8_t ring = static_cast<std::uint8_t>(text[1] - '1');
    switch (text[0]) {
        case 'A': sensor = ring; return true;
        case 'B': sensor = static_cast<std::uint8_t>(8U + ring); return true;
        case 'D': sensor = static_cast<std::uint8_t>(17U + ring); return true;
        case 'E': sensor = static_cast<std::uint8_t>(25U + ring); return true;
        default: return false;
    }
}

inline void append_char(char*& cursor, char value) noexcept {
    *cursor++ = value;
}

inline void append_text(char*& cursor, std::string_view value) noexcept {
    std::memcpy(cursor, value.data(), value.size());
    cursor += value.size();
}

inline bool append_u32(char*& cursor, char* end, std::uint32_t value) noexcept {
    const auto result = std::to_chars(cursor, end, value);
    if (result.ec != std::errc{}) {
        return false;
    }
    cursor = result.ptr;
    return true;
}

inline void append_touch_name(char*& cursor, std::uint8_t sensor) noexcept {
    if (sensor < 8U) {
        *cursor++ = 'A';
        *cursor++ = static_cast<char>('1' + sensor);
    } else if (sensor < 16U) {
        *cursor++ = 'B';
        *cursor++ = static_cast<char>('1' + sensor - 8U);
    } else if (sensor == 16U) {
        *cursor++ = 'C';
    } else if (sensor < 25U) {
        *cursor++ = 'D';
        *cursor++ = static_cast<char>('1' + sensor - 17U);
    } else {
        *cursor++ = 'E';
        *cursor++ = static_cast<char>('1' + sensor - 25U);
    }
}

[[nodiscard]] bool set_event_word_fast(CompactEvent& event) {
    std::array<char, 96> buffer{};
    char* cursor = buffer.data();
    char* const end = buffer.data() + buffer.size();

    if (!append_u32(cursor, end, event.tick)) {
        return false;
    }
    append_char(cursor, ' ');

    switch (event.type) {
        case EventType::Tap:
            append_text(cursor, "TAP B");
            append_char(cursor, static_cast<char>('1' + event.lane));
            append_text(cursor, " 128");
            break;
        case EventType::Break:
            append_text(cursor, "BREAK B");
            append_char(cursor, static_cast<char>('1' + event.lane));
            append_text(cursor, " 255 3");
            break;
        case EventType::Hold:
            append_text(cursor, "HOLD B");
            append_char(cursor, static_cast<char>('1' + event.lane));
            append_char(cursor, ' ');
            if (!append_u32(cursor, end, event.end_tick)) {
                return false;
            }
            append_text(cursor, " 128");
            break;
        case EventType::Slide:
            append_text(cursor, "SLIDE B");
            append_char(cursor, static_cast<char>('1' + event.lane));
            append_char(cursor, ' ');
            if (!append_u32(cursor, end, event.end_tick)) {
                return false;
            }
            append_text(cursor, " B");
            append_char(cursor, static_cast<char>('1' + event.target_lane));
            append_text(cursor, " 0");
            break;
        case EventType::TouchTap:
            append_text(cursor, "TOUCH ");
            append_touch_name(cursor, event.touch_sensor);
            append_text(cursor, " 192");
            break;
        case EventType::TouchHold:
            append_text(cursor, "TOUCH_HOLD ");
            append_touch_name(cursor, event.touch_sensor);
            append_char(cursor, ' ');
            if (!append_u32(cursor, end, event.end_tick)) {
                return false;
            }
            append_text(cursor, " 192");
            break;
    }

    event.word = event_word(std::string_view(
        buffer.data(), static_cast<std::size_t>(cursor - buffer.data())));
    return true;
}

inline void increment_count_fast(Counts& counts, EventType type) noexcept {
    switch (type) {
        case EventType::Tap:
            ++counts.tap;
            break;
        case EventType::Break:
            ++counts.break_count;
            break;
        case EventType::Hold:
            ++counts.hold;
            break;
        case EventType::Slide:
            ++counts.slide;
            break;
        case EventType::TouchTap:
        case EventType::TouchHold:
            ++counts.touch;
            break;
    }
    ++counts.total;
}

[[nodiscard]] bool parse_duration_fast(std::string_view token,
                                       std::uint64_t& duration) noexcept {
    const std::size_t colon = token.find(':');
    if (colon == std::string_view::npos ||
        token.find(':', colon + 1U) != std::string_view::npos) {
        return false;
    }
    std::uint32_t denominator = 0;
    std::uint32_t numerator = 0;
    if (!parse_uint_fast(token.substr(0U, colon), denominator) ||
        !parse_uint_fast(token.substr(colon + 1U), numerator) ||
        denominator == 0U || numerator == 0U ||
        std::gcd(denominator, numerator) != 1U) {
        return false;
    }
    const std::uint64_t scaled =
        static_cast<std::uint64_t>(kTicksPerWhole) * numerator;
    if (scaled % denominator != 0U) {
        return false;
    }
    duration = scaled / denominator;
    return true;
}

[[nodiscard]] bool parse_note_fast(FastParser& parser,
                                   std::string_view note,
                                   std::uint32_t byte_offset) {
    if (note.empty() || parser.tick > kMaxTimelineTick) {
        return false;
    }

    CompactEvent event;
    event.tick = static_cast<std::uint32_t>(parser.tick);
    event.location.byte_offset = byte_offset;

    const bool button = note.front() >= '1' && note.front() <= '8';
    std::size_t sensor_length = 0U;
    if (button) {
        event.lane = static_cast<std::uint8_t>(note.front() - '1');
        sensor_length = 1U;
    } else if (note.front() == 'C') {
        if (!touch_sensor_fast(note.substr(0U, 1U), event.touch_sensor)) {
            return false;
        }
        sensor_length = 1U;
    } else if (note.front() == 'A' || note.front() == 'B' ||
               note.front() == 'D' || note.front() == 'E') {
        if (note.size() < 2U ||
            !touch_sensor_fast(note.substr(0U, 2U), event.touch_sensor)) {
            return false;
        }
        sensor_length = 2U;
    } else {
        return false;
    }

    if (note.size() == sensor_length) {
        event.type = button ? EventType::Tap : EventType::TouchTap;
        event.strength = button ? 128U : 192U;
    } else if (button && note.size() == 2U && note[1] == 'b') {
        event.type = EventType::Break;
        event.strength = 255U;
    } else {
        const std::size_t open = note.find('[', sensor_length);
        if (open == std::string_view::npos || note.back() != ']' ||
            note.find('[', open + 1U) != std::string_view::npos ||
            note.find(']', open + 1U) != note.size() - 1U) {
            return false;
        }
        const std::string_view head =
            note.substr(sensor_length, open - sensor_length);
        std::uint64_t duration = 0;
        if (!parse_duration_fast(
                note.substr(open + 1U, note.size() - open - 2U), duration)) {
            return false;
        }
        const std::uint64_t end_tick = parser.tick + duration;
        if (end_tick > kMaxTimelineTick) {
            return false;
        }
        event.end_tick = static_cast<std::uint32_t>(end_tick);
        parser.chart.max_tick =
            std::max(parser.chart.max_tick, event.end_tick);

        if (head == "h") {
            event.type = button ? EventType::Hold : EventType::TouchHold;
            event.strength = button ? 128U : 192U;
        } else if (button && head.size() == 2U && head.front() == '-' &&
                   head.back() >= '1' && head.back() <= '8') {
            event.type = EventType::Slide;
            event.target_lane =
                static_cast<std::uint8_t>(head.back() - '1');
            if (event.target_lane == event.lane) {
                return false;
            }
        } else {
            return false;
        }
    }

    const bool touch =
        event.type == EventType::TouchTap ||
        event.type == EventType::TouchHold;
    const std::uint8_t source_index =
        touch ? event.touch_sensor : event.lane;
    SourceState& source_state =
        touch ? parser.touch_sources[source_index]
              : parser.button_sources[source_index];

    if ((source_state.have_tick && source_state.tick == event.tick) ||
        source_state.active_end > event.tick) {
        return false;
    }
    source_state.have_tick = true;
    source_state.tick = event.tick;
    if (event.type == EventType::Hold ||
        event.type == EventType::Slide ||
        event.type == EventType::TouchHold) {
        source_state.active_end =
            std::max(source_state.active_end, event.end_tick);
    }

    increment_count_fast(parser.chart.counts, event.type);
    ++parser.parsed_events;
    if (parser.parsed_events > kMaxEvents) {
        return false;
    }
    if (!set_event_word_fast(event)) {
        return false;
    }
    parser.chart.events.push_back(std::move(event));
    return true;
}

[[nodiscard]] bool parse_controls_fast(FastParser& parser,
                                       std::string_view cell,
                                       std::size_t& cursor) noexcept {
    while (cursor < cell.size() &&
           (cell[cursor] == '(' || cell[cursor] == '{')) {
        const char opening = cell[cursor];
        const char closing = opening == '(' ? ')' : '}';
        const std::size_t close = cell.find(closing, cursor + 1U);
        if (close == std::string_view::npos) {
            return false;
        }
        const std::string_view value =
            cell.substr(cursor + 1U, close - cursor - 1U);

        if (opening == '(') {
            std::uint32_t bpm = 0;
            if (!parse_bpm_fast(value, bpm)) {
                return false;
            }
            const bool conventional =
                parser.initial_bpm_control &&
                parser.initial_division_control &&
                parser.tick == 0U &&
                bpm == parser.chart.metadata.whole_bpm_milli &&
                cell.substr(close + 1U, 3U) == "{4}";
            if (bpm == parser.bpm_milli && !conventional) {
                ++parser.warning_count;
            }
            if (bpm < kMinBpmMilli || bpm > kMaxBpmMilli) {
                return false;
            }
            parser.bpm_milli = bpm;
            parser.initial_bpm_control = false;
            parser.conventional_division_pending = conventional;
        } else {
            std::uint32_t division = 0;
            if (!parse_uint_fast(value, division) ||
                division == 0U ||
                kTicksPerWhole % division != 0U) {
                return false;
            }
            const bool conventional =
                parser.conventional_division_pending && division == 4U;
            if (division == parser.division && !conventional) {
                ++parser.warning_count;
            }
            parser.division = division;
            parser.initial_division_control = false;
            parser.conventional_division_pending = false;
        }
        cursor = close + 1U;
    }
    return true;
}

[[nodiscard]] bool parse_cell_fast(FastParser& parser,
                                   std::string_view cell,
                                   bool has_comma,
                                   std::uint32_t cell_byte_offset) {
    if (cell == "E") {
        if (has_comma || parser.ended) {
            return false;
        }
        parser.ended = true;
        return true;
    }
    if (parser.ended) {
        return false;
    }

    std::size_t cursor = 0U;
    if (!parse_controls_fast(parser, cell, cursor)) {
        return false;
    }

    std::size_t note_count = 0U;
    if (cursor < cell.size()) {
        std::size_t note_begin = cursor;
        while (note_begin <= cell.size()) {
            const std::size_t slash = cell.find('/', note_begin);
            const std::size_t note_end =
                slash == std::string_view::npos ? cell.size() : slash;
            if (!parse_note_fast(
                    parser,
                    cell.substr(note_begin, note_end - note_begin),
                    static_cast<std::uint32_t>(
                        cell_byte_offset + note_begin))) {
                return false;
            }
            ++note_count;
            if (slash == std::string_view::npos) {
                break;
            }
            note_begin = slash + 1U;
        }
    }

    if (note_count >= kDenseEachThreshold) {
        ++parser.warning_count;
    }

    parser.tick += kTicksPerWhole / parser.division;
    if (parser.tick > kMaxTimelineTick) {
        return false;
    }
    parser.chart.max_tick =
        std::max(parser.chart.max_tick,
                 static_cast<std::uint32_t>(parser.tick));
    return true;
}

[[nodiscard]] bool parse_body_line_fast(
    FastParser& parser,
    std::string_view line,
    std::uint32_t line_byte_offset,
    bool last_line) {
    if (!last_line && (line.empty() || line.back() != ',')) {
        return false;
    }

    std::size_t begin = 0U;
    for (std::size_t i = 0U; i < line.size(); ++i) {
        if (line[i] == ' ') {
            return false;
        }
        if (line[i] != ',') {
            continue;
        }
        if (!parse_cell_fast(
                parser, line.substr(begin, i - begin), true,
                static_cast<std::uint32_t>(line_byte_offset + begin))) {
            return false;
        }
        begin = i + 1U;
    }

    if (begin < line.size() || (last_line && begin == 0U)) {
        if (!parse_cell_fast(
                parser, line.substr(begin), false,
                static_cast<std::uint32_t>(line_byte_offset + begin))) {
            return false;
        }
    } else if (last_line) {
        return false;
    }
    return true;
}

[[nodiscard]] bool fast_parse_valid(std::string_view text,
                                    FastChart& output) {
    if (text.empty() || text.size() > kMaxChartBytes ||
        text.back() != '\n') {
        return false;
    }

    std::array<std::string_view, 7> metadata{};
    std::size_t position = 0U;
    std::uint32_t line_number = 1U;

    for (std::size_t i = 0U; i < metadata.size(); ++i, ++line_number) {
        const std::size_t newline = text.find('\n', position);
        if (newline == std::string_view::npos) {
            return false;
        }
        const std::string_view line =
            text.substr(position, newline - position);
        if (line.empty() || line.size() > kMaxMetadataLineBytes ||
            line.back() == ' ') {
            return false;
        }
        for (const unsigned char byte : line) {
            if (byte < 0x20U || byte > 0x7eU) {
                return false;
            }
        }
        metadata[i] = line;
        position = newline + 1U;
    }

    if (position >= text.size() ||
        !metadata_value_fast(metadata[0], "&title=") ||
        !metadata_value_fast(metadata[1], "&artist=") ||
        !metadata_value_fast(metadata[2], "&des=") ||
        metadata[3] != "&first=0" ||
        !metadata[4].starts_with("&wholebpm=") ||
        !metadata_value_fast(metadata[5], "&lv_1=") ||
        metadata[6] != "&inote_1=") {
        return false;
    }

    std::uint32_t whole_bpm = 0U;
    if (!parse_bpm_fast(metadata[4].substr(10U), whole_bpm) ||
        whole_bpm < kMinBpmMilli || whole_bpm > kMaxBpmMilli) {
        return false;
    }

    FastParser parser;
    parser.chart.metadata.whole_bpm_milli = whole_bpm;
    parser.bpm_milli = whole_bpm;
    parser.chart.events.reserve(
        std::min<std::size_t>(kMaxEvents, text.size() / 12U + 8U));

    while (position < text.size()) {
        const std::size_t newline = text.find('\n', position);
        if (newline == std::string_view::npos) {
            return false;
        }
        const std::string_view line =
            text.substr(position, newline - position);
        if (line.empty() || line.size() > kMaxChartLineBytes ||
            line.back() == ' ') {
            return false;
        }
        for (const unsigned char byte : line) {
            if (byte < 0x20U || byte > 0x7eU) {
                return false;
            }
        }

        const bool last_line = newline + 1U == text.size();
        if (!parse_body_line_fast(
                parser, line, static_cast<std::uint32_t>(position),
                last_line)) {
            return false;
        }
        position = newline + 1U;
        ++line_number;
    }

    if (!parser.ended) {
        return false;
    }

    if (parser.chart.events.empty() &&
        parser.chart.counts.total == 0U) {
        ++parser.warning_count;
    }

    parser.chart.counts.warning = parser.warning_count;
    parser.chart.warning_count = parser.warning_count;

    std::sort(
        parser.chart.events.begin(), parser.chart.events.end(),
        [](const CompactEvent& lhs, const CompactEvent& rhs) {
            return std::tuple(
                       lhs.tick, static_cast<std::uint8_t>(lhs.type),
                       lhs.location.byte_offset) <
                   std::tuple(
                       rhs.tick, static_cast<std::uint8_t>(rhs.type),
                       rhs.location.byte_offset);
        });

    output = std::move(parser.chart);
    return true;
}


constexpr std::uint64_t kStateMixMultiplier = 0xd6e8feb86659fd93ULL;
constexpr std::uint64_t kStateMixIncrement = 0xa5a3564e27f8862fULL;

inline void add_button_energy_fast(
    State& state, std::uint8_t lane, std::uint32_t value) noexcept {
    if (lane >= state.energy.size()) {
        return;
    }
    state.energy[lane] = static_cast<std::uint16_t>(
        std::min<std::uint32_t>(
            65535U,
            static_cast<std::uint32_t>(state.energy[lane]) + value));
}

inline void add_touch_energy_fast(
    State& state, std::uint8_t sensor, std::uint32_t value) noexcept {
    if (sensor >= state.touch_energy.size()) {
        return;
    }
    state.touch_energy[sensor] = static_cast<std::uint16_t>(
        std::min<std::uint32_t>(
            65535U,
            static_cast<std::uint32_t>(state.touch_energy[sensor]) + value));
}

inline void apply_compact_event(
    State& state, const CompactEvent& event) noexcept {
    switch (event.type) {
        case EventType::Tap:
            if (event.lane < state.energy.size()) {
                add_button_energy_fast(state, event.lane, event.strength);
                state.taps[event.lane] = LatestTap{true, event.strength};
            }
            break;
        case EventType::Hold:
            if (event.lane < state.energy.size()) {
                add_button_energy_fast(
                    state, event.lane,
                    std::max<std::uint32_t>(1U, event.strength / 4U));
                state.holds[event.lane] =
                    ActiveHold{true, event.end_tick, event.strength};
            }
            break;
        case EventType::Slide:
            if (event.lane < state.energy.size() &&
                event.target_lane < state.energy.size()) {
                const std::uint32_t distance =
                    event.target_lane > event.lane
                        ? event.target_lane - event.lane
                        : event.lane - event.target_lane;
                add_button_energy_fast(
                    state, event.lane, 8U * (distance + 1U));
                state.slides[event.lane] = ActiveSlide{
                    true, event.end_tick, event.target_lane, event.path_id};
            }
            break;
        case EventType::Break:
            if (event.lane < state.energy.size()) {
                add_button_energy_fast(
                    state, event.lane, 2U * event.strength);
                state.breaks[event.lane] =
                    LatestBreak{true, event.strength, 3};
                state.break_remaining = 3;
            }
            break;
        case EventType::TouchTap:
            if (event.touch_sensor < state.touch_energy.size()) {
                add_touch_energy_fast(
                    state, event.touch_sensor, event.strength);
                state.touches[event.touch_sensor] =
                    LatestTouch{true, false, event.strength};
            }
            break;
        case EventType::TouchHold:
            if (event.touch_sensor < state.touch_energy.size()) {
                add_touch_energy_fast(
                    state, event.touch_sensor,
                    std::max<std::uint32_t>(1U, event.strength / 4U));
                state.touch_holds[event.touch_sensor] =
                    ActiveTouchHold{true, event.end_tick, event.strength};
                state.touches[event.touch_sensor] =
                    LatestTouch{true, true, event.strength};
            }
            break;
    }

    state.phase =
        std::rotl(state.phase ^ event.word, 13) *
            kStateMixMultiplier +
        kStateMixIncrement;
}

[[nodiscard]] EndpointStates evolve_compact_states(
    std::uint64_t chart_id, const FastChart& chart) {
    EndpointStates endpoints;
    State current;

    const std::uint64_t bpm_q16 =
        (static_cast<std::uint64_t>(chart.metadata.whole_bpm_milli) *
             65536ULL +
         500ULL) /
        1000ULL;
    current.phase =
        0x6a09e667f3bcc909ULL ^ chart_id ^ (bpm_q16 << 17U) ^
        static_cast<std::uint64_t>(kTouchSensorCount);
    current.warning_count = chart.warning_count;
    if (chart.warning_count != 0U) {
        current.latest_warning_severity = 1U;
    }

    std::size_t next_event = 0U;
    while (next_event < chart.events.size() &&
           chart.events[next_event].tick == 0U) {
        apply_compact_event(current, chart.events[next_event]);
        ++next_event;
    }
    endpoints.frame_begin = current;

    for (std::size_t frame = 0U; frame + 1U < kFrameCount; ++frame) {
        State next = current;
        const std::uint32_t target_tick =
            frame_tick(chart.max_tick, frame + 1U);

        for (std::size_t lane = 0U; lane < next.energy.size(); ++lane) {
            next.energy[lane] = static_cast<std::uint16_t>(
                (7U * static_cast<std::uint32_t>(next.energy[lane])) / 8U);

            if (current.holds[lane].present) {
                add_button_energy_fast(
                    next, static_cast<std::uint8_t>(lane),
                    std::max<std::uint32_t>(
                        1U, current.holds[lane].strength / 8U));
            }
            if (current.slides[lane].present) {
                const std::uint8_t target =
                    current.slides[lane].target_lane;
                const std::uint32_t distance =
                    target > lane
                        ? target - static_cast<std::uint32_t>(lane)
                        : static_cast<std::uint32_t>(lane) - target;
                add_button_energy_fast(
                    next, static_cast<std::uint8_t>(lane),
                    4U * (distance + 1U));
            }
        }

        for (std::size_t sensor = 0U;
             sensor < next.touch_energy.size(); ++sensor) {
            next.touch_energy[sensor] = static_cast<std::uint16_t>(
                (7U * static_cast<std::uint32_t>(
                          next.touch_energy[sensor])) /
                8U);
            if (current.touch_holds[sensor].present) {
                add_touch_energy_fast(
                    next, static_cast<std::uint8_t>(sensor),
                    std::max<std::uint32_t>(
                        1U,
                        current.touch_holds[sensor].strength / 8U));
            }
        }

        if (next.break_remaining > 0U) {
            --next.break_remaining;
        }
        for (LatestBreak& note_break : next.breaks) {
            if (note_break.flash_frames > 0U) {
                --note_break.flash_frames;
            }
        }
        for (std::size_t lane = 0U;
             lane < next.holds.size(); ++lane) {
            if (next.holds[lane].present &&
                next.holds[lane].end_tick <= target_tick) {
                next.holds[lane] = ActiveHold{};
            }
            if (next.slides[lane].present &&
                next.slides[lane].end_tick <= target_tick) {
                next.slides[lane] = ActiveSlide{};
            }
        }
        for (ActiveTouchHold& hold : next.touch_holds) {
            if (hold.present && hold.end_tick <= target_tick) {
                hold = ActiveTouchHold{};
            }
        }

        next.phase =
            std::rotl(
                next.phase ^ static_cast<std::uint64_t>(frame + 1U),
                7) *
                kStateMixMultiplier +
            kStateMixIncrement;

        while (next_event < chart.events.size() &&
               chart.events[next_event].tick <= target_tick) {
            apply_compact_event(next, chart.events[next_event]);
            ++next_event;
        }

        current = std::move(next);
        if (frame + 1U == kLastSampleFrame) {
            endpoints.frame_end = current;
        }
    }

    return endpoints;
}

}  // namespace

#if defined(__GNUC__)
__attribute__((hot))
#endif
void process_chart(
    std::uint64_t chart_id,
    std::string_view chart_text,
    std::span<std::uint8_t, kEncodedChartBytes> output) {
    FastChart fast_chart;
    ParsedChart parsed_chart;
    const bool fast = fast_parse_valid(chart_text, fast_chart);
    if (!fast) [[unlikely]] {
        parsed_chart = parse_maidata(chart_text);
    }

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

    const SemanticError* error =
        fast ? nullptr : select_semantic_error(parsed_chart);

    if (error != nullptr) [[unlikely]] {
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

        begin_digest = fast_frame_sha256(std::span<const std::uint8_t>(
            begin_payload.data(), begin_payload.size()));
        end_digest = fast_frame_sha256(std::span<const std::uint8_t>(
            end_payload.data(), end_payload.size()));
        result_digest = error_result_digest(
            chart_id,
            std::span<const std::uint8_t, 16>(error_fields),
            begin_digest, end_digest);
    } else [[likely]] {
        EndpointStates states;
        ParsedChart render_stub;
        const ParsedChart* render_chart = nullptr;

        if (fast) [[likely]] {
            counts = fast_chart.counts;
            states = evolve_compact_states(chart_id, fast_chart);
            render_stub.max_tick = fast_chart.max_tick;
            render_chart = &render_stub;
        } else {
            counts = parsed_chart.counts;
            states = evolve_endpoint_states(chart_id, parsed_chart);
            render_chart = &parsed_chart;
        }

        render_output_rgba_into(
            chart_id, *render_chart, states.frame_begin,
            kFirstSampleFrame, begin_payload);
        render_output_rgba_into(
            chart_id, *render_chart, states.frame_end,
            kLastSampleFrame, end_payload);

        std::array<std::uint8_t, kSerializedStateBytes> serialized_state{};
        serialize_state_into(states.frame_end, serialized_state);

        begin_digest = fast_frame_sha256(std::span<const std::uint8_t>(
            begin_payload.data(), begin_payload.size()));
        end_digest = fast_frame_sha256(std::span<const std::uint8_t>(
            end_payload.data(), end_payload.size()));
        result_digest = valid_result_digest(
            chart_id, counts, serialized_state,
            begin_digest, end_digest);
    }

    write_frame_header(
        output.data(), chart_id, FrameRole::Begin,
        status, begin_digest);
    write_frame_header(
        output.data() + kOutputFrameFileBytes,
        chart_id, FrameRole::End, status, end_digest);
    write_result(
        output.data() + 2U * kOutputFrameFileBytes,
        chart_id, status, counts, result_error_fields,
        begin_digest, end_digest, result_digest);
}

}  // namespace maimoe::kernel
