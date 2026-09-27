#pragma once

#include "kernel_operator.h"

namespace RsmSolution {

constexpr uint32_t kScoreTile = 4096;
constexpr uint32_t kMaxD = 256;
constexpr uint32_t kAlignment = 8;
// Reserve enough UB space to keep common complete segments resident. 16384
// FP32 x elements is 64 KiB (plus 32 KiB for the FP16 source), which still
// leaves room for score/reduction and moment scratch on 910B3.
constexpr uint32_t kXTileElements = 16384;
constexpr uint32_t kVarianceGroupRows = 8;

class ComputeCore {
public:
    __aicore__ inline void Init(
        GM_ADDR score, GM_ADDR x, GM_ADDR offsets, GM_ADDR mean, GM_ADDR rstd,
        GM_ADDR logsumexp, uint32_t n, uint32_t d, float epsilon)
    {
        d_ = d;
        epsilon_ = epsilon;
        score_gm_.SetGlobalBuffer(reinterpret_cast<__gm__ half *>(score), n);
        x_gm_.SetGlobalBuffer(reinterpret_cast<__gm__ half *>(x), n * d);
        offsets_gm_.SetGlobalBuffer(reinterpret_cast<__gm__ int32_t *>(offsets));
        mean_gm_.SetGlobalBuffer(reinterpret_cast<__gm__ float *>(mean));
        rstd_gm_.SetGlobalBuffer(reinterpret_cast<__gm__ float *>(rstd));
        logsumexp_gm_.SetGlobalBuffer(reinterpret_cast<__gm__ float *>(logsumexp));

        pipe_.InitBuffer(score_half_buffer_, kScoreTile * sizeof(half));
        pipe_.InitBuffer(score_float_buffer_, kScoreTile * sizeof(float));
        pipe_.InitBuffer(reduce_output_buffer_, 64 * sizeof(float));
        pipe_.InitBuffer(reduce_work_buffer_, kScoreTile * sizeof(float));
        pipe_.InitBuffer(x_half_buffer_, kXTileElements * sizeof(half));
        pipe_.InitBuffer(x_float_buffer_, kXTileElements * sizeof(float));
        pipe_.InitBuffer(mean_buffer_, kMaxD * sizeof(float));
        pipe_.InitBuffer(m2_buffer_, kMaxD * sizeof(float));
        pipe_.InitBuffer(mean_residual_buffer_, kMaxD * sizeof(float));
        pipe_.InitBuffer(variance_chunk_buffer_, kMaxD * sizeof(float));
        pipe_.InitBuffer(mean_correction_buffer_, kMaxD * sizeof(float));
        pipe_.InitBuffer(delta_buffer_, kMaxD * sizeof(float));
        pipe_.InitBuffer(temp_buffer_, kMaxD * sizeof(float));
    }

    __aicore__ inline void ProcessSegment(uint32_t segment)
    {
        const uint32_t begin = static_cast<uint32_t>(offsets_gm_.GetValue(segment));
        const uint32_t end = static_cast<uint32_t>(offsets_gm_.GetValue(segment + 1));
        const uint32_t rows = end - begin;

        if (rows == 1) {
            ProcessSingle(segment, begin);
            return;
        }

        // The common short/medium path keeps both weights and the complete x
        // segment in UB. x is fetched/cast once, then all numerically stable
        // mean-refinement and centered-variance passes operate locally.
        if (rows <= kScoreTile && rows * d_ <= kXTileElements) {
            const float maximum = PrepareWeights(begin, rows);
            LoadRows(begin, rows);
            const float normalizer = AccumulateLocalMoments(rows);
            FinalizeMoments(normalizer);
            StoreResults(segment, maximum, normalizer);
            return;
        }

        // If x does not fit but all scores do, still keep the normalized
        // weights resident. This removes two score loads, casts and Exp passes.
        if (rows <= kScoreTile) {
            const float maximum = PrepareWeights(begin, rows);
            const float normalizer = AccumulateCachedWeightMoments(begin, rows);
            FinalizeMoments(normalizer);
            StoreResults(segment, maximum, normalizer);
            return;
        }

        const float maximum = SegmentMaximum(begin, end);
        const float normalizer = AccumulateMoments(begin, end, maximum);
        FinalizeMoments(normalizer);
        StoreResults(segment, maximum, normalizer);
    }

private:
    __aicore__ inline uint32_t Minimum(uint32_t lhs, uint32_t rhs) const
    {
        return lhs < rhs ? lhs : rhs;
    }

    __aicore__ inline uint32_t RowsPerXTile() const
    {
        // Keep streaming tile boundaries aligned with the baseline's 8-row
        // variance accumulation groups so the FP32 summation order remains
        // conservative on adversarial precision cases.
        const uint32_t rows = kXTileElements / d_;
        return (rows / kVarianceGroupRows) * kVarianceGroupRows;
    }

    __aicore__ inline void LoadScores(
        uint32_t begin, uint32_t count, float maximum, bool exponentiate)
    {
        auto score_half = score_half_buffer_.Get<half>();
        const uint16_t bytes = static_cast<uint16_t>(count * sizeof(half));
        AscendC::DataCopyPad(score_half, score_gm_[begin], {1, bytes, 0, 0}, {});
        AscendC::PipeBarrier<PIPE_ALL>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Cast(score_float, score_half, AscendC::RoundMode::CAST_NONE, count);
        if (exponentiate) {
            AscendC::Adds(score_float, score_float, -maximum, count);
            AscendC::Exp(score_float, score_float, count);
        }
        AscendC::PipeBarrier<PIPE_V>();
    }

    __aicore__ inline void LoadRows(uint32_t begin, uint32_t rows)
    {
        const uint32_t elements = rows * d_;
        const uint16_t bytes = static_cast<uint16_t>(elements * sizeof(half));
        auto x_half = x_half_buffer_.Get<half>();
        AscendC::DataCopyPad(x_half, x_gm_[begin * d_], {1, bytes, 0, 0}, {});
        AscendC::PipeBarrier<PIPE_ALL>();
        AscendC::Cast(
            x_float_buffer_.Get<float>(), x_half,
            AscendC::RoundMode::CAST_NONE, elements);
        AscendC::PipeBarrier<PIPE_V>();
    }

    __aicore__ inline float PrepareWeights(uint32_t begin, uint32_t count)
    {
        LoadScores(begin, count, 0.0f, false);
        auto score_float = score_float_buffer_.Get<float>();
        auto reduce_output = reduce_output_buffer_.Get<float>();
        auto reduce_work = reduce_work_buffer_.Get<float>();
        AscendC::ReduceMax(reduce_output, score_float, reduce_work, count, false);
        AscendC::PipeBarrier<PIPE_V>();
        const float maximum = reduce_output.GetValue(0);
        AscendC::Adds(score_float, score_float, -maximum, count);
        AscendC::Exp(score_float, score_float, count);
        AscendC::PipeBarrier<PIPE_V>();
        return maximum;
    }

    __aicore__ inline float SegmentMaximum(uint32_t begin, uint32_t end)
    {
        auto score_float = score_float_buffer_.Get<float>();
        auto reduce_output = reduce_output_buffer_.Get<float>();
        auto reduce_work = reduce_work_buffer_.Get<float>();
        float maximum = -65504.0f;
        for (uint32_t base = begin; base < end; base += kScoreTile) {
            const uint32_t count = Minimum(kScoreTile, end - base);
            LoadScores(base, count, 0.0f, false);
            AscendC::ReduceMax(reduce_output, score_float, reduce_work, count, false);
            AscendC::PipeBarrier<PIPE_V>();
            const float tile_maximum = reduce_output.GetValue(0);
            maximum = maximum > tile_maximum ? maximum : tile_maximum;
        }
        return maximum;
    }

    __aicore__ inline void ProcessSingle(uint32_t segment, uint32_t row)
    {
        LoadScores(row, 1, 0.0f, false);
        LoadRows(row, 1);

        auto mean = mean_buffer_.Get<float>();
        auto rstd = m2_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        AscendC::DataCopy(mean, x_float, d_);
        AscendC::Duplicate(rstd, epsilon_, d_);
        AscendC::Rsqrt(rstd, rstd, d_);
        AscendC::PipeBarrier<PIPE_ALL>();

        const uint16_t bytes = static_cast<uint16_t>(d_ * sizeof(float));
        AscendC::DataCopyPad(mean_gm_[segment * d_], mean, {1, bytes, 0, 0});
        AscendC::DataCopyPad(rstd_gm_[segment * d_], rstd, {1, bytes, 0, 0});

        auto temp = temp_buffer_.Get<float>();
        temp.SetValue(0, score_float_buffer_.Get<float>().GetValue(0));
        AscendC::PipeBarrier<PIPE_ALL>();
        AscendC::DataCopyPad(logsumexp_gm_[segment], temp, {1, sizeof(float), 0, 0});
        AscendC::PipeBarrier<PIPE_ALL>();
    }

    // Complete-segment UB path. The arithmetic intentionally mirrors the
    // baseline compensated mean + residual refinement + centered variance,
    // but all x accesses after the first load are local.
    __aicore__ inline float AccumulateLocalMoments(uint32_t rows)
    {
        auto sum = mean_buffer_.Get<float>();
        auto correction = mean_correction_buffer_.Get<float>();
        auto term = temp_buffer_.Get<float>();
        auto delta = delta_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(sum, 0.0f, d_);
        AscendC::Duplicate(correction, 0.0f, d_);
        AscendC::PipeBarrier<PIPE_V>();

        float weight_total = 0.0f;
        float weight_correction = 0.0f;
        for (uint32_t row = 0; row < rows; ++row) {
            const float weight = score_float.GetValue(row);
            const float y = weight - weight_correction;
            const float next = weight_total + y;
            weight_correction = (next - weight_total) - y;
            weight_total = next;

            AscendC::Muls(term, x_float[row * d_], weight, d_);
            AscendC::Sub(delta, term, correction, d_);
            AscendC::Add(term, sum, delta, d_);
            AscendC::Sub(correction, term, sum, d_);
            AscendC::Sub(correction, correction, delta, d_);
            AscendC::DataCopy(sum, term, d_);
        }
        AscendC::PipeBarrier<PIPE_V>();
        AscendC::Muls(sum, sum, 1.0f / weight_total, d_);
        AscendC::PipeBarrier<PIPE_V>();

        RefineLocalMean(rows, weight_total);
        AccumulateLocalVariance(rows);
        return weight_total;
    }

    __aicore__ inline void RefineLocalMean(uint32_t rows, float normalizer)
    {
        auto mean = mean_buffer_.Get<float>();
        auto residual = mean_residual_buffer_.Get<float>();
        auto correction = mean_correction_buffer_.Get<float>();
        auto term = temp_buffer_.Get<float>();
        auto delta = delta_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(residual, 0.0f, d_);
        AscendC::Duplicate(correction, 0.0f, d_);
        for (uint32_t row = 0; row < rows; ++row) {
            AscendC::Sub(term, x_float[row * d_], mean, d_);
            AscendC::Muls(term, term, score_float.GetValue(row), d_);
            AscendC::Sub(delta, term, correction, d_);
            AscendC::Add(term, residual, delta, d_);
            AscendC::Sub(correction, term, residual, d_);
            AscendC::Sub(correction, correction, delta, d_);
            AscendC::DataCopy(residual, term, d_);
        }
        AscendC::Muls(residual, residual, 1.0f / normalizer, d_);
        AscendC::Add(mean, mean, residual, d_);
        AscendC::PipeBarrier<PIPE_V>();
    }

    __aicore__ inline void AccumulateLocalVariance(uint32_t rows)
    {
        auto mean = mean_buffer_.Get<float>();
        auto m2 = m2_buffer_.Get<float>();
        auto temp = temp_buffer_.Get<float>();
        auto chunk = variance_chunk_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(m2, 0.0f, d_);
        for (uint32_t row_base = 0; row_base < rows; row_base += kVarianceGroupRows) {
            const uint32_t group_rows = Minimum(kVarianceGroupRows, rows - row_base);
            AscendC::Duplicate(chunk, 0.0f, d_);
            for (uint32_t row = 0; row < group_rows; ++row) {
                const uint32_t index = row_base + row;
                AscendC::Sub(temp, x_float[index * d_], mean, d_);
                AscendC::Mul(temp, temp, temp, d_);
                AscendC::Muls(temp, temp, score_float.GetValue(index), d_);
                AscendC::Add(chunk, chunk, temp, d_);
            }
            AscendC::Add(m2, m2, chunk, d_);
        }
        AscendC::PipeBarrier<PIPE_V>();
    }

    // Score-resident path for segments whose x matrix does not fit in UB.
    __aicore__ inline float AccumulateCachedWeightMoments(
        uint32_t begin, uint32_t rows)
    {
        auto sum = mean_buffer_.Get<float>();
        auto correction = mean_correction_buffer_.Get<float>();
        auto term = temp_buffer_.Get<float>();
        auto delta = delta_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(sum, 0.0f, d_);
        AscendC::Duplicate(correction, 0.0f, d_);
        AscendC::PipeBarrier<PIPE_V>();

        float weight_total = 0.0f;
        float weight_correction = 0.0f;
        const uint32_t rows_per_tile = RowsPerXTile();
        for (uint32_t row_base = 0; row_base < rows; row_base += rows_per_tile) {
            const uint32_t tile_rows = Minimum(rows_per_tile, rows - row_base);
            LoadRows(begin + row_base, tile_rows);
            for (uint32_t row = 0; row < tile_rows; ++row) {
                const float weight = score_float.GetValue(row_base + row);
                const float y = weight - weight_correction;
                const float next = weight_total + y;
                weight_correction = (next - weight_total) - y;
                weight_total = next;

                AscendC::Muls(term, x_float[row * d_], weight, d_);
                AscendC::Sub(delta, term, correction, d_);
                AscendC::Add(term, sum, delta, d_);
                AscendC::Sub(correction, term, sum, d_);
                AscendC::Sub(correction, correction, delta, d_);
                AscendC::DataCopy(sum, term, d_);
            }
            AscendC::PipeBarrier<PIPE_V>();
        }

        AscendC::Muls(sum, sum, 1.0f / weight_total, d_);
        AscendC::PipeBarrier<PIPE_V>();
        RefineMeanCachedWeights(begin, rows, weight_total);
        AccumulateCenteredVarianceCachedWeights(begin, rows);
        return weight_total;
    }

    __aicore__ inline void RefineMeanCachedWeights(
        uint32_t begin, uint32_t rows, float normalizer)
    {
        auto mean = mean_buffer_.Get<float>();
        auto residual = mean_residual_buffer_.Get<float>();
        auto correction = mean_correction_buffer_.Get<float>();
        auto term = temp_buffer_.Get<float>();
        auto delta = delta_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(residual, 0.0f, d_);
        AscendC::Duplicate(correction, 0.0f, d_);

        const uint32_t rows_per_tile = RowsPerXTile();
        for (uint32_t row_base = 0; row_base < rows; row_base += rows_per_tile) {
            const uint32_t tile_rows = Minimum(rows_per_tile, rows - row_base);
            LoadRows(begin + row_base, tile_rows);
            for (uint32_t row = 0; row < tile_rows; ++row) {
                AscendC::Sub(term, x_float[row * d_], mean, d_);
                AscendC::Muls(term, term, score_float.GetValue(row_base + row), d_);
                AscendC::Sub(delta, term, correction, d_);
                AscendC::Add(term, residual, delta, d_);
                AscendC::Sub(correction, term, residual, d_);
                AscendC::Sub(correction, correction, delta, d_);
                AscendC::DataCopy(residual, term, d_);
            }
            AscendC::PipeBarrier<PIPE_V>();
        }

        AscendC::Muls(residual, residual, 1.0f / normalizer, d_);
        AscendC::Add(mean, mean, residual, d_);
        AscendC::PipeBarrier<PIPE_V>();
    }

    __aicore__ inline void AccumulateCenteredVarianceCachedWeights(
        uint32_t begin, uint32_t rows)
    {
        auto mean = mean_buffer_.Get<float>();
        auto m2 = m2_buffer_.Get<float>();
        auto temp = temp_buffer_.Get<float>();
        auto chunk = variance_chunk_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(m2, 0.0f, d_);

        const uint32_t rows_per_tile = RowsPerXTile();
        for (uint32_t row_base = 0; row_base < rows; row_base += rows_per_tile) {
            const uint32_t tile_rows = Minimum(rows_per_tile, rows - row_base);
            LoadRows(begin + row_base, tile_rows);
            for (uint32_t group_base = 0; group_base < tile_rows;
                 group_base += kVarianceGroupRows) {
                const uint32_t group_rows =
                    Minimum(kVarianceGroupRows, tile_rows - group_base);
                AscendC::Duplicate(chunk, 0.0f, d_);
                for (uint32_t row = 0; row < group_rows; ++row) {
                    const uint32_t local_row = group_base + row;
                    AscendC::Sub(temp, x_float[local_row * d_], mean, d_);
                    AscendC::Mul(temp, temp, temp, d_);
                    AscendC::Muls(temp, temp,
                        score_float.GetValue(row_base + local_row), d_);
                    AscendC::Add(chunk, chunk, temp, d_);
                }
                AscendC::Add(m2, m2, chunk, d_);
            }
            AscendC::PipeBarrier<PIPE_V>();
        }
    }

    // Generic streaming path for segments longer than the score buffer.
    __aicore__ inline float AccumulateMoments(
        uint32_t begin, uint32_t end, float maximum)
    {
        auto sum = mean_buffer_.Get<float>();
        auto correction = mean_correction_buffer_.Get<float>();
        auto term = temp_buffer_.Get<float>();
        auto delta = delta_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(sum, 0.0f, d_);
        AscendC::Duplicate(correction, 0.0f, d_);
        AscendC::PipeBarrier<PIPE_V>();

        float weight_total = 0.0f;
        float weight_correction = 0.0f;
        const uint32_t rows_per_tile = RowsPerXTile();
        for (uint32_t base = begin; base < end; base += kScoreTile) {
            const uint32_t count = Minimum(kScoreTile, end - base);
            LoadScores(base, count, maximum, true);
            for (uint32_t row_base = 0; row_base < count; row_base += rows_per_tile) {
                const uint32_t rows = Minimum(rows_per_tile, count - row_base);
                LoadRows(base + row_base, rows);
                for (uint32_t row = 0; row < rows; ++row) {
                    const float weight = score_float.GetValue(row_base + row);
                    const float y = weight - weight_correction;
                    const float next = weight_total + y;
                    weight_correction = (next - weight_total) - y;
                    weight_total = next;

                    AscendC::Muls(term, x_float[row * d_], weight, d_);
                    AscendC::Sub(delta, term, correction, d_);
                    AscendC::Add(term, sum, delta, d_);
                    AscendC::Sub(correction, term, sum, d_);
                    AscendC::Sub(correction, correction, delta, d_);
                    AscendC::DataCopy(sum, term, d_);
                }
                AscendC::PipeBarrier<PIPE_V>();
            }
        }
        AscendC::PipeBarrier<PIPE_V>();
        AscendC::Muls(sum, sum, 1.0f / weight_total, d_);
        AscendC::PipeBarrier<PIPE_V>();
        RefineMean(begin, end, maximum, weight_total);
        AccumulateCenteredVariance(begin, end, maximum);
        return weight_total;
    }

    __aicore__ inline void RefineMean(
        uint32_t begin, uint32_t end, float maximum, float normalizer)
    {
        auto mean = mean_buffer_.Get<float>();
        auto residual = mean_residual_buffer_.Get<float>();
        auto correction = mean_correction_buffer_.Get<float>();
        auto term = temp_buffer_.Get<float>();
        auto delta = delta_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(residual, 0.0f, d_);
        AscendC::Duplicate(correction, 0.0f, d_);

        const uint32_t rows_per_tile = RowsPerXTile();
        for (uint32_t base = begin; base < end; base += kScoreTile) {
            const uint32_t count = Minimum(kScoreTile, end - base);
            LoadScores(base, count, maximum, true);
            for (uint32_t row_base = 0; row_base < count; row_base += rows_per_tile) {
                const uint32_t rows = Minimum(rows_per_tile, count - row_base);
                LoadRows(base + row_base, rows);
                for (uint32_t row = 0; row < rows; ++row) {
                    AscendC::Sub(term, x_float[row * d_], mean, d_);
                    AscendC::Muls(term, term, score_float.GetValue(row_base + row), d_);
                    AscendC::Sub(delta, term, correction, d_);
                    AscendC::Add(term, residual, delta, d_);
                    AscendC::Sub(correction, term, residual, d_);
                    AscendC::Sub(correction, correction, delta, d_);
                    AscendC::DataCopy(residual, term, d_);
                }
                AscendC::PipeBarrier<PIPE_V>();
            }
        }
        AscendC::Muls(residual, residual, 1.0f / normalizer, d_);
        AscendC::Add(mean, mean, residual, d_);
        AscendC::PipeBarrier<PIPE_V>();
    }

    __aicore__ inline void AccumulateCenteredVariance(
        uint32_t begin, uint32_t end, float maximum)
    {
        auto mean = mean_buffer_.Get<float>();
        auto m2 = m2_buffer_.Get<float>();
        auto temp = temp_buffer_.Get<float>();
        auto chunk = variance_chunk_buffer_.Get<float>();
        auto x_float = x_float_buffer_.Get<float>();
        auto score_float = score_float_buffer_.Get<float>();
        AscendC::Duplicate(m2, 0.0f, d_);

        const uint32_t rows_per_tile = RowsPerXTile();
        for (uint32_t base = begin; base < end; base += kScoreTile) {
            const uint32_t count = Minimum(kScoreTile, end - base);
            LoadScores(base, count, maximum, true);
            for (uint32_t row_base = 0; row_base < count; row_base += rows_per_tile) {
                const uint32_t rows = Minimum(rows_per_tile, count - row_base);
                LoadRows(base + row_base, rows);
                for (uint32_t group_base = 0; group_base < rows;
                     group_base += kVarianceGroupRows) {
                    const uint32_t group_rows =
                        Minimum(kVarianceGroupRows, rows - group_base);
                    AscendC::Duplicate(chunk, 0.0f, d_);
                    for (uint32_t row = 0; row < group_rows; ++row) {
                        const uint32_t local_row = group_base + row;
                        AscendC::Sub(temp, x_float[local_row * d_], mean, d_);
                        AscendC::Mul(temp, temp, temp, d_);
                        AscendC::Muls(temp, temp,
                            score_float.GetValue(row_base + local_row), d_);
                        AscendC::Add(chunk, chunk, temp, d_);
                    }
                    AscendC::Add(m2, m2, chunk, d_);
                }
                AscendC::PipeBarrier<PIPE_V>();
            }
        }
    }

    __aicore__ inline void FinalizeMoments(float normalizer)
    {
        auto m2 = m2_buffer_.Get<float>();
        AscendC::Muls(m2, m2, 1.0f / normalizer, d_);
        AscendC::Maxs(m2, m2, 0.0f, d_);
        AscendC::Adds(m2, m2, epsilon_, d_);
        AscendC::Rsqrt(m2, m2, d_);
        AscendC::PipeBarrier<PIPE_ALL>();
    }

    __aicore__ inline void StoreResults(
        uint32_t segment, float maximum, float normalizer)
    {
        auto mean = mean_buffer_.Get<float>();
        auto m2 = m2_buffer_.Get<float>();
        auto temp = temp_buffer_.Get<float>();
        const uint16_t bytes = static_cast<uint16_t>(d_ * sizeof(float));
        AscendC::DataCopyPad(mean_gm_[segment * d_], mean, {1, bytes, 0, 0});
        AscendC::DataCopyPad(rstd_gm_[segment * d_], m2, {1, bytes, 0, 0});

        auto scalar = reduce_output_buffer_.Get<float>();
        AscendC::Duplicate(scalar, normalizer, kAlignment);
        AscendC::Ln(scalar, scalar, kAlignment);
        AscendC::PipeBarrier<PIPE_V>();
        temp.SetValue(0, scalar.GetValue(0) + maximum);
        AscendC::PipeBarrier<PIPE_ALL>();
        AscendC::DataCopyPad(
            logsumexp_gm_[segment], temp, {1, sizeof(float), 0, 0});
        AscendC::PipeBarrier<PIPE_ALL>();
    }

    AscendC::TPipe pipe_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> score_half_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> score_float_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> reduce_output_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> reduce_work_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> x_half_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> x_float_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> mean_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> m2_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> mean_residual_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> variance_chunk_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> mean_correction_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> delta_buffer_;
    AscendC::TBuf<AscendC::TPosition::VECCALC> temp_buffer_;
    AscendC::GlobalTensor<half> score_gm_;
    AscendC::GlobalTensor<half> x_gm_;
    AscendC::GlobalTensor<int32_t> offsets_gm_;
    AscendC::GlobalTensor<float> mean_gm_;
    AscendC::GlobalTensor<float> rstd_gm_;
    AscendC::GlobalTensor<float> logsumexp_gm_;
    uint32_t d_ = 0;
    float epsilon_ = 0.0f;
};

}  // namespace RsmSolution
