#include "rsm_solution_host.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace {
constexpr uint32_t kVectorCores = 40;
constexpr uint32_t kShardCols = 16;

uint64_t SegmentCost(uint32_t rows, uint32_t d)
{
    return static_cast<uint64_t>(rows) * (static_cast<uint64_t>(d) + 8U)
        + static_cast<uint64_t>(d) * 4U + 64U;
}
}

void ConfigureSolutionLaunch(uint32_t n, uint32_t d, uint32_t s, float epsilon,
    const int32_t *offsets, RsmSolutionTiling *tiling,
    uint32_t *block_dim, uint32_t *tiling_key)
{
    tiling->n = n;
    tiling->d = d;
    tiling->s = s;
    tiling->epsilon = epsilon;
    tiling->mode = 0;
    tiling->block_count = 0;
    tiling->shard_cols = 0;
    tiling->shards_per_segment = 0;
    tiling->task_count = 0;

    for (uint32_t i = 0; i < kRsmMaxBlocks; ++i) {
        tiling->segment_begin[i] = 0;
        tiling->segment_end[i] = 0;
    }

    // The two public families with pathological segment-level parallelism are
    // both low-S and wide-D. Column sharding is generic (shape-driven only):
    // every shard independently recomputes the small score reduction, while
    // expensive x/moment work is split across vector cores with disjoint output.
    if (s <= 4 && d >= 192) {
        tiling->mode = 1;
        tiling->shard_cols = kShardCols;
        tiling->shards_per_segment = (d + kShardCols - 1U) / kShardCols;
        tiling->task_count = s * tiling->shards_per_segment;
        tiling->block_count = std::min(kVectorCores, tiling->task_count);
        *block_dim = tiling->block_count;
        *tiling_key = 3;
        return;
    }

    const uint32_t blocks = std::min(s, kVectorCores);
    tiling->block_count = blocks;
    *block_dim = blocks;
    *tiling_key = 2;

    std::vector<uint64_t> prefix(static_cast<size_t>(s) + 1U, 0U);
    for (uint32_t segment = 0; segment < s; ++segment) {
        const uint32_t rows = static_cast<uint32_t>(
            offsets[segment + 1] - offsets[segment]);
        prefix[segment + 1] = prefix[segment] + SegmentCost(rows, d);
    }

    uint32_t start = 0;
    for (uint32_t block = 0; block < blocks; ++block) {
        tiling->segment_begin[block] = start;
        const uint32_t remaining_blocks = blocks - block;
        if (remaining_blocks == 1) {
            tiling->segment_end[block] = s;
            break;
        }

        const uint64_t remaining_cost = prefix[s] - prefix[start];
        const uint64_t target = remaining_cost / remaining_blocks;
        const uint32_t latest_end = s - (remaining_blocks - 1U);

        uint32_t end = start + 1U;
        while (end < latest_end && prefix[end] - prefix[start] < target) {
            ++end;
        }

        if (end > start + 1U) {
            const uint64_t after = prefix[end] - prefix[start];
            const uint64_t before = prefix[end - 1U] - prefix[start];
            const uint64_t after_error = after > target ? after - target : target - after;
            const uint64_t before_error = before > target ? before - target : target - before;
            if (before_error <= after_error) {
                --end;
            }
        }

        tiling->segment_end[block] = end;
        start = end;
    }
}
