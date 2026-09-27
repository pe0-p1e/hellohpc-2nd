#pragma once

#include <cstdint>

constexpr uint32_t kRsmMaxBlocks = 40;

struct RsmSolutionTiling {
    uint32_t n;
    uint32_t d;
    uint32_t s;
    float epsilon;

    // mode 0: host-balanced whole segments
    // mode 1: low-S, wide-D column sharding
    uint32_t mode;
    uint32_t block_count;

    uint32_t shard_cols;
    uint32_t shards_per_segment;
    uint32_t task_count;

    uint32_t segment_begin[40];
    uint32_t segment_end[40];
};
