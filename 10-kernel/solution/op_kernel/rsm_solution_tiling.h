#pragma once

#include <cstdint>

constexpr uint32_t kRsmMaxBlocks = 40;

struct RsmSolutionTiling {
    uint32_t n;
    uint32_t d;
    uint32_t s;
    float epsilon;
    uint32_t block_count;
    uint32_t segment_begin[kRsmMaxBlocks];
    uint32_t segment_end[kRsmMaxBlocks];
};
