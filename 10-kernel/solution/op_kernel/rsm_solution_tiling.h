#pragma once

#include <cstdint>

// Keep a named constant for host-side planning. The generated AscendC host
// stub copies only the tiling struct declaration, so struct array extents must
// be literal constants rather than referring to this symbol.
constexpr uint32_t kRsmMaxBlocks = 40;

struct RsmSolutionTiling {
    uint32_t n;
    uint32_t d;
    uint32_t s;
    float epsilon;
    uint32_t block_count;
    uint32_t segment_begin[40];
    uint32_t segment_end[40];
};
