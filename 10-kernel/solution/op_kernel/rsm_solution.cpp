#include "rsm_solution_tiling.h"
#include "rsm_solution_core.h"

extern "C" __global__ __aicore__ void rsm_solution(
    GM_ADDR score, GM_ADDR x, GM_ADDR offsets, GM_ADDR mean, GM_ADDR rstd,
    GM_ADDR logsumexp, GM_ADDR workspace, RsmSolutionTiling tiling)
{
    (void)workspace;
    const uint32_t block = AscendC::GetBlockIdx();
    if (block >= tiling.block_count) {
        return;
    }

    RsmSolution::ComputeCore core;
    core.Init(score, x, offsets, mean, rstd, logsumexp,
        tiling.n, tiling.d, tiling.epsilon);

    for (uint32_t segment = tiling.segment_begin[block];
         segment < tiling.segment_end[block]; ++segment) {
        core.ProcessSegment(segment);
    }
}
