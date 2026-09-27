#!/usr/bin/env bash
set -euo pipefail

: "${HELLOHPC_MODEL_PATH:?HELLOHPC_MODEL_PATH is set by the platform}"
: "${HELLOHPC_SERVICE_BASE_URL:?HELLOHPC_SERVICE_BASE_URL is set by the platform}"
: "${HELLOHPC_MODEL_ID:?HELLOHPC_MODEL_ID is set by the platform}"
: "${HELLOHPC_SERVICE_LOG:?HELLOHPC_SERVICE_LOG is set by the platform}"
: "${HELLOHPC_ASSIGNED_NPU_COUNT:?HELLOHPC_ASSIGNED_NPU_COUNT is set by the platform}"
: "${HELLOHPC_CASE_ID:?HELLOHPC_CASE_ID is set by the platform}"
: "${HELLOHPC_STAGE:?HELLOHPC_STAGE is set by the platform}"

if [[ -f /usr/local/Ascend/ascend-toolkit/set_env.sh ]]; then
    # shellcheck disable=SC1091
    source /usr/local/Ascend/ascend-toolkit/set_env.sh
fi

read -r HOST PORT < <(
    python3 - "${HELLOHPC_SERVICE_BASE_URL}" <<'PY'
import sys
from urllib.parse import urlparse
u = urlparse(sys.argv[1])
if u.scheme != "http" or not u.hostname or u.port is None:
    raise SystemExit(f"invalid HELLOHPC_SERVICE_BASE_URL: {sys.argv[1]!r}")
print(u.hostname, u.port)
PY
)

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export PYTORCH_NPU_ALLOC_CONF=expandable_segments:True
export TASK_QUEUE_ENABLE=1
export HCCL_OP_EXPANSION_MODE=AIV
export OMP_PROC_BIND=false

if [[ -r /usr/lib/aarch64-linux-gnu/libjemalloc.so.2 ]]; then
    export LD_PRELOAD="/usr/lib/aarch64-linux-gnu/libjemalloc.so.2${LD_PRELOAD:+:${LD_PRELOAD}}"
fi

MODEL="${HELLOHPC_MODEL_PATH}"
MODEL_ID="${HELLOHPC_MODEL_ID}"
LOG="${HELLOHPC_SERVICE_LOG}"
NPU_COUNT="${HELLOHPC_ASSIGNED_NPU_COUNT}"

case "${HELLOHPC_CASE_ID}" in
    qwen-performance)
        [[ "${NPU_COUNT}" == "1" ]] || {
            echo "qwen-performance expects exactly 1 NPU, got ${NPU_COUNT}" >&2
            exit 2
        }

        CMD=(
            vllm serve "${MODEL}"
            --host "${HOST}"
            --port "${PORT}"
            --served-model-name "${MODEL_ID}"
            --tensor-parallel-size 1
            --distributed-executor-backend mp
            --dtype bfloat16
            --trust-remote-code
            --max-model-len 10240
            --max-num-batched-tokens 12288
            --max-num-seqs 1
            --gpu-memory-utilization 0.92
            --no-enable-prefix-caching
            --compilation-config '{"cudagraph_mode":"FULL_DECODE_ONLY","cudagraph_capture_sizes":[1]}'
        )
        ;;

    dpsk-stage2)
        [[ "${NPU_COUNT}" == "8" ]] || {
            echo "dpsk-stage2 expects exactly 8 NPUs, got ${NPU_COUNT}" >&2
            exit 2
        }

        export VLLM_PREFIX_CACHE_RETENTION_INTERVAL=4096
        export OMP_NUM_THREADS=16

        CMD=(
            vllm serve "${MODEL}"
            --host "${HOST}"
            --port "${PORT}"
            --served-model-name "${MODEL_ID}"
            --data-parallel-size 1
            --tensor-parallel-size 8
            --distributed-executor-backend mp
            --enable-expert-parallel
            --tokenizer-mode deepseek_v4
            --quantization ascend
            --max-model-len 102400
            --max-num-batched-tokens 8192
            --max-num-seqs 32
            --gpu-memory-utilization 0.93
            --enable-prefix-caching
            --block-size 32
            --attention_config.indexer_kv_dtype int8
            --model-loader-extra-config '{"enable_multithread_load":true,"num_threads":64}'
            --speculative-config '{"method":"dspark","num_speculative_tokens":7,"enforce_eager":true}'
            --compilation-config '{"cudagraph_mode":"FULL_DECODE_ONLY"}'
            --additional-config '{"ascend_compilation_config":{"enable_npugraph_ex":true,"enable_static_kernel":false},"enable_cpu_binding":true,"enable_dsa_cp":true,"multistream_overlap_shared_expert":true}'
        )
        ;;

    *)
        echo "unsupported HELLOHPC_CASE_ID=${HELLOHPC_CASE_ID}" >&2
        exit 2
        ;;
esac

mkdir -p "$(dirname "${LOG}")"
{
    echo "stage=${HELLOHPC_STAGE}"
    echo "case=${HELLOHPC_CASE_ID}"
    echo "model=${MODEL_ID}"
    echo "npu_count=${NPU_COUNT}"
    printf 'launch:'
    printf ' %q' "${CMD[@]}"
    printf '\n'
} >> "${LOG}"

nohup "${CMD[@]}" >> "${LOG}" 2>&1 </dev/null &
echo "vllm pid=$!" >> "${LOG}"
exit 0
