#!/usr/bin/env python3
"""Run independent collision searches, using idle cores to hedge stragglers."""

import os
import secrets
import subprocess
import sys
import time


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <tasks.txt>", file=sys.stderr)
        return 2

    binary = os.path.join(os.path.dirname(os.path.abspath(__file__)), "md5fastcoll")
    if not os.access(binary, os.X_OK):
        print(f"{binary} not found or not executable", file=sys.stderr)
        return 1

    with open(sys.argv[1], encoding="utf-8") as task_file:
        rows = [line.split() for line in task_file if line.strip()]
    if any(len(row) != 3 for row in rows):
        print("each task must have three paths", file=sys.stderr)
        return 1

    tasks = [dict(paths=row, attempts=[], done=False, started=None) for row in rows]
    workers = min(32, len(os.sched_getaffinity(0)))
    active = []
    next_task = 0
    completed = 0
    serial = 0

    def start(task):
        nonlocal serial
        serial += 1
        prefix, output1, output2 = task["paths"]
        temporary1 = f"{output1}.hedge-{os.getpid()}-{serial}"
        temporary2 = f"{output2}.hedge-{os.getpid()}-{serial}"
        command = [binary, "-q", "--seed1", str(secrets.randbits(32)),
                   "--seed2", str(secrets.randbits(32)), "-p", prefix,
                   "-o", temporary1, temporary2]
        process = subprocess.Popen(
            command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        attempt = dict(task=task, process=process, temporary=(temporary1, temporary2))
        task["attempts"].append(attempt)
        if task["started"] is None:
            task["started"] = time.monotonic()
        active.append(attempt)

    def stop(attempt):
        process = attempt["process"]
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=0.2)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        for path in attempt["temporary"]:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass

    try:
        while completed < len(tasks):
            while len(active) < workers and next_task < len(tasks):
                start(tasks[next_task])
                next_task += 1

            if next_task == len(tasks):
                while len(active) < workers:
                    eligible = [task for task in tasks if not task["done"]]
                    if not eligible:
                        break

                    def hedge_key(task):
                        live = sum(1 for attempt in active if attempt["task"] is task)
                        return (live, task["started"])

                    start(min(eligible, key=hedge_key))

            for attempt in list(active):
                if attempt not in active:
                    continue
                status = attempt["process"].poll()
                if status is None:
                    continue

                task = attempt["task"]
                if status != 0:
                    raise RuntimeError(
                        f"collision search failed on {task['paths'][0]}: {status}")

                output1, output2 = task["paths"][1:]
                temporary1, temporary2 = attempt["temporary"]
                os.replace(temporary1, output1)
                os.replace(temporary2, output2)
                task["done"] = True
                completed += 1

                for peer in task["attempts"]:
                    if peer["process"].poll() is None:
                        peer["process"].terminate()
                for peer in task["attempts"]:
                    stop(peer)
                    if peer in active:
                        active.remove(peer)

            if completed < len(tasks):
                time.sleep(0.002)
    except (OSError, RuntimeError) as error:
        print(error, file=sys.stderr)
        return 1
    finally:
        for attempt in active:
            stop(attempt)

    return 0


if __name__ == "__main__":
    sys.exit(main())
