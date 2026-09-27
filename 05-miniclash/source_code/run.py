#!/usr/bin/env python3
"""Low-overhead parallel/hedged scheduler for independent fastcoll jobs."""
import os
import secrets
import subprocess
import sys


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <tasks.txt>", file=sys.stderr)
        return 2

    root = os.path.dirname(os.path.abspath(__file__))
    binary = os.path.join(root, "md5fastcoll")
    if not os.access(binary, os.X_OK):
        print(f"{binary} not found or not executable", file=sys.stderr)
        return 1

    try:
        with open(sys.argv[1], encoding="utf-8") as f:
            rows = [line.split() for line in f if line.strip()]
    except OSError as e:
        print(e, file=sys.stderr)
        return 1
    if any(len(row) != 3 for row in rows):
        print("each task must have three paths", file=sys.stderr)
        return 1
    if not rows:
        return 0

    cpus = sorted(os.sched_getaffinity(0))
    workers = min(32, len(cpus))
    tasks = [dict(paths=row, attempts=[], done=False, order=i) for i, row in enumerate(rows)]
    active_by_pid = {}
    free_cpus = cpus[:workers]
    next_task = 0
    completed = 0
    serial = 0

    def start(task):
        nonlocal serial
        serial += 1
        cpu = free_cpus.pop()
        prefix, out1, out2 = task["paths"]
        tmp1 = f"{out1}.hedge-{os.getpid()}-{serial}"
        tmp2 = f"{out2}.hedge-{os.getpid()}-{serial}"
        cmd = [binary, "-q", "--seed1", str(secrets.randbits(32)),
               "--seed2", str(secrets.randbits(32)), "-p", prefix,
               "-o", tmp1, tmp2]
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            os.sched_setaffinity(p.pid, {cpu})
        except OSError:
            pass
        a = dict(task=task, process=p, temporary=(tmp1, tmp2), cpu=cpu)
        task["attempts"].append(a)
        active_by_pid[p.pid] = a

    def cleanup_file(path):
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass

    def retire(a, terminate=False):
        p = a["process"]
        if terminate and p.poll() is None:
            p.terminate()
        try:
            p.wait(timeout=0.05 if terminate else None)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
        active_by_pid.pop(p.pid, None)
        free_cpus.append(a["cpu"])
        for path in a["temporary"]:
            cleanup_file(path)
        try:
            a["task"]["attempts"].remove(a)
        except ValueError:
            pass

    def fill():
        nonlocal next_task
        # First guarantee one independent search for every not-yet-started task.
        while free_cpus and next_task < len(tasks):
            start(tasks[next_task])
            next_task += 1
        if next_task < len(tasks):
            return
        # Thereafter every spare core becomes a hedge. Balance copies across
        # unfinished tasks; unlike the old scheduler, allow the last straggler
        # to consume all available cores instead of stopping at eight copies.
        while free_cpus:
            eligible = [t for t in tasks if not t["done"]]
            if not eligible:
                break
            task = min(eligible, key=lambda t: (len(t["attempts"]), t["order"]))
            start(task)

    try:
        fill()
        while completed < len(tasks):
            # Block in the kernel until a worker exits: no 20 ms polling holes
            # and effectively zero supervisor CPU while all cores are busy.
            pid, status = os.waitpid(-1, 0)
            a = active_by_pid.get(pid)
            if a is None:
                continue
            p = a["process"]
            p.returncode = os.waitstatus_to_exitcode(status)
            task = a["task"]

            if task["done"]:
                retire(a, terminate=False)
                fill()
                continue
            if p.returncode != 0:
                raise RuntimeError(f"collision search failed on {task['paths'][0]}: {p.returncode}")

            out1, out2 = task["paths"][1:]
            tmp1, tmp2 = a["temporary"]
            os.replace(tmp1, out1)
            os.replace(tmp2, out2)
            task["done"] = True
            completed += 1

            peers = list(task["attempts"])
            # Release the winning slot without deleting the now-renamed files.
            active_by_pid.pop(p.pid, None)
            free_cpus.append(a["cpu"])
            task["attempts"].remove(a)
            for peer in peers:
                if peer is a:
                    continue
                if peer["process"].poll() is None:
                    peer["process"].terminate()
            for peer in peers:
                if peer is a:
                    continue
                retire(peer, terminate=True)
            fill()
    except (OSError, RuntimeError) as e:
        print(e, file=sys.stderr)
        return 1
    finally:
        for a in list(active_by_pid.values()):
            retire(a, terminate=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
