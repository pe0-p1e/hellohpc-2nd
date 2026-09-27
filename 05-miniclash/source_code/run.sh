#!/usr/bin/env bash
set -e
if [ "$#" -ne 1 ]; then
    echo "usage: $0 <tasks.txt>" >&2
    exit 2
fi
exec "$(dirname "$0")/md5fastcoll" --batch "$1"
