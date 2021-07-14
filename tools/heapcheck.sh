#!/bin/sh
set -e
LD_PRELOAD=tools/tcmalloc.so
export LD_PRELOAD
exec "$@"
