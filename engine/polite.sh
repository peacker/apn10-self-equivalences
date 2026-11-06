#!/bin/bash
# polite.sh CMD...: wait while the machine is busy (1-min load above LIMIT, default 200), then run CMD at nice 19
LIMIT=${LIMIT:-200}
while [ "$(cut -d' ' -f1 /proc/loadavg | cut -d. -f1)" -gt "$LIMIT" ]; do sleep 60; done
exec nice -n 19 "$@"
