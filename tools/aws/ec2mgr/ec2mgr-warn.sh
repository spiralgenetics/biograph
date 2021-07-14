#!/bin/bash
if [ ! -f /etc/cron.d/ec2mgr ]; then
    return
fi

if [ -f /tmp/do_not_stop ]; then
    cat <<-EOM

*** NOTE ***

Automatic shutdown has been overridden for this instance.

To automatically shutdown after one hour with no active SSH connections:

  $ rm /tmp/do_not_stop

EOM
    uptime

else
    cat <<-EOM

*** NOTE ***

This instance will automatically shut down after one hour if there are
no active SSH connections.

To prevent automatic shutdown:

  $ touch /tmp/do_not_stop

EOM
fi
