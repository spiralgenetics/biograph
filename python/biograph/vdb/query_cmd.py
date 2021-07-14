#!/usr/bin/env python3
'''
VDB: interactive query commands
'''
import argparse
import atexit
import readline
import sys

from pathlib import Path

import biograph.vdb.athena as athena
from biograph.tools.log import setup_logging

def cmd_interactive():
    ''' Interactive mode '''
    db = athena.connect()
    while True:
        try:
            for reply in db.query(input(f"vdb:{db.database}> ")):
                if isinstance(reply, tuple) and len(reply) == 1:
                    print(reply[0])
                else:
                    print(reply)
        except (KeyboardInterrupt, EOFError):
            raise SystemExit('')
        except (Exception, SystemExit) as e: # pylint: disable=broad-except
            print(e)

def main(clargs):
    ''' No subcommands for query '''
    parser = argparse.ArgumentParser(prog="query",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("sql", nargs="?", type=str, help="Run a single SQL query and exit. Leave blank to run interactive mode.")
    parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args(clargs)

    setup_logging(debug_mode=args.debug, simple=True)

    if args.sql:
        db = athena.connect()
        print(db.query(args.sql))
        exit(0)

    histfile = Path("~/.config/spiral/vdb.history").expanduser()
    try:
        readline.read_history_file(histfile)
        readline.set_history_length(1000)
    except FileNotFoundError:
        histfile.parent.mkdir(parents=True, exist_ok=True)

    atexit.register(readline.write_history_file, histfile)

    cmd_interactive()

# top level command
CMD = 'biograph vdb query'

if __name__ == '__main__':
    main(sys.argv[1:])
