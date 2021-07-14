'''
Unified logging methods, featuring several convenience features:

 * setup_logging() is called automatically the first time a log message is sent

 * setup_logging() can be called again at any time to change the log format

 * The default log format is pcmp-style TIME: [LEVEL] MSG, and can be
   overridden by passing log_format, or set simple=True for plain log messages
   with no timestamp

 * Call setup_logging(debug_mode=True) to enable enhanced debugging, with
   stack depth, timestamp, and calling function name

 * `import biograph.tools.log as log` for log.setup_logging, log.info, log.warning, etc.
   or `from biograph.tools.log import *` to import everything

 * All logging methods take any number of arguments and convert them to str on the fly

 * The log() function logging level is logging.INFO. Call warn(), error(),
   crit(), etc. to log at a specific level.

'''

import inspect
import logging
import sys
import warnings

DEFAULT_LOG_LEVEL = logging.INFO
LOGGING_IS_SET_UP = False

__all__ = ["setup_logging", "debug", "info", "warning", "error", "critical", "warn", "crit", "log"]

def setup_logging(debug_mode=False, stream=sys.stderr, log_format="%(asctime)s [%(levelname)s] %(message)s", simple=False):
    '''
    If setup_logging() hasn't been called when the first message is logged, it
    is called automatically with defaults.

    debug_mode: log time, stack depth, and calling function name
    stream: any filehandle, STDERR by default
    log_format: pcmp-style logging, override with your preferred format
    simple: set to True for simple one-line, no timestamp logging
    '''
    level = logging.DEBUG if debug_mode else DEFAULT_LOG_LEVEL

    if simple:
        if debug_mode:
            log_format = "%(asctime)s : %(message)s"
        else:
            log_format = "%(message)s"

    # in python 3.8 or later we can just pass force=True to basicConfig. In the meantime:
    global LOGGING_IS_SET_UP # pylint: disable=global-statement
    if LOGGING_IS_SET_UP:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

    logging.basicConfig(stream=stream, level=level, format=log_format)

    # pylint:disable=unused-argument
    def sendWarningsToLog(message, category, filename, lineno, *args, **kwargs):
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno, category.__name__, message)

    warnings.showwarning = sendWarningsToLog

    LOGGING_IS_SET_UP = True

def debug(*args):
    ''' Auto debug logging with caller name and visual stack depth '''
    # Note that default setup doesn't show the debug log.
    # If you want debugging, turn it on explicitly.
    if not LOGGING_IS_SET_UP:
        setup_logging()

    logging.debug(f"{'>' * len(inspect.stack())} {inspect.stack()[1].function}: {' '.join([str(a) for a in args])}")

def info(*args):
    ''' info logging '''
    if not LOGGING_IS_SET_UP:
        setup_logging()

    logging.info(' '.join([str(a) for a in args]))

def warning(*args):
    ''' warning logging '''
    if not LOGGING_IS_SET_UP:
        setup_logging()

    logging.warning(' '.join([str(a) for a in args]))

def error(*args):
    ''' error logging '''
    if not LOGGING_IS_SET_UP:
        setup_logging()

    logging.error(' '.join([str(a) for a in args]))

def critical(*args):
    ''' critical logging '''
    if not LOGGING_IS_SET_UP:
        setup_logging()

    logging.critical(' '.join([str(a) for a in args]))

# short aliases
warn = warning
crit = critical
log = info
