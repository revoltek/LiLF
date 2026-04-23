# Logging utilities for the LiLF pipeline.
# Provides a coloured console handler and a Logger class that sets up
# per-run log files and convenience symlinks.

import os
import sys
import logging
import time

import warnings
from astropy.wcs import FITSFixedWarning

# Suppress noisy FITS WCS warnings that are irrelevant to pipeline logic.
warnings.simplefilter('ignore', FITSFixedWarning)

class _ColorStreamHandler(logging.StreamHandler):
    """A StreamHandler that wraps log messages in ANSI colour escape codes
    based on severity level."""

    # ANSI reset and colour codes.
    DEFAULT = '\x1b[0m'
    RED     = '\x1b[31m'
    GREEN   = '\x1b[32m'
    YELLOW  = '\x1b[33m'
    CYAN    = '\x1b[36m'

    # Map log levels to colours.
    CRITICAL = RED
    ERROR    = RED
    WARNING  = YELLOW
    INFO     = GREEN
    DEBUG    = CYAN

    @classmethod
    def _get_color(cls, level):
        """Return the ANSI colour code appropriate for the given log level."""
        if level >= logging.CRITICAL:  return cls.CRITICAL
        elif level >= logging.ERROR:   return cls.ERROR
        elif level >= logging.WARNING: return cls.WARNING
        elif level >= logging.INFO:    return cls.INFO
        elif level >= logging.DEBUG:   return cls.DEBUG
        else:                          return cls.DEFAULT

    def format(self, record):
        """Wrap the message in colour codes without permanently mutating the
        shared LogRecord (which would cause codes to stack when multiple
        handlers or repeated calls process the same record)."""
        color = self._get_color(record.levelno)
        # Save original message; record.msg may be any type, cast to str.
        original_msg = record.msg
        record.msg = color + str(record.msg) + self.DEFAULT
        formatted = logging.StreamHandler.format(self, record)
        # Restore the original so other handlers see the unmodified record.
        record.msg = original_msg
        return formatted

class Logger():
    """Sets up the LiLF logging infrastructure for a pipeline run.

    Creates a timestamped log directory and log file, configures the
    'LiLF' logger with both a file handler (DEBUG) and a coloured console
    handler (INFO), and maintains 'latest' symlinks for convenience.
    """

    def __init__(self, pipename):
        # Silence the root logger to prevent duplicate output from any
        # third-party library that calls logging.basicConfig() or adds
        # handlers to the root logger.
        root = logging.getLogger()
        root.propagate = False
        # Close existing handlers before discarding them to avoid file
        # descriptor leaks.
        for h in root.handlers[:]:
            h.close()
        root.handlers = []

        # Build timestamped names for this run's log file and log directory.
        timestamp = time.strftime('%Y-%m-%d_%H:%M:%S', time.localtime())
        self.logfile = pipename + '_' + timestamp + '.logger'
        self.log_dir = 'logs_' + pipename + '_' + timestamp
        os.makedirs(self.log_dir)

        self.set_logger(self.logfile)

        # Create 'latest' symlinks so users can always find the most recent
        # logs without knowing the timestamp.  A temp-then-rename pattern is
        # used so the replacement is atomic and never leaves a broken link.
        os.symlink(self.log_dir, f'logs_{pipename}.tmp')
        os.rename(f'logs_{pipename}.tmp', f'logs_{pipename}')
        os.symlink(self.logfile, pipename + '.logger.tmp')
        os.rename(pipename + '.logger.tmp', pipename + '.logger')

    def set_logger(self, logfile):
        """Configure the 'LiLF' named logger with a file handler and a
        coloured console handler."""

        logger = logging.getLogger("LiLF")
        # Accept all levels here; individual handlers filter further.
        logger.setLevel(logging.DEBUG)
        # Prevent messages from propagating to the root logger.
        # Without this, every record is handled by our handlers AND then
        # forwarded to the root logger which — when it has no handlers —
        # falls back to Python's lastResort handler, producing the bare
        # "INFO: message" duplicate lines.
        logger.propagate = False

        # Remove any handlers attached in a previous Logger() call
        # (e.g. PiLL's top-level logger) so we don't accumulate duplicates.
        for h in logger.handlers[:]:
            h.close()
        logger.handlers = []

        # File handler: capture everything down to DEBUG for post-mortem analysis.
        handlerFile = logging.FileHandler(logfile)
        handlerFile.setLevel(logging.DEBUG)

        # Console handler: show INFO and above with colour coding.
        handlerConsole = _ColorStreamHandler(stream=sys.stdout)
        handlerConsole.setLevel(logging.INFO)

        # Shared formatter: timestamp, level, and message.
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
        handlerFile.setFormatter(formatter)
        handlerConsole.setFormatter(formatter)

        logger.addHandler(handlerFile)
        logger.addHandler(handlerConsole)

        logger.info('Logging initialised in %s (file: %s)' % (os.getcwd(), logfile))


# Module-level logger used by all LiLF library modules.
# Note: Logger() must be instantiated by the pipeline entry point before
# this logger emits any output, otherwise messages go to the root logger.
logger = logging.getLogger("LiLF")