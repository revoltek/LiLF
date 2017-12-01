import logging
import lib_util

def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[1].levelno
        if (levelno >= 50):
            color = "\x1b[31m" # red
        elif (levelno >= 40):
            color = "\x1b[31m" # red
        elif (levelno >= 30):
            color = "\x1b[33m" # yellow
        elif (levelno >= 20):
            color = "\x1b[32m" # green 
        elif (levelno >= 10):
            color = "\x1b[35m" # pink
        else:
            color = "\x1b[0m" # normal
        args[1].msg = color + args[1].msg + "\x1b[0m"  # normal
        return fn(*args)
    return new


def set_logger(filename = "pipeline.logging"):
    # hopefully kill other loggers
    logger = logging.getLogger()
    logger.propagate = False
    for l in logger.handlers: l.setLevel("CRITICAL")

    logger = logging.getLogger("PiLL")
    logger.setLevel(logging.DEBUG)
    logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)
    
    # create file handler which logs even debug messages
    lib_util.check_rm(filename)
    handlerFile = logging.FileHandler(filename)
    handlerFile.setLevel(logging.DEBUG)
    
    # create console handler with a higher log level
    handlerConsole = logging.StreamHandler()
    handlerConsole.setLevel(logging.INFO)
    
    # create formatter and add it to the handlers
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    handlerFile.setFormatter(formatter)
    handlerConsole.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.addHandler(handlerFile)
    logger.addHandler(handlerConsole)

# this is used by all libraries for logging
logger = logging.getLogger("PiLL")
