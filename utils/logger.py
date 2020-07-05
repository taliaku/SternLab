import os
import logging

def _logger_already_exists(logger, log_file):
    if logger.hasHandlers(): # If logger exists, just return the existing logger.
        if logger.handlers[1].baseFilename != log_file:
            if log_file is not None:
                logger.warn(f"Logger already exists! Sticking with log file: {logger.handlers[1].baseFilename}")
        return_value = True
    else:
        return_value = False
    return return_value

def _create_new_logger(logger, log_file):

    logger.setLevel(logging.DEBUG)
    # create console handler and set level to info
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create file handler and set level to debug
    if not os.path.exists(log_file):
        os.mknod(log_file)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # add formatter
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add ch & fh to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.info(f'Log started! Outputing to: {log_file}')
    return logger

def pipeline_logger(logger_name, log_folder=None):
    # create logger
    if log_folder is not None:
        log_file = os.path.join(log_folder, '.log')
    else:
        log_file = None
    logger = logging.getLogger(logger_name)
    if not _logger_already_exists(logger, log_file):
        if log_file is None:
            raise ValueError("First instance of logger must be initiated with an output file!")
        logger = _create_new_logger(logger, log_file)
    return logger
