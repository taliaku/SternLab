import logging

def pipeline_logger(log_file='.log'):
    # create logger
    logger_name = 'PipelineLogger'
    logger = logging.getLogger(logger_name)
    if logging.getLogger(logger_name).hasHandlers(): # If logger exists, just return the existing logger.
        logger.warn(f"Logger already exists! Sticking with log file: {logger.handlers[1].baseFilename}")
        return logger
    logger.setLevel(logging.DEBUG)

    # create console handler and set level to info
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # create file handler and set level to debug
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
    
    return logger
