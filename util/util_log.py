import sys
import logging
from imp import reload

def main(debug:bool = False):

    reload(logging)
    if debug:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO
    logging.basicConfig(
        level=logging_level,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("run.log"),
        ],
    )
    logging.debug(f"Logging level: {logging_level}")
