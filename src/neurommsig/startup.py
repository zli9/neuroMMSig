# -*- coding: utf-8 -*-

import logging
import os
from pathlib import Path

# Default Directory Paths
home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".neurommsig")
os.makedirs(PROJECT_DIR, exist_ok=True)

# Logging Configuration
LOG_FILE_PATH = os.path.join(PROJECT_DIR, "neurommsig.log")
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", filename=LOG_FILE_PATH
)
