# -*- coding: utf-8 -*-

import logging
import os
from pathlib import Path

# Default Directory Paths
home_dir = Path.home()
PROJECT_DIR = home_dir.joinpath(".mechanrich")
os.makedirs(PROJECT_DIR, exist_ok=True)

# Logging Configuration
LOG_FILE_PATH = os.path.join(PROJECT_DIR, "mechanrich.log")
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", filename=LOG_FILE_PATH
)
