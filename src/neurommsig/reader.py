import logging

import pandas as pd

from src.neurommsig.constants import (
    GEO_FILE,
    GEO_FILE_COLS,
    MAPPING_FILE,
    MAPPING_FILE_COLS,
    PATHWAY_FILE,
    PATHWAY_FILE_COLS,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class DataReader:
    def __init__(self):
        self.__gene_file = GEO_FILE
        self.__mapping_file = MAPPING_FILE
        self.__pathway_file = PATHWAY_FILE

        try:
            self._gene_data = self.__read_geo_data()
            self._pathway_data = self.__read_pathway_data()
            self._mapping_data = self.__read_mapping_data()
            logger.info("Reading data completed.")
        except Exception as e:
            logger.error(e)
            print(e)

    def __read_geo_data(self) -> pd.DataFrame:
        cols = list(GEO_FILE_COLS.keys())
        df = pd.read_csv(self.__gene_file, sep="\t")[cols]
        df.columns = [GEO_FILE_COLS[col] for col in df.columns]
        return df

    def __read_pathway_data(self) -> pd.DataFrame:
        cols = list(PATHWAY_FILE_COLS.keys())
        df = pd.read_csv(self.__pathway_file, sep="\t", header=None)[cols]
        df.columns = [PATHWAY_FILE_COLS[col] for col in df.columns]
        return df

    def __read_mapping_data(self) -> pd.DataFrame:
        cols = list(MAPPING_FILE_COLS.keys())
        df = pd.read_csv(self.__mapping_file, sep="\t", index_col=0)[cols]
        df.columns = [MAPPING_FILE_COLS[col] for col in df.columns]
        return df

    def get_gene_set(self) -> pd.DataFrame:
        """Return gene set from a GEO experiment."""
        return self._gene_data

    def get_pathway(self) -> pd.DataFrame:
        """Return a pathway."""
        return (
            self._pathway_data.merge(self._mapping_data, how="left", on=["source", "target"])[
                ["source", "target", "relation"]
            ]
            .drop_duplicates()
            .reset_index(drop=True)
        )
