import logging

from src.neurommsig.constants import FC_THRED, P_THRED
from src.neurommsig.reader import DataReader

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class PreProcessing:
    """Pre-processing dataframe of user input."""

    def __init__(self):
        reader = DataReader()
        self._gene_set = reader.get_gene_set()
        self._upregu_genes = None
        self._downregu_genes = None

        self.__ext_regu_genes()

    def __ext_regu_genes(self):
        """Extract up- and down-regulated genes"""
        self._upregu_genes = self._gene_set.gene_symbol[
            (self._gene_set["p_value"] < P_THRED) & (self._gene_set["log_fold_change"] > FC_THRED)
        ].tolist()
        self._downregu_genes = self._gene_set.gene_symbol[
            (self._gene_set["p_value"] < P_THRED) & (self._gene_set["log_fold_change"] < -FC_THRED)
        ].tolist()

    def get_upregu_genes(self):
        """Return up-regulated genes."""
        return self._upregu_genes

    def get_downregu_genes(self):
        """Return down-regulated genes."""
        return self._downregu_genes
