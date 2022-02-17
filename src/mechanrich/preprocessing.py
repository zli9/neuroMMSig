import logging

from src.mechanrich.constants import FC_THRED, P_THRED
from src.mechanrich.reader import DataReader

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class PreProcessing:
    """Pre-processing dataframe of user input."""

    def __init__(self):
        reader = DataReader()
        self.gene_set = reader.get_gene_set()
        self.upregu_genes = None
        self.downregu_genes = None

        self.__ext_regu_genes()

    def __ext_regu_genes(self):
        """Extract up- and down-regulated genes"""
        self.upregu_genes = self.gene_set.gene_symbol[
            (self.gene_set["p_value"] < P_THRED) & (self.gene_set["log_fold_change"] > FC_THRED)
        ].tolist()
        self.downregu_genes = self.gene_set.gene_symbol[
            (self.gene_set["p_value"] < P_THRED) & (self.gene_set["log_fold_change"] < -FC_THRED)
        ].tolist()

    def get_upregu_genes(self):
        """Return up-regulated genes."""
        return self.upregu_genes

    def get_downregu_genes(self):
        """Return down-regulated genes."""
        return self.downregu_genes
