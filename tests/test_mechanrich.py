import os

import pytest

from src.mechanrich.mechanrich import RCR, RCRstat, Graph


class TestMechanrich:
    """Tests for mechanrich.py."""

    def testGraph(self):
        """Test for Graph class."""
        # test plot_full_network function
        graph = Graph()
        output_path = "test_data/fake_network.pdf"
        graph.plot_full_network(output=output_path)

        assert os.path.exists(output_path) is True
        os.remove(output_path)

        wrong_output_path = "test_data/fake_network.ext"
        with pytest.raises(ValueError):
            graph.plot_full_network(output=wrong_output_path)

    def testRCRstat(self):
        """Test for RCRstat class."""
        # test get_stat function
        stat = RCRstat()
        output_path = "test_data/fake_network_stat.txt"
        stat.get_stat(output_path)

        assert os.path.exists(output_path) is True
        os.remove(output_path)
