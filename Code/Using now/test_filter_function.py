import unittest
from collections import OrderedDict
from utilss import sortORFs, sortORFs_by_strand
from Filter import tile_filtering, storf_type_score



# --------------------------------------------------------
# Class: MockOptions
# Purpose: Simulate argparse.Namespace to pass options to tile_filtering
# --------------------------------------------------------
class MockOptions:
    """
    Mock argparse.Namespace to simulate CLI options passed into tile_filtering().
    """
    def __init__(self, priority_strategy='length', overlap_nt=50, storf_order='start_pos'):
        self.priority_strategy = priority_strategy
        self.overlap_nt = overlap_nt
        self.storf_order = storf_order


# --------------------------------------------------------
# Class: TestTileFiltering
# Purpose: Contains all unit tests for tile_filtering() behaviour
# --------------------------------------------------------
class TestTileFiltering(unittest.TestCase):

    def setUp(self):
        """
        This method runs before each test. It provides example StORFs to filter.
        Each key is a position string ("start,stop"), and the value is:
        [sequence, frame, strand, length, type, id]
        """
        self.test_data = OrderedDict({
            "100,300": ["ATG...", 1, '+', 201, 'StORF', 0],
            "150,320": ["ATG...", 1, '+', 171, 'StORF', 1],
            "310,500": ["ATG...", 1, '+', 191, 'Con-StORF', 2],
            "490,700": ["ATG...", 1, '+', 211, 'StORF', 3],
            "650,850": ["ATG...", 1, '+', 201, 'Con-StORF', 4],
            "900,1100": ["ATG...", 1, '+', 201, 'StORF', 5]
        })

    def test_length_priority_default(self):
        """
        Test that longer ORFs are prioritised when overlapping,
        and shorter ones are removed.
        """
        options = MockOptions(priority_strategy='length', overlap_nt=50)
        result = tile_filtering(self.test_data, options)

        # These should remain
        expected_keys = ["100,300", "310,500", "490,700", "900,1100"]
        for key in expected_keys:
            self.assertIn(key, result)

        # This one is shorter and overlaps — should be gone
        self.assertNotIn("150,320", result)
    
    
    def test_storf_type_priority(self):
        """
         Test that Con-StORFs are prioritised over regular StORFs.
     """
        options = MockOptions(priority_strategy='storf_type', overlap_nt=50)
        result = tile_filtering(self.test_data, options)

         # Make sure result is not empty
        self.assertTrue(len(result) > 0)

    # Con-StORFs should be in the filtered result
        kept_types = [v[4] for v in result.values()]
        self.assertIn('Con-StORF', kept_types)
        con_count = kept_types.count('Con-StORF')
        storf_count = kept_types.count('StORF')
        self.assertGreaterEqual(con_count, 1)
        
    
    def test_overlap_threshold(self):
        """
        Test that a stricter (lower) overlap threshold removes more ORFs.
        """
        options = MockOptions(priority_strategy='length', overlap_nt=10)
        result = tile_filtering(self.test_data, options)

        # Expect fewer ORFs to survive with tight filtering
        self.assertLess(len(result), len(self.test_data))

        
    def test_storf_order_sorting(self):
        """
        Test that sorting by strand maintains order based on internal IDs.
        """
        options = MockOptions(priority_strategy='length', storf_order='strand') 
        result = tile_filtering(self.test_data, options)

        # Extract final IDs
        storf_indices = [v[5] for v in result.values()]
        self.assertEqual(storf_indices, sorted(storf_indices))


# ------------------------------------------
# Run all tests
# ------------------------------------------
if __name__ == '__main__':
    unittest.main()