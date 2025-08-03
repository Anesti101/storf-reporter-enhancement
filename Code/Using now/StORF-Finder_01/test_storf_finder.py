import unittest
import StORF_Finder
from collections import OrderedDict
from types import SimpleNamespace
from your_module_name import tile_filtering  # Replace with actual module name

class TestTileFiltering(unittest.TestCase):

    def setUp(self):
        # Simulated input: OrderedDict of StORFs with positions as keys
        # Format: [sequence, frame, strand, length, type, storf_index]
        self.storfs = OrderedDict({
            "100,200": ["ATG...TAA", "1", "+", 100, "StORF", 1],
            "180,250": ["GTG...TAG", "1", "+", 70, "Con-StORF", 2],
            "210,300": ["TTG...TGA", "1", "+", 90, "StORF", 3],
            "240,320": ["ATG...TAA", "1", "+", 80, "Con-StORF", 4]
        })

        # Common options object with all flags needed
        self.default_options = SimpleNamespace(
            overlap_nt=20,
            storf_order='start_pos',
            priority_strategy='length'
        )

        self.priority_options = SimpleNamespace(
            overlap_nt=20,
            storf_order='start_pos',
            priority_strategy='storf_type'
        )

    def test_length_priority_removes_shorter_overlapping(self):
        result = tile_filtering(self.storfs, self.default_options)
        keys = list(result.keys())

        # Should retain 100-200 and remove 180-250 due to overlap and length
        self.assertIn("100,200", keys)
        self.assertNotIn("180,250", keys)

    def test_storf_type_priority_keeps_con_storf_first(self):
        result = tile_filtering(self.storfs, self.priority_options)
        keys = list(result.keys())

        # Should retain 180-250 over 100-200 due to Con-StORF priority
        self.assertIn("180,250", keys)
        self.assertNotIn("100,200", keys)

    def test_final_ordering_by_start_pos(self):
        result = tile_filtering(self.storfs, self.default_options)
        starts = [int(k.split(',')[0]) for k in result.keys()]
        self.assertEqual(starts, sorted(starts))

    def test_overlapping_threshold_exclusion(self):
        # Set a strict threshold to allow everything
        loose_options = SimpleNamespace(
            overlap_nt=100,
            storf_order='start_pos',
            priority_strategy='length'
        )
        result = tile_filtering(self.storfs, loose_options)
        self.assertEqual(len(result), len(self.storfs))  # All should pass

if __name__ == '__main__':
    unittest.main()
