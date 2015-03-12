from unittest import TestCase
from alignment import *
import unittest

class TestAlignment(TestCase):
    def setUp(self):
        pass

    def test_segnment_align(self):
        with open('data/6664499.original') as oh, open('data/6664499.altered') as ah:
            aligner = SegmentAlignment()
            original_text = oh.read()
            altered_text = ah.read()
            aligned_gold, aligned_altered = aligner.align(original_text, altered_text,
                                                          segment_half=True, base_alignment='Hirschberg')
            print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')


if __name__ == '__main__':
    unittest.main()