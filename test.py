from unittest import TestCase
from .alignment import *
import unittest
import os

class TestAlignment(TestCase):
    def setUp(self):
        self.path = os.path.dirname(__file__)
        pass

    def test_segnment_align(self):
        # found the bug of range index, [:end+1] is OK, [end+1]
        # throws IndexError
        of = os.path.join(self.path, 'data/6664499.original')
        af = os.path.join(self.path, 'data/6664499.altered')
        with open(of) as oh, open(af) as ah:
            aligner = SegmentAlignment()
            original_text = oh.read()
            altered_text = ah.read()
            aligned_gold, aligned_altered = aligner.align(original_text, altered_text,
                                                          segment_half=True, base_alignment='Hirschberg')
            print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')


    def test_segment_align_2(self):
        # found the bug of extra space in original text
        # looks like a bug in brat
        # try to fix it by using correct original text
        of = os.path.join(self.path, 'data/16741954.original')
        af = os.path.join(self.path, 'data/16741954.altered')
        with open(of) as oh, open(af) as ah:
            aligner = SegmentAlignment()
            original_text = oh.read()
            altered_text = ah.read()
            aligned_gold, aligned_altered = aligner.align(original_text, altered_text,
                                                          segment_half=True, base_alignment='Hirschberg')
            alter2gold = aligner.map_alignment(aligned_gold, aligned_altered)
            print(alter2gold[670], alter2gold[674])
            print(original_text[alter2gold[670]:alter2gold[674]])
            # print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')

if __name__ == '__main__':
    unittest.main()