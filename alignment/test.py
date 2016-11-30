from __future__ import unicode_literals, print_function
from unittest import TestCase
import unittest
import os
import codecs
from alignment import *
from align import *


class TestAlignment(TestCase):
    def setUp(self):
        self.path = os.path.dirname(__file__)

    def test_segnment_align(self):
        # found the bug of range index, [:end+1] is OK, [end+1]
        # throws IndexError
        of = os.path.join(self.path, 'data/6664499.original')
        af = os.path.join(self.path, 'data/6664499.altered')
        with open(of) as oh, open(af) as ah:
            aligner = SegmentAlignment()
            original_text = oh.read()
            altered_text = ah.read()
            aligned_gold, aligned_altered = aligner.align(list(original_text),
                                                          list(altered_text),
                                                          segment_half=True,
                                                          base_alignment='Hirschberg')
            score = aligner.score(aligned_gold, aligned_altered)
            print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')
            print(score, score/float(len(original_text)))

    def test_segment_align_2(self):
        # found the bug of extra space in original text
        # looks like a bug in brat
        # try to fix it by using correct original text
        of = os.path.join(self.path, 'data/24742516.original')
        af = os.path.join(self.path, 'data/24742516.altered')
        with open(of) as oh, open(af) as ah:
            aligner = SegmentAlignment()
            original_text = oh.read()
            altered_text = ah.read()
            aligned_gold, aligned_altered = aligner.align(list(original_text),
                                                          list(altered_text),
                                                          segment_half=True,
                                                          base_alignment='Hirschberg')
            alter2gold = aligner.map_alignment(aligned_gold, aligned_altered)
            # print(alter2gold[1515], alter2gold[1525])
            # print(original_text[alter2gold[1515]:alter2gold[1525]])
            # print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')

    def test_segment_align_3(self):
        # Problematic case, not fiexed yet.
        # found the bug of extra space in original text
        # looks like a bug in brat
        # try to fix it by using correct original text
        of = os.path.join(self.path, 'data/26662996.original')
        af = os.path.join(self.path, 'data/26662996.altered')
        with open(of) as oh, open(af) as ah:
            aligner = SegmentAlignment()
            original_text = oh.read()
            altered_text = ah.read()
            aligned_gold, aligned_altered = aligner.align(list(original_text),
                                                          list(altered_text),
                                                          segment_half=True,
                                                          base_alignment='Hirschberg')
            alter2gold = aligner.map_alignment(aligned_gold, aligned_altered)
            score = aligner.score(aligned_gold, aligned_altered)
            # print(alter2gold[1515], alter2gold[1525])
            # print(original_text[alter2gold[1515]:alter2gold[1525]])
            print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')
            print(score, score/float(len(original_text)))

    def test_segment_align_4(self):
        # Problematic case, not fiexed yet.
        # found the bug of extra space in original text
        # looks like a bug in brat
        # try to fix it by using correct original text
        of = os.path.join(self.path, 'data/27600506.original')
        af = os.path.join(self.path, 'data/27600506.altered')
        with codecs.open(of, encoding='utf8') as oh, codecs.open(af, encoding='utf8') as ah:
            aligner = SegmentAlignment()
            original_text = oh.read().lower()
            altered_text = ah.read().lower()
            aligned_gold, aligned_altered = aligner.align(list(original_text),
                                                          list(altered_text),
                                                          segment_half=True,
                                                          base_alignment='Hirschberg')
            alter2gold = aligner.map_alignment(aligned_gold, aligned_altered)
            score = aligner.score(aligned_gold, aligned_altered)
            # print(alter2gold[1515], alter2gold[1525])
            # print(original_text[alter2gold[1515]:alter2gold[1525]])
            print(''.join(aligned_gold), ''.join(aligned_altered), sep='\n')
            print(score, score/float(len(original_text)))


class TestSet(TestCase):
    def test(self):
        import os

        h = Hirschberg()
        s = SegmentAlignment()

        for root, _, files in os.walk('data/raw'):
            for f in files:
                # if f != 'PMID-2355960.txt':
                # continue
                print(f)
                raw_text_file = open(os.path.join(root, f), 'r')
                raw_text = raw_text_file.read()
                raw_text_file.close()

                altered_text_file = open(os.path.join('data/altered', f), 'r')
                altered_text = altered_text_file.read()
                altered_text_file.close()

                # golden global alignment with Hirschberg
                ha, hb = h.align(list(raw_text), list(altered_text))

                nsa, nsb = s.align(list(raw_text), list(altered_text),
                                   segment_half=True,
                                   base_alignment='Needleman')
                hsa, hsb = s.align(list(raw_text), list(altered_text),
                                   segment_half=True,
                                   base_alignment='Hirschberg')

                print('%22s' % 'golden-semiglobal', '%6s' % (ha == nsa),
                      '%6s' % (hb == nsb),
                      '%.4f' % (len(nsa) / len(ha)),
                      '%.4f' % (len(nsb) / len(hb)),
                      '\n%22s' % 'golden-segment-half', '%6s' % (ha == hsa),
                      '%6s' % (hb == hsb),
                      '%.4f' % (len(hsa) / len(ha)),
                      '%.4f' % (len(hsb) / len(hb)))
                print()

                res = open('data/aligned/' + f, 'w')
                res.write(''.join(ha) + '\n' + ''.join(hb) + '\n\n')
                res.write(''.join(nsa) + '\n' + ''.join(nsb) + '\n\n\n')

                # reverse
                # ha, hb = h.align(list(altered_text), list(raw_text))
                #
                # nsa, nsb = s.align(list(altered_text), list(raw_text), segment_half=True, base_alignment='Needleman')
                # hsa, hsb = s.align(list(altered_text), list(raw_text), segment_half=True, base_alignment='Hirschberg')
                #
                # print('%22s' % 'golden-semiglobal', '%6s' % (ha == nsa), '%6s' % (hb == nsb),
                #       '%.4f' % (len(nsa) / len(ha)), '%.4f' % (len(nsb) / len(hb)),
                #       '\n%22s' % 'golden-segment-half', '%6s' % (ha == hsa), '%6s' % (hb == hsb),
                #       '%.4f' % (len(hsa) / len(ha)), '%.4f' % (len(hsb) / len(hb)))
                # print()
                #
                # res.write(''.join(hb) + '\n' + ''.join(ha) + '\n\n')
                # res.write(''.join(nsb) + '\n' + ''.join(nsa))

                res.close()


class TestFunction(TestCase):
    def test_functions(self):
        seqa = list('12345678')
        seqb = list('123478908')

        n = Needleman()
        a, b = n.align(seqa, seqb)
        print(a)
        print(b)

        h = Hirschberg()
        a, b = h.align(seqa, seqb)
        print(a)
        print(b)

        seqa = list('TI - Transcription factor AP-2 activity is modulated')
        seqb = list('Transcription factor AP-2 activity is modulated')
        s = SegmentAlignment()
        a, b = s.align(seqa, seqb, semi_global=True)
        print(a)
        print(b)

        # semi-global
        seqa = list('CGTACGTGAGTGA')
        seqb = list('CGATTA')
        seqa = list('TI - Transcription factor AP-2 activity is modulated')
        seqb = list('Transcription factor AP-2 activity is modulated')
        # ['C', 'G', '|', 'T', '|', 'A', 'C', 'G', 'T', 'G', 'A', 'G', 'T', 'G', 'A']
        # ['C', 'G', 'A', 'T', 'T', 'A', '|', '|', '|', '|', '|', '|', '|', '|', '|']
        n = Needleman()
        a, b = n.align(seqa, seqb)
        print(a)
        print(b)


class TestAlignEntity(TestCase):
    def test_align_entity(self):
        original_text = 'I have a book.'
        altered_text = ' I  have a book.'
        entities = [{'charStart': 4, 'charEnd': 7}]
        align_entity(original_text, altered_text, entities)
        print(entities)

if __name__ == '__main__':
    unittest.main()
