# -*- coding: utf-8 -*-
from __future__ import unicode_literals, print_function
import sys
import traceback
from alignment import SegmentAlignment


def align_entity(original_text, altered_text, entities):
    """ align annotation with original text
    """
    gold_text = list(original_text)
    altered_text = list(altered_text)
    aligner = SegmentAlignment()

    # base_alginment = Hirschberg, segment_half = True, segment = 50, diff = 50
    aligned_gold, aligned_altered = aligner.align(
        gold_text, altered_text, segment_half=True, base_alignment='Hirschberg')
        
    alter2gold = aligner.map_alignment(aligned_gold, aligned_altered)

    for entity in entities:
        start = int(entity.get('charStart'))
        end = int(entity.get('charEnd'))

        try:
            entity['charStart'] = alter2gold[start]
            if end >= len(alter2gold):
                # end is an index in a range, so it could
                # equal to the length of length of the altered string
                entity['charEnd'] = alter2gold[-1]
            elif end > 0 and alter2gold[end] - alter2gold[end - 1] > 1:
                entity['charEnd'] = alter2gold[end - 1] + 1
            else:
                entity['charEnd'] = alter2gold[end]
        except IndexError:
            traceback.print_exc()
            print(len(alter2gold), start, end, sep="\t", file=sys.stderr)
        entity['entityText'] = original_text[entity.get('charStart'):entity.get('charEnd')+1]
