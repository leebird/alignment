# -*- coding: utf-8 -*-
from __future__ import unicode_literals, print_function
import traceback
from alignment import SegmentAlignment


def align_entity(doc, original_text):
    """ align annotation with original text
    :param doc: a document matching document.proto.
    :type doc: dict
    :param original_text: the original text
    :type original_text: str
    """
    aligner = SegmentAlignment()
    altered_text = list(doc.get('text'))
    original_text = list(original_text)
    
    # base_alginment = Hirschberg, segment_half = True, segment = 50, diff = 50
    aligned_gold, aligned_altered = aligner.align(
        original_text, altered_text, segment_half=True, base_alignment='Hirschberg')
        
    original_text = ''.join(original_text)
    alter2gold = aligner.map_alignment(aligned_gold, aligned_altered)

    for entity_id, entity in doc.get('entity').items():
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
            print(doc.get('docId'), len(alter2gold), start, end, sep="\t")
        entity['entityText'] = original_text[entity.get('charStart'):entity.get('charEnd')+1]

    doc['text'] = original_text
