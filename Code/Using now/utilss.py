from collections import OrderedDict


def sortORFs(orf_dict):
    """
    Sort ORFs by their start position.

    Parameters:
        orf_dict (OrderedDict): Dictionary with position keys (e.g. "100,300")

    Returns:
        OrderedDict: Sorted by start position (ascending)
    """
    return OrderedDict(
        sorted(orf_dict.items(), key=lambda x: int(x[0].split(',')[0]))
    )


def sortORFs_by_strand(orf_dict):
    """
    Optional: Sort ORFs by strand order based on ID.

    Parameters:
        orf_dict (OrderedDict): ORF dictionary

    Returns:
        OrderedDict: ORFs sorted by strand order (based on internal ID)
    """
    final = OrderedDict()
    ids = sorted([orf[-1] for orf in orf_dict.values()])
    for i in ids:
        for key, value in orf_dict.items():
            if value[-1] == i:
                final[key] = value
                break
    return final
