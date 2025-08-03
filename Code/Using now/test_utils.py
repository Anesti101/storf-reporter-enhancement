
from collections import OrderedDict
from utilss import sortORFs, sortORFs_by_strand

print("test_utils.py started!!!!!!!")


test_orfs = OrderedDict({
    "300,500": ["ATG...", 1, '+', 201],
    "100,200": ["ATG...", 1, '+', 101],
    "250,400": ["ATG...", 1, '+', 151]
})

sorted_orfs = sortORFs(test_orfs)

print("Sorted ORFs by start position:")
for key, value in sorted_orfs.items():
    print(f"{key} -> {value}")
