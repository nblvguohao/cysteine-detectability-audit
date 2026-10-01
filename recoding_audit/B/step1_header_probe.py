import csv, sys, os

import os  # v3.1.1: paths made repository-relative (were absolute Windows paths of the authoring machine)
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
SRC = os.path.join(_ROOT, 'supplemental', 'Supplemental_Data_9_survey_coding_table.csv')

with open(SRC, newline="", encoding="utf-8-sig") as f:
    r = csv.reader(f)
    header = next(r)
    n = sum(1 for _ in r)

print("ROWCOUNT_EXCL_HEADER:", n)
print("NCOL:", len(header))
for i, h in enumerate(header):
    print(i, repr(h))
