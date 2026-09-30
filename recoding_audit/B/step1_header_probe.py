import csv, sys, os

SRC = r"C:\Users\admin\Desktop\MCP\supplemental\Supplemental_Data_9_survey_coding_table.csv"

with open(SRC, newline="", encoding="utf-8-sig") as f:
    r = csv.reader(f)
    header = next(r)
    n = sum(1 for _ in r)

print("ROWCOUNT_EXCL_HEADER:", n)
print("NCOL:", len(header))
for i, h in enumerate(header):
    print(i, repr(h))
