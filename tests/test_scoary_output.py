# Script for verifying Scoary Test output

# Pseudocode:
# For all test scenarios
# - Open results file
# - Verify top hit

import sys
import os
import csv

reference = ["TetRCG","","A fictitious gene known to cause resistance against tetracycline",
             29,8,3,60,90.625,88.2352941176,72.5,1.08621066108E-014,6.45209132679E-011,
             6.45209132679E-011,25,25,1,5.96046447754E-008,1.54972076416E-006]


for Test in ["1","2","3","4"]:
    with open(os.getcwd() + "/Test" + Test + "/Tetracycline_resistance.results.csv" ,"rU") as resfile:
        tab = csv.reader(resfile, delimiter=",")
        for i in range(2):
            if i == 0:
                tab.next()
                continue
            if i == 1:
                data = tab.next()
                assert data[0] == reference[0]
                assert data[1] == reference[1]
                assert data[2] == reference[2]
                assert int(data[3]) == reference[3]
                assert int(data[4]) == reference[4]
                assert int(data[5]) == reference[5]
                assert int(data[6]) == reference[6]
                assert abs(float(data[7]) - reference[7]) <= 0.01
                assert abs(float(data[8]) - reference[8]) <= 0.01
                assert abs(float(data[9]) - reference[9]) <= 0.1
                assert abs(float(data[10]) - reference[10]) <= 1E-15
                assert abs(float(data[11]) - reference[11]) <= 1E-12
                assert abs(float(data[12]) - reference[12]) <= 1E-12
                if not Test == "2":
                    assert int(data[13]) == reference[13]
                    assert int(data[14]) == reference[14]
                    assert int(data[15]) == reference[15]
                    assert abs(float(data[16]) - reference[16]) <= 1E-9
                    assert abs(float(data[17]) - reference[17]) <= 1E-7
