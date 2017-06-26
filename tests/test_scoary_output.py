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

referencevcf = ["NC_000962", "4013","0","T","C","9999","0","TYPE=snp","GT","False","0","1",
                "1","1"]




for Test in ["1","2","4"]:
    with open(os.getcwd() + "/Test" + Test + "/Tetracycline_resistance.results.csv" ,"rU") as resfile:
        tab = csv.reader(resfile, delimiter=",")
        for i in range(2):
            if i == 0:
                next(tab)
            if i == 1:
                data = next(tab)
                try:
                    assert data[0] == reference[0]
                except AssertionError:
                    print("Not equal at Test %s col 0: %s %s" % (Test, data[0], reference[0]))
                    sys.exit(-1)
                try:
                    assert data[1] == reference[1]
                except AssertionError:
                    print("Not equal at Test %s col 1: %s %s" % (Test, data[1], reference[1]))
                    sys.exit(-1)
                try:
                    assert data[2] == reference[2]
                except AssertionError:
                    print("Not equal at Test %s col 2: %s %s" % (Test, data[2], reference[2]))
                    sys.exit(-1)
                try:
                    assert int(data[3]) == reference[3]
                except AssertionError:
                    print("Not equal at Test %s col 3: %s %s" % (Test, data[3], str(reference[3])))
                    sys.exit(-1)
                try:
                    assert int(data[4]) == reference[4]
                except AssertionError:
                    print("Not equal at Test %s col 4: %s %s" % (Test, data[4], str(reference[4])))
                    sys.exit(-1)
                try:
                    assert int(data[5]) == reference[5]
                except AssertionError:
                    print("Not equal at Test %s col 5: %s %s" % (Test, data[5], str(reference[5])))
                    sys.exit(-1)
                try:
                    assert int(data[6]) == reference[6]
                except AssertionError:
                    print("Not equal at Test %s col 6: %s %s" % (Test, data[6], str(reference[6])))
                    sys.exit(-1)
                try:
                    assert abs(float(data[7]) - reference[7]) <= 0.01
                except AssertionError:
                    print("Not equal at Test %s col 7: %s %s" % (Test, data[7], str(reference[7])))
                    sys.exit(-1)
                try:
                    assert abs(float(data[8]) - reference[8]) <= 0.01
                except AssertionError:
                    print("Not equal at Test %s col 8: %s %s" % (Test, data[8], str(reference[8])))
                    sys.exit(-1)
                try:
                    assert abs(float(data[9]) - reference[9]) <= 0.1
                except AssertionError:
                    print("Not equal at Test %s col 9: %s %s" % (Test, data[9], str(reference[9])))
                    sys.exit(-1)
                try:
                    assert abs(float(data[10]) - reference[10]) <= 1E-15
                except AssertionError:
                    print("Not equal at Test %s col 10: %s %s" % (Test, data[10], str(reference[10])))
                    sys.exit(-1)
                try:
                    assert abs(float(data[11]) - reference[11]) <= 1E-12
                except AssertionError:
                    print("Not equal at Test %s col 11: %s %s" % (Test, data[11], str(reference[11])))
                    sys.exit(-1)
                try:
                    assert abs(float(data[12]) - reference[12]) <= 1E-12
                except AssertionError:
                    print("Not equal at Test %s col 12: %s %s" % (Test, data[12], str(reference[12])))
                    sys.exit(-1)

                if not Test == "2":
                    try:
                        assert int(data[13]) == reference[13]
                    except AssertionError:
                        print("Not equal at Test %s col 13: %s %s" % (Test, data[13], str(reference[13])))
                        sys.exit(-1)
                    try:
                        assert int(data[14]) == reference[14]
                    except AssertionError:
                        print("Not equal at Test %s col 14: %s %s" % (Test, data[14], str(reference[14])))
                        sys.exit(-1)
                    try:
                        assert int(data[15]) == reference[15]
                    except AssertionError:
                        print("Not equal at Test %s col 15: %s %s" % (Test, data[15], str(reference[15])))
                        sys.exit(-1)
                    try:
                        assert abs(float(data[16]) - reference[16]) <= 1E-9
                    except AssertionError:
                        print("Not equal at Test %s col 16: %s %s" % (Test, data[16], str(reference[16])))
                        sys.exit(-1)
                    try:
                        assert abs(float(data[17]) - reference[17]) <= 1E-7
                    except AssertionError:
                        print("Not equal at Test %s col 17: %s %s" % (Test, data[17], str(reference[17])))
                        sys.exit(-1)

with open(os.getcwd() + "/mutations_presence_absence.csv" ,"rU") as vcfresfile:
    tab = csv.reader(vcfresfile, delimiter=",")
    for i in range(2):
        if i == 0:
            next(tab)
        if i == 1:
            data = next(tab)
            try:
                assert data == referencevcf
            except AssertionError:
                print("Got: %s" % ",".join(data))
                print("Expected: %s" % ",".join(referencevcf))
                print("VCF conversion did not produce the expected output")
                sys.exit(-1)

