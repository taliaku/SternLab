#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file")


    (options, args) = parser.parse_args()
    file = options.file
    if "SIV" in file or "check" in file:
        print("SIV")
        return
    base = file.split("/")[-1].split("-translated")[0]
    directory = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHC_by_seq_250/{base}"
    #directory_temp = f"/sternadi/nobackup/volume1/talia_temp/netMHC_by_seq/{base}"
    if os.path.isfile(f"{directory}.all_alleles.csv"):
        print(f"csv exists")
        return
    if not os.path.isdir(directory):
        os.makedirs(directory)
    #if not os.path.isdir(directory_temp):
    #    os.makedirs(directory_temp)
    results = glob.glob(f"{directory}/*")
    print(len(results))
    if len(results) == 249:
        print("finished")
        return
    #alleles = open("/sternadi/home/volume1/taliakustin/SternLab/alleles.txt", "r").read()
    #alleles = alleles.split("\n")[:-1]
    alleles = ['HLA-A01:01', 'HLA-A02:01', 'HLA-A03:01', 'HLA-A11:01', 'HLA-A23:01', 'HLA-A24:02', 'HLA-A25:01', 'HLA-A26:01', 'HLA-A29:01', 'HLA-A30:01', 'HLA-A31:01', 'HLA-A32:01', 'HLA-A33:01', 'HLA-A34:01', 'HLA-A36:01', 'HLA-A43:01', 'HLA-A66:01', 'HLA-A68:01', 'HLA-A69:01', 'HLA-A74:01', 'HLA-A80:01', 'HLA-B07:02', 'HLA-B08:01', 'HLA-B13:01', 'HLA-B14:01', 'HLA-B15:01', 'HLA-B18:01', 'HLA-B27:01', 'HLA-B35:01', 'HLA-B37:01', 'HLA-B38:01', 'HLA-B39:01', 'HLA-B40:01', 'HLA-B41:01', 'HLA-B42:01', 'HLA-B44:02', 'HLA-B45:01', 'HLA-B46:01', 'HLA-B47:01', 'HLA-B48:01', 'HLA-B49:01', 'HLA-B50:01', 'HLA-B51:01', 'HLA-B52:01', 'HLA-B53:01', 'HLA-B54:01', 'HLA-B55:01', 'HLA-B56:01', 'HLA-B57:01', 'HLA-B58:01', 'HLA-B59:01', 'HLA-B67:01', 'HLA-B73:01', 'HLA-B78:01', 'HLA-B81:01', 'HLA-B82:01', 'HLA-B83:01', 'HLA-C01:02', 'HLA-C02:02', 'HLA-C03:01', 'HLA-C04:01', 'HLA-C05:01', 'HLA-C06:02', 'HLA-C07:01', 'HLA-C08:01', 'HLA-C12:02', 'HLA-C14:02', 'HLA-C15:02', 'HLA-C16:01', 'HLA-C17:01', 'HLA-C18:01', 'Mamu-A01', 'Mamu-A1:00101', 'Mamu-A2:0101', 'Mamu-A3:1301', 'Mamu-A4:0101', 'Mamu-A5:30010', 'Mamu-A6:0101', 'Mamu-A7:0101', 'Mamu-AG:01', 'Mamu-B:00101', 'Patr-A0101', 'Patr-A0201', 'Patr-A0301', 'Patr-A0302', 'Patr-A0401', 'Patr-A0402', 'Patr-A0404', 'Patr-A0501', 'Patr-A0601', 'Patr-A0602', 'Patr-A0701', 'Patr-A0801', 'Patr-A0802', 'Patr-A0803', 'Patr-A0901', 'Patr-A0902', 'Patr-A1001', 'Patr-A1101', 'Patr-A1201', 'Patr-A1301', 'Patr-A1401', 'Patr-A1501', 'Patr-A1502', 'Patr-A1601', 'Patr-A1701', 'Patr-A1702', 'Patr-A1703', 'Patr-A1801', 'Patr-A2301', 'Patr-A2401', 'Patr-B0101', 'Patr-B0102', 'Patr-B0201', 'Patr-B0203', 'Patr-B0301', 'Patr-B0302', 'Patr-B0401', 'Patr-B0402', 'Patr-B0501', 'Patr-B0502', 'Patr-B0601', 'Patr-B0701', 'Patr-B0702', 'Patr-B0801', 'Patr-B0802', 'Patr-B0901', 'Patr-B1001', 'Patr-B1101', 'Patr-B1102', 'Patr-B1202', 'Patr-B1301', 'Patr-B1401', 'Patr-B1601', 'Patr-B1602', 'Patr-B1701', 'Patr-B1702', 'Patr-B1703', 'Patr-B1801', 'Patr-B1901', 'Patr-B2001', 'Patr-B2101', 'Patr-B2201', 'Patr-B2202', 'Patr-B2301', 'Patr-B2302', 'Patr-B2303', 'Patr-B2401', 'Patr-B2402', 'Patr-B2501', 'Patr-B2601', 'Patr-B2701', 'Patr-B2801', 'Patr-B2901', 'Patr-B3001', 'Patr-B3501', 'Patr-B3601', 'Patr-B3701', 'Patr-C0201', 'Patr-C0202', 'Patr-C0203', 'Patr-C0204', 'Patr-C0205', 'Patr-C0206', 'Patr-C0301', 'Patr-C0302', 'Patr-C0303', 'Patr-C0304', 'Patr-C0401', 'Patr-C0501', 'Patr-C0502', 'Patr-C0601', 'Patr-C0701', 'Patr-C0801', 'Patr-C0901', 'Patr-C0902', 'Patr-C0903', 'Patr-C0904', 'Patr-C0905', 'Patr-C1001', 'Patr-C1101', 'Patr-C1201', 'Patr-C1301', 'Patr-C1302', 'Patr-C1501', 'Patr-C1601', 'SLA-1:0101', 'SLA-2:0101', 'SLA-3:0101', 'SLA-6:0101', 'SLA-1-CHANGDA', 'SLA-1-HB01', 'SLA-1-HB02', 'SLA-1-HB03', 'SLA-1-HB04', 'SLA-1-LWH', 'SLA-1-TPK', 'SLA-1-YC', 'SLA-1-YDL01', 'SLA-1-YTH', 'SLA-2-YDL02', 'SLA-3-CDY', 'SLA-3-HB01', 'SLA-3-LWH', 'SLA-3-TPK', 'SLA-3-YC', 'SLA-3-YDL', 'SLA-3-YDY01', 'SLA-3-YDY02', 'SLA-3-YTH', 'H-2-Db', 'H-2-Dd', 'H-2-Kb', 'H-2-Kd', 'H-2-Kk', 'H-2-Ld', 'H-2-Qa2', 'H-2-Qa1', 'H2-Db', 'H2-Dd', 'H2-Kb', 'H2-Kd', 'H2-Kk', 'H2-Ld', 'H2-Qa2', 'H2-Qa1', 'Gogo-B0101', 'HLA-G01:01', 'HLA-E01:01', 'BoLA-N:00101', 'BoLA-NC1:00101', 'BoLA-NC2:00101', 'BoLA-NC3:00101', 'BoLA-NC4:00101', 'BoLA-HD6', 'BoLA-JSP.1', 'BoLA-T2c', 'BoLA-T2b', 'BoLA-T2a', 'BoLA-T7', 'BoLA-D18.4', 'BoLA-AW10', 'BoLA-T5', 'BoLA-1:00901', 'BoLA-2:00501', 'BoLA-3:00101', 'BoLA-4:02401', 'BoLA-5:00301', 'BoLA-6:01301']
    len_alleles = len(alleles)
    for a in (alleles):
        output = f"{directory}/{base}.{a}.csv"
        #output_temp = f"{directory_temp}/{base}.{a}.csv"
        if os.path.isfile(output):
            print(f"output exists {output}")
            continue
        os.system(f"/sternadi/home/volume1/taliakustin/software/netMHCpan-4.0/netMHCpan {file} -t 0.1 -l 9 -a {a} -xls -xlsfile {output}")








if __name__ == "__main__":
    main()

