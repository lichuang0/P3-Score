# -*- coding:utf-8 -*-
# python3
"""
Purpose: Predicting protein-ligand binding energy
usage: python3 P2-score.py --p protein_file --l ligand_file (--ac float(from 2 to 3))
author: Li Chuang
data: 2022-3-25
"""
from math import sqrt
import os
import click


class Get_info_tools(object):
    """
    Calculate the characteristic information of protein-ligand
    """

    # element's van der Waals radius
    dict_vdw_r = {"C.3": 1.70, "C": 1.77, "C.2": 1.70, "C.1": 1.77, "C.ar": 1.77, "C.cat": 1.70, "N.4": 1.55,
                  "N.3": 1.55, "N": 1.55, "O": 1.50, "S": 1.80,
                  "N.2": 1.55, "N.1": 1.60, "N.ar": 1.60, "N.am": 1.60, "N.pl3": 1.60, "O.3": 1.52, "O.2": 1.50,
                  "O.co2": 1.50, "S.3": 1.80, "S.2": 1.80, "S.o": 1.80, "S.o2": 1.80, "P.3": 1.88, "H.spc": 1.15,
                  "O.spc": 1.43, "H.t3p": 1.15, "O.t3p": 1.43, "H": 1.10, "Na": 2.38, "K": 2.52, "Mg": 2.00,
                  "Ca": 2.27, "Li": 2.14, "Al": 1.92, "Br": 1.85, "Cl": 1.75, "F": 1.50, "I": 1.97, "Si": 2.05,
                  "Se": 1.92, "Fe": 2.04, "Cu": 1.96, "Zn": 2.01, "Sn": 2.23, "Mo": 2.17, "Mn": 2.05, "Cr.oh": 2.06,
                  "Cr.th": 2.06, "Co.oh": 2.00}

    # Classification of amino acids
    """"
    Nonpolar amino acids are labeled as A1 ：ALA VAL LEU ILE PRO PHE TRP MET
    Polar and neutral amino acids are labeled as B1 ：GLY SER THR CYS TYR ASN GLN
                 Basic amino acids are labeled as B2 ：LYS ARG HIS
                 Acidic amino acids are labeled as B3 ：ASP GLU
    """
    dict_amino_acid = {"ALA": "A1", "VAL": "A1", "LEU": "A1", "ILE": "A1", "PRO": "A1", "PHE": "A1", "TRP": "A1",
                       "MET": "A1", "GLY": "B1", "SER": "B1", "THR": "B1", "CYS": "B1", "TYR": "B1", "ASN": "B1",
                       "GLN": "B1", "LYS": "B2", "ARG": "B2", "HIS": "B2", "ASP": "B3", "GLU": "B3", "HOH": "H2O"
                       }

    # Identify nitrogen and oxygen atoms
    list_N = ["N.4", "N.3", "N.2", "N.1", "N.ar", "N.am", "N.pl3"]
    list_O = ["O.3", "O.2", "O.co2", "O.spc", "O.t3p"]

    def __init__(self, file_name1, file_name2, num):
        """
        initialization
        :param file_name1: The name of the ligand file
        :param file_name2: The name of the protein file
        :param num: Accuracy of volume calculation
        """

        self.atom = []
        self.atom_polar = []
        self.ring = []
        self.ligand_name = []
        self.ligand_info, self.protein_info = [], []
        self.volume_list = []
        self.amino_acid = []
        self.xscore = []
        self.accuracy = float(num)
        self.name1 = file_name1
        self.name2 = file_name2
        self.get_file_info()
        self.get_xscore(self.name1, self.name2)

    def get_file_info(self):
        """
        Read the file information including ligand and proteins and save the necessary contents
        """

        # Get information on protein files
        with open(f"{self.name2}", 'r') as fi:
            while True:
                content = fi.readline().strip()
                if content == "":
                    break

                # Read information includes: element type, amino acid belonging to, number of amino acid, 3D coordinates
                if content[0:4] == "ATOM" or content[0:6] == "HETATM":
                    temp_info = [content[12:16].strip(), content[16:20].strip(), content[22:27].strip(),
                                 float(content[30:38].strip()), float(content[38:46].strip()),
                                 float(content[46:54].strip())]
                    self.protein_info.append(temp_info)
            fi.close()

        # Get information on ligand files
        with open(f"{self.name1}", 'r') as fi:
            temp_info = []
            read_state, exit_state, temp_num = 0, 0, 0
            read_state_list = ["@<TRIPOS>BOND\n", "@<TRIPOS>MOLECULE\n", "@<TRIPOS>ATOM\n"]
            while True:
                content = fi.readline()

                # Categorize what you want to read into various states
                if content in read_state_list:
                    read_state = read_state_list.index(content)

                    # Use the "molecular" label as the cut-off point when a complete ligand is read
                    # Launch tool for calculating ligand information
                    if read_state == 1:
                        if temp_info != []:
                            self.volume_list.append(self.get_volume(temp_info))
                            self.atom_polar.append(self.get_polar(temp_info))
                            self.amino_acid.append(self.get_ami(temp_info))
                        temp_info = []
                        temp_num = 0
                        exit_state = 0
                    continue

                # Get the name, number of atoms and number of rings of the ligand
                if read_state == 1:
                    temp_num += 1
                    if temp_num == 1:
                        self.ligand_name.append(content.strip())
                    elif temp_num == 2:
                        temp = content.strip().split()[0:2]
                        self.atom.append(int(temp[0]))
                        self.ring.append(int(temp[1]) - int(temp[0]) + 1)
                        continue

                # Read the 3D coordinates and atom type of ligands
                elif read_state == 2:

                    # Catch exception and set two read methods
                    try:
                        temp_info.append(
                            [float(content[16:26].strip()), float(content[26:36].strip()),
                             float(content[36:46].strip()),
                             content[46:53].strip()])
                    except:
                        list_content = content.strip().split()
                        temp_info.append(
                            [float(list_content[2]), float(list_content[3]), float(list_content[4]),
                             list_content[5]])

                # When the content cannot be read for twenty consecutive times, no longer read
                if content == "":
                    exit_state += 1
                    if exit_state == 20:
                        if temp_info == []:
                            self.volume_list.append(0)
                            self.atom_polar.append([0, 0])
                            self.amino_acid.append([0, 0, 0, 0, 0])
                            break
                        self.volume_list.append(self.get_volume(temp_info))
                        self.atom_polar.append(self.get_polar(temp_info))
                        self.amino_acid.append(self.get_ami(temp_info))
                        break
            fi.close()

    def get_volume(self, ligand):
        """
        Calculate the volume of the ligand
        :param ligand: List of ligand molecule information
        :return: The result of the calculation of the volume
        """

        # Save the three-dimensional coordinate information and the corresponding atom type respectively
        a, b, c, d_atom = [], [], [], []
        for i in ligand:
            d_atom.append(self.dict_vdw_r.get(i[3], 0))
        for i in ligand:
            a.append([i[0], d_atom[ligand.index(i)]])
            b.append([i[1], d_atom[ligand.index(i)]])
            c.append([i[2], d_atom[ligand.index(i)]])

        # Sort 3D coordinates by size respectively
        a1 = sorted(a, key=lambda a: a[0], reverse=False)
        b1 = sorted(b, key=lambda b: b[0], reverse=False)
        c1 = sorted(c, key=lambda c: c[0], reverse=False)

        # Calculate the length, width and height of the ligand
        len_x = round(a1[-1][0] - a1[0][0] + a1[-1][1] + a1[0][1], 4)
        len_y = round(b1[-1][0] - b1[0][0] + b1[-1][1] + b1[0][1], 4)
        len_z = round(c1[-1][0] - c1[0][0] + c1[-1][1] + c1[0][1], 4)

        # Computational Ligand Center
        self.center = [a1[0][0] - a1[0][1] + len_x * 0.5, b1[0][0] - b1[0][1] + len_y * 0.5,
                       c1[0][0] - c1[0][1] + len_z * 0.5, max([len_x, len_y, len_z]) * 0.5]

        # Determining the Accuracy of Calculations
        count1 = int(len_x * self.accuracy // 1)
        count2 = int(len_y * self.accuracy // 1)
        count3 = int(len_z * self.accuracy // 1)

        # Place the ligand in a 3D point set and count the points that fall inside the ligand
        num = 0
        count1 = int(len_x * self.accuracy // 1)
        count2 = int(len_y * self.accuracy // 1)
        count3 = int(len_z * self.accuracy // 1)
        d = [a1[0][0] - a1[0][1] - 1, b1[0][0] - b1[0][1] - 1, c1[0][0] - c1[0][1] - 1]
        for j1 in range(count1 + 1):
            d[0] += ((len_x + 2) / count1)
            for j2 in range(count2 + 1):
                d[1] += (len_y + 2) / count2
                for j3 in range(count3 + 1):
                    d[2] += (len_z + 2) / count3
                    for i in ligand:
                        if (i[0] - d[0]) ** 2 + (i[1] - d[1]) ** 2 + (i[2] - d[2]) ** 2 <= d_atom[ligand.index(i)] ** 2:
                            num += 1
                            break
                d[2] = c1[0][0] - c1[0][1]
            d[1] = b1[0][0] - b1[0][1]
        # Calculate volume
        volume_t = (len_x + 2) * (len_y + 2) * (len_z + 2)
        volume = round(num * volume_t / (count1 * count2 * count3), 3)

        return volume

    def get_polar(self, ligand):
        """
        Count the number of nitrogen and oxygen atoms
        :param ligand: List of ligand molecule information
        :return: List of nitrogen and oxygen numbers
        """
        num_N, num_O = 0, 0
        for i in ligand:
            if i[3] in self.list_N:
                num_N += 1
            elif i[3] in self.list_O:
                num_O += 1

        return [num_N, num_O]

    def get_ami(self, ligand):
        """
        Count and classify amino acids around ligands
        :param ligand: List of ligand molecule information
        :return: A list that records the number of four types of amino acids
        """
        # Obtain protein information around ligands to simplify calculations
        center_protein_info = []
        for k in self.protein_info:
            d = sqrt(
                (float(k[3]) - self.center[0]) ** 2 + (float(k[4]) - self.center[1]) ** 2 + (
                        float(k[5]) - self.center[2]) ** 2)
            if d < self.center[3] + 3:
                center_protein_info.append(k)

        # Record the 2.5 angstroms of amino acids around the ligand
        d_atom = []
        for i in ligand:
            d_atom.append(self.dict_vdw_r.get(i[3], 0))
        count, di = 0, {}
        for k in center_protein_info:
            for i in ligand:
                d = sqrt(
                    (float(k[3]) - float(i[0])) ** 2 + (float(k[4]) - float(i[1])) ** 2 + (
                            float(k[5]) - float(i[2])) ** 2)
                if d <= 2.5 + d_atom[ligand.index(i)]:
                    count += 1

                    # Prevent amino acid duplication records
                    di[k[2]] = k[1]
                    break

        # Count the number of various amino acids according to the recorded labels
        list_ami = [0, 0, 0, 0, 0]
        for k, v in di.items():
            a = self.dict_amino_acid.get(v, v)
            if a == "A1":
                list_ami[0] += 1
            elif a == "B1":
                list_ami[1] += 1
            elif a == "B2":
                list_ami[2] += 1
            elif a == "B3":
                list_ami[3] += 1
            elif a == "H2O":
                pass
            else:
                list_ami[4] += 1

        return list_ami

    def get_xscore(self, name1, name2):
        """
        Invoke xscore calculation and record features
        :param name1: The name of the ligand file
        :param name2: The name of the protein file
        :return: A list of records HMscore features
        """

        # Call it directly and check if it is working properly
        temp_path = os.getcwd()
        if os.path.exists(f"{temp_path}/xscore.log"):
            os.remove(f'{temp_path}/xscore.log')
        try:
            os.system(f"xscore -score {name2} {name1}")
            if os.path.exists(f"{temp_path}/xscore.log") == False:
                raise 1 / 0
        except:
            print("pleas input: python P3-Score_predict.py --help")
            print("Please make sure your Xscore is running.")

        # Read the features of HMscore
        else:
            with open('xscore.log', 'r') as fl:
                count = 0
                row = [1, 2, 4, 5, 6]
                while True:
                    temp = []
                    content = fl.readline().strip().split()
                    if content == []:
                        count += 1
                        if count == 10:
                            break
                        continue
                    if content[0] == "Total":
                        for i in row:
                            temp.append(float(content[i]))
                        self.xscore.append(temp)
                        count = 0

        return self.xscore


class calc(object):
    """
    Provides the calculation of the scoring function
    """
    row = [0, 1, 2, 3, 4, 5, 8, 10, 11, 12, 13, 14, 15, 16, 17, 24, 25, 26, 27, 28, 31, 33, 38, 40, 41, 46, 47, 50,
           61, 65, 66, 67, 68, 72, 74, 79, 86, 97, 100, 101, 103, 107, 108, 110, 111, 112, 115, 117, 119, 121, 128,
           132, 136, 142, 143, 146, 148, 155, 170, 171, 172, 177, 182, 187, 192, 193, 199, 204, 205, 213, 216, 217,
           219, 233, 236, 243, 258, 264, 270, 274, 278, 283, 289, 291, 295, 300, 304, 311, 312, 321, 336, 345, 346,
           348, 354, 358, 359, 360, 361, 372, 375, 376, 384, 385, 399, 403, 413, 427, 433, 435, 458, 478, 479, 480,
           481, 482, 485, 489, 492, 497, 506, 512, 521, 524, 537, 541, 542, 550, 554, 557, 563, 566, 570, 572, 573,
           575, 584, 589, 590, 606, 609, 611, 623, 637, 639, 651, 655, 662, 670, 671, 676]
    scale0 = [227.55524632769988, 2.973228393473277, 2.0126773360994417, 129.8330860295337, 4.716392837993995,
              151.50695268073721, 3.9224218409893914, 2.950142076415521, 1.6614690589883119, 1.472453623203242,
              0.4882252141187505, 296675.8439772615, 2390.585187665551, 1844.5247606647386, 126551.3094327611,
              2781.9925309489317, 1200.2875208040216, 1236.9258170988335, 234.5213735988161, 34.848395928370635,
              38.28959896573556, 9.908179313924798, 15.157885561767685, 3.1399184173117693, 13.483313120870669,
              11.083426240261288, 14.053004473222394, 5.630752038827854, 430.45563705487973, 3080.322669159909,
              16.758951023934355, 38.24044782306181, 55.55321832654137, 17.540759707131766, 119895.36813987223,
              1660.872858871134, 17.723170667216927, 1.505586271857734, 34.54530005320871, 19.0859676245226,
              3.074722028186388, 15.168746823178392, 2.79303047189354, 12.74564108451648, 11.426149161016449,
              2.673848858265198, 1.8490945014514202, 1.504661188139918, 337087085.34892535, 1898412.207337988,
              3629291.810627584, 158097.42547301843, 38041.15412350887, 25488.76082780236, 11978.258191779978,
              12882.741362166946, 27383.387475175983, 4764.707067539868, 3151660.860778027, 16602.493300107413,
              38029.23300803305, 17682.50323776452, 2088428.413472216, 85752.29091970452, 12689.70548282853,
              5000.737185630874, 22788.134660142245, 31859.780924320174, 31943.39541565111, 1595.0785494253412,
              10518.378756567166, 1672.0125274474592, 4195.299746673912, 325.28685635203124, 31.6595226589258,
              94.13719747034624, 2847.482400295553, 497.79685456422715, 909269.7764155063, 13400.548200964617,
              960.9263382086389, 92.08209189171593, 173.1215184528277, 81.62428759447226, 290.3434346393887,
              289.94943860858274, 17.141769579683107, 13.81196428429597, 53.02325948657239, 91.72734933992824,
              2822.275151812391, 71.50807438854775, 89.02412113751599, 766302.3128652817, 2990.3689415692647,
              47.57130476186433, 57.51867035356248, 101.74935312625712, 60.848408879205785, 141.8257067690038,
              41.585121972204846, 47.24760746026051, 42.49196546813215, 49.05464651925387, 762222.4110814757,
              18963.919804205874, 598.5427333391859, 4392.3774656091455, 8412.96537484996, 2330.592429156345,
              122.10235672885506, 85.08967381360814, 110.72095304127103, 198.64174073350821, 169.75195524146616,
              138.6179844046895, 6.8696045599333475, 386.3481025375999, 13.888196762274807, 183.67948603482253,
              161.8255969694386, 9.653872763918685, 463928.2365132873, 5953.249754497957, 574.95880618938,
              9526.664892671424, 7463.851039868166, 6342.460412876236, 2376.4686690119074, 408.57357248607127,
              77.99926364280881, 3.883247248946046, 94.15015204908993, 32.984845447234825, 4.328000320838915,
              139.54798040763157, 7.587530567862162, 32.115606622121255, 16.78428983496099, 74.80566582195964,
              203.3576732822364, 93.19774428879137, 1413.5764825137435, 22.40417729428934, 56.737323470894786,
              109.22148157601208, 10.30571517329792, 23.355889631818858, 31.283356349891193, 9.900689466784234,
              6.946949396305995]
    mean0 = [563.3200582100416, 3.9995149163230623, 2.2696361872422983, 171.9067911714767, 5.038927965073976,
             319.1441314576769, 4.9975745816153285, 5.459131700218288, 2.1765704584040746, 1.851564394858113,
             0.22580645161290322, 369110.87811302475, 2411.8634804753865, 1538.6111869997574, 114956.53021343648,
             3323.6922144069854, 1252.3133883094838, 1174.374193548384, 102.45590589376654, 24.836206645646264,
             25.349248120300807, 9.529105020616063, 10.753965559058953, 1.130244967256852, 9.202118481688144,
             7.85629396070822, 10.267501819063789, 4.363870967741932, 350.82789231142374, 2140.9271633519297,
             13.531772980839195, 20.659228716953674, 34.80875576036866, 11.367329614358477, 124807.33335447781,
             1860.1889219015306, 19.708949793839437, 0.5556633519282076, 29.756973077855932, 13.1040504487024,
             1.1695367450885277, 13.04487024011642, 1.2131942760126122, 12.067911714770798, 10.692456948823672,
             1.0904681057482415, 0.7613388309483385, 0.4921173902498181, 273765279.98989236, 1176401.236311038,
             3020734.1879674965, 55241.341380063146, 18809.429282076133, 16290.902626728137, 6641.46239631336,
             6923.476831414035, 12157.433113509594, 2985.706025709434, 1675024.4062378625, 10457.76517099199,
             16271.67815910742, 8788.523878243981, 1334924.306302936, 30371.79026808148, 10215.338782439969,
             3577.6517099199586, 13551.245282561218, 23116.307324763577, 20076.624011642005, 601.5072277467868,
             7245.341207858341, 564.4970652437544, 2627.3635217074934, 170.82125151588687, 8.050443851564392,
             43.59910308028131, 1424.412437545475, 217.8332524860539, 566047.8540014391, 9672.858692117392,
             313.71858763036664, 61.32794081979135, 100.02420567547902, 37.38120300751886, 172.37863206403065,
             209.5282803783656, 5.836623817608537, 4.064807179238416, 28.037181663836996, 44.56167278680574,
             1326.6257674023773, 30.996592287169594, 36.49887460586955, 433833.6814629463, 1798.211552255639,
             27.508273102110092, 33.66864176570464, 64.66860053359214, 39.60736114479736, 68.63499636187228,
             21.784292990540788, 22.882325976230923, 24.305212224108622, 25.405748241571608, 422018.00029347645,
             5271.970004850836, 147.87365995634246, 2990.2519282076178, 4943.090565122483, 1239.1602231384916,
             34.205748241571676, 51.967620664564635, 58.8254911472229, 97.0158864904196, 109.05760368663594,
             83.84210526315789, 1.8828522920203734, 144.60938636914867, 3.2667960223138492, 82.8649041959738,
             73.9140189182634, 2.475624545234053, 285101.19398070674, 4045.67643099685, 172.61789328159094,
             4914.601451612912, 4165.38503080282, 4292.399211253935, 1460.8416606839692, 128.5282813485326,
             55.6871210283774, 1.0434149890856173, 60.55687606112054, 20.93233082706767, 1.2318699975745817,
             99.12563667232598, 2.4319670143099685, 17.641280620907107, 10.533349502789232, 36.676206645646374,
             137.1290322580645, 51.88552025224351, 453.93524132912927, 6.615086102352656, 26.63667232597623,
             75.29177783167596, 4.107203492602474, 6.793111811787534, 14.645646373999515, 3.14649527043415,
             1.7125879214164443]
    coef = [1.1441663782499432, 0.6362632563970597, 0.33253447912590384, -0.18331956732261662, 0.32611170095566766,
            0.21370181260429325, -0.22620778244896292, 0.01748181869439621, 0.14748592890410492, 0.013837389850036769,
            0.2284668633423691, 0.30094884843060427, -0.34522835263792523, 0.19535978463960862, 0.4784858981554066,
            -0.037227232420745174, -0.1709838878603382, 0.26100767195978863, -0.24471936484221868, -0.0411826799569834,
            -0.33557623945614046, 0.20556589162253705, -0.2034684186312401, 0.3354167880118972, -0.24637761212732867,
            0.20894505878650657, 0.5177536306901687, 0.17867247231744113, -0.056550480542314326, -0.34220016044803697,
            -0.35759795974338227, -0.7561972422887501, -0.03923279437949557, -0.28828176640057557, -0.16810059209701209,
            -0.2457057081107467, 0.33475447658651525, 0.24020531697399522, 0.05812973242990665, -0.21574986656403236,
            -0.4444867559025378, 0.1798627080179703, 0.3257855317069273, -0.14224135716212, 0.18236229625179967,
            0.5488065919350236, 0.15648922285203348, -0.2985823268522644, -0.23209465259466497, -0.23400493641989356,
            -0.34602532153462967, -0.3024298582344142, 0.08568464720672193, -0.15487644296087474, 0.09245386770235073,
            -0.3079387382437739, -0.21289176740691554, -0.026202414557244164, -0.18077482380563795, 0.2481902543044735,
            0.49791030283794085, 0.4733864185797321, -0.1266317220860125, 0.3623946261861203, 0.19299131215823448,
            0.16273826946940367, -0.1709328305033089, 0.17948218886859266, -0.13909485283748962, -0.3311184518395256,
            -0.011555420400063817, -0.15552547575805323, -0.21286340241931345, 0.394750693007009, -0.26892137725519444,
            -0.15568293975676814, -0.20289903577712573, 0.06461356990603107, 0.21267687039845123, 0.17333207362528094,
            0.23803409599011882, -0.35198261612962684, -0.5099845148381258, 0.2342643844683067, -0.41655426104968374,
            0.23503975555916493, -0.2963956285579127, 0.16778712576240448, 0.17428585770150512, 0.4005703106463118,
            -0.3400071763981737, -0.08684607815372794, 0.3330292358292457, 0.017130454666645468, -0.27501765825975344,
            0.08308075784670002, -0.1902949948715159, -0.2603409229578681, -0.08210247207272314, -0.12134044482440598,
            -0.17174508179634543, 0.42705354592343203, 0.2585788961713932, -0.1896116240867361, 0.15089675951330955,
            -0.2684555328237243, 0.15560550133832685, 0.1759904382599455, -0.40913202430899603, -0.1077637481815873,
            -0.050932256856665706, -0.06946997360000177, 0.2220797691961497, 0.10161301470311854, 0.25288434173999225,
            -0.14126628663968294, -0.0773046080507367, 0.16343507568291898, -0.23499606284633812, 0.35034408135430906,
            -0.15392774371317683, 0.1433633023730116, -0.04873045770896722, 0.06960522278792992, 0.10392744519097917,
            0.63033080636774, -0.35266328578514616, 0.14256599497178377, -0.04056136955217704, 0.1483629938091759,
            0.07343473843474249, -0.14359719969936524, -0.09667304649129936, -0.09605202635556305, -0.24545620796243228,
            -0.35702824905195396, 0.4698884097131349, -0.20172520267926403, 0.09474665131798554, -0.04915413846172816,
            0.5812404219989707, 0.10352389578923392, 0.1762901752874345, 0.21422338283614653, -0.4885168791704469,
            -0.36325720696925035, 0.2666496877427053, -0.2908737215608939, 0.36613674115710576, -0.25534977122701547,
            -0.1812676960508424, 6.388217317487271]

    def __int__(self):
        pass

    @classmethod
    def get_predict(cls, features):
        """
        Predicting binding energy from features
        :param features: a list of record features
        :return: predicted binding energy
        """

        # Extend and trim features
        num = len(features)
        expand_features = []
        expand_features.extend(features[0:4])

        count = 3
        for i in range(4, num):
            count += 1
            if count in cls.row:
                expand_features.append(features[i])
        for i in range(num):
            for j in range(i, num):
                count += 1
                if count in cls.row:
                    expand_features.append(features[i] * features[j])
        for i in range(num):
            for j in range(i, num):
                for k in range(j, num):
                    count += 1
                    if count in cls.row:
                        expand_features.append(features[i] * features[j] * features[k])

        # Calculate the corresponding data
        predict = 0
        for i in range(len(expand_features)):
            predict += ((expand_features[i] - cls.mean0[i]) / cls.scale0[i]) * float(cls.coef[i])
        predict += float(cls.coef[-1])

        return round(predict, 4)


def integration(*args):
    """
    Provides integration with lists
    :param args: One-dimensional and two-dimensional lists
    :return: Consolidated list
    """

    integ_list = []
    for k in range(len(args[0])):
        temp = []
        for j in args:

            # Determine if a list is one-dimensional or two-dimensional
            if type(j[0]) is list:
                temp.extend(j[k])
            else:
                temp.append(j[k])
        integ_list.append(temp)

    return integ_list


def write_file(file_name, target, pridict_info):
    """
    Provide file write function
    :param file_name: the name of the file
    :param target: A list of record labels
    :param pridict_info: A two-dimensional list of recorded data
    """
    with open(f"{file_name}", "w")as fi:
        fi.write(",".join("%s" % i for i in target))
        fi.write("\n")
        for k in pridict_info:
            fi.write(",".join("%s" % i for i in k))
            fi.write("\n")
        fi.close()


def show(ligand_name, predict_list):
    """
    Provide printed content
    :param ligand_name: A list of recorded ligand names
    :param predict_list: a list of recorded binding energies
    :return:
    """

    if len(predict_list) == 1:
        print(f"****************     fix     ****************")
        print(f"Predict -log(Kd) = {predict_list[0]}\n")
        print(f"Predict binding energy = {round(predict_list[0] * (-1.3634), 4)}\n")
        print(f"The more information in the file of 'predict_info.txt'")
        print(f"**********************************************")
    else:
        print(f"****************     fix     ****************")
        print(f"id\t\tpredict(Pkd)\tbind_energy")
        for i in range(len(predict_list)):
            print("{:<15}{:<15}{:<8}".format(ligand_name[i], predict_list[i], round(predict_list[i] * (-1.3634), 4)))
        print(f"The more information in the file of 'predict_info.txt'")
        print(f"***********************************************")


@click.command()
@click.option("--protein_file", '--p', nargs=1, required=True, help="protein.pdb")
@click.option("--ligand_file", '--l', nargs=1, required=True, help="ligand.mol2")
@click.option('--ac', nargs=1, default=3, help="float or int")
def reception_and_display(protein_file, ligand_file, ac):
    """
    Provide all working procedures, call functions here
    :param protein_file: The name of the ligand file
    :param ligand_file: The name of the protein file
    :param ac: float form 2 to 3
    """

    # Extract features
    tool = Get_info_tools(ligand_file, protein_file, ac)

    # integration features
    predict_info = integration(tool.xscore, tool.volume_list, tool.ring, tool.atom_polar,
                               tool.amino_acid)

    # predicted binding energy
    predict_list = []
    for i in predict_info:
        predict_list.append(calc.get_predict(i))

    # write out data
    target = ["ligand_name", "VDW", "HB", "HM", "HS","RT", "volume", "num_ring", "num_N", "num_O", "ami_A1",
              "ami_B1", "ami_B2", "ami_B3","metal", "predict_PKi/Pkd"]
    new_predict_info = integration(tool.ligand_name, predict_info, predict_list)
    write_file(f"predict_info.txt", target, new_predict_info)

    # show the output
    show(tool.ligand_name, predict_list)


if __name__ == '__main__':
    reception_and_display()
