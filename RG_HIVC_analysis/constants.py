gag_ET86_interval = (170, 1684)
pol_ET86_interval = (1456, 4488)
IN_ET86_interval = (3619, 4483)
PR_ET86_interval = ()
RT_ET86_interval = (1939, 3619)
RT_short_ET86_interval = (1939, 2629) # 230 aa's
env_ET86_interval = (5637, 8196)
ET86_length= 9031

gag_HXB2_interval = (790, 2292)
pol_HXB2_interval = (2085, 5096)
env_HXB2_interval = (6225, 8795)

orig_samples= ['100888_S14','130945_S2','504185_S29','504190_S34','504194_S38','504198_S42','504202_S46','504206_S50','504210_S54','504215_S59','504221_S65','504226_S70','504231_S75','79504_S23','83476_S27','84864_S47','TASPX100926_S26','X145364_S75','X83354_S92','105090_S50','504181_S25','504186_S30','504191_S35','504195_S39','504199_S43','504203_S47','504207_S51','504211_S55','504217_S61','504223_S67','504227_S71','504233_S77','81065_S64','84298_S2','87415_S75','TASPX119494_S74','X145656_S76','X84335_S20','105094_S45','504182_S26','504188_S32','504192_S36','504196_S40','504200_S44','504204_S48','504208_S52','504212_S56','504218_S62','504224_S68','504228_S72','504234_S78','83351_S15','84434_S32','TASPX100702_S59','X100748_S73','X160138_S81','X84434_S3','105257_S39','504184_S28','504189_S33','504193_S37','504197_S41','504201_S45','504205_S49','504209_S53','504214_S58','504220_S64','504225_S69','504230_S74','504235_S79','83456_S87','84785_S56','TASPX100711_S29','X145364-R_S95','X83322_S89','X84994_S90']
orig_excluded_samples = ('X84335_S20', '504214_S58', '504184_S28', '504186_S30', '84864_S47', '504206_S50', '504190_S34', '504191_S35','504192_S36', '504198_S42', 'X84434_S3', 'X145364-R_S95')
orig_high_quality_patients = ('13003', '15664', '16207', '22097', '22763', '22828', '26892', '29447', '31254', '47939')
orig_high_quality_samples = ('504228_S72','504212_S56','504234_S78','504225_S69','504197_S41','504207_S51','504210_S54','504189_S33','504218_S62','TASPX100711_S29','504195_S39','504202_S46','504193_S37','504217_S61','504205_S49','504233_S77','504194_S38','504223_S67','504220_S64','504231_S75','504215_S59','504192_S36','504201_S45','504230_S74','504209_S53','504208_S52','504198_S42','504204_S48')
orig_first_samples= ('81065_S64', '504234_S78', '87415_S75', '504230_S74', '504203_S47', '84298_S2', '84434_S32', '84785_S56', 'X83322_S89', '504233_S77', '83476_S27', 'X84994_S90', '504235_S79', '83456_S87', 'X83354_S92', '105090_S50', '83351_S15', '105257_S39',)

# TODO- can be updated one sort the missing\mixed sample names
control_patients = ['78292_S1', 'X105350_S13', 'X108456_S11', 'X117807_S26', 'X122087_S18', 'X128277_S27', 'X135623_S20', 'X161610_S28', 'X165309_S33', 'X83744_S2', 'X96026_S5', 'X101530_S9', 'X105354_S14', 'X111406_S10', 'X122054_S15', 'X122107_S17', 'X130823_S22', 'X160392_S24', 'X165276_S32', 'X178870_S29', 'X87707_S25', 'X96051_S7', 'X102662_S8', 'X108262_S21', 'X111489_S12', 'X122061_S16', 'X122107_S27', 'X132700_S84', 'X160433_S25', 'X165290_S34', 'X83348_S4', 'X88469_S6', '78486_S9', '79649_S49', '81088_S63', 'X105130_S8', 'X105416_S169', 'X112003_S13', 'X117597_S55', 'X117659_S88', 'X119113_S320', 'X122207_S32', 'X160394_S25', 'X84355_S11', 'X84505_S40', 'X84772_S12', 'X88273_S273', 'X94731_S49', '79578', '81003_S48', 'X101706', 'X105289_S56', 'X105500_S302', 'X117547_S16', 'X117607_S11', 'X117890_S46', 'X121017_S309', 'X135606_S69', 'X83676_S242', 'X84469_S265', 'X84641_S47', 'X84990_S248', 'X88373_S267', 'X96465_S54']
control_excluded_patients = ['18355', '11484', '27416', '4012', '6402', '4887', '21369', '21134', '5609']
control_excluded_samples  = ['X117597','X105416','X84505','X160394','X105130','79578','X101706','X102662','X122107','X160392','X165273','X132700','X87707']
# control_chosen_patients = ['24277','6773','26755','4956','7965','8992','22992','18901','1689','13694','4845','15687','8516','8670','21859','14201','324','7878','4179','12643']
control_chosen_patients = ['24277', '6773', '26755', '4956', '7965', '8992', '22992', '1689', '13694', '4845', '15687', '8670', '14201', '324', '7878', '4179']

control_bl_samples = ['X112003','X121017','X117597','X84355','X84469','X117607','X117659','X160394','X84772','X88373','81003','X105416','X105500','X83676','X84641','X105289','X117547','X135606','79578','X94731','X101706','X105130','X96465','X84990','X84505','79649','78486','X119113','X88273','X122207','X117890','81088']
control_fl_samples = ['X117807','X135623','X165309','X96026','X122107','X161610','X165273','X165290','X83348','X101530','X102662','X105350','X105354','X128277','X130823','X160433','X165276','X178870','X83744','X96051','X122054','X122087','X122061','X111406','X88469','X108456','78292','X132700','X111489','X108262','X160392','X87707']

# Adequte positions for MR analysis according to ZN (synonymous unconserved positions - filtering details below)
# TODO- produce this array per patient
orig_high_syn_unconserved_pos_for_mr= [174, 193, 198, 201, 211, 213, 228, 229, 240, 250, 252, 259, 261, 267, 268, 273, 285, 315, 318, 324, 333, 336, 345, 349, 351, 352, 354, 358, 363, 370, 372, 381, 387, 391, 396, 405, 408, 414, 420, 423, 436, 438, 439, 441, 447, 451, 462, 468, 474, 477, 478, 495, 499, 501, 507, 508, 510, 513, 516, 519, 522, 525, 528, 531, 534, 537, 539, 540, 543, 558, 561, 567, 570, 574, 576, 594, 597, 600, 601, 603, 609, 615, 618, 636, 642, 645, 648, 651, 654, 660, 663, 685, 705, 711, 712, 717, 724, 732, 735, 759, 771, 774, 783, 795, 801, 804, 805, 816, 819, 822, 831, 840, 844, 846, 849, 852, 855, 867, 904, 906, 918, 921, 930, 936, 942, 945, 964, 969, 984, 1002, 1026, 1032, 1038, 1047, 1056, 1071, 1080, 1092, 1098, 1119, 1140, 1155, 1176, 1185, 1215, 1245, 1248, 1269, 1275, 1278, 1281, 1284, 1288, 1296, 1302, 1305, 1308, 1311, 1314, 1315, 1317, 1326, 1329, 1335, 1338, 1356, 1359, 1360, 1365, 1371, 1377, 1380, 1386, 1392, 1413, 1422, 1434, 1440, 1443, 1446, 1694, 1697, 1698, 1703, 1706, 1712, 1718, 1748, 1751, 1754, 1757, 1763, 1766, 1787, 1799, 1805, 1817, 1823, 1830, 1832, 1847, 1868, 1869, 1874, 1877, 1880, 1901, 1904, 1907, 1908, 1919, 1925, 1928, 1931, 1946, 1949, 1952, 1958, 1970, 1974, 1976, 2000, 2003, 2018, 2027, 2051, 2057, 2060, 2066, 2111, 2126, 2138, 2141, 2147, 2150, 2174, 2180, 2198, 2207, 2210, 2219, 2228, 2237, 2249, 2258, 2264, 2267, 2276, 2288, 2294, 2297, 2309, 2312, 2315, 2318, 2339, 2348, 2363, 2369, 2384, 2387, 2393, 2405, 2420, 2423, 2441, 2447, 2456, 2460, 2465, 2471, 2472, 2477, 2480, 2489, 2501, 2516, 2528, 2549, 2552, 2559, 2561, 2564, 2565, 2567, 2570, 2573, 2582, 2585, 2594, 2597, 2615, 2624, 2663, 2675, 2678, 2684, 2687, 2690, 2714, 2741, 2753, 2759, 2762, 2765, 2771, 2774, 2775, 2777, 2795, 2805, 2813, 2825, 2849, 2855, 2868, 2888, 2900, 2912, 2915, 2921, 2925, 2936, 2939, 2942, 2945, 2948, 2966, 2972, 2978, 2981, 2984, 2985, 2987, 2990, 2996, 2999, 3014, 3017, 3059, 3062, 3068, 3069, 3071, 3095, 3098, 3111, 3119, 3122, 3125, 3128, 3137, 3138, 3152, 3179, 3185, 3191, 3197, 3204, 3221, 3224, 3225, 3227, 3230, 3236, 3242, 3257, 3284, 3287, 3290, 3294, 3296, 3299, 3308, 3314, 3320, 3329, 3332, 3335, 3341, 3344, 3345, 3368, 3374, 3377, 3380, 3386, 3389, 3390, 3395, 3398, 3401, 3411, 3422, 3434, 3440, 3452, 3455, 3458, 3467, 3476, 3488, 3489, 3494, 3497, 3500, 3509, 3512, 3519, 3521, 3527, 3530, 3533, 3539, 3557, 3587, 3596, 3602, 3608, 3611, 3617, 3620, 3632, 3638, 3641, 3650, 3653, 3662, 3695, 3702, 3704, 3710, 3734, 3752, 3753, 3761, 3770, 3785, 3791, 3794, 3797, 3812, 3821, 3822, 3830, 3836, 3840, 3842, 3851, 3866, 3872, 3881, 3884, 3887, 3905, 3908, 3917, 3919, 3920, 3924, 3950, 3956, 3962, 3974, 3977, 3986, 3998, 4001, 4022, 4028, 4040, 4079, 4103, 4109, 4112, 4121, 4130, 4199, 4214, 4226, 4236, 4241, 4250, 4257, 4265, 4316, 4319, 4340, 4343, 4361, 4370, 4379, 4382, 4388, 4394, 4400, 4406, 4409, 4421, 4430, 4501, 4509, 4515, 4524, 4530, 4537, 4539, 4540, 4548, 4551, 4566, 4569, 4602, 4606, 4614, 4617, 4644, 4645, 4647, 4650, 4662, 4665, 4671, 4702, 4704, 4710, 4713, 4716, 4728, 4734, 4735, 4743, 4746, 4752, 4758, 4761, 4782, 4794, 4801, 4812, 4821, 4833, 4851, 4863, 4864, 4872, 4873, 4888, 4890, 4896, 4902, 4911, 4930, 4932, 4935, 4944, 4948, 5012, 5013, 5015, 5016, 5021, 5024, 5027, 5030, 5033, 5048, 5058, 5060, 5066, 5072, 5073, 5075, 5078, 5081, 5084, 5087, 5091, 5093, 5096, 5102, 5105, 5108, 5114, 5120, 5126, 5127, 5129, 5130, 5138, 5139, 5147, 5148, 5150, 5151, 5159, 5177, 5180, 5183, 5195, 5198, 5199, 5201, 5202, 5204, 5207, 5241, 5242, 5247, 5256, 5277, 5280, 5286, 5292, 5295, 5298, 5301, 5307, 5310, 5316, 5319, 5323, 5391, 5397, 5399, 5406, 5418, 5421, 5422, 5424, 5427, 5461, 5463, 5466, 5467, 5469, 5470, 5472, 5473, 5475, 5476, 5478, 5479, 5484, 5490, 5491, 5505, 5506, 5508, 5509, 5512, 5518, 5520, 5521, 5523, 5524, 5535, 5538, 5544, 5550, 5553, 5554, 5565, 5569, 5571, 5572, 5583, 5592, 5593, 5595, 5604, 5613, 5719, 5722, 5723, 5725, 5728, 5731, 5732, 5734, 5746, 5749, 5764, 5773, 5779, 5786, 5791, 5818, 5819, 5821, 5824, 5827, 5836, 5842, 5857, 5872, 5882, 5884, 5887, 5888, 5890, 5893, 5908, 5914, 5921, 5929, 5935, 5938, 5941, 5950, 5959, 5963, 5971, 5977, 5983, 5995, 5996, 6013, 6016, 6017, 6022, 6025, 6026, 6027, 6028, 6029, 6031, 6032, 6034, 6037, 6040, 6041, 6043, 6046, 6047, 6049, 6052, 6053, 6055, 6058, 6061, 6062, 6064, 6065, 6066, 6067, 6070, 6073, 6076, 6077, 6079, 6080, 6082, 6085, 6086, 6088, 6090, 6094, 6097, 6109, 6118, 6119, 6124, 6130, 6131, 6133, 6134, 6136, 6137, 6139, 6142, 6148, 6149, 6151, 6160, 6163, 6167, 6172, 6178, 6181, 6182, 6183, 6184, 6187, 6190, 6193, 6196, 6214, 6226, 6227, 6235, 6253, 6259, 6262, 6274, 6286, 6292, 6304, 6313, 6316, 6317, 6319, 6331, 6343, 6346, 6352, 6355, 6367, 6376, 6382, 6388, 6391, 6401, 6403, 6404, 6415, 6418, 6419, 6430, 6433, 6439, 6442, 6445, 6455, 6457, 6463, 6469, 6472, 6484, 6487, 6490, 6493, 6494, 6496, 6499, 6505, 6509, 6511, 6514, 6517, 6520, 6523, 6538, 6539, 6544, 6547, 6550, 6571, 6577, 6581, 6582, 6583, 6586, 6598, 6613, 6619, 6625, 6626, 6628, 6631, 6634, 6640, 6643, 6646, 6647, 6652, 6653, 6655, 6659, 6661, 6664, 6667, 6671, 6673, 6676, 6679, 6682, 6685, 6686, 6688, 6691, 6694, 6698, 6700, 6704, 6706, 6709, 6712, 6721, 6724, 6725, 6739, 6745, 6755, 6763, 6775, 6785, 6787, 6788, 6790, 6796, 6799, 6802, 6803, 6804, 6805, 6806, 6808, 6809, 6810, 6811, 6812, 6813, 6814, 6816, 6817, 6820, 6821, 6823, 6824, 6826, 6827, 6829, 6832, 6833, 6834, 6835, 6837, 6838, 6841, 6844, 6847, 6850, 6853, 6856, 6862, 6874, 6883, 6884, 6886, 6893, 6901, 6907, 6913, 6922, 6925, 6937, 6943, 6953, 6956, 6958, 6959, 6961, 6965, 6967, 6973, 6976, 6977, 6979, 6981, 6982, 6983, 6984, 6985, 6987, 6988, 6991, 6994, 6995, 6997, 7000, 7009, 7024, 7030, 7033, 7054, 7060, 7066, 7078, 7084, 7088, 7093, 7096, 7099, 7100, 7102, 7129, 7141, 7142, 7381, 7382, 7384, 7390, 7420, 7423, 7435, 7447, 7454, 7455, 7456, 7459, 7462, 7471, 7480, 7489, 7495, 7507, 7510, 7519, 7528, 7529, 7531, 7534, 7540, 7543, 7546, 7561, 7562, 7564, 7573, 7576, 7580, 7582, 7594, 7600, 7603, 7604, 7606, 7621, 7630, 7634, 7636, 7654, 7657, 7660, 7675, 7681, 7684, 7696, 7702, 7703, 7705, 7714, 7717, 7741, 7742, 7744, 7747, 7750, 8002, 8006, 8008, 8011, 8014, 8015, 8018, 8020, 8023, 8032, 8033, 8035, 8038, 8039, 8059, 8060, 8063, 8071, 8072, 8092, 8104, 8107, 8108, 8110, 8113, 8114, 8116, 8117, 8123, 8125, 8131, 8140, 8141, 8143, 8146, 8147, 8152, 8155, 8167, 8170, 8171, 8173, 8182, 8204, 8213, 8214, 8216, 8217, 8219, 8222, 8225, 8228, 8234, 8237, 8240, 8246, 8250, 8252, 8255, 8256, 8262, 8263, 8264, 8266, 8276, 8285, 8291, 8294, 8296, 8300, 8306, 8309, 8315, 8318, 8319, 8321, 8324, 8327, 8330, 8333, 8336, 8337, 8339, 8342, 8345, 8351, 8354, 8360, 8364, 8366, 8372, 8375, 8378, 8381, 8384, 8387, 8390, 8396, 8405, 8408, 8411, 8414, 8417, 8420, 8421, 8426, 8435, 8438, 8441, 8447, 8450, 8456, 8489, 8493, 8499, 8501, 8504, 8505, 8508, 8509, 8510, 8513, 8516, 8519, 8522, 8531, 8537, 8552, 8564, 8570, 8582, 8600, 8603, 8604, 8615, 8621, 8627, 8628, 8630, 8633, 8639, 8642, 8646, 8648, 8651, 8654, 8660, 8663, 8666, 8669, 8675, 8678, 8681, 8684, 8685, 8687, 8688, 8690, 8696, 8699, 8705, 8708, 8714, 8717, 8723, 8726, 8727, 8735, 8736, 8738, 8739, 8741, 8756, 8759, 8762, 8774, 8775, 8777, 8783, 8787, 8789, 8798, 8801, 8807, 8810, 8813, 8815]
# print(len(orig_high_syn_unconserved_pos_for_mr))

# approximate stats for these positions-
#  per pos-
#     zero_proteins_filter: 416
#     two_or_more_proteins_filter: 697
#     excluded_protein_filter: 0
#     gaps_filter: 36
#     rna_filter:
#     entropy_filter: 5826   (above 0.3 entropy according to Subtype-C entropy analysis (not conserved))
#  per base-
#     derived_filter: 1762 # TODO- needs to be verified, after running improved make_ref_from_con
#     non_syn_filter:
#  per timepoint-
#     coverage:

