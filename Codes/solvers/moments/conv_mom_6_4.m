function cent = conv_mom_6_4(m100000,m010000,m001000,m000100,m000010,m000001,m200000,m110000,m020000,m101000,m011000,m002000,m100100,m010100,m001100,m000200,m100010,m010010,m001010,m000110,m000020,m100001,m010001,m001001,m000101,m000011,m000002,m300000,m210000,m120000,m030000,m201000,m111000,m021000,m102000,m012000,m003000,m200100,m110100,m020100,m101100,m011100,m002100,m100200,m010200,m001200,m000300,m200010,m110010,m020010,m101010,m011010,m002010,m100110,m010110,m001110,m000210,m100020,m010020,m001020,m000120,m000030,m200001,m110001,m020001,m101001,m011001,m002001,m100101,m010101,m001101,m000201,m100011,m010011,m001011,m000111,m000021,m100002,m010002,m001002,m000102,m000012,m000003,m400000,m310000,m220000,m130000,m040000,m301000,m211000,m121000,m031000,m202000,m112000,m022000,m103000,m013000,m004000,m300100,m210100,m120100,m030100,m201100,m111100,m021100,m102100,m012100,m003100,m200200,m110200,m020200,m101200,m011200,m002200,m100300,m010300,m001300,m000400,m300010,m210010,m120010,m030010,m201010,m111010,m021010,m102010,m012010,m003010,m200110,m110110,m020110,m101110,m011110,m002110,m100210,m010210,m001210,m000310,m200020,m110020,m020020,m101020,m011020,m002020,m100120,m010120,m001120,m000220,m100030,m010030,m001030,m000130,m000040,m300001,m210001,m120001,m030001,m201001,m111001,m021001,m102001,m012001,m003001,m200101,m110101,m020101,m101101,m011101,m002101,m100201,m010201,m001201,m000301,m200011,m110011,m020011,m101011,m011011,m002011,m100111,m010111,m001111,m000211,m100021,m010021,m001021,m000121,m000031,m200002,m110002,m020002,m101002,m011002,m002002,m100102,m010102,m001102,m000202,m100012,m010012,m001012,m000112,m000022,m100003,m010003,m001003,m000103,m000013,m000004)
%CONV_MOM_6_4
%    CENT = CONV_MOM_6_4(M100000,M010000,M001000,M000100,M000010,M000001,M200000,M110000,M020000,M101000,M011000,M002000,M100100,M010100,M001100,M000200,M100010,M010010,M001010,M000110,M000020,M100001,M010001,M001001,M000101,M000011,M000002,M300000,M210000,M120000,M030000,M201000,M111000,M021000,M102000,M012000,M003000,M200100,M110100,M020100,M101100,M011100,M002100,M100200,M010200,M001200,M000300,M200010,M110010,M020010,M101010,M011010,M002010,M100110,M010110,M001110,M000210,M100020,M010020,M001020,M000120,M000030,M200001,M110001,M020001,M101001,M011001,M002001,M100101,M010101,M001101,M000201,M100011,M010011,M001011,M000111,M000021,M100002,M010002,M001002,M000102,M000012,M000003,M400000,M310000,M220000,M130000,M040000,M301000,M211000,M121000,M031000,M202000,M112000,M022000,M103000,M013000,M004000,M300100,M210100,M120100,M030100,M201100,M111100,M021100,M102100,M012100,M003100,M200200,M110200,M020200,M101200,M011200,M002200,M100300,M010300,M001300,M000400,M300010,M210010,M120010,M030010,M201010,M111010,M021010,M102010,M012010,M003010,M200110,M110110,M020110,M101110,M011110,M002110,M100210,M010210,M001210,M000310,M200020,M110020,M020020,M101020,M011020,M002020,M100120,M010120,M001120,M000220,M100030,M010030,M001030,M000130,M000040,M300001,M210001,M120001,M030001,M201001,M111001,M021001,M102001,M012001,M003001,M200101,M110101,M020101,M101101,M011101,M002101,M100201,M010201,M001201,M000301,M200011,M110011,M020011,M101011,M011011,M002011,M100111,M010111,M001111,M000211,M100021,M010021,M001021,M000121,M000031,M200002,M110002,M020002,M101002,M011002,M002002,M100102,M010102,M001102,M000202,M100012,M010012,M001012,M000112,M000022,M100003,M010003,M001003,M000103,M000013,M000004)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    09-Apr-2017 13:06:09

t2 = m100000.^2;
t3 = m010000.^2;
t4 = m001000.^2;
t5 = m000100.^2;
t6 = m000010.^2;
t7 = m000001.^2;
cent = [0.0;0.0;0.0;0.0;0.0;0.0;m200000-t2;m110000-m010000.*m100000;m020000-t3;m101000-m001000.*m100000;m011000-m001000.*m010000;m002000-t4;m100100-m000100.*m100000;m010100-m000100.*m010000;m001100-m000100.*m001000;m000200-t5;m100010-m000010.*m100000;m010010-m000010.*m010000;m001010-m000010.*m001000;m000110-m000010.*m000100;m000020-t6;m100001-m000001.*m100000;m010001-m000001.*m010000;m001001-m000001.*m001000;m000101-m000001.*m000100;m000011-m000001.*m000010;m000002-t7;m300000-m100000.*m200000.*3.0+m100000.*t2.*2.0;m210000-m010000.*m200000-m100000.*m110000.*2.0+m010000.*t2.*2.0;m120000-m010000.*m110000.*2.0-m020000.*m100000+m100000.*t3.*2.0;m030000-m010000.*m020000.*3.0+m010000.*t3.*2.0;m201000-m001000.*m200000-m100000.*m101000.*2.0+m001000.*t2.*2.0;m111000-m001000.*m110000-m010000.*m101000-m011000.*m100000+m001000.*m010000.*m100000.*2.0;m021000-m001000.*m020000-m010000.*m011000.*2.0+m001000.*t3.*2.0;m102000-m001000.*m101000.*2.0-m002000.*m100000+m100000.*t4.*2.0;m012000-m001000.*m011000.*2.0-m002000.*m010000+m010000.*t4.*2.0;m003000-m001000.*m002000.*3.0+m001000.*t4.*2.0;m200100-m000100.*m200000-m100000.*m100100.*2.0+m000100.*t2.*2.0;m110100-m000100.*m110000-m010000.*m100100-m010100.*m100000+m000100.*m010000.*m100000.*2.0;m020100-m000100.*m020000-m010000.*m010100.*2.0+m000100.*t3.*2.0;m101100-m000100.*m101000-m001000.*m100100-m001100.*m100000+m000100.*m001000.*m100000.*2.0;m011100-m000100.*m011000-m001000.*m010100-m001100.*m010000+m000100.*m001000.*m010000.*2.0;m002100-m000100.*m002000-m001000.*m001100.*2.0+m000100.*t4.*2.0;m100200-m000100.*m100100.*2.0-m000200.*m100000+m100000.*t5.*2.0;m010200-m000100.*m010100.*2.0-m000200.*m010000+m010000.*t5.*2.0;m001200-m000100.*m001100.*2.0-m000200.*m001000+m001000.*t5.*2.0;m000300-m000100.*m000200.*3.0+m000100.*t5.*2.0;m200010-m000010.*m200000-m100000.*m100010.*2.0+m000010.*t2.*2.0;m110010-m000010.*m110000-m010000.*m100010-m010010.*m100000+m000010.*m010000.*m100000.*2.0;m020010-m000010.*m020000-m010000.*m010010.*2.0+m000010.*t3.*2.0;m101010-m000010.*m101000-m001000.*m100010-m001010.*m100000+m000010.*m001000.*m100000.*2.0;m011010-m000010.*m011000-m001000.*m010010-m001010.*m010000+m000010.*m001000.*m010000.*2.0;m002010-m000010.*m002000-m001000.*m001010.*2.0+m000010.*t4.*2.0;m100110-m000010.*m100100-m000100.*m100010-m000110.*m100000+m000010.*m000100.*m100000.*2.0;m010110-m000010.*m010100-m000100.*m010010-m000110.*m010000+m000010.*m000100.*m010000.*2.0;m001110-m000010.*m001100-m000100.*m001010-m000110.*m001000+m000010.*m000100.*m001000.*2.0;m000210-m000010.*m000200-m000100.*m000110.*2.0+m000010.*t5.*2.0;m100020-m000010.*m100010.*2.0-m000020.*m100000+m100000.*t6.*2.0;m010020-m000010.*m010010.*2.0-m000020.*m010000+m010000.*t6.*2.0;m001020-m000010.*m001010.*2.0-m000020.*m001000+m001000.*t6.*2.0;m000120-m000010.*m000110.*2.0-m000020.*m000100+m000100.*t6.*2.0;m000030-m000010.*m000020.*3.0+m000010.*t6.*2.0;m200001-m000001.*m200000-m100000.*m100001.*2.0+m000001.*t2.*2.0;m110001-m000001.*m110000-m010000.*m100001-m010001.*m100000+m000001.*m010000.*m100000.*2.0;m020001-m000001.*m020000-m010000.*m010001.*2.0+m000001.*t3.*2.0;m101001-m000001.*m101000-m001000.*m100001-m001001.*m100000+m000001.*m001000.*m100000.*2.0;m011001-m000001.*m011000-m001000.*m010001-m001001.*m010000+m000001.*m001000.*m010000.*2.0;m002001-m000001.*m002000-m001000.*m001001.*2.0+m000001.*t4.*2.0;m100101-m000001.*m100100-m000100.*m100001-m000101.*m100000+m000001.*m000100.*m100000.*2.0;m010101-m000001.*m010100-m000100.*m010001-m000101.*m010000+m000001.*m000100.*m010000.*2.0;m001101-m000001.*m001100-m000100.*m001001-m000101.*m001000+m000001.*m000100.*m001000.*2.0;m000201-m000001.*m000200-m000100.*m000101.*2.0+m000001.*t5.*2.0;m100011-m000001.*m100010-m000010.*m100001-m000011.*m100000+m000001.*m000010.*m100000.*2.0;m010011-m000001.*m010010-m000010.*m010001-m000011.*m010000+m000001.*m000010.*m010000.*2.0;m001011-m000001.*m001010-m000010.*m001001-m000011.*m001000+m000001.*m000010.*m001000.*2.0;m000111-m000001.*m000110-m000010.*m000101-m000011.*m000100+m000001.*m000010.*m000100.*2.0;m000021-m000001.*m000020-m000010.*m000011.*2.0+m000001.*t6.*2.0;m100002-m000001.*m100001.*2.0-m000002.*m100000+m100000.*t7.*2.0;m010002-m000001.*m010001.*2.0-m000002.*m010000+m010000.*t7.*2.0;m001002-m000001.*m001001.*2.0-m000002.*m001000+m001000.*t7.*2.0;m000102-m000001.*m000101.*2.0-m000002.*m000100+m000100.*t7.*2.0;m000012-m000001.*m000011.*2.0-m000002.*m000010+m000010.*t7.*2.0;m000003-m000001.*m000002.*3.0+m000001.*t7.*2.0;m400000-m100000.*m300000.*4.0+m200000.*t2.*6.0-t2.^2.*3.0;m310000-m010000.*m300000-m100000.*m210000.*3.0+m110000.*t2.*3.0+m010000.*m100000.*m200000.*3.0-m010000.*m100000.*t2.*3.0;m220000-m010000.*m210000.*2.0-m100000.*m120000.*2.0+m020000.*t2+m200000.*t3-t2.*t3.*3.0+m010000.*m100000.*m110000.*4.0;m130000-m010000.*m120000.*3.0-m030000.*m100000+m110000.*t3.*3.0+m010000.*m020000.*m100000.*3.0-m010000.*m100000.*t3.*3.0;m040000-m010000.*m030000.*4.0+m020000.*t3.*6.0-t3.^2.*3.0;m301000-m001000.*m300000-m100000.*m201000.*3.0+m101000.*t2.*3.0+m001000.*m100000.*m200000.*3.0-m001000.*m100000.*t2.*3.0;m211000-m001000.*m210000-m010000.*m201000-m100000.*m111000.*2.0+m011000.*t2+m001000.*m010000.*m200000+m001000.*m100000.*m110000.*2.0+m010000.*m100000.*m101000.*2.0-m001000.*m010000.*t2.*3.0;m121000-m001000.*m120000-m010000.*m111000.*2.0-m021000.*m100000+m101000.*t3+m001000.*m010000.*m110000.*2.0+m001000.*m020000.*m100000+m010000.*m011000.*m100000.*2.0-m001000.*m100000.*t3.*3.0;m031000-m001000.*m030000-m010000.*m021000.*3.0+m011000.*t3.*3.0+m001000.*m010000.*m020000.*3.0-m001000.*m010000.*t3.*3.0;m202000-m001000.*m201000.*2.0-m100000.*m102000.*2.0+m002000.*t2+m200000.*t4-t2.*t4.*3.0+m001000.*m100000.*m101000.*4.0;m112000-m001000.*m111000.*2.0-m010000.*m102000-m012000.*m100000+m110000.*t4+m001000.*m010000.*m101000.*2.0+m001000.*m011000.*m100000.*2.0+m002000.*m010000.*m100000-m010000.*m100000.*t4.*3.0;m022000-m001000.*m021000.*2.0-m010000.*m012000.*2.0+m002000.*t3+m020000.*t4-t3.*t4.*3.0+m001000.*m010000.*m011000.*4.0;m103000-m001000.*m102000.*3.0-m003000.*m100000+m101000.*t4.*3.0+m001000.*m002000.*m100000.*3.0-m001000.*m100000.*t4.*3.0;m013000-m001000.*m012000.*3.0-m003000.*m010000+m011000.*t4.*3.0+m001000.*m002000.*m010000.*3.0-m001000.*m010000.*t4.*3.0;m004000-m001000.*m003000.*4.0+m002000.*t4.*6.0-t4.^2.*3.0;m300100-m000100.*m300000-m100000.*m200100.*3.0+m100100.*t2.*3.0+m000100.*m100000.*m200000.*3.0-m000100.*m100000.*t2.*3.0;m210100-m000100.*m210000-m010000.*m200100-m100000.*m110100.*2.0+m010100.*t2+m000100.*m010000.*m200000+m000100.*m100000.*m110000.*2.0+m010000.*m100000.*m100100.*2.0-m000100.*m010000.*t2.*3.0;m120100-m000100.*m120000-m010000.*m110100.*2.0-m020100.*m100000+m100100.*t3+m000100.*m010000.*m110000.*2.0+m000100.*m020000.*m100000+m010000.*m010100.*m100000.*2.0-m000100.*m100000.*t3.*3.0;m030100-m000100.*m030000-m010000.*m020100.*3.0+m010100.*t3.*3.0+m000100.*m010000.*m020000.*3.0-m000100.*m010000.*t3.*3.0;m201100-m000100.*m201000-m001000.*m200100-m100000.*m101100.*2.0+m001100.*t2+m000100.*m001000.*m200000+m000100.*m100000.*m101000.*2.0+m001000.*m100000.*m100100.*2.0-m000100.*m001000.*t2.*3.0;m111100-m000100.*m111000-m001000.*m110100-m010000.*m101100-m011100.*m100000+m000100.*m001000.*m110000+m000100.*m010000.*m101000+m000100.*m011000.*m100000+m001000.*m010000.*m100100+m001000.*m010100.*m100000+m001100.*m010000.*m100000-m000100.*m001000.*m010000.*m100000.*3.0;m021100-m000100.*m021000-m001000.*m020100-m010000.*m011100.*2.0+m001100.*t3+m000100.*m001000.*m020000+m000100.*m010000.*m011000.*2.0+m001000.*m010000.*m010100.*2.0-m000100.*m001000.*t3.*3.0;m102100-m000100.*m102000-m001000.*m101100.*2.0-m002100.*m100000+m100100.*t4+m000100.*m001000.*m101000.*2.0+m000100.*m002000.*m100000+m001000.*m001100.*m100000.*2.0-m000100.*m100000.*t4.*3.0;m012100-m000100.*m012000-m001000.*m011100.*2.0-m002100.*m010000+m010100.*t4+m000100.*m001000.*m011000.*2.0+m000100.*m002000.*m010000+m001000.*m001100.*m010000.*2.0-m000100.*m010000.*t4.*3.0;m003100-m000100.*m003000-m001000.*m002100.*3.0+m001100.*t4.*3.0+m000100.*m001000.*m002000.*3.0-m000100.*m001000.*t4.*3.0;m200200-m000100.*m200100.*2.0-m100000.*m100200.*2.0+m000200.*t2+m200000.*t5-t2.*t5.*3.0+m000100.*m100000.*m100100.*4.0;m110200-m000100.*m110100.*2.0-m010000.*m100200-m010200.*m100000+m110000.*t5+m000100.*m010000.*m100100.*2.0+m000100.*m010100.*m100000.*2.0+m000200.*m010000.*m100000-m010000.*m100000.*t5.*3.0;m020200-m000100.*m020100.*2.0-m010000.*m010200.*2.0+m000200.*t3+m020000.*t5-t3.*t5.*3.0+m000100.*m010000.*m010100.*4.0;m101200-m000100.*m101100.*2.0-m001000.*m100200-m001200.*m100000+m101000.*t5+m000100.*m001000.*m100100.*2.0+m000100.*m001100.*m100000.*2.0+m000200.*m001000.*m100000-m001000.*m100000.*t5.*3.0;m011200-m000100.*m011100.*2.0-m001000.*m010200-m001200.*m010000+m011000.*t5+m000100.*m001000.*m010100.*2.0+m000100.*m001100.*m010000.*2.0+m000200.*m001000.*m010000-m001000.*m010000.*t5.*3.0;m002200-m000100.*m002100.*2.0-m001000.*m001200.*2.0+m000200.*t4+m002000.*t5-t4.*t5.*3.0+m000100.*m001000.*m001100.*4.0;m100300-m000100.*m100200.*3.0-m000300.*m100000+m100100.*t5.*3.0+m000100.*m000200.*m100000.*3.0-m000100.*m100000.*t5.*3.0;m010300-m000100.*m010200.*3.0-m000300.*m010000+m010100.*t5.*3.0+m000100.*m000200.*m010000.*3.0-m000100.*m010000.*t5.*3.0;m001300-m000100.*m001200.*3.0-m000300.*m001000+m001100.*t5.*3.0+m000100.*m000200.*m001000.*3.0-m000100.*m001000.*t5.*3.0;m000400-m000100.*m000300.*4.0+m000200.*t5.*6.0-t5.^2.*3.0;m300010-m000010.*m300000-m100000.*m200010.*3.0+m100010.*t2.*3.0+m000010.*m100000.*m200000.*3.0-m000010.*m100000.*t2.*3.0;m210010-m000010.*m210000-m010000.*m200010-m100000.*m110010.*2.0+m010010.*t2+m000010.*m010000.*m200000+m000010.*m100000.*m110000.*2.0+m010000.*m100000.*m100010.*2.0-m000010.*m010000.*t2.*3.0;m120010-m000010.*m120000-m010000.*m110010.*2.0-m020010.*m100000+m100010.*t3+m000010.*m010000.*m110000.*2.0+m000010.*m020000.*m100000+m010000.*m010010.*m100000.*2.0-m000010.*m100000.*t3.*3.0;m030010-m000010.*m030000-m010000.*m020010.*3.0+m010010.*t3.*3.0+m000010.*m010000.*m020000.*3.0-m000010.*m010000.*t3.*3.0;m201010-m000010.*m201000-m001000.*m200010-m100000.*m101010.*2.0+m001010.*t2+m000010.*m001000.*m200000+m000010.*m100000.*m101000.*2.0+m001000.*m100000.*m100010.*2.0-m000010.*m001000.*t2.*3.0;m111010-m000010.*m111000-m001000.*m110010-m010000.*m101010-m011010.*m100000+m000010.*m001000.*m110000+m000010.*m010000.*m101000+m000010.*m011000.*m100000+m001000.*m010000.*m100010+m001000.*m010010.*m100000+m001010.*m010000.*m100000-m000010.*m001000.*m010000.*m100000.*3.0;m021010-m000010.*m021000-m001000.*m020010-m010000.*m011010.*2.0+m001010.*t3+m000010.*m001000.*m020000+m000010.*m010000.*m011000.*2.0+m001000.*m010000.*m010010.*2.0-m000010.*m001000.*t3.*3.0;m102010-m000010.*m102000-m001000.*m101010.*2.0-m002010.*m100000+m100010.*t4+m000010.*m001000.*m101000.*2.0+m000010.*m002000.*m100000+m001000.*m001010.*m100000.*2.0-m000010.*m100000.*t4.*3.0;m012010-m000010.*m012000-m001000.*m011010.*2.0-m002010.*m010000+m010010.*t4+m000010.*m001000.*m011000.*2.0+m000010.*m002000.*m010000+m001000.*m001010.*m010000.*2.0-m000010.*m010000.*t4.*3.0;m003010-m000010.*m003000-m001000.*m002010.*3.0+m001010.*t4.*3.0+m000010.*m001000.*m002000.*3.0-m000010.*m001000.*t4.*3.0;m200110-m000010.*m200100-m000100.*m200010-m100000.*m100110.*2.0+m000110.*t2+m000010.*m000100.*m200000+m000010.*m100000.*m100100.*2.0+m000100.*m100000.*m100010.*2.0-m000010.*m000100.*t2.*3.0;m110110-m000010.*m110100-m000100.*m110010-m010000.*m100110-m010110.*m100000+m000010.*m000100.*m110000+m000010.*m010000.*m100100+m000010.*m010100.*m100000+m000100.*m010000.*m100010+m000100.*m010010.*m100000+m000110.*m010000.*m100000-m000010.*m000100.*m010000.*m100000.*3.0;m020110-m000010.*m020100-m000100.*m020010-m010000.*m010110.*2.0+m000110.*t3+m000010.*m000100.*m020000+m000010.*m010000.*m010100.*2.0+m000100.*m010000.*m010010.*2.0-m000010.*m000100.*t3.*3.0;m101110-m000010.*m101100-m000100.*m101010-m001000.*m100110-m001110.*m100000+m000010.*m000100.*m101000+m000010.*m001000.*m100100+m000010.*m001100.*m100000+m000100.*m001000.*m100010+m000100.*m001010.*m100000+m000110.*m001000.*m100000-m000010.*m000100.*m001000.*m100000.*3.0;m011110-m000010.*m011100-m000100.*m011010-m001000.*m010110-m001110.*m010000+m000010.*m000100.*m011000+m000010.*m001000.*m010100+m000010.*m001100.*m010000+m000100.*m001000.*m010010+m000100.*m001010.*m010000+m000110.*m001000.*m010000-m000010.*m000100.*m001000.*m010000.*3.0;m002110-m000010.*m002100-m000100.*m002010-m001000.*m001110.*2.0+m000110.*t4+m000010.*m000100.*m002000+m000010.*m001000.*m001100.*2.0+m000100.*m001000.*m001010.*2.0-m000010.*m000100.*t4.*3.0;m100210-m000010.*m100200-m000100.*m100110.*2.0-m000210.*m100000+m100010.*t5+m000010.*m000100.*m100100.*2.0+m000010.*m000200.*m100000+m000100.*m000110.*m100000.*2.0-m000010.*m100000.*t5.*3.0;m010210-m000010.*m010200-m000100.*m010110.*2.0-m000210.*m010000+m010010.*t5+m000010.*m000100.*m010100.*2.0+m000010.*m000200.*m010000+m000100.*m000110.*m010000.*2.0-m000010.*m010000.*t5.*3.0;m001210-m000010.*m001200-m000100.*m001110.*2.0-m000210.*m001000+m001010.*t5+m000010.*m000100.*m001100.*2.0+m000010.*m000200.*m001000+m000100.*m000110.*m001000.*2.0-m000010.*m001000.*t5.*3.0;m000310-m000010.*m000300-m000100.*m000210.*3.0+m000110.*t5.*3.0+m000010.*m000100.*m000200.*3.0-m000010.*m000100.*t5.*3.0;m200020-m000010.*m200010.*2.0-m100000.*m100020.*2.0+m000020.*t2+m200000.*t6-t2.*t6.*3.0+m000010.*m100000.*m100010.*4.0;m110020-m000010.*m110010.*2.0-m010000.*m100020-m010020.*m100000+m110000.*t6+m000010.*m010000.*m100010.*2.0+m000010.*m010010.*m100000.*2.0+m000020.*m010000.*m100000-m010000.*m100000.*t6.*3.0;m020020-m000010.*m020010.*2.0-m010000.*m010020.*2.0+m000020.*t3+m020000.*t6-t3.*t6.*3.0+m000010.*m010000.*m010010.*4.0;m101020-m000010.*m101010.*2.0-m001000.*m100020-m001020.*m100000+m101000.*t6+m000010.*m001000.*m100010.*2.0+m000010.*m001010.*m100000.*2.0+m000020.*m001000.*m100000-m001000.*m100000.*t6.*3.0;m011020-m000010.*m011010.*2.0-m001000.*m010020-m001020.*m010000+m011000.*t6+m000010.*m001000.*m010010.*2.0+m000010.*m001010.*m010000.*2.0+m000020.*m001000.*m010000-m001000.*m010000.*t6.*3.0;m002020-m000010.*m002010.*2.0-m001000.*m001020.*2.0+m000020.*t4+m002000.*t6-t4.*t6.*3.0+m000010.*m001000.*m001010.*4.0;m100120-m000010.*m100110.*2.0-m000100.*m100020-m000120.*m100000+m100100.*t6+m000010.*m000100.*m100010.*2.0+m000010.*m000110.*m100000.*2.0+m000020.*m000100.*m100000-m000100.*m100000.*t6.*3.0;m010120-m000010.*m010110.*2.0-m000100.*m010020-m000120.*m010000+m010100.*t6+m000010.*m000100.*m010010.*2.0+m000010.*m000110.*m010000.*2.0+m000020.*m000100.*m010000-m000100.*m010000.*t6.*3.0;m001120-m000010.*m001110.*2.0-m000100.*m001020-m000120.*m001000+m001100.*t6+m000010.*m000100.*m001010.*2.0+m000010.*m000110.*m001000.*2.0+m000020.*m000100.*m001000-m000100.*m001000.*t6.*3.0;m000220-m000010.*m000210.*2.0-m000100.*m000120.*2.0+m000020.*t5+m000200.*t6-t5.*t6.*3.0+m000010.*m000100.*m000110.*4.0;m100030-m000010.*m100020.*3.0-m000030.*m100000+m100010.*t6.*3.0+m000010.*m000020.*m100000.*3.0-m000010.*m100000.*t6.*3.0;m010030-m000010.*m010020.*3.0-m000030.*m010000+m010010.*t6.*3.0+m000010.*m000020.*m010000.*3.0-m000010.*m010000.*t6.*3.0;m001030-m000010.*m001020.*3.0-m000030.*m001000+m001010.*t6.*3.0+m000010.*m000020.*m001000.*3.0-m000010.*m001000.*t6.*3.0;m000130-m000010.*m000120.*3.0-m000030.*m000100+m000110.*t6.*3.0+m000010.*m000020.*m000100.*3.0-m000010.*m000100.*t6.*3.0;m000040-m000010.*m000030.*4.0+m000020.*t6.*6.0-t6.^2.*3.0;m300001-m000001.*m300000-m100000.*m200001.*3.0+m100001.*t2.*3.0+m000001.*m100000.*m200000.*3.0-m000001.*m100000.*t2.*3.0;m210001-m000001.*m210000-m010000.*m200001-m100000.*m110001.*2.0+m010001.*t2+m000001.*m010000.*m200000+m000001.*m100000.*m110000.*2.0+m010000.*m100000.*m100001.*2.0-m000001.*m010000.*t2.*3.0;m120001-m000001.*m120000-m010000.*m110001.*2.0-m020001.*m100000+m100001.*t3+m000001.*m010000.*m110000.*2.0+m000001.*m020000.*m100000+m010000.*m010001.*m100000.*2.0-m000001.*m100000.*t3.*3.0;m030001-m000001.*m030000-m010000.*m020001.*3.0+m010001.*t3.*3.0+m000001.*m010000.*m020000.*3.0-m000001.*m010000.*t3.*3.0;m201001-m000001.*m201000-m001000.*m200001-m100000.*m101001.*2.0+m001001.*t2+m000001.*m001000.*m200000+m000001.*m100000.*m101000.*2.0+m001000.*m100000.*m100001.*2.0-m000001.*m001000.*t2.*3.0;m111001-m000001.*m111000-m001000.*m110001-m010000.*m101001-m011001.*m100000+m000001.*m001000.*m110000+m000001.*m010000.*m101000+m000001.*m011000.*m100000+m001000.*m010000.*m100001+m001000.*m010001.*m100000+m001001.*m010000.*m100000-m000001.*m001000.*m010000.*m100000.*3.0;m021001-m000001.*m021000-m001000.*m020001-m010000.*m011001.*2.0+m001001.*t3+m000001.*m001000.*m020000+m000001.*m010000.*m011000.*2.0+m001000.*m010000.*m010001.*2.0-m000001.*m001000.*t3.*3.0;m102001-m000001.*m102000-m001000.*m101001.*2.0-m002001.*m100000+m100001.*t4+m000001.*m001000.*m101000.*2.0+m000001.*m002000.*m100000+m001000.*m001001.*m100000.*2.0-m000001.*m100000.*t4.*3.0;m012001-m000001.*m012000-m001000.*m011001.*2.0-m002001.*m010000+m010001.*t4+m000001.*m001000.*m011000.*2.0+m000001.*m002000.*m010000+m001000.*m001001.*m010000.*2.0-m000001.*m010000.*t4.*3.0;m003001-m000001.*m003000-m001000.*m002001.*3.0+m001001.*t4.*3.0+m000001.*m001000.*m002000.*3.0-m000001.*m001000.*t4.*3.0;m200101-m000001.*m200100-m000100.*m200001-m100000.*m100101.*2.0+m000101.*t2+m000001.*m000100.*m200000+m000001.*m100000.*m100100.*2.0+m000100.*m100000.*m100001.*2.0-m000001.*m000100.*t2.*3.0;m110101-m000001.*m110100-m000100.*m110001-m010000.*m100101-m010101.*m100000+m000001.*m000100.*m110000+m000001.*m010000.*m100100+m000001.*m010100.*m100000+m000100.*m010000.*m100001+m000100.*m010001.*m100000+m000101.*m010000.*m100000-m000001.*m000100.*m010000.*m100000.*3.0;m020101-m000001.*m020100-m000100.*m020001-m010000.*m010101.*2.0+m000101.*t3+m000001.*m000100.*m020000+m000001.*m010000.*m010100.*2.0+m000100.*m010000.*m010001.*2.0-m000001.*m000100.*t3.*3.0;m101101-m000001.*m101100-m000100.*m101001-m001000.*m100101-m001101.*m100000+m000001.*m000100.*m101000+m000001.*m001000.*m100100+m000001.*m001100.*m100000+m000100.*m001000.*m100001+m000100.*m001001.*m100000+m000101.*m001000.*m100000-m000001.*m000100.*m001000.*m100000.*3.0;m011101-m000001.*m011100-m000100.*m011001-m001000.*m010101-m001101.*m010000+m000001.*m000100.*m011000+m000001.*m001000.*m010100+m000001.*m001100.*m010000+m000100.*m001000.*m010001+m000100.*m001001.*m010000+m000101.*m001000.*m010000-m000001.*m000100.*m001000.*m010000.*3.0;m002101-m000001.*m002100-m000100.*m002001-m001000.*m001101.*2.0+m000101.*t4+m000001.*m000100.*m002000+m000001.*m001000.*m001100.*2.0+m000100.*m001000.*m001001.*2.0-m000001.*m000100.*t4.*3.0;m100201-m000001.*m100200-m000100.*m100101.*2.0-m000201.*m100000+m100001.*t5+m000001.*m000100.*m100100.*2.0+m000001.*m000200.*m100000+m000100.*m000101.*m100000.*2.0-m000001.*m100000.*t5.*3.0;m010201-m000001.*m010200-m000100.*m010101.*2.0-m000201.*m010000+m010001.*t5+m000001.*m000100.*m010100.*2.0+m000001.*m000200.*m010000+m000100.*m000101.*m010000.*2.0-m000001.*m010000.*t5.*3.0;m001201-m000001.*m001200-m000100.*m001101.*2.0-m000201.*m001000+m001001.*t5+m000001.*m000100.*m001100.*2.0+m000001.*m000200.*m001000+m000100.*m000101.*m001000.*2.0-m000001.*m001000.*t5.*3.0;m000301-m000001.*m000300-m000100.*m000201.*3.0+m000101.*t5.*3.0+m000001.*m000100.*m000200.*3.0-m000001.*m000100.*t5.*3.0;m200011-m000001.*m200010-m000010.*m200001-m100000.*m100011.*2.0+m000011.*t2+m000001.*m000010.*m200000+m000001.*m100000.*m100010.*2.0+m000010.*m100000.*m100001.*2.0-m000001.*m000010.*t2.*3.0;m110011-m000001.*m110010-m000010.*m110001-m010000.*m100011-m010011.*m100000+m000001.*m000010.*m110000+m000001.*m010000.*m100010+m000001.*m010010.*m100000+m000010.*m010000.*m100001+m000010.*m010001.*m100000+m000011.*m010000.*m100000-m000001.*m000010.*m010000.*m100000.*3.0;m020011-m000001.*m020010-m000010.*m020001-m010000.*m010011.*2.0+m000011.*t3+m000001.*m000010.*m020000+m000001.*m010000.*m010010.*2.0+m000010.*m010000.*m010001.*2.0-m000001.*m000010.*t3.*3.0;m101011-m000001.*m101010-m000010.*m101001-m001000.*m100011-m001011.*m100000+m000001.*m000010.*m101000+m000001.*m001000.*m100010+m000001.*m001010.*m100000+m000010.*m001000.*m100001+m000010.*m001001.*m100000+m000011.*m001000.*m100000-m000001.*m000010.*m001000.*m100000.*3.0;m011011-m000001.*m011010-m000010.*m011001-m001000.*m010011-m001011.*m010000+m000001.*m000010.*m011000+m000001.*m001000.*m010010+m000001.*m001010.*m010000+m000010.*m001000.*m010001+m000010.*m001001.*m010000+m000011.*m001000.*m010000-m000001.*m000010.*m001000.*m010000.*3.0;m002011-m000001.*m002010-m000010.*m002001-m001000.*m001011.*2.0+m000011.*t4+m000001.*m000010.*m002000+m000001.*m001000.*m001010.*2.0+m000010.*m001000.*m001001.*2.0-m000001.*m000010.*t4.*3.0;m100111-m000001.*m100110-m000010.*m100101-m000100.*m100011-m000111.*m100000+m000001.*m000010.*m100100+m000001.*m000100.*m100010+m000001.*m000110.*m100000+m000010.*m000100.*m100001+m000010.*m000101.*m100000+m000011.*m000100.*m100000-m000001.*m000010.*m000100.*m100000.*3.0;m010111-m000001.*m010110-m000010.*m010101-m000100.*m010011-m000111.*m010000+m000001.*m000010.*m010100+m000001.*m000100.*m010010+m000001.*m000110.*m010000+m000010.*m000100.*m010001+m000010.*m000101.*m010000+m000011.*m000100.*m010000-m000001.*m000010.*m000100.*m010000.*3.0;m001111-m000001.*m001110-m000010.*m001101-m000100.*m001011-m000111.*m001000+m000001.*m000010.*m001100+m000001.*m000100.*m001010+m000001.*m000110.*m001000+m000010.*m000100.*m001001+m000010.*m000101.*m001000+m000011.*m000100.*m001000-m000001.*m000010.*m000100.*m001000.*3.0;m000211-m000001.*m000210-m000010.*m000201-m000100.*m000111.*2.0+m000011.*t5+m000001.*m000010.*m000200+m000001.*m000100.*m000110.*2.0+m000010.*m000100.*m000101.*2.0-m000001.*m000010.*t5.*3.0;m100021-m000001.*m100020-m000010.*m100011.*2.0-m000021.*m100000+m100001.*t6+m000001.*m000010.*m100010.*2.0+m000001.*m000020.*m100000+m000010.*m000011.*m100000.*2.0-m000001.*m100000.*t6.*3.0;m010021-m000001.*m010020-m000010.*m010011.*2.0-m000021.*m010000+m010001.*t6+m000001.*m000010.*m010010.*2.0+m000001.*m000020.*m010000+m000010.*m000011.*m010000.*2.0-m000001.*m010000.*t6.*3.0;m001021-m000001.*m001020-m000010.*m001011.*2.0-m000021.*m001000+m001001.*t6+m000001.*m000010.*m001010.*2.0+m000001.*m000020.*m001000+m000010.*m000011.*m001000.*2.0-m000001.*m001000.*t6.*3.0;m000121-m000001.*m000120-m000010.*m000111.*2.0-m000021.*m000100+m000101.*t6+m000001.*m000010.*m000110.*2.0+m000001.*m000020.*m000100+m000010.*m000011.*m000100.*2.0-m000001.*m000100.*t6.*3.0;m000031-m000001.*m000030-m000010.*m000021.*3.0+m000011.*t6.*3.0+m000001.*m000010.*m000020.*3.0-m000001.*m000010.*t6.*3.0;m200002-m000001.*m200001.*2.0-m100000.*m100002.*2.0+m000002.*t2+m200000.*t7-t2.*t7.*3.0+m000001.*m100000.*m100001.*4.0;m110002-m000001.*m110001.*2.0-m010000.*m100002-m010002.*m100000+m110000.*t7+m000001.*m010000.*m100001.*2.0+m000001.*m010001.*m100000.*2.0+m000002.*m010000.*m100000-m010000.*m100000.*t7.*3.0;m020002-m000001.*m020001.*2.0-m010000.*m010002.*2.0+m000002.*t3+m020000.*t7-t3.*t7.*3.0+m000001.*m010000.*m010001.*4.0;m101002-m000001.*m101001.*2.0-m001000.*m100002-m001002.*m100000+m101000.*t7+m000001.*m001000.*m100001.*2.0+m000001.*m001001.*m100000.*2.0+m000002.*m001000.*m100000-m001000.*m100000.*t7.*3.0;m011002-m000001.*m011001.*2.0-m001000.*m010002-m001002.*m010000+m011000.*t7+m000001.*m001000.*m010001.*2.0+m000001.*m001001.*m010000.*2.0+m000002.*m001000.*m010000-m001000.*m010000.*t7.*3.0;m002002-m000001.*m002001.*2.0-m001000.*m001002.*2.0+m000002.*t4+m002000.*t7-t4.*t7.*3.0+m000001.*m001000.*m001001.*4.0;m100102-m000001.*m100101.*2.0-m000100.*m100002-m000102.*m100000+m100100.*t7+m000001.*m000100.*m100001.*2.0+m000001.*m000101.*m100000.*2.0+m000002.*m000100.*m100000-m000100.*m100000.*t7.*3.0;m010102-m000001.*m010101.*2.0-m000100.*m010002-m000102.*m010000+m010100.*t7+m000001.*m000100.*m010001.*2.0+m000001.*m000101.*m010000.*2.0+m000002.*m000100.*m010000-m000100.*m010000.*t7.*3.0;m001102-m000001.*m001101.*2.0-m000100.*m001002-m000102.*m001000+m001100.*t7+m000001.*m000100.*m001001.*2.0+m000001.*m000101.*m001000.*2.0+m000002.*m000100.*m001000-m000100.*m001000.*t7.*3.0;m000202-m000001.*m000201.*2.0-m000100.*m000102.*2.0+m000002.*t5+m000200.*t7-t5.*t7.*3.0+m000001.*m000100.*m000101.*4.0;m100012-m000001.*m100011.*2.0-m000010.*m100002-m000012.*m100000+m100010.*t7+m000001.*m000010.*m100001.*2.0+m000001.*m000011.*m100000.*2.0+m000002.*m000010.*m100000-m000010.*m100000.*t7.*3.0;m010012-m000001.*m010011.*2.0-m000010.*m010002-m000012.*m010000+m010010.*t7+m000001.*m000010.*m010001.*2.0+m000001.*m000011.*m010000.*2.0+m000002.*m000010.*m010000-m000010.*m010000.*t7.*3.0;m001012-m000001.*m001011.*2.0-m000010.*m001002-m000012.*m001000+m001010.*t7+m000001.*m000010.*m001001.*2.0+m000001.*m000011.*m001000.*2.0+m000002.*m000010.*m001000-m000010.*m001000.*t7.*3.0;m000112-m000001.*m000111.*2.0-m000010.*m000102-m000012.*m000100+m000110.*t7+m000001.*m000010.*m000101.*2.0+m000001.*m000011.*m000100.*2.0+m000002.*m000010.*m000100-m000010.*m000100.*t7.*3.0;m000022-m000001.*m000021.*2.0-m000010.*m000012.*2.0+m000002.*t6+m000020.*t7-t6.*t7.*3.0+m000001.*m000010.*m000011.*4.0;m100003-m000001.*m100002.*3.0-m000003.*m100000+m100001.*t7.*3.0+m000001.*m000002.*m100000.*3.0-m000001.*m100000.*t7.*3.0;m010003-m000001.*m010002.*3.0-m000003.*m010000+m010001.*t7.*3.0+m000001.*m000002.*m010000.*3.0-m000001.*m010000.*t7.*3.0;m001003-m000001.*m001002.*3.0-m000003.*m001000+m001001.*t7.*3.0+m000001.*m000002.*m001000.*3.0-m000001.*m001000.*t7.*3.0;m000103-m000001.*m000102.*3.0-m000003.*m000100+m000101.*t7.*3.0+m000001.*m000002.*m000100.*3.0-m000001.*m000100.*t7.*3.0;m000013-m000001.*m000012.*3.0-m000003.*m000010+m000011.*t7.*3.0+m000001.*m000002.*m000010.*3.0-m000001.*m000010.*t7.*3.0;m000004-m000001.*m000003.*4.0+m000002.*t7.*6.0-t7.^2.*3.0];
