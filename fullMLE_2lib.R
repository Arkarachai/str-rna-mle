#!/usr/bin/env Rscript

######################################################
## estimate ML of RDD,RT error rates/expansion bias ##
######################################################


#######################
## argument handle  ###
args<-commandArgs(TRUE)

if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("
      The R Script to estimate rates and expansion/contraction of RNA-DNA difference and RT reverse transcriptase
 
      Arguments:
      --filename=string - input file name [MANDATORY]
      --bin_size=integer [2] - number of cDNA and RNA molecule that will be used in Likelihood estimation [2-50]
        Use larger integer for loci with high number of mapped RNA-seq read (high expression level).
        Use bin_size = 2 for loci with 3-5 mapped RNA-seq reads
        Use bin_size = 3 for loci with 4-7 mapped RNA-seq reads
        Use bin_size = 5 for loci with 6-16 mapped RNA-seq reads
        Use bin_size = 40 for loci with 49-102 mapped RNA-seq reads
        The larger the bin size, the longer the running time
      --N_core=integer [2]   - number of core for parallel processing among loci [2-12]
      --lower_RDD_RT=float [1e-9] - upper limit of RDD and RT error rates
      --upper_RDD_RT=float [0.5] - lower limit of RDD and RT error rates
      --step_RDD_RT=float [0.00001] - step of change in MLE to infer RDD and RT error rates
      --step_RDD_RT_expansion=float [0.01] - step of change in MLE to infer RDD and RT error expansion probabilities
      --help                 - print this text

      Example:
      ./fullMLE_2lib.R --filename=inputfile.txt --bin_size=3 --N_core=4 > outputfile.txt \n\n")
 
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
#argsL <- as.list(argsDF$V2)
names(argsL) <- argsDF$V1
 
if(is.null(argsL$filename)) {
  print('input filename is missing')
  q()
} 
if(is.null(argsL$bin_size)) {
  argsL$bin_size=2
}
argsL$bin_size=as.integer(argsL$bin_size)
if(is.null(argsL$N_core)) {
  argsL$N_core=2
}
argsL$N_core=as.integer(argsL$N_core)
if(is.null(argsL$lower_RDD_RT)) {
  argsL$lower_RDD_RT=1e-9
}
argsL$lower_RDD_RT=as.numeric(argsL$lower_RDD_RT)
if(is.null(argsL$upper_RDD_RT)) {
  argsL$upper_RDD_RT=0.5
}
argsL$upper_RDD_RT=as.numeric(argsL$upper_RDD_RT)
if(is.null(argsL$step_RDD_RT)) {
  argsL$step_RDD_RT=0.00001
}
argsL$step_RDD_RT=as.numeric(argsL$step_RDD_RT)
if(is.null(argsL$step_RDD_RT_expansion)) {
  argsL$step_RDD_RT_expansion=0.01
}
argsL$step_RDD_RT_expansion=as.numeric(argsL$step_RDD_RT_expansion)

####################
## import library ##
####################

library(combinat)
library(snow)
 
#######################
## Special function  ##
#######################

SeqErrorSwtich <- function(motif_temp,oriAllele,deriveAllele){
nameindex=paste(motif_temp,oriAllele,deriveAllele,sep="_")
switch(nameindex,
        A_10_10=0.979752747841,
        A_10_9=0.0134874552739,
        A_10_11=0.00648126325713,
        A_10_8=0.000171405309279,
        A_10_12=8.57026546397e-05,
        A_10_7=2.14256636599e-05,
        A_11_11=0.956766146349,
        A_11_10=0.0292166276629,
        A_11_12=0.0124007816835,
        A_11_9=0.0013028058578,
        A_11_13=0.000217134309634,
        A_11_8=7.23781032112e-05,
        A_11_14=2.41260344037e-05,
        A_12_12=0.92692761605,
        A_12_11=0.0484854445319,
        A_12_13=0.0205055074744,
        A_12_10=0.003589693155,
        A_12_14=0.000393391030685,
        A_12_9=4.91738788356e-05,
        A_12_8=4.91738788356e-05,
        A_13_13=0.88634981319,
        A_13_12=0.0768963419932,
        A_13_14=0.0278043270484,
        A_13_11=0.00721174732818,
        A_13_15=0.00104266226431,
        A_13_10=0.00069510817621,
        A_14_14=0.85197291618,
        A_14_13=0.0954938127481,
        A_14_15=0.0352556619192,
        A_14_12=0.0144758346953,
        A_14_11=0.00140088722858,
        A_14_16=0.00116740602382,
        A_14_10=0.000233481204763,
        A_15_15=0.788640595903,
        A_15_14=0.130353817505,
        A_15_16=0.0558659217877,
        A_15_13=0.0186219739292,
        A_15_17=0.00372439478585,
        A_15_12=0.00279329608939,
        A_16_16=0.671586715867,
        A_16_15=0.221402214022,
        A_16_17=0.0516605166052,
        A_16_14=0.0442804428044,
        A_16_13=0.00369003690037,
        A_16_12=0.00369003690037,
        A_16_18=0.00369003690037,
        A_17_17=0.578947368421,
        A_17_16=0.263157894737,
        A_17_15=0.105263157895,
        A_17_18=0.0526315789474,
        A_18_18=0.666666666667,
        A_18_17=0.333333333333,
        A_5_5=0.999995620676,
        A_5_6=4.37932389808e-06,
        A_6_6=0.999905994121,
        A_6_5=6.9436160572e-05,
        A_6_7=2.45697183563e-05,
        A_7_7=0.999633738347,
        A_7_6=0.000285173874064,
        A_7_8=7.56211870522e-05,
        A_7_5=5.4665918351e-06,
        A_8_8=0.998270905829,
        A_8_7=0.0010455616314,
        A_8_9=0.000670023991182,
        A_8_6=8.10512892558e-06,
        A_8_10=5.40341928372e-06,
        A_9_9=0.993543188655,
        A_9_8=0.00379780833691,
        A_9_10=0.00261065749893,
        A_9_7=2.68586162441e-05,
        A_9_11=2.14868929953e-05,
        C_10_10=0.938775510204,
        C_10_9=0.0612244897959,
        C_5_5=0.999994833948,
        C_5_6=5.16605166052e-06,
        C_6_6=0.999917290007,
        C_6_5=7.23712438032e-05,
        C_6_7=1.03387491147e-05,
        C_7_7=0.999819616685,
        C_7_6=0.000180383314543,
        C_8_8=0.99798115747,
        C_8_9=0.00100942126514,
        C_8_7=0.00100942126514,
        C_9_9=0.986089644513,
        C_9_8=0.0123647604328,
        C_9_7=0.0015455950541,
        AC_10_10=0.99981145416,
        AC_10_8=0.000141409380156,
        AC_10_12=4.71364600519e-05,
        AC_11_11=0.999335901182,
        AC_11_9=0.000664098817904,
        AC_12_12=0.997359300215,
        AC_12_10=0.00247565604885,
        AC_12_14=0.00016504373659,
        AC_13_13=0.99197808648,
        AC_13_11=0.00782625709255,
        AC_13_15=0.000195656427314,
        AC_14_14=0.985659655832,
        AC_14_12=0.0140216698534,
        AC_14_10=0.00031867431485,
        AC_15_15=0.988888888889,
        AC_15_13=0.0107638888889,
        AC_15_17=0.000347222222222,
        AC_16_16=0.975478927203,
        AC_16_14=0.0229885057471,
        AC_16_12=0.00153256704981,
        AC_17_17=0.979607250755,
        AC_17_15=0.0203927492447,
        AC_18_18=0.965019011407,
        AC_18_16=0.0326996197719,
        AC_18_20=0.00152091254753,
        AC_18_14=0.000760456273764,
        AC_19_19=0.974175035868,
        AC_19_17=0.0243902439024,
        AC_19_21=0.00143472022956,
        AC_20_20=0.945747800587,
        AC_20_18=0.049853372434,
        AC_20_22=0.00293255131965,
        AC_20_16=0.00146627565982,
        AC_21_21=0.96015936255,
        AC_21_19=0.0371845949535,
        AC_21_17=0.00132802124834,
        AC_21_23=0.00132802124834,
        AC_22_22=0.939736346516,
        AC_22_20=0.0546139359699,
        AC_22_18=0.00564971751412,
        AC_23_23=0.942408376963,
        AC_23_21=0.0523560209424,
        AC_23_19=0.00349040139616,
        AC_23_25=0.00174520069808,
        AC_24_24=0.903755868545,
        AC_24_22=0.0892018779343,
        AC_24_26=0.00469483568075,
        AC_24_20=0.00234741784038,
        AC_25_25=0.94212962963,
        AC_25_23=0.0509259259259,
        AC_25_27=0.00462962962963,
        AC_25_21=0.00231481481481,
        AC_26_26=0.880136986301,
        AC_26_24=0.102739726027,
        AC_26_22=0.0102739726027,
        AC_26_28=0.00342465753425,
        AC_26_20=0.00342465753425,
        AC_27_27=0.911290322581,
        AC_27_25=0.0752688172043,
        AC_27_23=0.00806451612903,
        AC_27_29=0.00537634408602,
        AC_28_28=0.897777777778,
        AC_28_26=0.0755555555556,
        AC_28_30=0.0266666666667,
        AC_29_29=0.879365079365,
        AC_29_27=0.0920634920635,
        AC_29_31=0.0190476190476,
        AC_29_25=0.00952380952381,
        AC_30_30=0.87969924812,
        AC_30_28=0.0902255639098,
        AC_30_32=0.0225563909774,
        AC_30_18=0.00751879699248,
        AC_31_31=0.857142857143,
        AC_31_29=0.107142857143,
        AC_31_27=0.0238095238095,
        AC_31_33=0.0119047619048,
        AC_32_32=0.782945736434,
        AC_32_30=0.178294573643,
        AC_32_28=0.015503875969,
        AC_32_34=0.015503875969,
        AC_32_26=0.0077519379845,
        AC_33_33=0.848,
        AC_33_31=0.12,
        AC_33_35=0.024,
        AC_33_29=0.008,
        AC_34_34=0.825,
        AC_34_32=0.175,
        AC_35_35=0.807692307692,
        AC_35_33=0.153846153846,
        AC_35_31=0.0384615384615,
        AC_36_36=0.923076923077,
        AC_36_34=0.0769230769231,
        AC_37_37=0.666666666667,
        AC_37_35=0.2,
        AC_37_31=0.0666666666667,
        AC_37_39=0.0666666666667,
        AC_38_38=0.8,
        AC_38_36=0.2,
        AC_6_6=1,
        AC_7_7=1,
        AC_8_8=0.999985141711,
        AC_8_6=1.4858289068e-05,
        AC_9_9=0.999934150959,
        AC_9_7=4.9386780805e-05,
        AC_9_11=1.64622602683e-05,
        AG_10_10=0.999816311536,
        AG_10_8=0.000137766348273,
        AG_10_12=4.59221160911e-05,
        AG_11_11=0.999208120366,
        AG_11_9=0.000719890576632,
        AG_11_13=7.19890576632e-05,
        AG_12_12=0.996988485385,
        AG_12_10=0.00230292294066,
        AG_12_14=0.000708591674048,
        AG_13_13=0.996231434272,
        AG_13_11=0.00376856572822,
        AG_14_14=0.992158912703,
        AG_14_12=0.00784108729744,
        AG_15_15=0.99263803681,
        AG_15_13=0.00736196319018,
        AG_16_16=0.977293369664,
        AG_16_14=0.0217983651226,
        AG_16_12=0.000908265213442,
        AG_17_17=0.984234234234,
        AG_17_15=0.0135135135135,
        AG_17_19=0.00112612612613,
        AG_17_13=0.00112612612613,
        AG_18_18=0.962298025135,
        AG_18_16=0.0359066427289,
        AG_18_14=0.00179533213645,
        AG_19_19=0.957482993197,
        AG_19_17=0.0425170068027,
        AG_20_20=0.93488372093,
        AG_20_18=0.0651162790698,
        AG_21_21=0.962962962963,
        AG_21_19=0.037037037037,
        AG_22_22=0.943181818182,
        AG_22_20=0.0568181818182,
        AG_23_23=0.960784313725,
        AG_23_21=0.0326797385621,
        AG_23_25=0.00653594771242,
        AG_24_24=0.951923076923,
        AG_24_22=0.0384615384615,
        AG_24_18=0.00961538461538,
        AG_25_25=0.939393939394,
        AG_25_23=0.0454545454545,
        AG_25_27=0.0151515151515,
        AG_26_26=0.826086956522,
        AG_26_24=0.173913043478,
        AG_27_27=0.857142857143,
        AG_27_25=0.107142857143,
        AG_27_23=0.0357142857143,
        AG_28_28=0.875,
        AG_28_26=0.125,
        AG_29_29=0.666666666667,
        AG_29_27=0.277777777778,
        AG_29_31=0.0555555555556,
        AG_30_30=0.777777777778,
        AG_30_28=0.222222222222,
        AG_31_31=0.636363636364,
        AG_31_27=0.272727272727,
        AG_31_23=0.0909090909091,
        AG_32_32=0.8,
        AG_32_28=0.2,
        AG_6_6=1,
        AG_7_7=0.999998539801,
        AG_7_9=1.46019934641e-06,
        AG_8_8=1,
        AG_9_9=0.999907351034,
        AG_9_7=7.94133996876e-05,
        AG_9_11=1.32355666146e-05,
        AT_10_10=0.999395272128,
        AT_10_8=0.000384826827927,
        AT_10_12=0.00021990104453,
        AT_11_11=0.999220142602,
        AT_11_9=0.000557040998217,
        AT_11_13=0.000222816399287,
        AT_12_12=0.997958350347,
        AT_12_10=0.00163331972234,
        AT_12_14=0.000408329930584,
        AT_13_13=0.993578308955,
        AT_13_11=0.00606493043168,
        AT_13_15=0.000356760613628,
        AT_14_14=0.980475382003,
        AT_14_12=0.0169779286927,
        AT_14_16=0.00169779286927,
        AT_14_10=0.000848896434635,
        AT_15_15=0.976604278075,
        AT_15_13=0.0220588235294,
        AT_15_11=0.000668449197861,
        AT_15_17=0.000668449197861,
        AT_16_16=0.951298701299,
        AT_16_14=0.0454545454545,
        AT_16_18=0.0021645021645,
        AT_16_12=0.00108225108225,
        AT_17_17=0.9359375,
        AT_17_15=0.059375,
        AT_17_19=0.003125,
        AT_17_13=0.0015625,
        AT_18_18=0.908415841584,
        AT_18_16=0.0717821782178,
        AT_18_20=0.0173267326733,
        AT_18_14=0.00247524752475,
        AT_19_19=0.857692307692,
        AT_19_17=0.130769230769,
        AT_19_21=0.0115384615385,
        AT_20_20=0.850877192982,
        AT_20_18=0.122807017544,
        AT_20_16=0.0175438596491,
        AT_20_22=0.00877192982456,
        AT_21_21=0.759493670886,
        AT_21_19=0.227848101266,
        AT_21_17=0.0126582278481,
        AT_22_22=0.697368421053,
        AT_22_20=0.197368421053,
        AT_22_24=0.0657894736842,
        AT_22_18=0.0394736842105,
        AT_23_23=0.916666666667,
        AT_23_21=0.0833333333333,
        AT_24_24=0.636363636364,
        AT_24_20=0.181818181818,
        AT_24_22=0.181818181818,
        AT_6_6=0.99999940189,
        AT_6_8=5.98110091732e-07,
        AT_8_8=0.999964201698,
        AT_8_10=2.55702158126e-05,
        AT_8_6=1.0228086325e-05,
        AT_9_9=0.999942943,
        AT_9_7=5.70569999429e-05,
        AAC_12_12=0.999839897534,
        AAC_12_9=0.000160102465578,
        AAC_14_14=0.999763313609,
        AAC_14_11=0.000236686390533,
        AAC_15_15=0.997459779848,
        AAC_15_12=0.00254022015241,
        AAC_16_16=0.99709513435,
        AAC_16_13=0.00290486564996,
        AAC_17_17=0.998410174881,
        AAC_17_14=0.00158982511924,
        AAC_18_18=0.997885835095,
        AAC_18_15=0.00211416490486,
        AAC_19_19=0.997811816193,
        AAC_19_16=0.00218818380744,
        AAC_20_20=0.997968855789,
        AAC_20_17=0.00203114421124,
        AAC_21_21=0.996960486322,
        AAC_21_18=0.00303951367781,
        AAC_23_23=0.99814471243,
        AAC_23_26=0.00185528756957,
        AAC_26_26=0.995833333333,
        AAC_26_23=0.00416666666667,
        AAC_28_28=0.966666666667,
        AAC_28_25=0.0333333333333,
        AAC_31_31=0.974358974359,
        AAC_31_28=0.025641025641,
        AAC_32_32=0.969072164948,
        AAC_32_29=0.0309278350515,
        AAC_35_35=0.982142857143,
        AAC_35_32=0.0178571428571,
        AAG_12_12=0.998267898383,
        AAG_12_9=0.00173210161663,
        AAG_15_15=0.992753623188,
        AAG_15_12=0.00724637681159,
        AAT_12_12=0.999916680553,
        AAT_12_9=8.33194467589e-05,
        AAT_13_13=0.999721059972,
        AAT_13_10=0.000278940027894,
        AAT_14_14=0.999451854559,
        AAT_14_11=0.000548145441257,
        AAT_15_15=0.998812821528,
        AAT_15_12=0.0011871784725,
        AAT_16_16=0.999423298731,
        AAT_16_13=0.000576701268743,
        AAT_17_17=0.997739261492,
        AAT_17_14=0.00226073850791,
        AAT_18_18=0.997073170732,
        AAT_18_15=0.00292682926829,
        AAT_19_19=0.994059405941,
        AAT_19_16=0.00594059405941,
        AAT_20_20=0.996503496503,
        AAT_20_17=0.0034965034965,
        AAT_21_21=0.99730458221,
        AAT_21_18=0.00269541778976,
        AAT_23_23=0.982035928144,
        AAT_23_20=0.0179640718563,
        AAT_24_24=0.993055555556,
        AAT_24_21=0.00694444444444,
        AAT_25_25=0.991869918699,
        AAT_25_22=0.00813008130081,
        AAT_26_26=0.957446808511,
        AAT_26_23=0.0425531914894,
        AAT_27_27=0.969696969697,
        AAT_27_24=0.030303030303,
        AAT_29_29=0.984615384615,
        AAT_29_26=0.0153846153846,
        AAT_30_30=0.965517241379,
        AAT_30_24=0.0344827586207,
        AAT_32_32=0.9,
        AAT_32_29=0.1,
        AAT_9_9=0.999994419113,
        AAT_9_12=5.58088657964e-06,
        ACC_22_22=0.833333333333,
        ACC_22_19=0.166666666667,
        ACC_28_28=0.96,
        ACC_28_25=0.04,
        ACT_41_41=0.833333333333,
        ACT_41_35=0.166666666667,
        AGC_17_17=0.961538461538,
        AGC_17_14=0.0384615384615,
        ATC_34_34=0.916666666667,
        ATC_34_31=0.0833333333333,
        AAAC_18_18=0.998926462695,
        AAAC_18_14=0.00107353730542,
        AAAC_19_19=0.999429332319,
        AAAC_19_15=0.000380445120791,
        AAAC_19_23=0.000190222560396,
        AAAC_20_20=0.998905908096,
        AAAC_20_16=0.00109409190372,
        AAAC_22_22=0.996086105675,
        AAAC_22_18=0.00391389432485,
        AAAC_23_23=0.99733924612,
        AAAC_23_19=0.00221729490022,
        AAAC_23_15=0.000443458980044,
        AAAC_24_24=0.993957703927,
        AAAC_24_20=0.00604229607251,
        AAAC_25_25=0.995670995671,
        AAAC_25_21=0.004329004329,
        AAAC_27_27=0.996363636364,
        AAAC_27_23=0.00363636363636,
        AAAC_28_28=0.994897959184,
        AAAC_28_24=0.00510204081633,
        AAAC_31_31=0.993975903614,
        AAAC_31_27=0.00602409638554,
        AAAG_12_12=0.998136501281,
        AAAG_12_16=0.00100939513937,
        AAAG_12_20=0.000698812019567,
        AAAG_12_18=0.000155291559904,
        AAAG_14_14=0.997841920691,
        AAAG_14_13=0.00215807930941,
        AAAG_15_15=0.998188874515,
        AAAG_15_17=0.00155239327296,
        AAAG_15_13=0.00025873221216,
        AAAG_17_17=0.998668442077,
        AAAG_17_13=0.00133155792277,
        AAAG_18_18=0.984455958549,
        AAAG_18_20=0.0129533678756,
        AAAG_18_14=0.00259067357513,
        AAAG_19_19=0.999141630901,
        AAAG_19_15=0.000858369098712,
        AAAG_23_23=0.996894409938,
        AAAG_23_19=0.00310559006211,
        AAAG_26_26=0.896551724138,
        AAAG_26_13=0.103448275862,
        AAAT_16_16=0.999355670103,
        AAAT_16_12=0.000644329896907,
        AAAT_17_17=0.99900990099,
        AAAT_17_13=0.000990099009901,
        AAAT_19_19=0.999323638823,
        AAAT_19_15=0.000338180588434,
        AAAT_19_14=0.000338180588434,
        AAAT_21_21=0.997194950912,
        AAAT_21_17=0.00280504908836,
        AAAT_22_22=0.996015936255,
        AAAT_22_18=0.00398406374502,
        AAAT_23_23=0.996987951807,
        AAAT_23_19=0.00301204819277,
        AAAT_24_24=0.992207792208,
        AAAT_24_20=0.00779220779221,
        AAAT_26_26=0.994623655914,
        AAAT_26_22=0.00537634408602,
        AAAT_27_27=0.992932862191,
        AAAT_27_23=0.00706713780919,
        AAAT_28_28=0.988023952096,
        AAAT_28_24=0.0119760479042,
        AAAT_34_34=0.96,
        AAAT_34_30=0.04,
        AAAT_35_35=0.953488372093,
        AAAT_35_31=0.046511627907,
        AAAT_39_39=0.952380952381,
        AAAT_39_35=0.047619047619,
        AAGG_14_14=0.997838616715,
        AAGG_14_15=0.0021613832853,
        AAGG_16_16=0.993506493506,
        AAGG_16_14=0.00649350649351,
        AAGT_22_22=0.977272727273,
        AAGT_22_26=0.0227272727273,
        AAGT_23_23=0.969696969697,
        AAGT_23_19=0.030303030303,
        AATG_27_27=0.95652173913,
        AATG_27_23=0.0434782608696,
        ACAT_12_12=0.999178981938,
        ACAT_12_15=0.000821018062397,
        ACAT_19_19=0.995815899582,
        ACAT_19_15=0.00418410041841,
        ACAT_44_44=0.888888888889,
        ACAT_44_40=0.111111111111,
        AGAT_22_22=0.92,
        AGAT_22_14=0.08,
        AGAT_36_36=0.928571428571,
        AGAT_36_32=0.0714285714286,
        AGAT_39_39=0.909090909091,
        AGAT_39_35=0.0909090909091,
        1.0e-09)# dummy default rate1.0e-09
}

########################

stringToNumList <- function(sequence){
splat <- strsplit(as.character(sequence), ",")[[1]]
numbs <- as.numeric(splat)
return(numbs)
}

########################

ProbErrorChoice <- function(errorParameter,biasParameter,upstreamRepeat,downstreamRepeat,motifsize){
    if (downstreamRepeat == upstreamRepeat+motifsize) errorParameter*biasParameter
    else if (downstreamRepeat == upstreamRepeat) 1-errorParameter
    else if (downstreamRepeat == upstreamRepeat-motifsize) errorParameter*(1-biasParameter)
    else 0
}


########################
MultinormCoeff<-function(sequence){
    alllength=base::sum(sequence)
    return(factorial(alllength)/(base::prod(factorial(sequence))))
}

########################

jointProbLoci<-function(locusdataset){
    
    run1=locusdataset[[3]]
    run2=locusdataset[[4]]

    count_r1=table(run1)
    count_r2=table(run2)
    allcount=list(count_r1,count_r2)
    names(allcount)=1:2

    ###################
    # cDNA loop start #
    ###################
    prob_value_collector=0.0
    
    Seq123list_value_collector=c()
    
    for (C1_minus2 in 0:bin_size){
    for (C1_minus1 in 0:(bin_size-C1_minus2)){
    for (C1_equal in 0:(bin_size-C1_minus2-C1_minus1)){
    for (C1_plus1 in 0:(bin_size-C1_minus2-C1_minus1-C1_equal)){
        
        C1_plus2=bin_size-C1_minus2-C1_minus1-C1_equal-C1_plus1
        C1_count=c(C1_minus2,C1_minus1,C1_equal,C1_plus1,C1_plus2)

        Seq123_value_collector=1
        for (Seqrun in 1:1){

            Seq_CNA1_value_collector=c()
            for (Seqform in as.integer(names(allcount[[as.character(Seqrun)]]))){
                Seq_CNA1_value_collector=c(Seq_CNA1_value_collector,sum(C1_count[1]*SeqErrorSwtich(motifSTR,fiveform[1],Seqform),C1_count[2]*SeqErrorSwtich(motifSTR,fiveform[2],Seqform),C1_count[3]*SeqErrorSwtich(motifSTR,fiveform[3],Seqform),C1_count[4]*SeqErrorSwtich(motifSTR,fiveform[4],Seqform),C1_count[5]*SeqErrorSwtich(motifSTR,fiveform[5],Seqform))/bin_size)
            }
            Seq123_value_collector=Seq123_value_collector*dmultinom(c(allcount[[as.character(Seqrun)]],0),prob=c(Seq_CNA1_value_collector,max(0,1-sum(Seq_CNA1_value_collector))))
        }
        Seq123list_value_collector=c(Seq123list_value_collector,Seq123_value_collector)
    }}}}

    Seq456list_value_collector=c()
    for (C2_minus2 in 0:bin_size){
    for (C2_minus1 in 0:(bin_size-C2_minus2)){
    for (C2_equal in 0:(bin_size-C2_minus2-C2_minus1)){
    for (C2_plus1 in 0:(bin_size-C2_minus2-C2_minus1-C2_equal)){
        
        C2_plus2=bin_size-C2_minus2-C2_minus1-C2_equal-C2_plus1
        C2_count=c(C2_minus2,C2_minus1,C2_equal,C2_plus1,C2_plus2)
   
        Seq456_value_collector=1
        for (Seqrun in 2:2){
            Seq_CNA2_value_collector=c()
            for (Seqform in as.integer(names(allcount[[as.character(Seqrun)]]))){
                Seq_CNA2_value_collector=c(Seq_CNA2_value_collector,sum(C2_count[1]*SeqErrorSwtich(motifSTR,fiveform[1],Seqform),C2_count[2]*SeqErrorSwtich(motifSTR,fiveform[2],Seqform),C2_count[3]*SeqErrorSwtich(motifSTR,fiveform[3],Seqform),C2_count[4]*SeqErrorSwtich(motifSTR,fiveform[4],Seqform),C2_count[5]*SeqErrorSwtich(motifSTR,fiveform[5],Seqform))/bin_size)
            }
            Seq456_value_collector=Seq456_value_collector*dmultinom(c(allcount[[as.character(Seqrun)]],0),prob=c(Seq_CNA2_value_collector,max(0,1-sum(Seq_CNA2_value_collector))))            
        }
        Seq456list_value_collector=c(Seq456list_value_collector,Seq456_value_collector)
    }}}}

    prob_value_collector=t(Seq123list_value_collector) %*% Pmatrix %*%  Seq456list_value_collector

    return(log10(prob_value_collector))
}

########################

jointProbPDF <- function(Parameter_estimate,tempdataset,bin_size,N_core){   #DNAform,motifSTR,run1,run2,run3,run4,run5,run6)

    ptm <- proc.time()
    ep1=Parameter_estimate[1]	
    Q1=Parameter_estimate[2]
    ep2=Parameter_estimate[3]
    Q2=Parameter_estimate[4]
    
    DNAform=tempdataset[[1]][[1]]
    motifSTR=tempdataset[[1]][[2]]
    motifSTRclass=nchar(motifSTR)
    threeform=c(DNAform-motifSTRclass,DNAform,DNAform+motifSTRclass)
    fiveform=c(DNAform-(2*motifSTRclass),DNAform-motifSTRclass,DNAform,DNAform+motifSTRclass,DNAform+(2*motifSTRclass))


    pdimension=nsimplex(5,bin_size)
    Pmatrix=matrix(0,pdimension,pdimension)
    for (R_minus1 in 0:bin_size){
    for (R_equal in 0:(bin_size-R_minus1)){
        R_plus1=bin_size-R_minus1-R_equal
        r_count=c(R_minus1,R_equal,R_plus1) #r_count=vector of rna count by rna length        
        R_value_collector=dmultinom(r_count,prob=c(ProbErrorChoice(ep1,Q1,DNAform,threeform[1],motifSTRclass),ProbErrorChoice(ep1,Q1,DNAform,threeform[2],motifSTRclass),ProbErrorChoice(ep1,Q1,DNAform,threeform[3],motifSTRclass)))      
    
    C_value_collector=c()
    for (C1_minus2 in 0:bin_size){
    for (C1_minus1 in 0:(bin_size-C1_minus2)){
    for (C1_equal in 0:(bin_size-C1_minus2-C1_minus1)){
    for (C1_plus1 in 0:(bin_size-C1_minus2-C1_minus1-C1_equal)){
        
        C1_plus2=bin_size-C1_minus2-C1_minus1-C1_equal-C1_plus1
        C1_count=c(C1_minus2,C1_minus1,C1_equal,C1_plus1,C1_plus2)
        CNA1_R1_value_collector=c()
        for (C1 in fiveform){
            CNA1_R1_value_collector=c(CNA1_R1_value_collector,sum(r_count[1]*ProbErrorChoice(ep2,Q2,threeform[1],C1,motifSTRclass),r_count[2]*ProbErrorChoice(ep2,Q2,threeform[2],C1,motifSTRclass),r_count[3]*ProbErrorChoice(ep2,Q2,threeform[3],C1,motifSTRclass))/bin_size)
        }

        #CNA1_value_collector=dmultinom(c(C1_count,0),prob=c(CNA1_R1_value_collector,1-sum(CNA1_R1_value_collector)) )
        CNA1_value_collector=dmultinom(C1_count,prob=CNA1_R1_value_collector)
        C_value_collector=c(C_value_collector,CNA1_value_collector)
    
    }}}}

    temp_C_matrix=C_value_collector %*% t(C_value_collector) * R_value_collector
    #print(temp_C_matrix)
    Pmatrix=Pmatrix+temp_C_matrix
    }}
    #print(Pmatrix)
    #locicombine=0.0
    clus <- makeCluster(N_core,type="SOCK")
    clusterExport(clus,"jointProbLoci")
    clusterExport(clus,"SeqErrorSwtich")
    clusterExport(clus,c("motifSTRclass","motifSTR","DNAform","Pmatrix","bin_size","fiveform"),envir=environment())
    #clusterExport(clus,"bin_size")
    locicombinelist=parLapply(clus,tempdataset,jointProbLoci)
    #print(locicombinelist)
    locicombine=sum(unlist(locicombinelist))
    stopCluster(clus)
    
    #End of locicombine loop
    print(c('time',proc.time() - ptm))

    print(locicombine)
    print(c("ep1",ep1))
    print(c("Q1",Q1))
    print(c("ep2",ep2))
    print(c("Q2",Q2))
    return(locicombine)
                
}

##################
## Data read in ##
##################


    
print(c('filename','bin_size','N_core','lower_RDD_RT','upper_RDD_RT','step_RDD_RT','step_RDD_RT_expansion'))
print('filename')
print(argsL$filename)
print('bin_size')
print(argsL$bin_size)
print('N_core')
print(argsL$N_core)
print('lower_RDD_RT')
print(argsL$lower_RDD_RT)
print('upper_RDD_RT')
print(argsL$upper_RDD_RT)
print('step_RDD_RT')
print(argsL$step_RDD_RT)
print('step_RDD_RT_expansion')
print(argsL$step_RDD_RT_expansion)
#filename=args[1]
oridataset=read.table(argsL$filename,colClasses = "character")
columnname=(c('locus','DNA','motif','class','functional','seq1','seq2'))
colnames(oridataset)=columnname
transformdataset=list()
for (line in 1:nrow(oridataset)){
    individualdataset=list(as.numeric(oridataset$DNA[line]),oridataset$motif[line]
    ,stringToNumList(oridataset$seq1[line]),stringToNumList(oridataset$seq2[line]))
    
    transformdataset[[length(transformdataset)+1]] <- individualdataset
}

################
## optimizing ##
################
    tempep1=10**runif(1,-9,-0.3)
    tempQ1=runif(1,0,1)
    tempep2=10**runif(1,-9,-0.3)
    tempQ2=runif(1,0,1)
    initialparameter=c(tempep1,tempQ1,tempep2,tempQ2)
    print(c("initialparameter",initialparameter))
    optimresult=optim(initialparameter,jointProbPDF,NULL,method="L-BFGS-B",
    lower=c(argsL$lower_RDD_RT,0,argsL$lower_RDD_RT,0),
    upper=c(argsL$upper_RDD_RT,1,argsL$upper_RDD_RT,1),
    control=list(ndeps=c(argsL$step_RDD_RT,argsL$step_RDD_RT_expansion,argsL$step_RDD_RT,argsL$step_RDD_RT_expansion),
    fnscale=-1),
    tempdataset=transformdataset,bin_size=argsL$bin_size,N_core=argsL$N_core)

print('done')
