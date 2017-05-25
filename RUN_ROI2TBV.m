clear all; close all; clc
odir = 'C:\Users\Tobias Hauser\Desktop\BEN\2017_05_25\';
structMRI = [odir 'anon_sMP01857-0013-00001-000224-01_T1w.nii'];
functMRI = [odir 'fMP01859-0006-00001-000001-01.nii'];
ROI = [odir 'T1_RightVentricle_LeftEye2.nii'];
options.out_fname = 'test_fun2.roi';
ROI2TBV(structMRI,functMRI,ROI,options)