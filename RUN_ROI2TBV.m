clear all; close all; clc
odir = 'C:\Users\Tobias Hauser\Desktop\BEN\2017_05_25\';
structMRI = [odir 'sMP01857-0013-00001-000224-01_T1w.nii'];
functMRI = [odir 'fMP01859-0006-00001-000001-01.nii'];
ROI = [odir 'SNVTA_MP01857.nii'];
options.out_fname = 'SNVTA_MP01857.roi';
ROI2TBV(structMRI,functMRI,ROI,options)