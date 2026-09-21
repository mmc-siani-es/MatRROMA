clc; clear;
download_zenodo_data()%% function to download the required files from Zenodo
matrroma("Lvl_crossing_RR.txt",50,512,1,"COVD",50,'downsampling_rate',64)
