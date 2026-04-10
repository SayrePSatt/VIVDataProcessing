This contains multiple processing codes to conduct free decay tests, kinematic data processing, and creating contour plots from SPIV data

free_decay_fitter_newfilestructure.m: This file processes the freedecay tests to output the natural frequencies and damping ratios
    Dependencies: None
    Input: MM-DD-YYYY_HH-MM_freedecay_(spring config)k_(air or water)_0(test).0_00.00.dat
    Output: MM-DD-YYYY_freedecay_results_(spring config)k_(air or water).dat

structural_processing_code.m: This code processes all of the data from the encoder to create output plots
    Dependencies: plot_fn.m, plot_psd_fn_newstructure.m, plot_fn_prc.m, psdd3_sayre., retrievephase1.m, retrievephase2.m, 
        FivePointDiff.m, norm_PSD_calc.m, ave_bounds_newstructure.m, gov_2005_amplift.csv, Gov_Will_SphereA.csv, GovWill_UStar.csv,
        griffin_govwill.csv, sareen2018b_ampphase.csv, Sareen_FS_Ustar.csv, unorm_reference.csv, pumpFit_freq2velo.mat
    Input: MM-DD-YYYY_HH-MM_(distance ratio)D_(sphere diameter ratio)D_(spring setup)k_(VIV sphere diameter)mm_VIV_(U*)_(pump frequency).dat
    Output: .png, .pdf, and .fig files

SPIV_phase_averaging.m: This code takes in kinematic data, bins the data, and exports which image indices to export from davis
    Dependencies:FivePointDiff.m
    Input: MM-DD-YYYY_HH-MM_(distance ratio)D_(sphere diameter ratio)D_(spring setup)k_(VIV sphere diameter)mm_VIV_(U*)_(pump frequency).dat
    Output: 

imx_readandplot: Takes the .imx files exported from DaVIS and exports images of the contours
    Dependencies: num2roman.m
    Input: .imx files
    Output: .pdf, .png, .fig files

