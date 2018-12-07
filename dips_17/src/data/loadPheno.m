function [ pheno ] = loadPheno(id, col, opt)

%------------------------------------------------------------------------%
% Given a list of subject id (cell array), extract the specific phenotype 
% in the col-th column of phenotype/Phenotypic_V1_0b_preprocessed1.csv
% 
% Phenotype column legend
% 	1			NULL
% 	2			Unnamed: 0
% 	3			SUB_ID
% 	4			X
% 	5			subject
% 	6			SITE_ID
% 	7			FILE_ID
% 	8			DX_GROUP
% 	9			DSM_IV_TR
% 	10		AGE_AT_SCAN
% 	11		SEX
% 	12		HANDEDNESS_CATEGORY
% 	13		HANDEDNESS_SCORES
% 	14		FIQ
% 	15		VIQ
% 	16		PIQ
% 	17		FIQ_TEST_TYPE
% 	18		VIQ_TEST_TYPE
% 	19		PIQ_TEST_TYPE
% 	20		ADI_R_SOCIAL_TOTAL_A
% 	21		ADI_R_VERBAL_TOTAL_BV
% 	22		ADI_RRB_TOTAL_C
% 	23		ADI_R_ONSET_TOTAL_D
% 	24		ADI_R_RSRCH_RELIABLE
% 	25		ADOS_MODULE
% 	26		ADOS_TOTAL
% 	27		ADOS_COMM
% 	28		ADOS_SOCIAL
% 	29		ADOS_STEREO_BEHAV
% 	30		ADOS_RSRCH_RELIABLE
% 	31		ADOS_GOTHAM_SOCAFFECT
% 	32		ADOS_GOTHAM_RRB
% 	33		ADOS_GOTHAM_TOTAL
% 	34		ADOS_GOTHAM_SEVERITY
% 	35		SRS_VERSION
% 	36		SRS_RAW_TOTAL
% 	37		SRS_AWARENESS
% 	38		SRS_COGNITION
% 	39		SRS_COMMUNICATION
% 	40		SRS_MOTIVATION
% 	41		SRS_MANNERISMS
% 	42		SCQ_TOTAL
% 	43		AQ_TOTAL
% 	44		COMORBIDITY
% 	45		CURRENT_MED_STATUS
% 	46		MEDICATION_NAME
% 	47		OFF_STIMULANTS_AT_SCAN
% 	48		VINELAND_RECEPTIVE_V_SCALED
% 	49		VINELAND_EXPRESSIVE_V_SCALED
% 	50		VINELAND_WRITTEN_V_SCALED
% 	51		VINELAND_COMMUNICATION_STANDARD
% 	52		VINELAND_PERSONAL_V_SCALED
% 	53		VINELAND_DOMESTIC_V_SCALED
% 	54		VINELAND_COMMUNITY_V_SCALED
% 	55		VINELAND_DAILYLVNG_STANDARD
% 	56		VINELAND_INTERPERSONAL_V_SCALED
% 	57		VINELAND_PLAY_V_SCALED
% 	58		VINELAND_COPING_V_SCALED
% 	59		VINELAND_SOCIAL_STANDARD
% 	60		VINELAND_SUM_SCORES
% 	61		VINELAND_ABC_STANDARD
% 	62		VINELAND_INFORMANT
% 	63		WISC_IV_VCI
% 	64		WISC_IV_PRI
% 	65		WISC_IV_WMI
% 	66		WISC_IV_PSI
% 	67		WISC_IV_SIM_SCALED
% 	68		WISC_IV_VOCAB_SCALED
% 	69		WISC_IV_INFO_SCALED
% 	70		WISC_IV_BLK_DSN_SCALED
% 	71		WISC_IV_PIC_CON_SCALED
% 	72		WISC_IV_MATRIX_SCALED
% 	73		WISC_IV_DIGIT_SPAN_SCALED
% 	74		WISC_IV_LET_NUM_SCALED
% 	75		WISC_IV_CODING_SCALED
% 	76		WISC_IV_SYM_SCALED
% 	77		EYE_STATUS_AT_SCAN
% 	78		AGE_AT_MPRAGE
% 	79		BMI
% 	80		anat_cnr
% 	81		anat_efc
% 	82		anat_fber
% 	83		anat_fwhm
% 	84		anat_qi1
% 	85		anat_snr
% 	86		func_efc
% 	87		func_fber
% 	88		func_fwhm
% 	89		func_dvars
% 	90		func_outlier
% 	91		func_quality
% 	92		func_mean_fd
% 	93		func_num_fd
% 	94		func_perc_fd
% 	95		func_gsr
% 	96		qc_rater_1
% 	97		qc_notes_rater_1
% 	98		qc_anat_rater_2
% 	99		qc_anat_notes_rater_2
% 	100		qc_func_rater_2
% 	101		qc_func_notes_rater_2
% 	102		qc_anat_rater_3
% 	103		qc_anat_notes_rater_3
% 	104		qc_func_rater_3
% 	105		qc_func_notes_rater_3
% 	106		SUB_IN_SMP
%------------------------------------------------------------------------%


	if nargin < 3
		opt = [];
	end

	% default phenotype: diagnose group 1 or 2
	if nargin < 2
		col = 8;
	end

	% import phenotypic data of ABIDE dataset	
	data = importdata('phenotype/Phenotypic_V1_0b_preprocessed1.csv');
	% id = unique(id);
	
	% initialize
	index = zeros(1,length(id));
    pheno = zeros(1,length(id));
    
	% loop through all subjects
	for i = 1:length(id),

		index(i) = find(strcmp(data.textdata(:,3), id{i}(3:end)));
        pheno(i) = str2double(data.textdata{index(i), col});
		% disp(['subject ',num2str(i),': pheno ', data.textdata(index(i), col)]);

    end

    pheno(isnan(pheno)) = -9999;

    if isfield(opt, 'meta_dir')
    	dlmwrite([opt.meta_dir,'/phenotype_sample_',num2str(col),'.txt'],...
    				pheno, 'delimiter', ' ');
    else
		dlmwrite(['phenotype_sample_',num2str(col),'.txt'],...
					pheno, 'delimiter', ' ');
	end


end