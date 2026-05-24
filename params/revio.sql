SET VARIABLE min_q               = 28;
SET VARIABLE min_len             = 3500;
SET VARIABLE max_len             = 6000;
SET VARIABLE fw1                 = 'CAAGCAGAAGACGGCATACGAGAT';
SET VARIABLE fw2                 = 'AGRGTTYGATYMTGGCTCAG';
SET VARIABLE rv1                 = 'AATGATACGGCGACCACCGAGATC';
SET VARIABLE rv2                 = 'CGACATCGAGGTGCCAAAC';
SET VARIABLE fw1_rc              = 'ATCTCGTATGCCGTCTTCTGCTTG';
SET VARIABLE fw2_rc              = 'CTGAGCCATGATCAAACTCT';
SET VARIABLE rv1_rc              = 'GATCTCGGTGGTCGCCGTATCATT';
SET VARIABLE rv2_rc              = 'GTTTGGCACCTCGATGTCG';
SET VARIABLE umi_pair_pattern    = '^([ATCG]{3}[CT][AG]){3}[ATCG]{3}([ATCG]{3}[CT][AG]){3}[ATCG]{3}$';
SET VARIABLE umi_cluster_id      = 0.95;
SET VARIABLE umi_max_nm_per_half = 2;
SET VARIABLE ume_mean_max        = 3.0;
SET VARIABLE ume_sd_max          = 30.0;
SET VARIABLE ro_frac             = 0.05;
SET VARIABLE bin_cluster_ratio   = 10.0;
SET VARIABLE min_bin_size        = 5;
SET VARIABLE umi_coverage_min    = 5;
SET VARIABLE subcluster_id       = 0.998;
SET VARIABLE variant_min_support = 2;
SET VARIABLE snp_min_alt_depth   = 2;
-- positive_ref_path: path to 16S reference FASTA for sortmerna filter
-- download from ftp.microbio.me (88_otus.fasta, same as deblur)
-- SET VARIABLE positive_ref_path = '/path/to/88_otus.fasta';
