
module load python3/3.5.0

source venv/bin/activate

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_long/m160422_005059_42183_c101007230310000001823234509161690_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_long/m160422_005059_42183_c101007230310000001823234509161690_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_long/m160422_005059_42183_c101007230310000001823234509161690_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/26APR16_PacBio_HelaS3_unsyn_long_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_short/m160422_070719_42183_c101007230310000001823234509161691_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_short/m160422_070719_42183_c101007230310000001823234509161691_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_short/m160422_070719_42183_c101007230310000001823234509161691_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/26APR16_PacBio_HelaS3_unsyn_short_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/07JUL16_Mitotic_HeLa_3CPacbio_long/m160702_060221_42183_c101009742550000001823230910211634_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/07JUL16_Mitotic_HeLa_3CPacbio_long/m160702_060221_42183_c101009742550000001823230910211634_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/07JUL16_Mitotic_HeLa_3CPacbio_long/m160702_060221_42183_c101009742550000001823230910211634_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/07JUL16_Mitotic_HeLa_3CPacbio_long_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/07JUL16_Mitotic_HeLa_3CPacbio_short/m160701_234601_42183_c101009742550000001823230910211633_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/07JUL16_Mitotic_HeLa_3CPacbio_short/m160701_234601_42183_c101009742550000001823230910211633_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/07JUL16_Mitotic_HeLa_3CPacbio_short/m160701_234601_42183_c101009742550000001823230910211633_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/07JUL16_Mitotic_HeLa_3CPacbio_short_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_long/m160812_063853_42183_c101012312550000001823230910211625_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_long/m160812_063853_42183_c101012312550000001823230910211625_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_long/m160812_063853_42183_c101012312550000001823230910211625_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/15AUG16_PacBio_HelaCluster1_long_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_short/m160812_125513_42183_c101012312550000001823230910211626_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_short/m160812_125513_42183_c101012312550000001823230910211626_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_short/m160812_125513_42183_c101012312550000001823230910211626_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/15AUG16_PacBio_HelaCluster1_short_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_MitoticHela_long/m160811_180639_42183_c101012312550000001823230910211623_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_MitoticHela_long/m160811_180639_42183_c101012312550000001823230910211623_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_MitoticHela_long/m160811_180639_42183_c101012312550000001823230910211623_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/15AUG16_PacBio_MitoticHela_long_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_MitoticHela_short/m160812_002233_42183_c101012312550000001823230910211624_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_MitoticHela_short/m160812_002233_42183_c101012312550000001823230910211624_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_MitoticHela_short/m160812_002233_42183_c101012312550000001823230910211624_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/15AUG16_PacBio_MitoticHela_short_digested.fastq"


bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_long_resequence/m170127_073632_42183_c101122332550000001823252305221794_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_long_resequence/m170127_073632_42183_c101122332550000001823252305221794_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_long_resequence/m170127_073632_42183_c101122332550000001823252305221794_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/26APR16_PacBio_HelaS3_unsyn_long_resequence_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_short_resequence/m170127_140122_42183_c101122332550000001823252305221795_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_short_resequence/m170127_140122_42183_c101122332550000001823252305221795_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/26APR16_PacBio_HelaS3_unsyn_short_resequence/m170127_140122_42183_c101122332550000001823252305221795_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/26APR16_PacBio_HelaS3_unsyn_short_resequence_digested.fastq"


bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_long_resequence/m170127_202630_42183_c101122332550000001823252305221796_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_long_resequence/m170127_202630_42183_c101122332550000001823252305221796_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_long_resequence/m170127_202630_42183_c101122332550000001823252305221796_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/15AUG16_PacBio_HelaCluster1_long_resequence_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_short_resequence/m170128_025137_42183_c101122332550000001823252305221797_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_short_resequence/m170128_025137_42183_c101122332550000001823252305221797_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/15AUG16_PacBio_HelaCluster1_short_resequence/m170128_025137_42183_c101122332550000001823252305221797_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/15AUG16_PacBio_HelaCluster1_short_resequence_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/21SEP_PacBio_long/m170921_220322_42183_c101240962550000001823288412071700_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/21SEP_PacBio_long/m170921_220322_42183_c101240962550000001823288412071700_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/21SEP_PacBio_long/m170921_220322_42183_c101240962550000001823288412071700_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/21SEP_PacBio_long_digested.fastq"

bsub -q short -R "rusage[mem=8000]" "source venv/bin/activate
python3 digest_roi.py '/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/21SEP_PacBio_short/m170922_041942_42183_c101240962550000001823288412071701_s1_p0.1.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/21SEP_PacBio_short/m170922_041942_42183_c101240962550000001823288412071701_s1_p0.2.ccs.fastq::/nl/umw_job_dekker/users/fc85w/pacbio_aligner/roi_files/21SEP_PacBio_short/m170922_041942_42183_c101240962550000001823288412071701_s1_p0.3.ccs.fastq' /nl/umw_job_dekker/users/fc85w/pacbio_aligner/digested_reads/21SEP_PacBio_short_digested.fastq"



