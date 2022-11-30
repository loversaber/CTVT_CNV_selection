####Step0. Get the matrix (meta_matrix.csv) with sampling information of samples
###Run the script:
python3 merge_matrix.py

################################################
##################END of Step0##################
################################################

####Step1. Copy bam files in project 1630 (447) and 1531 (453), total 900 bam files.
####1.1 Script for copying bam files:
####447 bam files in project 1630:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files/project_1630$ cat cp_bam.sh
for i in `ls /nfs/irods-cgp-*/intproj/1630/sample/*/*.bam`
do
	name=`echo $i | cut -d "/" -f 7`
	echo $name
	ls -alth $i
	cp $i ./
done
####End of the script####

####453 bam files in project 1531:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files/project_1531$ cat cp_bam.sh
for i in `ls /nfs/irods-cgp*/intproj/1531/sample/*/*.bam`
do
	name=`echo $i | cut -d "/" -f 7`
	echo $name
	ls -alth $i
	cp $i ./
done
####End of the script####

####1.2 Check 4 samples without data (There were 451 samples for the project 1630 on the canapps website, but only 447 were found on Farm5. So I checked which four samples are missing)
####Script:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files$ cat find4_missed_samples.py
import re
import pandas as pd

files451="project_1630_451samples"
files447="project_1630_447samples"

l451=pd.read_csv(files451,header=None)[0].to_list()

l447=[]
with open(files447,"r") as f:
    for i in f:
        i1=i.strip().split("/")[-1].split(".v1.sample.dupmarked.bam")[0]
        i2=i.strip().split("/")[-2]
        if i1!=i2:
            print(i1,i2,i1==i2)
        l447.append(i1)

l4=set(l451)-set(l447)
print(l4)
####End of the script####

####Run the script:
(base) ql4@farm5-head1:~/lustre_ql4/work/ctvt_mito/bam_files$ python3 find4_missed_samples.py
{'1607H-Dog', '1562B-Dog', '1740T-Dog', '1663T-Dog'}


####1.3 There are some duplicate samples in two projects. Therefore, their names should first be modified to avoid errors in the SNP calling. (add suffix like 'p1'(1531) or 'p2'(1630) to sample name)

####Script for checking duplicated samples:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files$ cat check_duplicate_samples.py
import re
import pandas as pd

files451="project_1630_451samples"
files447="project_1531_453samples"

l451=pd.read_csv(files451,header=None)[0].to_list()
deleted_4samples=['1740T-Dog', '1607H-Dog', '1562B-Dog', '1663T-Dog']
l_1630=[i for i in l451 if i not in deleted_4samples]

l_1531=pd.read_csv(files447,header=None)[0].to_list()
l4=[i for i in l_1630 if i in l_1531]

print(len(l4),len(l_1531),len(l_1630))

import glob

all1630=glob.glob("/nfs/irods-cgp-*/intproj/1630/sample/*/*.bam")
all1630=[i.split("/")[-1] for i in all1630]
print(len(all1630))
print(all1630[:5])
all1531=glob.glob("/nfs/irods-cgp-*/intproj/1531/sample/*/*.bam")
all1531=[i.split("/")[-1] for i in all1531]
print(len(all1531))
print(all1531[:5])
a=[i for i in all1630 if i in all1531]
print(a[:5])
print(len(a))

####End of the script####

####1.4 Script for replacing name:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files$ cat change_name_for_overlapped_samples.sh
samtools=/software/CGP/modules/installs/samtools/samtools-1.14/bin/samtools


for i in `ls ../$1/*.bam`
do
	a=${i%%.bam}
	b=${a##*/}
	d=`echo $b|cut -d"." -f1`
	c=$d"_"$2
	echo $c
	new_bam_file=$c".new.bam"
	cmd="$samtools addreplacerg -@ 20 -r "ID:$c" -r "SM:$c" -r "LB:$c" -r "PL:ILLUMINA" -o $new_bam_file $i"
	echo $cmd
	$cmd
	$samtools index -@ 20 -b $new_bam_file
	#ls -alth $i
done
####End of the script####

####1.4 Run the script for changing name of samples (add suffix 'p1', 'p2' to sample name) in these two projects:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files/all_bam_files$ cat generate_new_bam_files_p1.sh
script=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/bam_files

bash $script/change_name_for_overlapped_samples.sh project_1531 p1

####End of the script####

(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/bam_files/all_bam_files$ cat generate_new_bam_files_p2.sh
script=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/bam_files

bash $script/change_name_for_overlapped_samples.sh project_1630 p2

####End of the script####

####1.5 All 900 bam files for SNP calling are stored at "/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/bam_files/bam_files_900"

################################################
##################END of Step1##################
################################################

####Step2. Run Somatypus to do SNP calling (Very large memory needed when submitting the job on Farm5, max memory is about ~53G)
(base) ql4@farm5-head1:~/lustre_ql4/work/ctvt_mito/run_somatypus$ bsub -o run_somatypus900.o.%J -e run_somatypus900.e.%J -G team267-grp -q long -R "select[mem>100000]  rusage[mem=100000]" -M100000 -n8 -R"span[hosts=1]" bash cmd_somatypus.sh

####2.1 Script for running Somatypus:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/run_somatypus$ cat cmd_somatypus.sh
export MODULEPATH=/software/CGP/modules/modulefiles:$MODULEPATH
module load /software/CGP/modules/modulefiles/vcftools/0.1.16
#VCFTOOLS_PATH=/software/CGP/modules/installs/vcftools/vcftools-0.1.16/bin
VCFTOOLS=$(which vcftools)
VCFTOOLS_PATH=$(dirname $VCFTOOLS)

#Define a specialised PATH for working with legacy platypus
PLATYPUS_PATH=/nfs/dog_n_devil/adrian/software/somatypus/src:/nfs/dog_n_devil/adrian/software/somatypus/utils:/nfs/dog_n_devil/adrian/software/platypus/bin:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/opt/renci/icommands//bin:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/software/badger/bin:/software/oracle-ic-11.2:/software/bin:/software/CGP/bin:$VCFTOOLS_PATH
PLATYPUS_LDFLAGS=-L/nfs/dog_n_devil/adrian/software/platypus/htslib-1.2.1
PLATYPUS_LD_LIBRARY_PATH=/nfs/dog_n_devil/adrian/software/platypus/htslib-1.2.1

#Run somatypus
somatypus=/nfs/users/nfs_k/kg8/lustre_ms/projects/ctvt_horizontal_transfer/scripts/somatypus/src/somatypus

PATH=$PLATYPUS_PATH \
LDFLAGS=$PLATYPUS_LD_FLAGS \
LD_LIBRARY_PATH=$PLATYPUS_LD_LIBRARY_PATH \

python3=/nfs/users/nfs_q/ql4/anaconda3/bin/python3
script_path=/nfs/users/nfs_q/ql4/script/cmd_index_error_files
new_samtools_path=/software/CGP/modules/installs/samtools/samtools-1.14/bin
old_samtools_path=/nfs/users/nfs_q/ql4/lustre_ql4/tools/samtools-1.2/samtools

bam_files_path=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/bam_files/bam_files_900

somatypus_opt_path=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result

ref_path=/nfs/users/nfs_q/ql4/lustre_ql4/data/reference_sequences/CanFam3.1/Canis_familiaris.CanFam3.1.70.dna.toplevel.fa
mito_region=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/test_run_somatypus/mitochondrial_region.txt

$somatypus -i $bam_files_path -g $ref_path -o $somatypus_opt_path -r $mito_region -c 8 -p '--rmsmqThreshold=20 --qdThreshold=5 --maxReads=50000000 --bufferSize=2500'

####End of the script####

####2.2 The final VCF file (Somatypus_SNVs_final.vcf) is at:
/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result

################################################
##################END of Step2##################
################################################


####Step3. Construct the RAxML tree
####According to the substitution VCF file and MT reference sequence, I made a fasta file containing 900 samples' consensus sequences to be used directly as the multiple sequence alignment (MSA) file. Considering that MSA of 900 sequences takes a long time and only substitutions were used here, I thought it was not necessary to do MSA in software. 

####3.1 Get the fasta file (RaxML_MSA_v1.fa), which contains consensus sequences of all 900 samples:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result$ time python3 calc_VAF_simple_method_v2.1.py Somatypus_SNVs_final.vcf &>log_calc_VAF_simple_method_v2.1

####This script also generates files needed for making VAF plot, like VAF_matrix_added_category_v2.1.txt (clade-definiting substitutions are marked with corresponding colors in this file), Tumor_no_host.txt, Tumor_has_host.txt, the file numberID_samples_distribution.txt which contains tumour and corresponding host information, and files containing all the substitutions for each sample(they are in the fold "/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/VAF_files" with suffix ".tmp_vaf").


####3.2 Get the MT reference genome from CanFam3.1:
samtools faidx /nfs/users/nfs_q/ql4/lustre_ql4/data/reference_sequences/CanFam3.1/Canis_familiaris.CanFam3.1.70.dna.toplevel.fa MT > ref.fa

####3.3 Merge the MT reference sequence with two coyote sequences:
(base) ql4@farm5-head1:~/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result$ cat ../sequence_DQ480509.fasta ../sequence_DQ480510.fasta ../ref.fa > ref_with_coyote2.fa

####3.4 Run the MSA of these 3 sequences:
(base) ql4@farm5-head1:~/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result$ /nfs/users/nfs_q/ql4/anaconda3/bin/mafft --auto --preservecase ref_with_coyote2.fa > ref_with_coyote2_MSA.fa

####3.5 Using the differences between the MT ref and 2 coyote sequences (from ref_with_coyote2_MSA.fa) to adjust the 900 consensus sequences (inserting gaps) in RaxML_MSA_v1.fa. Finally, the file RaxML_MSA_v2.3.fa was generated, which contains 903 sequences (includes MT ref and 2 coyote sequences):
(base) ql4@farm5-head1:~/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result$ python3 get_new_ref2.py &>log_get_new_ref2

####3.6 Run RAxML:
####The file RaxML_MSA_v2.4.fa was used here. I have made some modifications to the sequence ID (sample names) in this new fasta file. There should be no blank in the sequence ID in fasta file, otherwise the sample names in the tree file output by RAxML will not contain the information after the first blank.

####Submit the job:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result$ bsub -o Raxml.o.%J -e Raxml.e.%J -G team267-grp -q basement -R "select[mem>25000]  rusage[mem=25000]" -M25000 -n8 -R"span[hosts=1]" bash run_raxml3.sh

####Script:
(base) ql4@farm5-head2:~/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result$ cat run_raxml3.sh
fa_dir=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result

/nfs/users/nfs_q/ql4/anaconda3/bin/raxmlHPC-PTHREADS-AVX2 -s $fa_dir/RaxML_MSA_v2.4.fa -f a -n raxml3 -m GTRGAMMA -c 4 -T 8 -x 12345 -p 12345 -# autoMRE
####End of the script

####3.7 Using the output file RAxML_bipartitions.raxml3 to display the tree
####Rotate the tree in iTOL (online version) to make the 2 coyote samples as the outgroup. Finally, I downloaded the tree saved as 2022-06-14-RAxML_Tree.pdf.

################################################
##################END of Step3##################
################################################

####Step4. VAF plot (based on file VAF_matrix_added_category_v2.1.txt)
####Scripts plot_All_T_without_H_v2.R and plot_All_T_with_H_v2.R are used to get VAF plots (VAF_T_without_H_v2.pdf, VAF_T_with_H_v2.pdf(top: Host, bottom: Tumour)).

Definition of true substitutions (indicated by black dots):

Host: 
"present": the number of reads for alternative allele (nAlt) >=3 and variant allele frequency (VAF)>=0.1 

Tumour with host:
1) If not present in matched host then "present": nAlt >=3 and VAF>=0.1
2) If present in matched host then "present": nAlt >=3 and VAF>=0.9

Tumour without host:
present: nAlt >=3 and VAF>=0.5


################################################
##################END of Step4##################
################################################


All scripts and all crucial output files can be found in each sub-folder or on Farm5.
