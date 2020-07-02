#! /usr/local/python_anaconda/bin/python3.4

from utils import pbs_jobs
import os
from os import path
from seqFileTools import convert_fasta_to_phylip, get_longest_sequence_name_in_fasta
from file_utilities import set_filenames_for_pbs_runs, check_filename, check_dirname, make_dir
from beast_utilities import configure_log_path

def baseml_runner(ctl, alias = "bml", cmdname="baseml"):
    """
    run baseml program from PAML on cluster
    :param ctl: ctl file path
    :param alias: job name (default: bml)
    :return: job id
    """
    ctl = check_filename(ctl)
    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias); tnum = 1; gmem = 2
    cmds = "echo %s \n/sternadi/home/volume1/taliakustin/software/paml4.8/bin/baseml %s" %(ctl, ctl)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def codeml_runner(ctl, alias = "cml", queue="adistzachi", cmdname="codeml.txt"):
    """
    run baseml program from PAML on cluster
    :param ctl: ctl file path
    :param alias: job name (default: bml)
    :return: job id
    """
    ctl = check_filename(ctl)
    base = os.path.split(ctl)[0]
    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias); tnum = 1; gmem = 2
    cmds = "cd %s\necho %s \n/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s" %(base, ctl, ctl)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds, queue=queue)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def script_runner(cmds, alias = "script", load_python=False, gmem=2, queue="adistzachi", run_after_job=None, cmdname = "script", toRun = True):
    """
    run script on cluster
    :param cmds: script running line
    :param alias: job name (default: script)
    :return: job id
    """
    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias); tnum=1; gmem=gmem
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds, load_python=load_python, run_after_job=run_after_job)
    if toRun:
        job_id = pbs_jobs.submit(cmdfile)
        return job_id


def array_script_runner(cmds, jnum, alias = "script", load_python=False, gmem=1, queue="adistzachi", run_after_job=None, toRun=False):

    """
    run script on cluster as a pbs array
    :param cmds: script running line, should include $PBS_ARRAY_INDEX
    :param alias: job name (default: script)
    :param jnum: number of jobs in the pbs array
    :return: job id
    """

    cmdfile = pbs_jobs.assign_cmdfile_path("script", alias); gmem=gmem
    print(cmdfile, alias, queue, jnum, gmem, cmds)
    pbs_jobs.create_array_pbs_cmd(cmdfile, jnum=jnum, alias=alias, queue=queue, gmem=gmem, cmds=cmds, load_python=load_python, run_after_job=run_after_job)
    if toRun:
        job_id = pbs_jobs.submit(cmdfile)
        return job_id



def phyml_runner(alignment, alias = "phyml", phylip=True):
    """
    run phyml on cluster (converts tpo phylip if the flag phylip==False)
    :param alignment: alignment file path
    :param alias: job name (default: phyml)
    :param phylip: True if phylip file, False if fasta file
    :return: job id
    """
    alignment = check_filename(alignment)
    if phylip == False:
        alignment = convert_fasta_to_phylip(alignment)
    cmdfile = pbs_jobs.assign_cmdfile_path("phyml", alias); tnum = 1; gmem = 2
    cmds = "/sternadi/home/volume1/shared/tools/PhyML/PhyML_3.0_linux64 -i %s -b 0 -o n" % alignment
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def phyml_aa_runner(alignment, alias = "phyml", phylip=True):
    """
    run phyml on aa alignment on cluster
    :param alignment: alignment file path
    :param alias: job name (default: phyml)
    :param phylip: True if phylip file, False if fasta file
    :return: job id
    """
    alignment = check_filename(alignment)
    if phylip == False:
        alignment = convert_fasta_to_phylip(alignment)
    cmdfile = pbs_jobs.assign_cmdfile_path("phyml", alias); tnum = 1; gmem = 2
    cmds = "/sternadi/home/volume1/shared/tools/PhyML/PhyML_3.0_linux64 -i %s -d aa -q -b 0" % alignment
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def fastml_runner(alignment, tree, outdir = None, log_file = None, alias = "fastml", get_cmdfile_dir = True,
                  additional_params=None, optBranchLen=False,
                  fastml_path="/sternadi/home/volume1/shared/tools/phylogenyCode/programs/fastml/fastml", queue="adistzachi"):
    """
    run fastml from phylogenyCode on cluster
    :param alignment: alignment file path
    :param tree: tree file path
    :param alias: job name (default: fastml)
    :param outdir: output directory for results (default: None - saves in the alignment's dir)
    :param optBranchLen: to optimize branch length? default - true
    :return: job id
    """
    alignment = check_filename(alignment)
    tree = check_filename(tree)
    if outdir == None:
        outdir = os.path.dirname(alignment)
    else:
        outdir = check_dirname(outdir)
    basename = os.path.basename(alignment).split(".")[0].split("_aln")[0]
    newick_tree = outdir + "/" + basename + ".tree.newick.txt"
    ancestor_tree = outdir + "/" + basename + ".tree.ancestor.txt"
    joint_seqs = outdir + "/" + basename + ".seq.joint.txt"
    marginal_seqs = outdir + "/" + basename + ".seq.marginal.txt"
    joint_prob = outdir + "/" + basename + ".prob.joint.txt"
    marginal_prob = outdir + "/" + basename + ".prob.marginal.txt"
    if log_file == None:
        log_file = outdir + "/" + basename + ".fastml.log"
    cmdfile = pbs_jobs.assign_cmdfile_path("fastml.txt", alias); tnum = 1; gmem = 1
    cmds = "%s -s %s -t %s -mn -x %s " \
           "-y %s -j %s -k %s -d %s -e %s -qf -R %s" % (fastml_path, alignment, tree, newick_tree, ancestor_tree, joint_seqs,
                                                 marginal_seqs, joint_prob, marginal_prob, log_file)
    if not optBranchLen:
        cmds += " -b"
    if additional_params != None:
        cmds += " %s" % additional_params
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds, queue=queue)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def mafft_runner(sequence, alignment = None, alias = "mafft"):
    """
    run mafft on cluster
    :param sequence: sequence file (fasta format)
    :param alignment: alignment output file (default: None)
    :param alias: job name (default: mafft)
    :return: job id
    """
    sequence = check_filename(sequence)
    if alignment == None:
        alignment = sequence.split(".fasta")[0] + ".aln"
    alignment = check_filename(alignment, Truefile=False)
    cmds = "/sternadi/home/volume1/taliakustin/software/mafft-7.300-with-extensions/scripts/mafft %s > %s"\
           % (sequence, alignment)
    cmdfile = pbs_jobs.assign_cmdfile_path("mafft", alias); tnum = 1; gmem = 1
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def prank_runner(sequence, alignment=None, alias = "prank"):
    """
    run prank on cluster
    :param sequence: sequence file path (fasta format)
    :param alignment: alignment output file (default: None)
    :param alias: job name (default: prank)
    :return: job id
    """
    if alignment == None:
        alignment = sequence.split(".fasta")[0] + ".aln"
    sequence = check_filename(sequence)
    alignment = check_filename(alignment, Truefile=False)
    cmds = "/powerapps/share/bin/prank -d=%s -o=%s -F" % (sequence, alignment)
    cmdfile = pbs_jobs.assign_cmdfile_path("prank_alignment", alias); tnum=1; gmem=1
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id



def prank_codon_runner(sequence, alignment=None, alias = "prank_codon", tree=None):
    """
    run prank codon on cluster
    :param sequence: sequence file path (fasta format)
    :param alignment: alignment output file (default: None)
    :param alias: job name (default: prank)
    :param tree: an option to add a tree
    :return: job id
    """
    if alignment == None:
        alignment = sequence.split(".fasta")[0] + ".codon_aln"
    sequence = check_filename(sequence)
    alignment = check_filename(alignment, Truefile=False)
    if tree == None:
        cmds = "/powerapps/share/bin/prank -d=%s -o=%s -F -codon" % (sequence, alignment)
    else:
        tree = check_filename(tree)
        cmds = "/powerapps/share/bin/prank -d=%s -t=%s -o=%s -F -codon" % (sequence, tree, alignment)
    cmdfile = pbs_jobs.assign_cmdfile_path("prank_alignment", alias); tnum=1; gmem=5
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def prank_runner_with_tree(sequence, tree, alignment=None, alias = "prank"):
    """
    run prank with tree on cluster
    :param sequence: sequence file path (fasta format)
    :param tree: tree file path
    :param alignment: alignment output file (default: None)
    :param alias: job name (default: prank)
    :return: job id
    """
    if alignment == None:
        alignment = sequence.split(".fasta")[0] + ".aln"
    sequence = check_filename(sequence)
    tree = check_filename(tree)
    alignment = check_filename(alignment, Truefile=False)
    cmds = "/powerapps/share/bin/prank -d=%s -t=%s -o=%s -F" % (sequence, tree, alignment)
    cmdfile = pbs_jobs.assign_cmdfile_path("prank_alignment_with_tree", alias); tnum=1; gmem=5
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def njTree_runner(alignment, tree=None, alias = "njTree"):
    """
    run neighbors-joining tree on cluster
    :param alignment: alignment file path
    :param tree: output tree path (default: None)
    :param alias: job name (default: njTree)
    :return: job id
    """
    if tree == None:
        tree = alignment.split(".")[0] + ".tree"
    alignment = check_filename(alignment)
    tree = check_filename(tree, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("njTree", alias); tnum=1; gmem=2
    cmds = "/sternadi/home/volume1/shared/tools/phylogenyCode/programs/treeUtil/njTreeJCdist -i %s -o %s -an"\
           % (alignment, tree)
    dir = "/sternadi/home/volume1/shared/tools/phylogenyCode/programs/treeUtil/"
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def njTree_codon_runner(alignment, tree=None, alias = "njCodonTree"):
    """
    run neighbors-joining tree on cluster
    :param alignment: alignment file path
    :param tree: output tree path (default: None)
    :param alias: job name (default: njTree)
    :return: job id
    """
    if tree == None:
        tree = alignment.split(".")[0] + ".codon_tree"
    alignment = check_filename(alignment)
    tree = check_filename(tree, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("njTree", alias); tnum=1; gmem=2
    cmds = "/sternadi/home/volume1/shared/tools/phylogenyCode/programs/treeUtil/njTreeJCdist -i %s -o %s -ac"\
           % (alignment, tree)
    dir = "/sternadi/home/volume1/shared/tools/phylogenyCode/programs/treeUtil/"
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def sampling_runner(alignment, amount, sampled_file=None, alias = "sampling", alphabet="an", random=False):
    """
    run sampling on cluster (doesn't sample random seqs)
    :param alignment: alignment file path
    :param amount: amount of sequences to sample
    :param sampled_file: output file (default: None)
    :param alias: job name (default: sampling)
    :param alphabet: type of alphabet to use - an - nucleutides, aa - amino acid, ac  codon
    :return: job id
    """
    alignment = check_filename(alignment)
    if sampled_file == None:
        sampled_file = alignment.split(".aln")[0] + "_sampled_%s.aln" % str(amount)
    if alphabet not in ["an", "aa", "ac"]:
        alphabet = "an"
        print("alphabet type is wrong - changed to default - nucleotides - an")
    output_file = check_filename(sampled_file, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("njTree", alias); tnum=1; gmem=5
    if random:
        cmds = "/sternadi/home/volume1/shared/tools/phylogenyCode/programs/sampling/sampling -i %s -n %s -o %s -%s -r" \
               % (alignment, amount, sampled_file, alphabet)
    else:
        cmds = "/sternadi/home/volume1/shared/tools/phylogenyCode/programs/sampling/sampling -i %s -n %s -o %s -%s"\
           % (alignment, amount, sampled_file, alphabet)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def gzip_runner(file, alias = "gzip"):
    """
    run gzip on cluster
    :param file: input file path
    :param alias: job name (default: gzip)
    :return: job id
    """
    file = check_filename(file)
    cmdfile = pbs_jobs.assign_cmdfile_path("gzip", alias); tnum=1; gmem=2
    cmds = "gzip %s" % file
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def cp_runner(file, dest_file, alias = "cp"):
    """
    run cp on cluster
    :param file: input file path
    :param dest_file: output file path
    :param alias: job name (default: cp)
    :return: job id
    """
    if not "*" in file:
        file = check_filename(file)
        dest_file = check_filename(dest_file, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("cp", alias); tnum=1; gmem=2
    cmds = "cp %s %s" % (file, dest_file)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def baliphy_runner(sequence, alias = "baliphy"):
    """
    run baliphy on cluster
    :param sequence: input sequence file path
    :param alias: job name (default: baliphy)
    :return: job id
    """
    sequence = check_filename(sequence)
    cmdfile = pbs_jobs.assign_cmdfile_path("baliphy_M8_cmd.txt", alias); tnum = 1; gmem = 2
    cmds = "/sternadi/home/volume1/taliakustin/software/bali-phy-2.3.7/bin/bali-phy"\
                                                            + " " + sequence\
                                                            + " -V "
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id



def pear_runner(forward, reverse, output, alias = "pear"):
    """
    run pear (pair-ended merger) on cluster
    :param forward: forward file path
    :param reverse: reverse file path
    :param output: output file path
    :param alias: job name (default pear)
    :return: job id
    """
    forward = check_filename(forward)
    reverse = check_filename(reverse)
    output = check_filename(output, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("pear", alias); tnum = 1; gmem = 2
    cmds = "/usr/local/bin/pear -f %s -r %s -o %s" % (forward, reverse, output)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def ufilter_runner(fastq, output, alias = "ufilter"):
    """
    run ufilter on cluster
    :param fastq: fastq file path
    :param output: output file path
    :param alias: job name (default: ufilter)
    :return: job id
    """
    fastq = check_filename(fastq)
    output = check_filename(output, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("ufilter", alias); tnum = 1; gmem = 2
    cmds = "/usr/local/bin/usearch -fastq_stripleft 5 -fastq_filter %s -fastqout %s" % (fastq, output)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def umerge_runner(forward_fastq, output, alias = "umerge"):
    """
    run umerge on cluster
    :param forward_fastq: forward fastq file path
    :param output: ouput file path
    :param alias: job name (default: umerge)
    :return: job id
    """
    forward_fastq = check_filename(forward_fastq)
    output = check_filename(output, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("umerge", alias); tnum = 1; gmem = 2
    cmds = "/usr/local/bin/usearch" \
           " -fastq_qmaxout 80 -fastq_qmax 80 -fastq_mergepairs %s -fastqout %s -report %s.report"\
           % (forward_fastq, output, output)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def blast_runner(seqfile, dbfile="/sternadi/home/volume1/shared/data/nt/nt", outfile=None, alias="blast", outfmt=6, task="megablast", max_target_seqs=50000):
    """
    run blast on cluster
    :param seqfile: sequence file path
    :param dbfile: db file (default: /sternadi/home/volume1/shared/data/nt/nt)
    :param outfile: output file path (default: None)
    :param alias: job name (blast)
    :return: job id
    """
    seqfile = check_filename(seqfile)
    if outfile != None:
        outfile = check_filename(outfile, Truefile=False)
    else:
        outfile = path.split(seqfile)[0] + "/blast_results.txt"
    cmdfile = pbs_jobs.assign_cmdfile_path("blast_cmd", alias); tnum = 1; gmem = 2
    cmds = f"/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/blastn " \
           f"-query {seqfile} -out {outfile} -db {dbfile} -outfmt  {outfmt} -task {task} -max_target_seqs {max_target_seqs}"
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def blast_output6_runner(seqfile, dbfile, outfile, alias = "blast", task="megablast"):
    """
    run blast on cluster - output as pipeline - format 6
    :param seqfile: sequence file path
    :param dbfile: db file
    :param outfile: output file path
    :param alias: job name (blast)
    :return: job id
    """
    seqfile = check_filename(seqfile)
    if outfile != None:
        outfile = check_filename(outfile, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("blast_cmd", alias); tnum = 1; gmem = 2
    cmds = "/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/blastn"\
                + " -query %s -task {task} -out %s -db %s -outfmt '6 sseqid qseqid qstart qend qstrand sstart send sstrand length btop' " \
                  "-num_alignments 100 -dust no -soft_masking F -perc_identity 85 -evalue 1e-7"\
                  % (seqfile, outfile, dbfile)

    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def blastx_output6_runner(seqfile, outfile, dbfile="/sternadi/home/volume1/shared/data/nr/nr", alias = "blast"):
    """
    run blast on cluster - output as pipeline - format 6
    :param seqfile: sequence file path
    :param dbfile: db file
    :param outfile: output file path
    :param alias: job name (blast)
    :return: job id
    """
    seqfile = check_filename(seqfile)
    if outfile != None:
        outfile = check_filename(outfile, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("blast_cmd", alias); tnum = 1; gmem = 2
    cmds = "/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin/blastx"\
                + " -query %s -out %s -db %s -outfmt 6" % (seqfile, outfile, dbfile)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def bowtie2_runner(bowtie_index_path, fastq_file, sam_output, alias="bowtie2"):
    """
    run bowtie2 - very fast local flag is on
    :param bowtie_index_path: bowtie index file path (output of bowtie2-build)
    :param fastq_file: fastq file path
    :param sam_output: output file for sam file
    :param alias: job name (bowtie2)
    :return: job id
    """
    bowtie_index_path = check_filename(bowtie_index_path, Truefile=False)
    fastq_file = check_filename(fastq_file)
    sam_output = check_filename(sam_output, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("bowtie2", alias); tnum = 1; gmem = 2
    cmds = "/usr/local/bin/bowtie2"\
           + " --very-fast-local -x  %s %s -S %s" % (bowtie_index_path, fastq_file, sam_output)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def bowtie2_build_runner(input_file, output_db_name=None, alias="bowtie2-build"):
    """

    :param input_file:
    :param output_db_name:
    :param alias:
    :return:
    """
    input_file = check_filename(input_file, Truefile=False)
    if output_db_name == "None":
        output_db_name = input_file.split(".fasta")[0].split(".fna")[0]
    else:
        output_db_name = check_filename(output_db_name, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("bowtie2-build", alias); tnum = 1; gmem = 2
    cmds = "/usr/local/bin/bowtie2-build %s %s" % (input_file, output_db_name)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def tophat2_runner(output_dir, bowtie_reference, fastq, alias="tophat2"):
    """
    tophat2 runner
    :param output_dir: output directory
    :param bowtie_reference: bowtie reference path
    :param fastq: fastq path
    :param alias: job name (tophat2)
    :return: job id
    """
    output_dir = check_dirname(output_dir, Truedir=False)
    bowtie_reference = check_filename(bowtie_reference, Truefile=False)
    fastq = check_filename(fastq)

    cmdfile = pbs_jobs.assign_cmdfile_path("tophat2", alias); tnum = 1; gmem = 2
    cmds = "/sternadi/home/volume1/taliakustin/software/tophat-2.1.1.Linux_x86_64/tophat2"\
           + " -o %s %s %s" % (output_dir, bowtie_reference, fastq)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds, load_python=False)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id



def r4s_runner(tree_file, seq_file, outfile, dirname, tree_outfile=None, unormelized_outfile=None, log_outfile=None, \
               ref_seq = None, n_categories = 4, alias = "r4s"):
    """
    run r4site on cluster
    :param tree_file: input tree file path
    :param seq_file: input sequence file path
    :param outfile: outfile path
    :param dirname: dirname for ouput files
    :param tree_outfile: output tree file path (default: None)
    :param unormelized_outfile: unormelized rated output file (default: None)
    :param log_outfile: output log file (default: None)
    :param alias: job name (default: r4s)
    :return: job id
    """
    tree_file = check_filename(tree_file)
    seq_file = check_filename(seq_file)
    dirname = check_dirname(dirname)

    if tree_outfile != None:
        tree_outfile = check_filename(tree_outfile, Truefile=False)
    else:
        tree_outfile = dirname + "/" + "out-tree"
    if unormelized_outfile != None:
        unormelized_outfile = check_filename(unormelized_outfile, Truefile=False)
    else:
        unormelized_outfile = dirname + "/out-unormelized"
    if log_outfile != None:
        log_outfile = check_filename(log_outfile, Truefile=False)
    else:
        log_outfile = dirname + "/out-log"


    cmdfile = pbs_jobs.assign_cmdfile_path("r4s_cmd.txt", alias); tnum = 1; gmem = 2
    ref_seq_parameter = " -a " + ref_seq if ref_seq is not None else ""
    if tree_file !=None:
        cmds = "/sternadi/home/volume1/shared/tools/rate4site"\
                                                            + " -t " + tree_file\
                                                            + " -s " + seq_file\
                                                            + " -o " + outfile\
                                                            + ref_seq_parameter \
                                                            + " -x " + tree_outfile\
                                                            + " -y " + unormelized_outfile\
                                                            + " -V 10"\
                                                            + " -l " + log_outfile\
                                                            + " -Mh -k " + n_categories
    else:
        cmds = "/sternadi/home/volume1/shared/tools/rate4site"\
                                                            + " -s " + seq_file\
                                                            + " -o " + outfile \
                                                            + ref_seq_parameter\
                                                            + " -x " + tree_outfile\
                                                            + " -y " + unormelized_outfile\
                                                            + " -V 10"\
                                                            + " -l " + log_outfile\
                                                            + " -Mh -k " + n_categories

    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds, queue="adistzachi")# added change
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def codeml_united_runner(clt1, clt2, rst1_name, rst2_name, alias ="cml"):
    """
    run codeml program from PAML on cluster - (runs both alternative and null model in one job).
    :param ctl1: ctl file path for null model
    :param ctl2: ctl file path for alternative model
    :param rst1_name: result file name for null model
    :param rst2_name: result file alternative for null model
    :param alias: job name (default: bml)
    :return: job id
    """

    clt1 = file_utilities.check_filename(clt1)
    clt2 = file_utilities.check_filename(clt2)
    base = os.path.split(clt1)[0]
    cmdfile = pbs_jobs.assign_cmdfile_path("codeml.txt", alias); tnum = 1; gmem = 2
    cmds = "cd %s\n" \
           "/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s\n" \
           "mv rst %s\n" \
           "/sternadi/home/volume1/taliakustin/software/paml4.8/bin/codeml %s\n" \
           "mv rst %s\n" % (base, clt1, rst1_name, clt2, rst2_name)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def selecton_runner(codon_aln, output_dir=None, tree=None, log=None, rate=None, output=None,
                    color=None, out_tree=None, query_seq = None, model="M8", alias="selecton", use_query_seq=False):
    codon_aln = check_filename(codon_aln)
    if output_dir == None:
        base = codon_aln.split(".")[0] + "_selecton"
    else:
        base = check_dirname(output_dir)
        base = base + codon_aln.split("/")[-1].split(".")[0] + "_selecton"
    log = set_filenames_for_pbs_runs(log, base, "log.txt")
    rate = set_filenames_for_pbs_runs(rate, base, "kaks.txt")
    output = set_filenames_for_pbs_runs(output, base, "output.txt")
    color = set_filenames_for_pbs_runs(color, base, "color.txt")
    out_tree = set_filenames_for_pbs_runs(out_tree, base, "output_tree.txt")

    if query_seq == None:
        query_seq = get_longest_sequence_name_in_fasta(codon_aln)

    if model == "M8":
        model = ""
    elif model == "M8a":
        model = "-w1 -Fw"
    elif model == "M7":
        model = "-p1 -Fp"

    if tree != None:
        tree = check_filename(tree)
        if use_query_seq == False:
            cmds = "selecton -i %s -u %s -l %s -r %s -o %s -c %s -t %s %s" \
                   % (codon_aln, tree, log, rate, output, color, out_tree, model)
        else:
            cmds = "selecton -i %s -u %s -l %s -r %s -o %s -c %s -t %s %s -q %s" \
                   % (codon_aln, tree, log, rate, output, color, out_tree, model, query_seq)
    else:
        if use_query_seq == False:
            cmds = "selecton -i %s -l %s -r %s -o %s -c %s -t %s %s" \
                   % (codon_aln, log, rate, output, color, out_tree, model)
        else:
            cmds = "selecton -i %s -l %s -r %s -o %s -c %s -t %s %s -q %s" \
                   % (codon_aln, log, rate, output, color, out_tree, model, query_seq)
    cmdfile = pbs_jobs.assign_cmdfile_path("selecton.txt", alias); tnum = 1; gmem = 2
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def pipeline_runner(input_dir, output_dir, ref_file, NGS_or_Cirseq, TYPE_OF_INPUT_FILE=None, start=None, end=None, gaps=None,
                    qscore=None, blast=None, rep=None, t=None, alias="pipeline"):
    input_dir = check_dirname(input_dir)
    output_dir = check_dirname(output_dir)
    ref_file = check_filename(ref_file)
    if NGS_or_Cirseq not in [1, 2]:
        raise Exception("NGS_or_Cirseq has to be 1 or 2")
    cmds = "python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py -i %s -o %s -r %s -NGS_or_Cirseq %i" \
           % (input_dir, output_dir, ref_file, NGS_or_Cirseq)
    if TYPE_OF_INPUT_FILE != None:
        cmds += " -t %s" % TYPE_OF_INPUT_FILE
    if start != None:
        cmds += " -s %i" % start
    if end != None:
        cmds += " -e %i" % end
    if gaps != None:
        cmds += " -g %s" % gaps
    if qscore != None:
        cmds += " -q %i" % qscore
    if blast != None:
        cmds += " -b %i" % blast
    if rep != None:
        cmds += " -rep %i" % int(rep)
    if t != None:
        cmds += " -t %s" % t


    print(cmds)
    cmdfile = pbs_jobs.assign_cmdfile_path("pipeline.txt", alias); tnum = 1; gmem = 2;
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds, load_python=True)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def fits_runner(inference_type, dataset_file, param_file,alias='FITS', posterior_file=None, summary_file=None, batch=None, queue="adistzachi"):
    """
    run fits currect version on cluster
    :param inference_type: the type of inference - fitness = 0, mutation rate =1, population size=2, simulate=3,
    :param dataset_file: dataset file. if batch != None should indicate $PBS_ARRAY_INDEX
    :param param_file: parameter file
    :param alias: job job_name. default is FITS.
    :param posterior_file: output posterior file. should be provided to all types except of simulate. if batch != None should indicate $PBS_ARRAY_INDEX
    :param summary_file: output summary file. should be provided to all types except of simulate. if batch != None should indicate $PBS_ARRAY_INDEX
    :param batch: the number of jobs in the array, if None run as a single job
    :return: sumbit a job\ job array to the cluster
    """

    if batch == None:
        dataset_file = check_filename(dataset_file)
    param_file = check_filename(param_file)

    if inference_type not in [0,1,2,3]:
        raise Exception('Inference type should be 0 (fitness), 1 (mutation rate), 2 (population size), or 3 (simulate)')


    if inference_type == 0: # fitness inference
        if posterior_file == None:
            posterior_file = os.path.join(os.path.dirname(dataset_file), 'posterior.txt')
        if summary_file == None:
            summary_file = os.path.join(os.path.dirname(dataset_file), 'summary.txt')
        cmds = 'module load gcc/gcc-7.3.0\n' +\
        '/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits_current_version -fitness ' \
        '{} {} {} {}'.format(param_file, dataset_file, posterior_file,summary_file)

        if batch == None:
            script_runner(cmds, alias, queue=queue)
        else:
            array_script_runner(cmds,batch,alias, queue=queue)

    elif inference_type == 1: # mutation rate inference
        if posterior_file == None:
            posterior_file = os.path.join(os.path.dirname(dataset_file), 'posterior.txt')
        if summary_file == None:
            summary_file = os.path.join(os.path.dirname(dataset_file), 'summary.txt')
        cmds = 'module load gcc/gcc-7.3.0\n' +\
        '/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits_current_version -mutation ' \
        '{} {} {} {}'.format(param_file, dataset_file, posterior_file,summary_file)

        if batch == None:
            script_runner(cmds, alias, queue=queue)
        else:
            array_script_runner(cmds,batch,alias, queue=queue)

    elif inference_type == 2:  # population size inference
        if posterior_file == None:
            posterior_file = os.path.join(os.path.dirname(dataset_file), 'posterior.txt')
        if summary_file == None:
            summary_file = os.path.join(os.path.dirname(dataset_file), 'summary.txt')
        cmds = 'module load gcc/gcc-7.3.0\n' + \
               '/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits_current_version -popsize ' \
               '{} {} {} {}'.format(param_file, dataset_file, posterior_file, summary_file)

        if batch == None:
            script_runner(cmds, alias, queue=queue)
        else:
            array_script_runner(cmds, batch, alias, queue=queue)

    else:  # simulations

        cmds = 'module load gcc/gcc-7.3.0\n' + \
               '/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits_current_version -simulate ' \
               '{} {}'.format(param_file, dataset_file)

        if batch == None:
            script_runner(cmds, alias, queue=queue)
        else:
            array_script_runner(cmds, batch, alias, queue=queue)



def dirSel_runner(dirSel_params, dirSel_path="/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/directionalSelection/directionalSelection",
                  alias = "dirSel", out_dir=None, queue="adistzachi", cmdname="dirSel_cmd.txt"):
    """
    run directional selection
    :param dirSel_params: params file
    :param dirSel_path: path of program
    :param alias: job name (default: dirSel)
    :return: job_id
    """
    dirSel_params = check_filename(dirSel_params)
    dirSel_path = check_filename(dirSel_path)
    if out_dir != None:
        cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, out_dir)
    else:
        cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias)
    tnum = 1; gmem = 20
    cmd = "%s %s" % (dirSel_path, dirSel_params)
    cmds = "echo %s \n%s" %(cmd, cmd)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds, queue=queue)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def netMHCpan_runner(input, output, path = "/sternadi/home/volume1/taliakustin/software/netMHCpan-4.0/netMHCpan", alias="netMHCpan", peptide_len="", allele=""):
    input = check_filename(input)
    output = check_filename(output, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path("netMHCpan_cmd.txt", alias); tnum = 1; gmem = 2
    cmd = "%s %s -t 0.2" % (path, input)
    if allele != "":
        cmd += " -a '%s'" % allele
    if peptide_len != "":
        cmd += " -l %i" % peptide_len
    cmd += "> %s\n echo DONE" % output
    cmds = "echo %s \n%s" % (cmd, cmd)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def merge_runner(fastq_r1, fastq_r2, output, alias = "merge"):
    """
    run merge on cluster with repeats = 60.
    :param fastq_r1: R1 fastq path
    :param fastq_r2: R2 fastq path
    :param output: output file path, will be gz or fastq according to fastq_r1 type
    :param alias: job name (default: dirSel)
    :return: job_id
    """
    fastq_r1 = check_filename(fastq_r1)
    fastq_r2 = check_filename(fastq_r2)
    cmdfile = pbs_jobs.assign_cmdfile_path("merge", alias); tnum = 1; gmem = 2
    cmds = "python /sternadi/home/volume1/shared/SternLab/scripts/merge_fastq_files.py" \
           " -f %s -e %s -o %s -r 60"\
           % (fastq_r1, fastq_r2, output)
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def threeseq_runner(alignment, alias = "3seq", identifier = None, start=0, end=-1):
    """
    run 3seq on cluster
    :param alignment: alignment file path, fasta or phylip
    :param alias: job name (default: 3seq)
    :param identifier: default None. run id, if not provided the program will automatically generate an id
    :param start: start position for recombination detection, default 0
    :param end: end position for recombination detection, default -1
    :return: job id
    """
    alignment = check_filename(alignment)

    cmdfile = pbs_jobs.assign_cmdfile_path("3seq", alias); tnum = 1; gmem = 2
    if identifier != None:
        if end > 0:
            if start >= end:
                raise Exception ("end should be larger than start")
            else:
                cmds = "/sternadi/home/volume1/shared/tools/3seq_build_170612/3seq -f {} -d -id {}" \
                       "-f{} -l{} -p".format(alignment, identifier, start, end)
        else:
            cmds = "/sternadi/home/volume1/shared/tools/3seq_build_170612/3seq -f {} -d -id {} -p".format(alignment, identifier)
    else: # no identifier
        if end > 0:
            if start >= end:
                raise Exception("end should be larger than start")
            else:
                cmds = "/sternadi/home/volume1/shared/tools/3seq_build_170612/3seq -f {} -d -id-auto" \
                       "-f{} -l{} -p".format(alignment, start, end)
        else:
            cmds = "/sternadi/home/volume1/shared/tools/3seq_build_170612/3seq -f {} -d -id-auto -p".format(alignment)

    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def r4s_runner_aln(tree_file, seq_file, outfile, dirname, tree_outfile=None, unormelized_outfile=None, log_outfile=None, \
               ref_seq = None, alias = "r4s"):
    """
    run r4site on cluster
    :param tree_file: input tree file path
    :param seq_file: input sequence file path
    :param outfile: outfile path
    :param dirname: dirname for ouput files
    :param tree_outfile: output tree file path (default: None)
    :param unormelized_outfile: unormelized rated output file (default: None)
    :param log_outfile: output log file (default: None)
    :param alias: job name (default: r4s)
    :return: job id
    """
    tree_file = check_filename(tree_file)
    seq_file = check_filename(seq_file)
    dirname = check_dirname(dirname)

    if tree_outfile != None:
        tree_outfile = check_filename(tree_outfile, Truefile=False)
    else:
        tree_outfile = dirname + "/" + "out-tree"
    if unormelized_outfile != None:
        unormelized_outfile = check_filename(unormelized_outfile, Truefile=False)
    else:
        unormelized_outfile = dirname + "/out-unormelized"
    if log_outfile != None:
        log_outfile = check_filename(log_outfile, Truefile=False)
    else:
        log_outfile = dirname + "/out-log"


    cmdfile = pbs_jobs.assign_cmdfile_path("r4s_cmd.txt", alias); tnum = 1; gmem = 2
    ref_seq_parameter = " -a " + ref_seq if ref_seq is not None else ""
    if tree_file !=None:
        cmds = "/sternadi/home/volume1/daniellem1/software/rate4site.3.2.source/sourceMar09/rate4site"\
                                                            + " -t " + tree_file\
                                                            + " -s " + seq_file\
                                                            + " -o " + outfile\
                                                            + ref_seq_parameter \
                                                            + " -x " + tree_outfile\
                                                            + " -y " + unormelized_outfile\
                                                            + " -V 10"\
                                                            + " -l " + log_outfile\

    else:
        cmds = "/sternadi/home/volume1/daniellem1/software/rate4site.3.2.source/sourceMar09/rate4site"\
                                                            + " -s " + seq_file\
                                                            + " -o " + outfile \
                                                            + ref_seq_parameter\
                                                            + " -x " + tree_outfile\
                                                            + " -y " + unormelized_outfile\
                                                            + " -V 10"\
                                                            + " -l " + log_outfile\

    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=tnum, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id


def beast_runner(xml_config, logpath=None, xml_config_out=None, multi_thread=True, gpu=False, alias="beast", ncpus=8, ngpus=4, gmem=2):
    """
    run beast on cluster - TODO - params definitions 
    :param seqfile: sequence file path
    :param dbfile: db file
    :param outfile: output file path
    :param alias: job name (blast)
    :return: job id
    """
    xml_config = check_filename(xml_config)
    beast_loc = "/sternadi/home/volume1/shared/tools/BEASTv1.10.4"

    if logpath != None:
        configure_log_path(xml_config, logpath, xml_config_out)
        if xml_config_out:
            xml_config = xml_config_out

    cmdfile = pbs_jobs.assign_cmdfile_path("beast_cmd", alias)
    if gpu:
        cmds = f"module load beagle-lib\n\
                cd {beast_loc}\n\
                ./bin/beast -beagle -beagle_GPU -beagle_cuda {xml_config}"
        pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds, ngpus=ngpus)
    elif multi_thread:
        cmds = f"module load beagle-lib\n\
                cd {beast_loc}\n\
                ./bin/beast -beagle {xml_config}"
        pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds, ncpus=ncpus)
    else:
        cmds = f"module load beagle-lib\n\
                cd {beast_loc}\n\
                ./bin/beast {xml_config}"        
        pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def phydyn_runner(xml_config, multi_thread=True, gpu=False, alias="PhyDyn", ncpus=8, ngpus=4, gmem=2, queue='adis'):
    """
    run beast on cluster - TODO - params definitions 
    :param seqfile: sequence file path
    :param dbfile: db file
    :param outfile: output file path
    :param alias: job name (blast)
    :return: job id
    """
    xml_config = check_filename(xml_config)
    beast_loc = "/sternadi/home/volume1/daniellem1/software/beast"

    cmdfile = pbs_jobs.assign_cmdfile_path("beast_cmd", alias)
    if gpu:
        cmds = f"module load beagle-lib\n\
                cd {beast_loc}\n\
                ./bin/beast -beagle -beagle_GPU -beagle_cuda {xml_config}"
        pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds, ngpus=ngpus, queue=queue)
    elif multi_thread:
        cmds = f"module load beagle-lib\n\
                cd {beast_loc}\n\
                ./bin/beast -beagle {xml_config}"
        pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds, ncpus=ncpus, queue=queue)
    else:
        cmds = f"module load beagle-lib\n\
                cd {beast_loc}\n\
                ./bin/beast {xml_config}"        
        pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, gmem=gmem, cmds=cmds, queue=queue)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id



def kraken2_runner(db_location, input_file, output_file, alias='kraken2', nthreads=None, ncpus=2, queue='dudulight'):
    """
    run kraken2 on cluster -
    :param db_location: db file path location
    :param input_file: input file - can be .fasta, .fastq or comprassed form
    :param output_file: output file path
    :param alias: job name (kraken2)
    :return: job id
    """
    input_file = check_filename(input_file)
    comprassed = ''
    threads = ''

    if '.fastq.gz':
        comprassed = '--gzip-compressed'
    if nthreads:
        thread = f'--threads {nthreads}'
    cmdfile = pbs_jobs.assign_cmdfile_path("kraken_cmd", alias)
    cmds = f"/davidb/local/software/kraken2/kraken2-2.0.8-beta/kraken2 --db {db_location} {threads} --use-names\
     {comprassed} {input_file} --output {output_file} --report {output_file.split('.')[0] + '_report.tab'}"
    pbs_jobs.create_pbs_cmd(cmdfile=cmdfile, alias=alias, ncpus=ncpus, cmds=cmds, queue=queue)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def grep_runner(input, text, output, alias = "grep",queue="adistzachi", cmdname = "grep"):
    input = check_filename(input)
    output = check_filename(output, Truefile=False)
    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias); tnum = 1; gmem = 1;
    cmds = f"grep {text} {input} > {output}"
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def make_blastDB(input, dbtype, alias="makeBlastDB", queue="adistzachi", cmdname = "makeBlastDB"):
    input = check_filename(input)
    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias);
    tnum = 1;
    gmem = 1;
    cmds = f"makeblastdb -in {input} -dbtype {dbtype}"
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)
    return job_id

def associvar_runner(input_dir, output_dir, start_position, end_position, job_num=100, alias="associvar", queue="adistzachi", cmdname = "associvar"):
    '''
    @input_dir: directory with blast results file/files (compatible with old pipeline)
    @out_dir: directory to write outputs into
    @start_position: start position to test associations, also requires all reads to span this
    @end_position: end position to test associations, also requires all reads to span this
    @job_num: how many jobs to split the association tests into (default 100)
    '''
    input_dir = check_dirname(input_dir)
    output_dir = check_dirname(output_dir, Truedir = False)
    make_dir(output_dir)
    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias);
    gmem = 5;
    cmds = f"python /sternadi/home/volume1/shared/tools/AssociVar/association_tests/parse_blasts.py -i {input_dir} -o {output_dir}\n"
    cmds += f"python /sternadi/home/volume1/shared/tools/AssociVar/association_tests/create_couples_index_file.py -s {start_position} -e {end_position} -j {job_num} -o {output_dir}/couples_index_file.csv"
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds)
    job_id = pbs_jobs.submit(cmdfile)

    # array job, waits for first job to finish
    associations_dir = check_dirname(f'{output_dir}/associations', Truedir=False)
    make_dir(associations_dir)
    cmds = f"python /sternadi/home/volume1/shared/tools/AssociVar/association_tests/association_test.py -b {output_dir}/blasts.csv -m {output_dir}/mutation.csv -i {output_dir}/couples_index_file.csv -o {associations_dir} -j $PBS_ARRAY_INDEX"
    array_job_id = array_script_runner(cmds, job_num, alias+'associations', queue=queue, gmem=gmem, toRun=True, run_after_job=job_id)
    # unify and normalize after array is done

    cmdfile = pbs_jobs.assign_cmdfile_path(cmdname, alias+'unify');
    gmem = 5;
    cmds = f"python /sternadi/home/volume1/shared/tools/AssociVar/association_tests/unify_association_results.py -i {associations_dir} -o {output_dir}/associations.csv\n"
    cmds += f"python /sternadi/home/volume1/shared/tools/AssociVar/association_tests/normalize_chi2.py -i {output_dir}/associations.csv -o {output_dir}/associations.ztest.csv\n"
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias+'unify', queue=queue, gmem=gmem, cmds=cmds, run_after_job=array_job_id)
    unify_job_id = pbs_jobs.submit(cmdfile)

    return (job_id, array_job_id, unify_job_id)