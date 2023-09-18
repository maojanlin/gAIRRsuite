# Wrap up python file for gAIRR-annotate
import subprocess
import sys
import os
import io
from contextlib import redirect_stdout
import argparse
from shutil import which

# make sure the package modules is in the path
sys.path.append(os.path.dirname(__file__))

import scripts.annotation_with_asm
import scripts.analysis.bed_generator
import scripts.analysis.collect_gene_from_bed



def get_index_name(path_name):
    dir_idx  = path_name.rfind('/')
    dir_name = path_name[:dir_idx] + '/bwa/'
    idx_name = path_name[:dir_idx] + '/bwa' + path_name[dir_idx:]
    return dir_name, idx_name


def get_max_thread():
    if sys.platform == "darwin":
        result = subprocess.run(["sysctl -n hw.ncpu"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    else:
        result = subprocess.run(["nproc"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    return int(result.stdout.strip())


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


def check_program_install(list_names):
    flag_violate = False
    for name in list_names:
        if is_tool(name) == False:
            print(name, "is a prerequisite program, please install it before running biastools")
            flag_violate = True
    if flag_violate:
        exit(2)




def main():
    parser = argparse.ArgumentParser(description="gAIRR-annotate pipeline of gAIRR-suite.")
    parser.add_argument('-wd', '--work_dir', help="Path to output directory ['target_annotation/'].", default="target_annotation/")
    parser.add_argument('-lc', '--locus', help="Target Locus [TR IG]", nargs='+', default=['TR', 'IG'])
    parser.add_argument('-id', '--sample_id', help="Sample ID ['sample']")
    parser.add_argument('-a1', '--assembly_1', help="Path to assembly haplotype 1.", required=True)
    parser.add_argument('-a2', '--assembly_2', help="Path to assembly haplotype 2.")
    parser.add_argument('-t', '--thread', help="Number of threads to use [max].", type=int)
    args = parser.parse_args()

    outer_dir     = args.work_dir
    target_locus  = args.locus
    person_name   = args.sample_id
    path_asm_1    = args.assembly_1
    path_asm_2    = args.assembly_2
    thread        = args.thread
    if thread == None:
        thread = get_max_thread()
    
    path_module   = os.path.dirname(__file__) + '/scripts/'
    path_material = os.path.dirname(__file__) + '/material/'
    
    check_program_install(["bwa"])

    ########## STEP 1: indexing the personal assemblies ##########
    print("[gAIRR-annotate] Indexing the personal assembly", path_asm_1)
    path_idx_dir, path_asm_1_idx = get_index_name(path_asm_1)
    subprocess.call("mkdir -p " + path_idx_dir, shell=True)
    if os.path.isfile(path_asm_1_idx + '.bwt') == False:
        command = "bwa index " + path_asm_1 + " -p " + path_asm_1_idx
        subprocess.call(command, shell=True)
    if path_asm_2 != None: # if path2 exist
        print("[gAIRR-annotate] Indexing the personal assembly", path_asm_2)
        path_idx_dir, path_asm_2_idx = get_index_name(path_asm_2)
        subprocess.call("mkdir -p " + path_idx_dir, shell=True)
        if os.path.isfile(path_asm_2_idx + '.bwt') == False:
            command = "bwa index " + path_asm_2 + " -p " + path_asm_2_idx
            subprocess.call(command, shell=True)


    ########## STEP 2: take the locus information and annotate ##########
    list_allele_names = []
    if "TR" in target_locus:
        list_allele_names += ["TCRV", "TCRJ", "TCRD_plusHep"]
    if "IG" in target_locus:
        list_allele_names += ["BCRV", "BCRJ", "BCRD_plusHep"]
    list_allele_names = [(allele_name, path_material + allele_name + '_alleles_parsed.fasta') for allele_name in list_allele_names]
    
    # make working directory
    subprocess.call("mkdir -p " + outer_dir, shell=True)
    subprocess.call("mkdir -p " + outer_dir+'/'+person_name, shell=True)
    for allele_name, path_allele in list_allele_names:
        if 'plusHep' in allele_name: # D genes, apply special bwa parameters
            command_prefix = "bwa mem -a -T 10 -t"
        else:
            command_prefix = "bwa mem -a -t"
        command = (command_prefix, str(thread), path_asm_1_idx, path_allele, '>', outer_dir+'/'+person_name+'/'+'bwa_'+person_name+'_'+allele_name+'_alleles_to_asm_1.sam')
        subprocess.call(' '.join(command), shell=True)
        if path_asm_2 != None:
            command = (command_prefix, str(thread), path_asm_2_idx, path_allele, '>', outer_dir+'/'+person_name+'/'+'bwa_'+person_name+'_'+allele_name+'_alleles_to_asm_2.sam')
            subprocess.call(' '.join(command), shell=True)

        print("[gAIRR-annotate] Parse the ${allele_name} alleles sam files to annotation report...")
        command = ['-foa',  outer_dir+'/'+person_name+'/annotation_'+person_name+'_'+allele_name+'.txt', \
                   '-foma', outer_dir+'/'+person_name+'/annotation_imperfect_'+person_name+'_'+allele_name+'.txt', \
                   '-fom',  outer_dir+'/'+person_name+'/novel_'+person_name+'_'+allele_name+'.fasta', \
                   '-fof',  outer_dir+'/'+person_name+'/flanking_'+person_name+'_'+allele_name+'.fasta', \
                   '-fos',  outer_dir+'/'+person_name+'/summary_'+person_name+'_'+allele_name+'.rpt', \
                   '-ext',  '200', \
                   '-fs1',  outer_dir+'/'+person_name+'/bwa_'+person_name+'_'+allele_name+'_alleles_to_asm_1.sam', \
                   '-fp1',  outer_dir+'/'+person_name+'/dict_'+person_name+'_asm_1.pickle', \
                   '-fasm1',path_asm_1]
        if path_asm_2 != None:
            command += [ \
                   '-fs2',  outer_dir+'/'+person_name+'/bwa_'+person_name+'_'+allele_name+'_alleles_to_asm_2.sam', \
                   '-fp2',  outer_dir+'/'+person_name+'/dict_'+person_name+'_asm_2.pickle', \
                   '-fasm2',path_asm_2]
        scripts.annotation_with_asm.main(command)
    
    ########## STEP 3: collecting the annotation results ##########
    path_report = [allele_info[0] for allele_info in list_allele_names]
    path_report = [outer_dir+'/'+person_name+'/annotation_imperfect_'+person_name+'_'+allele_name+'.txt' for allele_name in path_report]
    command = ['-rl'] + path_report + ['-o', outer_dir+'/'+person_name+'/group_genes']
    scripts.analysis.bed_generator.main(command)

    command = ['-b1', outer_dir+'/'+person_name+'/group_genes.1.bed', \
               '-b2', outer_dir+'/'+person_name+'/group_genes.2.bed', \
               '-fl', path_material+'/IGH_functional.txt']
    f = io.StringIO()
    with redirect_stdout(f):
        scripts.analysis.collect_gene_from_bed.main(command)
    out = f.getvalue()
    f = open(outer_dir+'/'+person_name+'/IGH_functional.rpt', 'w')
    f.write(out)

    print("[gAIRR-annotation]", person_name, "annotation Finished!")




if __name__ == "__main__":
    main()


