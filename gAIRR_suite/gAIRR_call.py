# Wrap up python file for gAIRR-call
import subprocess
import sys
import os
import argparse
from shutil import which




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




def main():
    parser = argparse.ArgumentParser(description="gAIRR-call pipeline of gAIRR-suite.")
    parser.add_argument('-wd', '--work_dir', help="Path to output directory ['target_call/'].", default="target_call/")
    parser.add_argument('-lc', '--locus', help="Target Locus [TRV TRJ TRD IGV IGJ IGD]", nargs='+', default=['TRV', 'TRJ', 'TRD', 'IGV', 'IGJ', 'IGD'])
    parser.add_argument('-id', '--sample_id', help="Sample ID ['sample']")
    parser.add_argument('-rd1', '--read1', help="Path to gAIRR-seq. Pair-end 1", required=True)
    parser.add_argument('-rd2', '--read2', help="Path to gAIRR-seq. Pair-end 2", required=True)
    parser.add_argument('--flanking', help="Option to do flanking sequence analysis.", action='store_true')
    parser.add_argument('--keep', help="Option specify to keep to sam files", action='store_true')
    parser.add_argument('-t', '--thread', help="Number of threads to use [max].", type=int)
    args = parser.parse_args()

    workspace     = args.work_dir
    target_locus  = args.locus
    person_name   = args.sample_id
    path_read1    = args.read1
    path_read2    = args.read2
    flag_flanking = args.flanking
    flag_keeping  = args.keep
    thread        = args.thread
    if thread == None:
        thread = get_max_thread()
    path_module   = os.path.dirname(__file__) + '/'
    path_material = path_module + '../example/material/'

    
    dict_locus_map = {'TRV':'TCRV', 'TRJ':'TCRJ', 'TRD':'TCRD_plusHep', 'IGV':'BCRV', 'IGJ':'BCRJ', 'IGD':'BCRD_plusHep'}
    for locus_name in target_locus:
        if dict_locus_map.get(locus_name) == None:
            print('Only TRV, TRJ, TRD, IGV, IGJ, IGD are allowed for target locus call.')
            parser.print_usage()
            exit(2)
    list_allele_names = [dict_locus_map[locus_name] for locus_name in target_locus]
    

    subprocess.call("mkdir -p " + workspace, shell=True)
    subprocess.call("mkdir -p " + workspace+'/'+person_name, shell=True)
    # call novel alleles
    for allele_name in list_allele_names:
        allele_path = path_material+allele_name+"_alleles_parsed.fasta"
        print("[gAIRR-call]", person_name, allele_name, "novel allele calling...")
        command = ' '.join(["bash", path_module+"novel_allele.sh", workspace+'/'+person_name, allele_name, allele_path, person_name, path_read1, path_read2, str(thread)])
        subprocess.call(command, shell=True)
        if flag_keeping == False:
            subprocess.call("rm " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'_novel/bwa_read_to_allele.sam' , shell=True)
        
        allele_path = workspace+'/'+person_name+'/'+person_name+'_'+allele_name+"_novel/"+allele_name+"_with_novel.fasta"
        print("[gAIRR-call]", person_name, allele_name, "allele calling...")
        command = ' '.join(["bash", path_module+"allele_calling.sh", workspace+'/'+person_name, allele_name, allele_path, person_name, path_read1, path_read2, str(thread)])
        subprocess.call(command, shell=True)
        if flag_keeping == False:
            subprocess.call("rm " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'/bwa_read_to_allele_all.sam' , shell=True)
        
        if flag_flanking:
            print("[gAIRR-call]", person_name, allele_name, "flanking sequence calling...")
            command = ' '.join(["bash", path_module+"flanking_sequence.sh", workspace+'/'+person_name, allele_name, allele_path, person_name, path_read1, path_read2, 'spades.py', str(thread)])
            subprocess.call(command, shell=True)
            if flag_keeping == False:
                subprocess.call("rm " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'_flanking/asm_check/*.sam' , shell=True)
                subprocess.call("rm " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'_flanking/haplotype_sam/bwa_reads_to_flanking.sam' , shell=True)
                subprocess.call("rm -rf " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'_flanking/asm_contigs/*/' , shell=True)
                subprocess.call("rm " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'_flanking/asm_contigs/*.fasta.*' , shell=True)
                subprocess.call("rm " + workspace+'/'+person_name+'/'+person_name+'_'+allele_name+'_flanking/group_allele_reads/*.fasta' , shell=True)

    




if __name__ == "__main__":
    main()


