IGHV_functional="example/material/IGHV_functional.txt"
IGKV_functional="example/material/IGKV_functional.txt"
IGLV_functional="example/material/IGLV_functional.txt"
TRAV_functional="example/material/TRAV_functional.txt"
TRBV_functional="example/material/TRBV_functional.txt"
TRGV_functional="example/material/TRGV_functional.txt"

# genotype of gAIRR-annotate
python3 scripts/collect_gene_from_bed.py -fl ${IGHV_functional} -b1 target_annotation/HG001-denovo/IGH_H1.bed -b2 target_annotation/HG001-denovo/IGH_H2.bed > RM_genotype/annotate_genotype_HG001-denovo.rpt
python3 scripts/collect_gene_from_bed.py -fl ${IGLV_functional} -b1 target_annotation/HG001-denovo/IGL_H1.bed -b2 target_annotation/HG001-denovo/IGL_H2.bed > RM_genotype/annotate_genotype_HG001-denovo.IGL.rpt
python3 scripts/collect_gene_from_bed.py -fl ${IGKV_functional} -b1 target_annotation/HG001-denovo/IGK_H1.bed -b2 target_annotation/HG001-denovo/IGK_H2.bed > RM_genotype/annotate_genotype_HG001-denovo.IGK.rpt
python3 scripts/collect_gene_from_bed.py -fl ${TRAV_functional} -b1 target_annotation/HG001-denovo/TRA_H1.bed -b2 target_annotation/HG001-denovo/TRA_H2.bed > RM_genotype/annotate_genotype_HG001-denovo.TRA.rpt
python3 scripts/collect_gene_from_bed.py -fl ${TRBV_functional} -b1 target_annotation/HG001-denovo/TRB_H1.bed -b2 target_annotation/HG001-denovo/TRB_H2.bed > RM_genotype/annotate_genotype_HG001-denovo.TRB.rpt
python3 scripts/collect_gene_from_bed.py -fl ${TRGV_functional} -b1 target_annotation/HG001-denovo/TRG_H1.bed -b2 target_annotation/HG001-denovo/TRG_H2.bed > RM_genotype/annotate_genotype_HG001-denovo.TRG.rpt
# HG002 IGH
python3 scripts/collect_gene_from_bed.py -fl ${IGHV_functional} -b1 target_annotation/HG002/IGH_HG002_H1.bed -b2 target_annotation/HG002/IGH_HG002_H2.bed   > RM_genotype/annotate_genotype_HG002.rpt
python3 scripts/collect_gene_from_bed.py -fl ${IGHV_functional} -b1 target_annotation/HG002-denovo/IGH_HG002_H1.bed -b2 target_annotation/HG002-denovo/IGH_HG002_H2.bed > RM_genotype/annotate_genotype_HG002-denovo.rpt
# HG002 TRA
python3 scripts/collect_gene_from_bed.py -fl ${TRAV_functional} -b1 target_annotation/HG002/TRA_H1.bed        -b2 target_annotation/HG002/TRA_H2.bed        > RM_genotype/annotate_genotype_HG002.TRA.rpt
python3 scripts/collect_gene_from_bed.py -fl ${TRAV_functional} -b1 target_annotation/HG002-denovo/TRA_H1.bed -b2 target_annotation/HG002-denovo/TRA_H2.bed > RM_genotype/annotate_genotype_HG002-denovo.TRA.rpt
# HG002 TRB
python3 scripts/collect_gene_from_bed.py -fl ${TRBV_functional} -b1 target_annotation/HG002/TRB_H1.bed        -b2 target_annotation/HG002/TRB_H2.bed        > RM_genotype/annotate_genotype_HG002.TRB.rpt
python3 scripts/collect_gene_from_bed.py -fl ${TRBV_functional} -b1 target_annotation/HG002-denovo/TRB_H1.bed -b2 target_annotation/HG002-denovo/TRB_H2.bed > RM_genotype/annotate_genotype_HG002-denovo.TRB.rpt
# HG002 TRG
python3 scripts/collect_gene_from_bed.py -fl ${TRGV_functional} -b1 target_annotation/HG002/TRG_H1.bed        -b2 target_annotation/HG002/TRG_H2.bed        > RM_genotype/annotate_genotype_HG002.TRG.rpt
python3 scripts/collect_gene_from_bed.py -fl ${TRGV_functional} -b1 target_annotation/HG002-denovo/TRG_H1.bed -b2 target_annotation/HG002-denovo/TRG_H2.bed > RM_genotype/annotate_genotype_HG002-denovo.TRG.rpt
# HG002 IGL
python3 scripts/collect_gene_from_bed.py -fl ${IGLV_functional} -b1 target_annotation/HG002/IGL_HG002_H1.bed        -b2 target_annotation/HG002/IGL_HG002_H2.bed        > RM_genotype/annotate_genotype_HG002.IGL.rpt
python3 scripts/collect_gene_from_bed.py -fl ${IGLV_functional} -b1 target_annotation/HG002-denovo/IGL_HG002_H1.bed -b2 target_annotation/HG002-denovo/IGL_HG002_H2.bed > RM_genotype/annotate_genotype_HG002-denovo.IGL.rpt
# HG002 IGK
python3 scripts/collect_gene_from_bed.py -fl ${IGKV_functional} -b1 target_annotation/HG002/IGK_HG002_H1.bed        -b2 target_annotation/HG002/IGK_HG002_H2.bed        > RM_genotype/annotate_genotype_HG002.IGK.rpt
python3 scripts/collect_gene_from_bed.py -fl ${IGKV_functional} -b1 target_annotation/HG002-denovo/IGK_HG002_H1.bed -b2 target_annotation/HG002-denovo/IGK_HG002_H2.bed > RM_genotype/annotate_genotype_HG002-denovo.IGK.rpt

# genotype of gAIRR-call
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG001/HG001_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG001.rpt 
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG002/HG002_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG002.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG003/HG003_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG003.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG004/HG004_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG004.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG005/HG005_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG005.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG006/HG006_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG006.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG007/HG007_BCRV/gAIRR-call_functional_report.rpt > RM_genotype/call_genotype_HG007.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/PL-B/PL-B_BCRV/gAIRR-call_functional.IGH.rpt      > RM_genotype/call_genotype_PL-B.IGH.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/PL-S/PL-S_BCRV/gAIRR-call_functional.IGH.rpt      > RM_genotype/call_genotype_PL-S.IGH.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG002-250-60X/HG002-250-60X_BCRV/gAIRR-call_functional.IGH.rpt > RM_genotype/call_genotype_HG002-250-60X.IGH.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGHV_functional} -rd target_call/HG002-novaseq-30x/HG002-novaseq-30x_BCRV/gAIRR-call_functional.IGH.rpt > RM_genotype/call_genotype_HG002-novaseq-30x.IGH.rpt

python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG001/HG001_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG001.TRA.rpt 
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG002/HG002_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG002.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG003/HG003_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG003.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG004/HG004_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG004.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG005/HG005_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG005.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG006/HG006_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG006.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG007/HG007_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG007.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/PL-B/PL-B_TCRV/gAIRR-call_functional.TRA.rpt   > RM_genotype/call_genotype_PL-B.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/PL-S/PL-S_TCRV/gAIRR-call_functional.TRA.rpt   > RM_genotype/call_genotype_PL-S.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG002-250-60X/HG002-250-60X_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG002-250-60X.TRA.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRAV_functional} -rd target_call/HG002-novaseq-30x/HG002-novaseq-30x_TCRV/gAIRR-call_functional.TRA.rpt > RM_genotype/call_genotype_HG002-novaseq-30x.TRA.rpt

python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG001/HG001_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG001.TRB.rpt 
python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG002/HG002_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG002.TRB.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG003/HG003_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG003.TRB.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG004/HG004_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG004.TRB.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG005/HG005_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG005.TRB.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG006/HG006_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG006.TRB.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRBV_functional} -rd target_call/HG007/HG007_TCRV/gAIRR-call_functional.TRB.rpt > RM_genotype/call_genotype_HG007.TRB.rpt

python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG001/HG001_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG001.TRG.rpt 
python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG002/HG002_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG002.TRG.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG003/HG003_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG003.TRG.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG004/HG004_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG004.TRG.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG005/HG005_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG005.TRG.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG006/HG006_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG006.TRG.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${TRGV_functional} -rd target_call/HG007/HG007_TCRV/gAIRR-call_functional.TRG.rpt > RM_genotype/call_genotype_HG007.TRG.rpt

python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG001/HG001_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG001.IGK.rpt 
python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG002/HG002_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG002.IGK.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG003/HG003_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG003.IGK.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG004/HG004_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG004.IGK.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG005/HG005_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG005.IGK.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG006/HG006_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG006.IGK.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGKV_functional} -rd target_call/HG007/HG007_BCRV/gAIRR-call_functional.IGK.rpt > RM_genotype/call_genotype_HG007.IGK.rpt

python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG001/HG001_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG001.IGL.rpt 
python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG002/HG002_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG002.IGL.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG003/HG003_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG003.IGL.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG004/HG004_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG004.IGL.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG005/HG005_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG005.IGL.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG006/HG006_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG006.IGL.rpt
python3 scripts/collect_gene_from_read_depth.py -fl ${IGLV_functional} -rd target_call/HG007/HG007_BCRV/gAIRR-call_functional.IGL.rpt > RM_genotype/call_genotype_HG007.IGL.rpt
