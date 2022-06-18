# Extract and call functional IGHV alleles
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG001/HG001_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG001/HG001_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG001/HG001_BCRV/read_depth_functional.rpt > target_call/HG001/HG001_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG002/HG002_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG002/HG002_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002/HG002_BCRV/read_depth_functional.rpt > target_call/HG002/HG002_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG003/HG003_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG003/HG003_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG003/HG003_BCRV/read_depth_functional.rpt > target_call/HG003/HG003_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG004/HG004_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG004/HG004_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG004/HG004_BCRV/read_depth_functional.rpt > target_call/HG004/HG004_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG005/HG005_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG005/HG005_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG005/HG005_BCRV/read_depth_functional.rpt > target_call/HG005/HG005_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG006/HG006_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG006/HG006_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG006/HG006_BCRV/read_depth_functional.rpt > target_call/HG006/HG006_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG007/HG007_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG007/HG007_BCRV/read_depth_functional.rpt
python3 scripts/calling_threshold.py -dp target_call/HG007/HG007_BCRV/read_depth_functional.rpt > target_call/HG007/HG007_BCRV/gAIRR-call_functional_report.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/PL-B/PL-B_BCRV/read_depth_calling_by_bwa.rpt > target_call/PL-B/PL-B_BCRV/read_depth_functional.IGH.rpt
python3 scripts/calling_threshold.py -dp target_call/PL-B/PL-B_BCRV/read_depth_functional.IGH.rpt > target_call/PL-B/PL-B_BCRV/gAIRR-call_functional.IGH.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/PL-S/PL-S_BCRV/read_depth_calling_by_bwa.rpt > target_call/PL-S/PL-S_BCRV/read_depth_functional.IGH.rpt
python3 scripts/calling_threshold.py -dp target_call/PL-S/PL-S_BCRV/read_depth_functional.IGH.rpt > target_call/PL-S/PL-S_BCRV/gAIRR-call_functional.IGH.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG002-250-60X/HG002-250-60X_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG002-250-60X/HG002-250-60X_BCRV/read_depth_functional.IGH.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002-250-60X/HG002-250-60X_BCRV/read_depth_functional.IGH.rpt > target_call/HG002-250-60X/HG002-250-60X_BCRV/gAIRR-call_functional.IGH.rpt
python3 scripts/extract_functional.py -l example/material/IGHV_functional.txt -r target_call/HG002-novaseq-30x/HG002-novaseq-30x_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG002-novaseq-30x/HG002-novaseq-30x_BCRV/read_depth_functional.IGH.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002-novaseq-30x/HG002-novaseq-30x_BCRV/read_depth_functional.IGH.rpt > target_call/HG002-novaseq-30x/HG002-novaseq-30x_BCRV/gAIRR-call_functional.IGH.rpt


# Extract and call functional TRV alleles
TRA_func="example/material/TRAV_functional.txt"
TRB_func="example/material/TRBV_functional.txt"
TRG_func="example/material/TRGV_functional.txt"
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG001/HG001_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG001/HG001_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG002/HG002_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG002/HG002_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG003/HG003_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG003/HG003_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG004/HG004_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG004/HG004_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG005/HG005_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG005/HG005_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG006/HG006_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG006/HG006_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG007/HG007_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG007/HG007_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/PL-B/PL-B_TCRV/read_depth_calling_by_bwa.rpt > target_call/PL-B/PL-B_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/PL-S/PL-S_TCRV/read_depth_calling_by_bwa.rpt > target_call/PL-S/PL-S_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG002-250-60X/HG002-250-60X_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG002-250-60X/HG002-250-60X_TCRV/read_depth_functional.TRA.rpt
python3 scripts/extract_functional.py -l ${TRA_func} -r target_call/HG002-novaseq-30x/HG002-novaseq-30x_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG002-novaseq-30x/HG002-novaseq-30x_TCRV/read_depth_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG001/HG001_TCRV/read_depth_functional.TRA.rpt > target_call/HG001/HG001_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002/HG002_TCRV/read_depth_functional.TRA.rpt > target_call/HG002/HG002_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG003/HG003_TCRV/read_depth_functional.TRA.rpt > target_call/HG003/HG003_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG004/HG004_TCRV/read_depth_functional.TRA.rpt > target_call/HG004/HG004_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG005/HG005_TCRV/read_depth_functional.TRA.rpt > target_call/HG005/HG005_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG006/HG006_TCRV/read_depth_functional.TRA.rpt > target_call/HG006/HG006_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG007/HG007_TCRV/read_depth_functional.TRA.rpt > target_call/HG007/HG007_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/PL-B/PL-B_TCRV/read_depth_functional.TRA.rpt > target_call/PL-B/PL-B_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/PL-S/PL-S_TCRV/read_depth_functional.TRA.rpt > target_call/PL-S/PL-S_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002-250-60X/HG002-250-60X_TCRV/read_depth_functional.TRA.rpt > target_call/HG002-250-60X/HG002-250-60X_TCRV/gAIRR-call_functional.TRA.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002-novaseq-30x/HG002-novaseq-30x_TCRV/read_depth_functional.TRA.rpt > target_call/HG002-novaseq-30x/HG002-novaseq-30x_TCRV/gAIRR-call_functional.TRA.rpt

python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG001/HG001_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG001/HG001_TCRV/read_depth_functional.TRB.rpt
python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG002/HG002_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG002/HG002_TCRV/read_depth_functional.TRB.rpt
python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG003/HG003_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG003/HG003_TCRV/read_depth_functional.TRB.rpt
python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG004/HG004_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG004/HG004_TCRV/read_depth_functional.TRB.rpt
python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG005/HG005_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG005/HG005_TCRV/read_depth_functional.TRB.rpt
python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG006/HG006_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG006/HG006_TCRV/read_depth_functional.TRB.rpt
python3 scripts/extract_functional.py -l ${TRB_func} -r target_call/HG007/HG007_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG007/HG007_TCRV/read_depth_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG001/HG001_TCRV/read_depth_functional.TRB.rpt > target_call/HG001/HG001_TCRV/gAIRR-call_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002/HG002_TCRV/read_depth_functional.TRB.rpt > target_call/HG002/HG002_TCRV/gAIRR-call_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG003/HG003_TCRV/read_depth_functional.TRB.rpt > target_call/HG003/HG003_TCRV/gAIRR-call_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG004/HG004_TCRV/read_depth_functional.TRB.rpt > target_call/HG004/HG004_TCRV/gAIRR-call_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG005/HG005_TCRV/read_depth_functional.TRB.rpt > target_call/HG005/HG005_TCRV/gAIRR-call_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG006/HG006_TCRV/read_depth_functional.TRB.rpt > target_call/HG006/HG006_TCRV/gAIRR-call_functional.TRB.rpt
python3 scripts/calling_threshold.py -dp target_call/HG007/HG007_TCRV/read_depth_functional.TRB.rpt > target_call/HG007/HG007_TCRV/gAIRR-call_functional.TRB.rpt

python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG001/HG001_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG001/HG001_TCRV/read_depth_functional.TRG.rpt
python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG002/HG002_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG002/HG002_TCRV/read_depth_functional.TRG.rpt
python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG003/HG003_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG003/HG003_TCRV/read_depth_functional.TRG.rpt
python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG004/HG004_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG004/HG004_TCRV/read_depth_functional.TRG.rpt
python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG005/HG005_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG005/HG005_TCRV/read_depth_functional.TRG.rpt
python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG006/HG006_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG006/HG006_TCRV/read_depth_functional.TRG.rpt
python3 scripts/extract_functional.py -l ${TRG_func} -r target_call/HG007/HG007_TCRV/read_depth_calling_by_bwa.rpt > target_call/HG007/HG007_TCRV/read_depth_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG001/HG001_TCRV/read_depth_functional.TRG.rpt > target_call/HG001/HG001_TCRV/gAIRR-call_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002/HG002_TCRV/read_depth_functional.TRG.rpt > target_call/HG002/HG002_TCRV/gAIRR-call_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG003/HG003_TCRV/read_depth_functional.TRG.rpt > target_call/HG003/HG003_TCRV/gAIRR-call_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG004/HG004_TCRV/read_depth_functional.TRG.rpt > target_call/HG004/HG004_TCRV/gAIRR-call_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG005/HG005_TCRV/read_depth_functional.TRG.rpt > target_call/HG005/HG005_TCRV/gAIRR-call_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG006/HG006_TCRV/read_depth_functional.TRG.rpt > target_call/HG006/HG006_TCRV/gAIRR-call_functional.TRG.rpt
python3 scripts/calling_threshold.py -dp target_call/HG007/HG007_TCRV/read_depth_functional.TRG.rpt > target_call/HG007/HG007_TCRV/gAIRR-call_functional.TRG.rpt

IGK_func="example/material/IGKV_functional.txt"
IGL_func="example/material/IGLV_functional.txt"
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG001/HG001_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG001/HG001_BCRV/read_depth_functional.IGK.rpt
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG002/HG002_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG002/HG002_BCRV/read_depth_functional.IGK.rpt
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG003/HG003_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG003/HG003_BCRV/read_depth_functional.IGK.rpt
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG004/HG004_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG004/HG004_BCRV/read_depth_functional.IGK.rpt
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG005/HG005_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG005/HG005_BCRV/read_depth_functional.IGK.rpt
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG006/HG006_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG006/HG006_BCRV/read_depth_functional.IGK.rpt
python3 scripts/extract_functional.py -l ${IGK_func} -r target_call/HG007/HG007_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG007/HG007_BCRV/read_depth_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG001/HG001_BCRV/read_depth_functional.IGK.rpt > target_call/HG001/HG001_BCRV/gAIRR-call_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002/HG002_BCRV/read_depth_functional.IGK.rpt > target_call/HG002/HG002_BCRV/gAIRR-call_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG003/HG003_BCRV/read_depth_functional.IGK.rpt > target_call/HG003/HG003_BCRV/gAIRR-call_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG004/HG004_BCRV/read_depth_functional.IGK.rpt > target_call/HG004/HG004_BCRV/gAIRR-call_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG005/HG005_BCRV/read_depth_functional.IGK.rpt > target_call/HG005/HG005_BCRV/gAIRR-call_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG006/HG006_BCRV/read_depth_functional.IGK.rpt > target_call/HG006/HG006_BCRV/gAIRR-call_functional.IGK.rpt
python3 scripts/calling_threshold.py -dp target_call/HG007/HG007_BCRV/read_depth_functional.IGK.rpt > target_call/HG007/HG007_BCRV/gAIRR-call_functional.IGK.rpt

python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG001/HG001_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG001/HG001_BCRV/read_depth_functional.IGL.rpt
python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG002/HG002_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG002/HG002_BCRV/read_depth_functional.IGL.rpt
python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG003/HG003_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG003/HG003_BCRV/read_depth_functional.IGL.rpt
python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG004/HG004_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG004/HG004_BCRV/read_depth_functional.IGL.rpt
python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG005/HG005_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG005/HG005_BCRV/read_depth_functional.IGL.rpt
python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG006/HG006_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG006/HG006_BCRV/read_depth_functional.IGL.rpt
python3 scripts/extract_functional.py -l ${IGL_func} -r target_call/HG007/HG007_BCRV/read_depth_calling_by_bwa.rpt > target_call/HG007/HG007_BCRV/read_depth_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG001/HG001_BCRV/read_depth_functional.IGL.rpt > target_call/HG001/HG001_BCRV/gAIRR-call_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG002/HG002_BCRV/read_depth_functional.IGL.rpt > target_call/HG002/HG002_BCRV/gAIRR-call_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG003/HG003_BCRV/read_depth_functional.IGL.rpt > target_call/HG003/HG003_BCRV/gAIRR-call_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG004/HG004_BCRV/read_depth_functional.IGL.rpt > target_call/HG004/HG004_BCRV/gAIRR-call_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG005/HG005_BCRV/read_depth_functional.IGL.rpt > target_call/HG005/HG005_BCRV/gAIRR-call_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG006/HG006_BCRV/read_depth_functional.IGL.rpt > target_call/HG006/HG006_BCRV/gAIRR-call_functional.IGL.rpt
python3 scripts/calling_threshold.py -dp target_call/HG007/HG007_BCRV/read_depth_functional.IGL.rpt > target_call/HG007/HG007_BCRV/gAIRR-call_functional.IGL.rpt

