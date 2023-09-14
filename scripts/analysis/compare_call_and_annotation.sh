IGHV_func="example/material/IGHV_functional.txt"
IGKV_func="example/material/IGKV_functional.txt"
IGLV_func="example/material/IGLV_functional.txt"
TRAV_func="example/material/TRAV_functional.txt"
TRBV_func="example/material/TRBV_functional.txt"
TRGV_func="example/material/TRGV_functional.txt"

collect="scripts/analysis/collect_gene_from_annotation.py"
compare="scripts/analysis/compare_call_and_annotation.py"

# IGHV
python3 ${collect} -r target_annotation/HG001-denovo/annotation_imperfect_HG001-denovo_BCRV.txt -l ${IGHV_func} > target_annotation/HG001-denovo/annotation_IGHV_functional.txt
python3 ${collect} -r target_annotation/HG002-denovo/annotation_imperfect_HG002-denovo_BCRV.txt -l ${IGHV_func} > target_annotation/HG002-denovo/annotation_IGHV_functional.txt
python3 ${collect} -r target_annotation/HG002/annotation_imperfect_HG002_BCRV.txt               -l ${IGHV_func} > target_annotation/HG002/annotation_IGHV_functional.txt
python3 ${collect} -r target_annotation/HG001-denovo/annotation_imperfect_HG001-denovo_BCRV.txt > target_annotation/HG001-denovo/annotation_IGHV.txt
python3 ${collect} -r target_annotation/HG002-denovo/annotation_imperfect_HG002-denovo_BCRV.txt > target_annotation/HG002-denovo/annotation_IGHV.txt
python3 ${collect} -r target_annotation/HG002/annotation_imperfect_HG002_BCRV.txt               > target_annotation/HG002/annotation_IGHV.txt

python3 ${compare} -an target_annotation/HG001-denovo/annotation_IGHV_curated.txt    -rc target_call/HG001/HG001_BCRV/read_depth_functional.rpt > target_call/HG001/HG001_BCRV/read_depth_annotation.rpt
python3 ${compare} -an target_annotation/HG002-denovo/annotation_IGHV_functional.txt -rc target_call/HG002/HG002_BCRV/read_depth_functional.rpt > target_call/HG002/HG002_BCRV/read_depth_annotation-denovo.rpt
python3 ${compare} -an target_annotation/HG002/annotation_IGHV_functional.txt        -rc target_call/HG002/HG002_BCRV/read_depth_functional.rpt > target_call/HG002/HG002_BCRV/read_depth_annotation.rpt

# IGKV
python3 ${collect} -r target_annotation/HG001-denovo/annotation_imperfect_HG001-denovo_BCRV.txt -l ${IGKV_func} > target_annotation/HG001-denovo/annotation_IGKV_functional.txt
python3 ${collect} -r target_annotation/HG002-denovo/annotation_imperfect_HG002-denovo_BCRV.txt -l ${IGKV_func} > target_annotation/HG002-denovo/annotation_IGKV_functional.txt
python3 ${collect} -r target_annotation/HG002/annotation_imperfect_HG002_BCRV.txt               -l ${IGKV_func} > target_annotation/HG002/annotation_IGKV_functional.txt

python3 ${compare} -an target_annotation/HG001-denovo/annotation_IGKV_functional.txt -rc target_call/HG001/HG001_BCRV/read_depth_functional.IGK.rpt > target_call/HG001/HG001_BCRV/read_depth_annotation.IGK.rpt
python3 ${compare} -an target_annotation/HG002-denovo/annotation_IGKV_functional.txt -rc target_call/HG002/HG002_BCRV/read_depth_functional.IGK.rpt > target_call/HG002/HG002_BCRV/read_depth_annotation-denovo.IGK.rpt
python3 ${compare} -an target_annotation/HG002/annotation_IGKV_functional.txt        -rc target_call/HG002/HG002_BCRV/read_depth_functional.IGK.rpt > target_call/HG002/HG002_BCRV/read_depth_annotation.IGK.rpt

#IGLV
python3 ${collect} -r target_annotation/HG001-denovo/annotation_imperfect_HG001-denovo_BCRV.txt -l ${IGLV_func} > target_annotation/HG001-denovo/annotation_IGLV_functional.txt
python3 ${collect} -r target_annotation/HG002-denovo/annotation_imperfect_HG002-denovo_BCRV.txt -l ${IGLV_func} > target_annotation/HG002-denovo/annotation_IGLV_functional.txt
python3 ${collect} -r target_annotation/HG002/annotation_imperfect_HG002_BCRV.txt               -l ${IGLV_func} > target_annotation/HG002/annotation_IGLV_functional.txt

python3 ${compare} -an target_annotation/HG001-denovo/annotation_IGLV_functional.txt -rc target_call/HG001/HG001_BCRV/read_depth_functional.IGL.rpt > target_call/HG001/HG001_BCRV/read_depth_annotation.IGL.rpt
python3 ${compare} -an target_annotation/HG002-denovo/annotation_IGLV_functional.txt -rc target_call/HG002/HG002_BCRV/read_depth_functional.IGL.rpt > target_call/HG002/HG002_BCRV/read_depth_annotation-denovo.IGL.rpt
python3 ${compare} -an target_annotation/HG002/annotation_IGLV_functional.txt        -rc target_call/HG002/HG002_BCRV/read_depth_functional.IGL.rpt > target_call/HG002/HG002_BCRV/read_depth_annotation.IGL.rpt
