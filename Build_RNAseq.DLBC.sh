

## build pipeline scripts

now=$(date +"%m-%d-%Y_%H:%M:%S")

## project info
project="DLBC"
SubmitRNAseqExe="Submit_RNAseq.$project.sh"
padding="example/"

## command
echo "START" `date` " Running build_rnaseq.py"
python3 SRC/Python/build_rnaseq.py \
	--projdir $PWD/${padding}$project \
	--metadata $PWD/${padding}$project.metadata.txt \
	--config $PWD/${padding}$project.pipeline.yaml \
	--systype cluster \
	--threads 4 \
	--log_file $PWD/Build_RNAseq.$project.$now.log

## submit pipeline master script
echo "START" `date` " Running ${padding}$project/Submit_RNAseq.$project.sh"
echo "bash ${padding}$project/Submit_RNAseq.$project.sh"

echo "END" `date`
