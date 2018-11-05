

## build pipeline scripts

now=$(date +"%m-%d-%Y_%H:%M:%S")

## project info
project="DLBC"
SubmitRNAseqExe="Submit_RNAseq.$project.sh"
padding="example/"

## command
echo "START" `date` " Running build_rnaseq.py"
python3 SRC/Python/build_rnaseq.py \
	--projdir $PWD \
	--metadata $PWD/${padding}$project.metadata.txt \
	--config $PWD/${padding}$project.pipeline.yaml \
	--systype cluster \
	--threads 4 \
	--log_file $PWD/Build_RNAseq.$project.$now.log

## submit pipeline master script
echo "START" `date` " Running Submit_$project.sh"
echo "bash Submit_$project.sh"

echo "END" `date`
