

## build pipeline scripts

now=$(date +"%m-%d-%Y_%H:%M:%S")

## project info
project="DLBC"
SubmitRNAseqExe="Submit_${PWD##*/}.sh"
padding="example/"

## command
echo "START" `date` " Running build_rnaseq.py"
python3 SRC/Python/build_rnaseq.py \
	--projdir $PWD \
	--metadata $PWD/${padding}$project.metadata.txt \
	--config $PWD/${padding}$project.pipeline.yaml \
	--systype cluster \
	--threads 8 \
	--log_file $PWD/Build_RNAseq.$project.$now.log

## submit pipeline master script
echo "START" `date` " Running $SubmitRNAseqExe"
echo "bash $SubmitRNAseqExe"

echo "END" `date`
