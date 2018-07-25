

## build pipeline scripts

now=$(date +"%m-%d-%Y_%H:%M:%S")

## project info
project=DLBC
SubmitRNAseqExe=Submit_RNAseq.DLBC.sh

## command
echo "START" `date`
echo "START" `date` " Running build_rnaseq.py"
python3 SRC/Python/build_rnaseq.py \
	--projdir $PWD/example/$project \
	--metadata $PWD/example/$project.metadata.txt \
	--config $PWD/example/$project.pipeline.yaml \
	--systype cluster \
	--threads 4 \
	--log_file $PWD/Build_RNAseq.$project.$now.log

## submit pipeline master script
echo "START" `date` " Running example/$project/Submit_RNAseq.$project.sh"
echo "bash example/$project/Submit_RNAseq.$project.sh"

echo "END" `date`
