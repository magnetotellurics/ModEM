#!/bin/bash
#
# ARGUMENTS: 
#		$1 - MODELS FULL PATH
#		$2 - DATA FULL PATH
#		$3 - EXECUTABLE
#		$4 - NUMBER OF CORES
#
MODELS_DIR=$1
DATA_DIR=$2
EXEC=$3
ncores=$4   #GET ENVIROMENT NUMBER OF CORES
#
emptyspace=""
#
# GENERATE MODEL DATA TEST OUTPUT FOLDER
now=$(date "+%Y_%m_%d_%H_%M_%S")
mkdir -p test_model_data_$now/
cd test_model_data_$now/
#
#
for MODEL in ../$MODELS_DIR/*.ws
do
	model_name=$MODEL
	model_name=${model_name/.ws/$emptyspace}
	model_name=${model_name/*\//$emptyspace}
	#
	for DATA in ../$DATA_DIR/*.dat
	do
		data_name=$DATA
		data_name=${data_name/.dat/$emptyspace}
		data_name=${data_name/*\//$emptyspace}
		#
		echo "### START Mod3DMT: $ncores CORES, FOR '$model_name' AND '$data_name' ###" | tee std_out.txt
		#
		# RUN MODEM-ON
		mpirun -n $ncores ../$EXEC -F ../${MODELS_DIR}${model_name}.ws ../${DATA_DIR}${data_name}.dat w${data_name}.dat | tee std_out.txt
		#
		echo "### FINISH Mod3DMT ###" | tee std_out.txt
		#
		# RUN COMPARE DATA SCRIPT
		bash ../scripts/compare_data.sh w${data_name}.dat ../${DATA_DIR}${data_name}.dat test_model_data_$now/
 	#
	done
done
#
cd ..
#
# MOVE MODEL_DATA_TEST OUTPUT FOLDER TO MAIN OUTPUT FOLDER
mv test_model_data_$now/ outputs/
#
#
# END OF SCRIPT
