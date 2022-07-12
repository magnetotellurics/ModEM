#!/bin/bash
#
# ARGUMENTS: 1 - OUTPUT FILE NAME
OUTPUT=$1
#
# STRING NOW
NOW=$(date "+%Y/%m/%d - %H:%M:%S")
#
#
THREAD_ID=0
#
THREAD_NAME="Master"
#
N_MASSIF_OUTS=0
MASSIF_OUTS=$(ls massif.out.*)
#
for MASSIF_OUT in $MASSIF_OUTS; do
	#
	N_MASSIF_OUTS=$((N_MASSIF_OUTS+1))
	#
done
#
for MASSIF_OUT in $MASSIF_OUTS; do
	#
	echo "$MASSIF_OUT => $THREAD_NAME"
	#
	echo "id Time(ms) $THREAD_NAME" > temp
	ms_print $MASSIF_OUT >> temp
	#
	# REMOVE LINES WITH SPECIFIC CHARACTERS
	sed -i '/|/d;/-/d;/:/d;/#/d;/MB/d' temp
	#
	# REMOVE EMPTY LINES
	sed -i '/^[[:space:]]*$/d' temp
	#
	# MERGE SPACES
	sed -i 's/  */ /g' temp
	#
	# REMOVE STARTING SPACES
	sed -i 's/^ //' temp
	#
	# REMOVE LINES STARTING WITH ZERO
	sed -i '/^0/d' temp
	#
	# ONLY TIME AND TOTAL MEM COLUMNS
	awk '{ print $2, $3 }' temp > $THREAD_NAME.txt
	rm temp
	#
	THREAD_ID=$((THREAD_ID+1))
	THREAD_NAME="Worker$THREAD_ID"
	#
	#
done
#
join -1 1 -2 1 -a1 -a2 -e"###" -o'0,1.2,2.2' Master.txt Worker1.txt > MemoryBrief.txt
#
# SPACES TO TABS
sed -i 's/[[:blank:]]/\t/g' MemoryBrief.txt
#
#
sed -i 's/###//' MemoryBrief.txt
#
exit 0
#
# END OF SCRIPT

