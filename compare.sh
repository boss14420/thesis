#!/bin/bash - 
#===============================================================================
#
#          FILE: blocktest.sh
# 
#         USAGE: ./blocktest.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: BOSS14420 (), 
#  ORGANIZATION: 
#       CREATED: 09/09/2014 08:34
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

PROG_HOME="$(dirname $0)"
PROG1_HOME="$(dirname $0)/GPU"
PROG2_HOME="$(dirname $0)/CPU"
PROG1="hydrodynamic-flow-cuda"
PROG2="hydrodynamic-flow2"
SRC1FILE="$PROG1".cu
SRC2FILE="$PROG2".cc
TIMEFILE="$PROG_HOME/time.csv"
GPUFILE="$PROG_HOME/gpu.txt"
CPUFILE="$PROG_HOME/cpu.txt"

export LC_NUMERIC="en_US.UTF-8"

printf "size,GPU,CPU,ratio\n" "$TIMEFILE" > "$TIMEFILE"
rm -rf "$GPUFILE" "$CPUFILE"

cd "$PROG1_HOME"
sed -rie "s/^#define cellPerThreadX .*/#define cellPerThreadX 2/" "$SRC1FILE"
sed -rie "s/^#define cellPerThreadY .*/#define cellPerThreadY 8/" "$SRC1FILE"
sed -rie "s/^#define threadPerBlockX .*/#define threadPerBlockX 16/" "$SRC1FILE"
sed -rie "s/^#define threadPerBlockY .*/#define threadPerBlockY 16/" "$SRC1FILE"
cd ..

for SZ in 128 256 512 1024 1536 2048 3072 4096
do
	printf "Grid size (%4s, %4s): \n" $SZ $SZ
	SZm1=$(( $SZ-1 ))

	# GPU
	printf "\tGPU: "
	cd "$PROG1_HOME"
	sed -rie "s/^#define STRIDE .*/#define STRIDE $SZ/" "$SRC1FILE"
	sed -rie "s/^#define WIDTH .*/#define WIDTH $SZm1/" "$SRC1FILE"
	sed -rie "s/^#define HEIGHT .*/#define HEIGHT $SZm1/" "$SRC1FILE"
	make 2>&1 >/dev/null

	t1=$(/usr/bin/time -f "%e" "./$PROG1" 2>&1 >/dev/null)
	printf "%.02f, " $t1
	sleep 1
	t2=$(/usr/bin/time -f "%e" "./$PROG1" 2>&1 >/dev/null)
	printf "%.02f, " $t2
	sleep 1
	t3=$(/usr/bin/time -f "%e" "./$PROG1" 2>&1 >/dev/null)
	TIME1=$(echo "($t1+$t2+$t3)/3" | bc -l)
	printf "%.02f,  avg: %.03f\n" $t3 $TIME1
	cd ..

	# CPU
	printf "\tCPU: "
	cd "$PROG2_HOME"
	sed -rie "s/^#define STRIDE .*/#define STRIDE $SZ/" "$SRC2FILE"
	sed -rie "s/^#define WIDTH .*/#define WIDTH $SZm1/" "$SRC2FILE"
	sed -rie "s/^#define HEIGHT .*/#define HEIGHT $SZm1/" "$SRC2FILE"
	make 2>&1 >/dev/null

	t1=$(/usr/bin/time -f "%e" "./$PROG2" 2>&1 >/dev/null)
	printf "%.02f, " $t1
	sleep 1
	t2=$(/usr/bin/time -f "%e" "./$PROG2" 2>&1 >/dev/null)
	printf "%.02f, " $t2
	sleep 1
	t3=$(/usr/bin/time -f "%e" "./$PROG2" 2>&1 >/dev/null)
	TIME2=$(echo "($t1+$t2+$t3)/3" | bc -l)
	printf "%.02f,  avg: %.03f\n" $t3 $TIME2
	cd ..

	ratio=$(echo "$TIME2 / $TIME1" | bc -l)
	printf "\tratio: %.03f\n" $ratio

	printf "%4s\t%.03f\n" $SZ $TIME2 >> "$CPUFILE"
	printf "%4s\t%.03f\n" $SZ $TIME1 >> "$GPUFILE"

	printf "%s,%.03f,%.03f,%.03f\n" $SZ $TIME1 $TIME2 $ratio >> "$TIMEFILE"
done
