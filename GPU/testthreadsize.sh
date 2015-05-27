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
PROG="$PROG_HOME/hydrodynamic-flow-cuda"
SRCFILE="$PROG.cu"
RESFILE="$PROG_HOME/threadsizetest.csv"

export LC_NUMERIC="en_US.UTF-8"

printf "M,N,time\n" > "$RESFILE"
sed -rie "s/^#define threadPerBlockX .*/#define threadPerBlockX 16/" "$SRCFILE"
sed -rie "s/^#define threadPerBlockY .*/#define threadPerBlockY 16/" "$SRCFILE"
sed -rie "s/^#define STRIDE .*/#define STRIDE 1024/" "$SRCFILE"
sed -rie "s/^#define WIDTH .*/#define WIDTH 1023/" "$SRCFILE"
sed -rie "s/^#define HEIGHT .*/#define HEIGHT 1023/" "$SRCFILE"

for M in 2 4 6 8 16 32
do
    for N in 1 2 4 6 8 16 32
    do
		if (( $M*$N > 1 ))
		then
			printf "Block size (%2s, %2s): " $M $N
			sed -rie "s/^#define cellPerThreadX .*/#define cellPerThreadX $M/" "$SRCFILE"
			sed -rie "s/^#define cellPerThreadY .*/#define cellPerThreadY $N/" "$SRCFILE"
			make 2>&1 >/dev/null

			t1=$(/usr/bin/time -f "%e" "$PROG" 2>&1 >/dev/null)
			sleep 1
			t2=$(/usr/bin/time -f "%e" "$PROG" 2>&1 >/dev/null)
			sleep 1
			t3=$(/usr/bin/time -f "%e" "$PROG" 2>&1 >/dev/null)
			TIME=$(echo "($t1+$t2+$t3)/3" | bc -l)
			printf "%.02f, %.02f, %.02f,  avg: %.03f\n" $t1 $t2 $t3 $TIME
			printf "%.02f,%.02f,%.03f\n" $M $N $TIME >> "$RESFILE"
		fi
    done
done
