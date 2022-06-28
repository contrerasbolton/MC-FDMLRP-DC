#!/bin/bash
# Script to run algorithms for MC-FDMLRP-DC

TYPE=$1
INSTANCE=$2
I=$3
J=$4
time=$5
memory=$6
OUTPUT="output/"
start=0
end=0
startSeed=0
endSeed=0
if [ "$TYPE" = "" ] && [ "$INSTANCE" = "" ] && [ "$I" = "" ] && [ "$J" = "" ] && [ "$time" = "" ] && [ "$memory" = "" ]
then
	echo "format: ./run algorithm instance [instance_i] [instance_j] [time] [memory]"
	exit
fi
if [ "$TYPE" = "M1-MTZ" ]
then
    algorithm=0
elif [ "$TYPE" = "M1-GG" ]
then
    algorithm=1
elif [ "$TYPE" = "M2-MTZ" ]
then
    algorithm=2
elif [ "$TYPE" = "M2-GG" ]
then
    algorithm=3
elif [ "$TYPE" = "M3-MTZ" ]
then
    algorithm=4
elif [ "$TYPE" = "M3-GG" ]
then
    algorithm=5
elif [ "$TYPE" = "MH" ]
then
    algorithm=6
    endSeed=9
elif [ "$TYPE" = "M1" ]
then
    algorithm=8
    endSeed=9
elif [ "$TYPE" = "M2" ]
then
    algorithm=7
    endSeed=9
else
    echo "algorithm does not exist"
    exit
fi

if [ "$INSTANCE" = "A" ]
then
    k1=1
    k2=2
    if [ "$I" = "" ]
    then
	start=1
	end=10
    elif [ "$I" != "" ] && [ "$J" != "" ]
    then
	start=$I
	end=$J
    else
	start=$I
	end=$I
    fi
elif [ "$INSTANCE" = "B" ]
then
    k1=1
    k2=3
    if [ "$I" = "" ]
    then
	start=1
	end=2
    elif [ "$I" != "" ] && [ "$J" != "" ]
    then
	start=$I
	end=$J
    else
	start=$I
	end=$I
    fi
else
    echo "instance does not exist"
    exit
fi

if [ "$time" != "" ]
then 
    time="-t $time"
fi

fixed="-a 40 -b 500000 -c 20000 -d 25000 -e 500 -f 0.136 -g 0.298 -h 0.26 -j 0.039 -k 0.268 -l 0.766 -m 5"
for i in `seq $start $end`;
do
    for k in `seq $k1 $k2`;
    do
        for seed in `seq $startSeed $endSeed`;
		do
		run="./drone -i $INSTANCE$i.$k.dat -o $algorithm -s $seed $fixed $time"
		echo "$seed $run > $OUTPUT$INSTANCE$i.$k"_"$TYPE.txt"
		$run > $OUTPUT$INSTANCE$i.$k"_"$TYPE"_"$seed".txt"
		done
	done
done;
