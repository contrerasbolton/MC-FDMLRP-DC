#!/bin/bash
# Script to run algorithms for MC-FDMLRP-DC

TYPE=$1
INSTANCE=$2
I=$3

OUTPUT="output/"
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
else
    echo "algorithm does not exist"
    exit
fi

if [ "$INSTANCE" = "A" ]
then
    if [ "$I" = "" ]
    then
	start=1
	end=15
    elif [ "$I" = "" ]
    then
	start=1
	end=15
    else
	start=$I
	end=$I
    fi
elif [ "$INSTANCE" = "B" ]
then
    if [ "$I" = "" ]
    then
	start=1
	end=5
    elif [ "$I" = "" ]
    then
	start=1
	end=5
    else
	start=$I
	end=$I
    fi
else
    echo "instance does not exist"
    exit
fi

for i in $(seq $star $end);
do
    run="./drone $INSTANCE$i.dat $algorithm"
    echo "$run > $OUTPUT$INSTANCE$i"_"$TYPE.txt"
    $run > $OUTPUT$INSTANCE$i"_"$TYPE.txt
done;
