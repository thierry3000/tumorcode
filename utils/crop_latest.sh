#~/bin/bash
N=0
for i in $(h5ls $1); do
	if [[ "$i" == "out"* ]]; then
		testarray[$N]="$i"
		echo "$N = $i" #to confirm

		let "N= $N + 1"
	fi
done
echo "laste out is: $N"
h5copy -i $1 -o cropped_$N.h5 -s parameters -d parameters
h5copy -i $1 -o cropped_$N.h5 -s field_ld -d field_ld
#h5copy -i $1 -o cropped_$N.h5 -s last_state -d last_state
h5copy -i $1 -o cropped_$N.h5 -s out0$N -d out0$N
