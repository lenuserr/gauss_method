#!/bin/bash

dir="."
log="test-run.log"

N=30
r=4

args=(
  "4 3 4 0 ${dir}/t.ttt"
  "4 3 4 0 ${dir}/s.txt"
  "6 3 6 0 ${dir}/g.txt"
  "4 3 4 0 ${dir}/a.txt"
  "4 3 4 0 ${dir}/a20.txt"
  "4 3 4 0 ${dir}/a40.txt"
  "4 3 4 0 ${dir}/b.txt"
  "6 3 6 0 ${dir}/c.txt"
  "6 3 6 0 ${dir}/d.txt"
  "6 3 6 0 ${dir}/e.txt"
  "4 3 4 0 ${dir}/f.txt"
)

start_time=$(date +%s)

#clear log
echo "########################################################################" | tee ${log}

# test on files
for arg_i in "${args[@]}"
do \
  fname=$(basename $(echo "${arg_i}" | cut -d ' ' -f 5))
  echo "============================== a.out ${fname} =========================" | tee -a ${log}
  echo "### Run ./a.out ${arg_i}" | tee -a ${log}
  result=$(valgrind ./a.out ${arg_i} 2>&1)
  echo "${result}" >> ${log}
  echo "${result}" | grep ': Task'
done

# test by formula
for (( s = 1; s <= 4; s++ ))
do \
  for (( n = 1; n <= $N; n++ ))
  do \
    for (( m = 3; m <= n; m += 3 ))
    do \
      echo "============================== a.out n=$n m=$m s=$s ====================" | tee -a ${log}
      echo "### Run ./a.out $n $m $r $s" | tee -a ${log}
      result=$(valgrind ./a.out $n $m $r $s 2>&1)
      echo "${result}" >> $log
      echo "${result}" | grep ': Task'
    done
  done
done

echo "########################################################################" | tee -a ${log}

end_time=$(date +%s)
echo "SCRIPT ELAPSED: $((end_time - start_time)) SEC"

