#!/bin/bash
clean_exit() {
    killall beacon-localize
    exit
}
trap 'clean_exit' SIGINT

evaluations=("evaluate-los" "evaluate-nlos" "evaluate-nlos2" "evaluate-nlos3")
sweeps=("by-measurements" "by-delay-sigma" "by-delay-mean")

for evaluation in "${evaluations[@]}"
do
    for sweep in "${sweeps[@]}"
    do
        flags="--$evaluation --$sweep"
        echo "Evaluating $flags..."
        file="evaluations/${evaluation}_$sweep.csv"

        bin/beacon-localize $flags --eval-run-i 0 > tmp0.txt &
        pid0=$!
        bin/beacon-localize $flags --eval-run-i 1 > tmp1.txt &
        pid1=$!
        bin/beacon-localize $flags --eval-run-i 2 > tmp2.txt &
        pid2=$!
        bin/beacon-localize $flags --eval-run-i 3 > tmp3.txt &
        pid3=$!
        bin/beacon-localize $flags --eval-run-i 4 > tmp4.txt &
        pid4=$!
        wait $pid0
        wait $pid1
        wait $pid2
        wait $pid3
        wait $pid4
        cat tmp0.txt | sed "s/^.*\r//" | awk '{ORS=" "; print $4, $9; for(i=10;i<=NF;++i) print $i; print "\n"}' | tr -d ',' > $file
        cat tmp1.txt | sed "s/^.*\r//" | awk '{ORS=" "; print $4, $9; for(i=10;i<=NF;++i) print $i; print "\n"}' | tr -d ',' >> $file
        cat tmp2.txt | sed "s/^.*\r//" | awk '{ORS=" "; print $4, $9; for(i=10;i<=NF;++i) print $i; print "\n"}' | tr -d ',' >> $file
        cat tmp3.txt | sed "s/^.*\r//" | awk '{ORS=" "; print $4, $9; for(i=10;i<=NF;++i) print $i; print "\n"}' | tr -d ',' >> $file
        cat tmp4.txt | sed "s/^.*\r//" | awk '{ORS=" "; print $4, $9; for(i=10;i<=NF;++i) print $i; print "\n"}' | tr -d ',' >> $file
        sed -i '/^ *$/d' $file
    done
done

methods=("" "--first-peak" "--max-peak" "--triangulation" "--fixed-delays" "--manual-delays")
for evaluation in "${evaluations[@]}"
do
    for method in "${methods[@]}"
    do
        flags="--$evaluation $method"
        echo "Evaluating $flags..."
        file="evaluations/${evaluation}_${method:2}.csv"

        bin/beacon-localize $flags --eval-run-i 0 > tmp0.txt &
        pid0=$!
        bin/beacon-localize $flags --eval-run-i 1 > tmp1.txt &
        pid1=$!
        bin/beacon-localize $flags --eval-run-i 2 > tmp2.txt &
        pid2=$!
        bin/beacon-localize $flags --eval-run-i 3 > tmp3.txt &
        pid3=$!
        bin/beacon-localize $flags --eval-run-i 4 > tmp4.txt &
        pid4=$!
        wait $pid0
        wait $pid1
        wait $pid2
        wait $pid3
        wait $pid4
        cat tmp0.txt | sed "s/^.*\r//" | awk '{print $6}' | tr -d ',' > $file
        cat tmp1.txt | sed "s/^.*\r//" | awk '{print $6}' | tr -d ',' >> $file
        cat tmp2.txt | sed "s/^.*\r//" | awk '{print $6}' | tr -d ',' >> $file
        cat tmp3.txt | sed "s/^.*\r//" | awk '{print $6}' | tr -d ',' >> $file
        cat tmp4.txt | sed "s/^.*\r//" | awk '{print $6}' | tr -d ',' >> $file
        sed -i '/^ *$/d' $file
    done
done

rm tmp0.txt
rm tmp1.txt
rm tmp2.txt
rm tmp3.txt
rm tmp4.txt
