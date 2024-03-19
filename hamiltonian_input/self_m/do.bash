rm all_*
python data.py
for ((j=1;j<=662;j=j+1))
    do
    cat input_"$j" >> all_input
    cat output_"$j" >> all_output
    done


