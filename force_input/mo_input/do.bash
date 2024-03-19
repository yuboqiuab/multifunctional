rm all_input
rm all_output
python data.py
for ((j=1;j<=22052;j=j+1))
    do
    cat input_"$j" >> all_input
    cat output_"$j" >> all_output
    done
