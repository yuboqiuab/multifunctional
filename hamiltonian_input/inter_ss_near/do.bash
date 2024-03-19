rm all_*
python data.py

for ((j=1;j<=662;j=j+1))
    do
    cat inputss1s_"$j"  >> all_inputss1s
    cat inputss2s_"$j"  >> all_inputss2s
    cat outputss_"$j"  >> all_outputss
    cat vectorss_"$j"  >> all_vectorss
    cat 25outputss_"$j"  >> all_25outputss
    done



cd ss1s
bash do.bash
cd ../


