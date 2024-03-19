rm all_*
python data.py

for ((j=1;j<=662;j=j+1))
    do
    cat inputmm1m_"$j"  >> all_inputmm1m
    cat inputmm2m_"$j"  >> all_inputmm2m
    cat outputmm_"$j"  >> all_outputmm
    cat vectormm_"$j"  >> all_vectormm

    cat inputms1m_"$j"  >> all_inputms1m
    cat inputms2s_"$j"  >> all_inputms2s
    cat outputms_"$j"  >> all_outputms
    cat vectorms_"$j"  >> all_vectorms

    cat inputsm1s_"$j"  >> all_inputsm1s
    cat inputsm2m_"$j"  >> all_inputsm2m
    cat outputsm_"$j"  >> all_outputsm
    cat vectorsm_"$j"  >> all_vectorsm


    done

cd mm1m
bash do.bash
cd ../


cd ms1s
bash do.bash
cd ../

cd sm1m
bash do.bash
cd ../
