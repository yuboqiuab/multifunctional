rm all_*

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

    cat inputss1s_"$j"  >> all_inputss1s
    cat inputss2s_"$j"  >> all_inputss2s
    cat outputss_"$j"  >> all_outputss
    cat vectorss_"$j"  >> all_vectorss

    done

