python data.py

rm -r mm2m
mkdir mm2m
cd mm2m
cp ../../momodel/*.tar .
cp ../../momodel/pmm2m.py .
cp ../../momodel/mm2m.py .
cp ../../momodel/utils.py .
python pmm2m.py
python mm2m.py
cd ../

rm -r mm1m
mkdir mm1m
cd mm1m
cp ../../momodel/*.tar .
cp ../../momodel/pmm1m.py .
cp ../../momodel/mm1m.py .
cp ../../momodel/utils.py .
cp ../../momodel/merge.py .
python pmm1m.py
python mm1m.py
cp reduce_inputmm1m 1
cp ../mm2m/reduce_inputmm2m 2
python merge.py
cp ../outputmm_1  all_output
cp ../vectormm_1  all_vector
cp ../indexmm_1  all_index
cd ../

rm -r sm2m
mkdir sm2m
cd sm2m
cp ../../momodel/*.tar .
cp ../../momodel/psm2m.py .
cp ../../momodel/sm2m.py .
cp ../../momodel/utils.py .
python psm2m.py
python sm2m.py
cd ../

rm -r sm1s
mkdir sm1s
cd sm1s
cp ../../smodel/*.tar .
cp ../../smodel/psm1s.py .
cp ../../smodel/sm1s.py .
cp ../../smodel/utils.py .
cp ../../smodel/merge.py .
python psm1s.py
python sm1s.py
cp reduce_inputsm1s 1
cp ../sm2m/reduce_inputsm2m 2
python merge.py
cp ../outputsm_1  all_output
cp ../vectorsm_1  all_vector
cp ../indexsm_1  all_index
cd ../




rm -r ms2s
mkdir ms2s
cd ms2s
cp ../../smodel/*.tar .
cp ../../smodel/pms2s.py .
cp ../../smodel/ms2s.py .
cp ../../smodel/utils.py .
python pms2s.py
python ms2s.py
cd ../

rm -r ms1m
mkdir ms1m
cd ms1m
cp ../../momodel/*.tar .
cp ../../momodel/pms1m.py .
cp ../../momodel/ms1m.py .
cp ../../momodel/utils.py .
cp ../../momodel/merge.py .
python pms1m.py
python ms1m.py
cp reduce_inputms1m 1
cp ../ms2s/reduce_inputms2s 2
python merge.py
cp ../outputms_1  all_output
cp ../vectorms_1  all_vector
cp ../indexms_1  all_index
cd ../


rm -r ss2s
mkdir ss2s
cd ss2s
cp ../../smodel/*.tar .
cp ../../smodel/pss2s.py .
cp ../../smodel/ss2s.py .
cp ../../smodel/utils.py .
python pss2s.py
python ss2s.py
cd ../

rm -r ss1s
mkdir ss1s
cd ss1s
cp ../../smodel/*.tar .
cp ../../smodel/pss1s.py .
cp ../../smodel/ss1s.py .
cp ../../smodel/utils.py .
cp ../../smodel/merge.py .
python pss1s.py
python ss1s.py
cp reduce_inputss1s 1
cp ../ss2s/reduce_inputss2s 2
python merge.py
cp ../outputss_1  all_output
cp ../vectorss_1  all_vector
cp ../indexss_1  all_index
cd ../


rm -r ss2sgap
mkdir ss2sgap
cd ss2sgap
cp ../../smodel/*.tar .
cp ../../smodel/pss2sgap.py .
cp ../../smodel/ss2sgap.py .
cp ../../smodel/utils.py .
python pss2sgap.py
python ss2sgap.py
cd ../

rm -r ss1sgap
mkdir ss1sgap
cd ss1sgap
cp ../../smodel/*.tar .
cp ../../smodel/pss1sgap.py .
cp ../../smodel/ss1sgap.py .
cp ../../smodel/utils.py .
cp ../../smodel/merge.py .
python pss1sgap.py
python ss1sgap.py
cp reduce_inputss1s 1
cp ../ss2sgap/reduce_inputss2s 2
python merge.py
cp ../outputssgap_1  all_output
cp ../vectorssgap_1  all_vector
cp ../indexssgap_1  all_index
cd ../



rm -r mm2mgap
mkdir mm2mgap
cd mm2mgap
cp ../../momodel/*.tar .
cp ../../momodel/pmm2mgap.py .
cp ../../momodel/mm2mgap.py .
cp ../../momodel/utils.py .
python pmm2mgap.py
python mm2mgap.py
cd ../

rm -r mm1mgap
mkdir mm1mgap
cd mm1mgap
cp ../../momodel/*.tar .
cp ../../momodel/pmm1mgap.py .
cp ../../momodel/mm1mgap.py .
cp ../../momodel/utils.py .
cp ../../momodel/merge.py .
python pmm1mgap.py
python mm1mgap.py
cp reduce_inputmm1m 1
cp ../mm2mgap/reduce_inputmm2m 2
python merge.py
cp ../outputmmgap_1  all_output
cp ../vectormmgap_1  all_vector
cp ../indexmmgap_1  all_index
cd ../



rm -r sm2mgap
mkdir sm2mgap
cd sm2mgap
cp ../../momodel/*.tar .
cp ../../momodel/psm2mgap.py .
cp ../../momodel/sm2mgap.py .
cp ../../momodel/utils.py .
python psm2mgap.py
python sm2mgap.py
cd ../

rm -r sm1sgap
mkdir sm1sgap
cd sm1sgap
cp ../../smodel/*.tar .
cp ../../smodel/psm1sgap.py .
cp ../../smodel/sm1sgap.py .
cp ../../smodel/utils.py .
cp ../../smodel/merge.py .
python psm1sgap.py
python sm1sgap.py
cp reduce_inputsm1s 1
cp ../sm2mgap/reduce_inputsm2m 2
python merge.py
cp ../outputsmgap_1  all_output
cp ../vectorsmgap_1  all_vector
cp ../indexsmgap_1  all_index
cd ../





rm -r ms2sgap
mkdir ms2sgap
cd ms2sgap
cp ../../smodel/*.tar .
cp ../../smodel/pms2sgap.py .
cp ../../smodel/ms2sgap.py .
cp ../../smodel/utils.py .
python pms2sgap.py
python ms2sgap.py
cd ../

rm -r ms1mgap
mkdir ms1mgap
cd ms1mgap
cp ../../momodel/*.tar .
cp ../../momodel/pms1mgap.py .
cp ../../momodel/ms1mgap.py .
cp ../../momodel/utils.py .
cp ../../momodel/merge.py .
python pms1mgap.py
python ms1mgap.py
cp reduce_inputms1m 1
cp ../ms2sgap/reduce_inputms2s 2
python merge.py
cp ../outputmsgap_1  all_output
cp ../vectormsgap_1  all_vector
cp ../indexmsgap_1  all_index
cd ../

