cd  intra_mm_near
cp  ../../intra/mm1m/all_* .
python 1.py
cd  ../

cd  intra_ms_near
cp  ../../intra/ms1m/all_* .
python 1.py
cd  ../

cd  intra_sm_near
cp  ../../intra/sm1s/all_* .
python 1.py
cd  ../

cd  intra_ss_near
cp  ../../intra/ss1s/all_* .
python 1.py
cd  ../

cd  inter_ss_near
cp  ../../intra/ss1sgap/all_* .
python 1.py
cd  ../

cd  intra_mm_long
cp  ../../long/mm1m/all_* .
python 1.py
cd  ../

cd  intra_ms_long
cp  ../../long/ms1m/all_* .
python 1.py
cd  ../

cd  intra_sm_long
cp  ../../long/sm1s/all_* .
python 1.py
cd  ../

cd  intra_ss_long
cp  ../../long/ss1s/all_* .
python 1.py
cd  ../

cd  inter_mm_long
cp  ../../long/mm1mgap/all_* .
python 1.py
cd  ../

cd  inter_ms_long
cp  ../../long/ms1m/all_* .
python 1.py
cd  ../

cd  inter_sm_long
cp  ../../long/sm1s/all_* .
python 1.py
cd  ../

cd  inter_ss_long
cp  ../../long/ss1s/all_* .
python 1.py
cd  ../

cd  self_m
bash do.bash
cd ../

cd  self_s
bash do.bash 
cd ../
