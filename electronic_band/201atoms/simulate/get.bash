rm merge
rm all_pt.dat

cd intra_ms_near
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd intra_sm_near
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd intra_mm_near
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd intra_ss_near
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd inter_ss_near
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../



cd self_m
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd self_s
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../


cd intra_ms_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd intra_sm_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd intra_mm_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd intra_ss_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd inter_ss_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd inter_sm_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd inter_ms_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../

cd inter_mm_long
paste all_index pred_target.dat > merge
cat pred_target.dat >> ../all_pt.dat
cat merge >> ../merge
cd ../
