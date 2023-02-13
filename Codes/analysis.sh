PATIENTS_PATH="/home/users/l/l/llovat/stroke_pilab/study/subjects/"

STATIC_REF_FILE="/home/users/l/l/llovat/Memoire/Atlas_Maps/"

NUM_FILES=18

FILES=("xtract_prob_Corticospinal_Tract_L" "xtract_prob_Corticospinal_Tract_R" "xtract_prob_Frontal_Aslant_Tract_L" "xtract_prob_Frontal_Aslant_Tract_R" "xtract_prob_Superior_Longitudinal_Fasciculus_1_L" "xtract_prob_Superior_Longitudinal_Fasciculus_2_L" "xtract_prob_Superior_Longitudinal_Fasciculus_1_R" "xtract_prob_Superior_Longitudinal_Fasciculus_2_R" "xtract_prob_Superior_Longitudinal_Fasciculus_3_L" "xtract_prob_Superior_Longitudinal_Fasciculus_3_R" "xtract_prob_Superior_Thalamic_Radiation_L" "xtract_prob_Superior_Thalamic_Radiation_R" "Cortico_Ponto_Cerebellum_Left" "Cortico_Ponto_Cerebellum_Right" "mni_prob_Frontal_Lobe" "mni_prob_Parietal_Lobe" "harvardoxford-subcortical_prob_Right_Thalamus" "harvardoxford-subcortical_prob_Left_Thalamus")

NUM_PATIENTS=47

PATIENTS=("20_11_02_E1" "20_11_02_E0" "20_11_02_E2" "20_01_01_E0" "20_02_01_E0" "20_08_02_E1" "20_10_01_E1" "20_10_02_E0" "20_06_02_E3" "20_05_01_E1" "20_06_02_E1" "20_05_01_E2" "20_01_01_E1" "20_05_01_E3" "20_03_01_E0" "20_02_02_E3" "20_08_02_E0" "20_02_02_E2" "20_05_02_E0" "20_08_02_E2" "20_06_01_E1" "20_06_02_E2" "20_07_01_E2" "20_06_01_E0" "20_10_02_E3" "20_10_02_E1" "20_08_01_E0" "20_10_01_E2" "20_05_02_E1" "20_02_01_E2" "20_08_02_E3" "20_07_01_E0" "20_03_01_E1" "20_11_02_E3" "20_06_01_E2" "20_08_01_E2" "20_10_02_E2" "20_01_01_E2" "20_02_02_E1" "20_02_02_E0" "20_02_01_E1" "20_07_01_E1" "20_05_01_E0" "20_06_02_E0" "20_05_02_E2" "20_08_01_E1" "20_03_01_E2")

python3 ./analysis_fixedPatient.py $PATIENTS_PATH $STATIC_REF_FILE $NUM_FILES ${FILES[@]} $NUM_PATIENTS ${PATIENTS[@]}
echo "Done !"