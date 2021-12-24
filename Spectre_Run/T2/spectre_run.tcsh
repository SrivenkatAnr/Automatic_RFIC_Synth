#tcsh
source ~/.cshrc
cd /home/ee18b038/cadence_project/PA_single_pt/T2/basic_tsmc_65_rcm
cp /home/ee18b038/Auto_Ckt_Synth_Codes/Automatic_RFIC_Synth/Netlists_Ref/circ.scs ./
spectre circ.scs =log circ_log.txt
exit