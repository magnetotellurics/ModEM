#!/bin/bash
#
# OO FWD
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_dip1d.txt -pd pred_csem_dip1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_ani_oo.dat
#
# ON FWD
#
cd ../../modem-on/src/
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat pred_csem_dip1d_iso_on.dat esoln control_dip1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat pred_csem_em1d_iso_on.dat esoln control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -F ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat pred_csem_em1d_ani_on.dat esoln control_em1d
#
# OO JMULT
#
cd ../../modem-oo/src/
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_dip1d.txt -jm jmhat_csem_dip1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -jm jmhat_csem_em1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -jm jmhat_csem_em1d_ani_oo.dat
#
# ON JMULT
#
cd ../../modem-on/src/
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat jmhat_csem_dip1d_iso_on.dat control_dip1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat jmhat_csem_em1d_iso_on.dat control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -M ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat jmhat_csem_em1d_ani_on.dat control_em1d
#
# OO JMULT_T
#
cd ../../modem-oo/src/
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_dip1d.txt -dm dsigma_csem_dip1d_iso_oo.mod
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -dm dsigma_csem_em1d_iso_oo.mod
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -dm dsigma_csem_em1d_ani_oo.mod
#
# ON JMULT_T
#
cd ../../modem-on/src/
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat dsigma_csem_dip1d_iso_on.mod control_dip1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat dsigma_csem_em1d_iso_on.mod control_em1d
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./Mod3DMT_STD -T ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat dsigma_csem_em1d_ani_on.mod control_em1d
#
# OO INV
#
cd ../../modem-oo/src/
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_dip1d.txt -ci inv_ctrl_template.txt -o csem_dip1d_iso_oo
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_template.txt -o csem_em1d_iso_oo
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey_DIF.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_template.txt -o csem_em1d_ani_oo
#
