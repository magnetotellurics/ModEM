#!/bin/bash
#
# OO FWD
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_dip1d.txt -pd pred_csem_dip1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -pd pred_csem_em1d_ani_oo.dat
#
# OO JMULT
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_dip1d.txt -jm jmhat_csem_dip1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -jm jmhat_csem_em1d_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -jm jmhat_csem_em1d_ani_oo.dat
#
# OO JMULT_T
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_dip1d.txt -dm dsigma_csem_dip1d_iso_oo.mod
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -dm dsigma_csem_em1d_iso_oo.mod
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -dm dsigma_csem_em1d_ani_oo.mod
#
# OO INV
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_dip1d.txt -ci inv_ctrl_template.txt -o csem_dip1d_iso_oo
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_template.txt -o csem_em1d_iso_oo
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Ey.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_template.txt -o csem_em1d_ani_oo
#
