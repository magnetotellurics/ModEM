#!/bin/bash
#
# OO FWD
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -pd pred_mt_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -f -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -pd pred_mt_ani_oo.dat
#
# OO JMULT
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -jm jmhat_mt_iso_oo.dat
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -j -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -pm ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -jm jmhat_mt_ani_oo.dat
#
# OO JMULT_T
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -dm dsigma_mt_iso_oo.mod
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -jt -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -dm dsigma_mt_ani_oo.mod
#
# OO INV NLCG
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_nlcg.txt -o nlcg_mt_iso_oo
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_nlcg.txt -o nlcg_mt_ani_oo
#
# OO INV DCG
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ISO.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_dcg.txt -o dcg_mt_iso_oo
#
mpirun -np 22 --hostfile ../../hostfile_22.txt ./ModEM_MPI -i -m ../../modem_inputs/ANI_CSEM/model_VTI_h10_v10_ANI_LOGE.mod -d ../../modem_inputs/ANI_CSEM/Data_Z.dat -cf fwd_ctrl_em1d.txt -ci inv_ctrl_dcg.txt -o dcg_mt_ani_oo
#
