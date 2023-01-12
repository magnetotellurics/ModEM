!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module GlobalVariables
    !
    use Grid3D_SG
    !
    use ModelParameterCell_SG
    !
    use ModelOperator_MF
    !
    use ForwardSolverIT_DC
    !
    use ModelCovarianceRec
    !
    use TransmitterArray
    !
    use ReceiverArray
    !
    use DataGroupTxArray
    !
    !> Global Variables
    class( Grid_t ), allocatable, target :: main_grid
    class( ModelParameter_t ), allocatable :: sigma0, pmodel
    class( ModelOperator_t ), allocatable :: model_operator
    !
    class( ForwardSolver_t ), allocatable, target :: forward_solver
    !
    class( ModelCovarianceRec_t ), allocatable :: model_cov
	!
    !> Program control variables
    character(50) :: outdir_name
    !
    character(:), allocatable :: control_file_name
    character(:), allocatable :: model_file_name
    character(:), allocatable :: pmodel_file_name
    character(:), allocatable :: data_file_name
    character(:), allocatable :: dsigma_file_name
    character(:), allocatable :: e_solution_file_name
    character(:), allocatable :: modem_job
    !
    !> Program control flags
    logical :: has_outdir_name
    !
    logical :: has_control_file
    logical :: has_model_file
    logical :: has_pmodel_file
    logical :: has_data_file
    logical :: has_e_solution_file
    logical :: verbosis
    !
end module GlobalVariables
!