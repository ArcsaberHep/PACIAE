    program rms_analysis
!Lei20230819    implicit none
!Lei20230819    real(kind=8) :: a, b, c, d, e, f, g, h, o, p, q, r, s, t, u, v, w, x, y, z
!Lei20230819    integer(kind=8) :: i, j, k, l, m, n
        implicit real(kind=8) (a-h, o-z)   !Lei20230821
        implicit integer(kind=8) (i-n)   !Lei20230821
        integer(kind=8) :: n_file, n_event_per_file, n_event_total
        character(len=50) :: name_file, prefix, suffix
        real(kind=8) :: distr_h(40,4,2), sum_distr_h(40,4,2), mult_h(2), sum_mult_h(2)   !Lei20230820 mult_h. sum_mult_h
!Lei20230819B---
        integer(kind=8) :: i_sim_mode, n_current_date_and_time(8)
        real(kind=8) :: Ncoll_in_coll_list, Npart_in_coll_list, Ncoll_max_in_coll_list,         &
                        Ncoll_call_PYTHIA, Ncoll_not_call_PYTHIA, Npart_real,                   &
                        Ncoll_over_Npart_Optical_Glauber_multi_string, Ncoll_Optical_Glauber,   &
                        Npart_call_PYTHIA, Ncoll_nn_call_PYTHIA, Ncoll_pp_call_PYTHIA,          &
                        Ncoll_np_call_PYTHIA, Ncoll_call_PYTHIA_2, Ncoll_lp_call_PYTHIA,        &
                        Npart_Optical_Glauber_proj, Npart_Optical_Glauber_targ,                 &
                        Overlap_function_Optical_Glauber, Ncoll_Optical_Glauber_2,              &
                        n_coll_parton_rescattering_success, n_coll_parton_rescattering_block,   &
                        n_coll_parton_rescattering_total, n_coll_parton_rescattering,           &
                        n_coll_hadron_rescattering_elastic, n_coll_hadron_rescattering_inelastic,   &
                        n_coll_hadron_rescattering_total, n_gluon_in_a_string_single_multi_string,  &
                        n_string_single_string, n_q_thrown, n_qbar_thrown,                      &
                        multiplicity_minus_partial, multiplicity_plus_partial, multiplicity_total_partial,  &
                        multiplicity_minus_full, multiplicity_plus_full, multiplicity_total_full,   &
                        multiplicity_quark_partial, multiplicity_qbar_partial,                  &
                        multiplicity_q_qbar_total_partial, multiplicity_gluon_partial,          &
                        multiplicity_quark_full, multiplicity_qbar_full,                        &
                        multiplicity_q_qbar_total_full, multiplicity_gluon_full,                &
                        kapa_effective
        real(kind=8) :: n_process_parton_rescattering(9),  sum_n_process_parton_rescattering(9),    &
                        n_process_hadron_rescattering(6,60), sum_n_process_hadron_rescattering(6,60),    &
                        multiplicity_specie_partial(21), sum_multiplicity_specie_partial(21),       &
                        multiplicity_specie_full(21), sum_multiplicity_specie_full(21),             &
                        multiplicity_parton_partial(14), sum_multiplicity_parton_partial(14),       &
                        multiplicity_parton_full(14), sum_multiplicity_parton_full(14)
        real(kind=8) :: x_coor(40,6)
        real(kind=8) :: distr_hadron_partial(21,40,6), sum_distr_hadron_partial(21,40,6),   &
                        distr_hadron_full(21,40,6), sum_distr_hadron_full(21,40,6),         &
                        distr_parton_partial(14,40,6), sum_distr_parton_partial(14,40,6),   &
                        distr_parton_full(14,40,6), sum_distr_parton_full(14,40,6)
        character(len=200) :: comment_line(400)
        character(len=26) :: name_specie(21), shebang_sign
        character(len=4) :: c_date_and_time(8)
        character(len=4) :: id_abscissa(6)
        data id_abscissa / "y", "pT", "eta", "mT", "mult", "pT" /
!Lei20230819E---


        open(98,file="input_rms_analysis.dat",status="unknown")
        read(98,*) n_file, n_event_per_file   !Lei20230820
        read(98,*) i_sim_mode, i_b_sampling_method   !Lei20230819
        close(98)
        n_event_total = n_file * n_event_per_file


! Initialization.
        sum_distr_h = 0D0
!Lei2230819B---
        comment_line = ""
        name_specie  = ""
        sum_Ncoll_in_coll_list = 0D0
        sum_Npart_in_coll_list = 0D0
        sum_Ncoll_max_in_coll_list = 0D0
        sum_Ncoll_call_PYTHIA = 0D0
        sum_Ncoll_not_call_PYTHIA = 0D0
        sum_Npart_real = 0D0
        sum_Ncoll_over_Npart_Optical_Glauber_multi_string = 0D0
        sum_Ncoll_optical_Glauber = 0D0
        sum_Npart_call_PYTHIA = 0D0
        sum_Ncoll_nn_call_PYTHIA = 0D0
        sum_Ncoll_pp_call_PYTHIA = 0D0
        sum_Ncoll_np_call_PYTHIA = 0D0
        sum_Ncoll_call_PYTHIA_2 = 0D0
        sum_Ncoll_lp_call_PYTHIA = 0D0
        sum_ave_b_param = 0D0
        sum_avneu = 0D0
        sum_Npart_Optical_Glauber_proj = 0D0
        sum_Npart_Optical_Glauber_targ = 0D0
        sum_Overlap_function_Optical_Glauber = 0D0
        sum_Ncoll_Optical_Glauber_2 = 0D0
        sum_E_gamma_1 = 0D0
        sum_E_gamma_2 = 0D0
        sum_E_gamma_3 = 0D0
        sum_E_gamma_4 = 0D0
        sum_n_coll_parton_rescattering_success = 0D0
        sum_n_coll_parton_rescattering_block = 0D0
        sum_n_coll_parton_rescattering_total = 0D0
        sum_n_process_parton_rescattering = 0D0
        sum_n_coll_hadron_rescattering_elastic = 0D0
        sum_n_coll_hadron_rescattering_inelastic = 0D0
        sum_n_coll_hadron_rescattering_total = 0D0
        sum_parj1_effective = 0D0
        sum_parj2_effective = 0D0
        sum_parj3_effective = 0D0
        sum_parj4_effective = 0D0
        sum_parj21_effective = 0D0
        sum_kapa_effective = 0D0
        sum_n_gluon_in_a_string_single_multi_string = 0D0
        sum_xi_factor_single_string = 0D0
        sum_pT_hardest_gluon_single_string = 0D0
        sum_n_string_single_string = 0D0
        sum_time_NN_collision = 0D0
        sum_time_parton_rescattering = 0D0
        sum_time_hadron_rescattering = 0D0
        sum_n_q_thrown = 0D0
        sum_n_qbar_thrown = 0D0
        sum_charge_thrown = 0D0
        sum_px_thrown = 0D0
        sum_py_thrown = 0D0
        sum_pz_thrown = 0D0
        sum_E_thrown = 0D0
        sum_multiplicity_minus_partial = 0D0
        sum_multiplicity_plus_partial = 0D0
        sum_multiplicity_minus_full = 0D0
        sum_multiplicity_plus_full = 0D0
        sum_multiplicity_specie_partial = 0D0
        sum_multiplicity_specie_full = 0D0
        sum_distr_hadron_partial = 0D0
        sum_distr_hadron_full = 0D0
        sum_distr_parton_partial = 0D0
        sum_distr_parton_full = 0D0
        sum_multiplicity_quark_partial = 0D0
        sum_multiplicity_qbar_partial = 0D0
        sum_multiplicity_gluon_partial = 0D0
        sum_multiplicity_quark_full = 0D0
        sum_multiplicity_qbar_full = 0D0
        sum_multiplicity_gluon_full = 0D0
        sum_n_process_hadron_rescattering = 0D0
!Lei2230819E---


        ! Cycle over all of the "rms_i.out" files.
        loop_file: do i_file=1,n_file,1
            ! Opens one of the "rms_i.out" files.
            write(suffix,*) i_file
            name_file = "rms_" // TRIM(ADJUSTL(suffix)) // ".out"
            open(99,file=name_file,status="old")
!Lei20230819B----------
            ! Reecords comment line.
            !1 #!***************************************************************************!#
            !2 #!*********************|    PACIAE  Analysis  Output    |********************!#
            !3 #!***************************************************************************!#
            !4
            !5 #! Now is   HH:MM:SS   DD/MM/YYYY
            !6 #! Seed (PYTHIA default=19780503) =   xxxxxxxxx
            !7
            !8 #!-----------------------------------------------------------------------------
            !9 #! parp81, parp82, bp, mstp82 =
            do j=1,9,1
                read(99,"(A200)") comment_line(j)
            end do
            !10
            read(99,*) PARP81, PARP82, b_param_current, MSTP82
            !11 #!-----------------------------------------------------------------------------
            !12 #! MC Glauber-like <N_coll>, <N_part> =
            do j=11,12,1
                read(99,"(A200)") comment_line(j)
            end do
            !13
            read(99,*) Ncoll_in_coll_list, Npart_in_coll_list
                sum_Ncoll_in_coll_list = sum_Ncoll_in_coll_list + Ncoll_in_coll_list
                sum_Npart_in_coll_list = sum_Npart_in_coll_list + Npart_in_coll_list
            !14 #! largest ave. # of NN collision pairs =
            read(99,"(A200)") comment_line(14)
            !15
            read(99,*) Ncoll_max_in_coll_list
                sum_Ncoll_max_in_coll_list = sum_Ncoll_max_in_coll_list + Ncoll_max_in_coll_list
            !16 #! ave. # of NN collision pairs calling pythia, not calling pythia =
            read(99,"(A200)") comment_line(16)
            !17
            read(99,*) Ncoll_call_PYTHIA, Ncoll_not_call_PYTHIA
                sum_Ncoll_call_PYTHIA = sum_Ncoll_call_PYTHIA + Ncoll_call_PYTHIA
                sum_Ncoll_not_call_PYTHIA = sum_Ncoll_not_call_PYTHIA + Ncoll_not_call_PYTHIA
            !18 #! ave. # of wounded nucleons in parini =
            read(99,"(A200)") comment_line(18)
            !19
            read(99,*) Npart_real
                sum_Npart_real = sum_Npart_real + Npart_real
            !20 #! colli. # suffered by projectile nucleon in target nucleus
            read(99,"(A200)") comment_line(20)
            !21
            read(99,*) Ncoll_over_Npart_Optical_Glauber_multi_string
                sum_Ncoll_over_Npart_Optical_Glauber_multi_string =   &
                    sum_Ncoll_over_Npart_Optical_Glauber_multi_string +   &
                    Ncoll_over_Npart_Optical_Glauber_multi_string

            !22 #! event averaged N_bin
            read(99,"(A200)") comment_line(22)
            !23
            read(99,*) Ncoll_Optical_Glauber
                sum_Ncoll_Optical_Glauber = sum_Ncoll_Optical_Glauber + Ncoll_Optical_Glauber
            !24 #! (Npart)mini-jet, Nnn, Npp=
            read(99,"(A200)") comment_line(24)
            !25
            read(99,*) Npart_call_PYTHIA, Ncoll_nn_call_PYTHIA, Ncoll_pp_call_PYTHIA
                sum_Npart_call_PYTHIA    = sum_Npart_call_PYTHIA    + Npart_call_PYTHIA
                sum_Ncoll_nn_call_PYTHIA = sum_Ncoll_nn_call_PYTHIA + Ncoll_nn_call_PYTHIA
                sum_Ncoll_pp_call_PYTHIA = sum_Ncoll_pp_call_PYTHIA + Ncoll_pp_call_PYTHIA
            !26 #! Nnp, Ntot, Nep=
            read(99,"(A200)") comment_line(26)
            !27
            read(99,*) Ncoll_np_call_PYTHIA, Ncoll_call_PYTHIA_2, Ncoll_lp_call_PYTHIA
                sum_Ncoll_np_call_PYTHIA = sum_Ncoll_np_call_PYTHIA + Ncoll_np_call_PYTHIA
                sum_Ncoll_call_PYTHIA_2  = sum_Ncoll_call_PYTHIA_2  + Ncoll_call_PYTHIA_2
                sum_Ncoll_lp_call_PYTHIA = sum_Ncoll_lp_call_PYTHIA + Ncoll_lp_call_PYTHIA
            !28-29
            ! For the case of fixed impact-parameter (b).
            if( i_b_sampling_method == 0 ) then
                read(99,*)
                read(99,*)
            ! For the case of random impact-parameter (b) sampling method.
            else if( i_b_sampling_method == 1 ) then
                !28 #! event averaged b, avneu, Npart_p, Npart_t, T_pt=
                read(99,"(A200)") comment_line(28)
                !29
                read(99,*) ave_b_param, avneu,  &
                           Npart_Optical_Glauber_proj, Npart_Optical_Glauber_targ,  &
                           Overlap_function_Optical_Glauber
                    sum_ave_b_param = sum_ave_b_param + ave_b_param
                    sum_avneu = sum_avneu + avneu   ! TODO(Lei20230819): ???
                    sum_Npart_Optical_Glauber_proj = sum_Npart_Optical_Glauber_proj +   &
                                                         Npart_Optical_Glauber_proj
                    sum_Npart_Optical_Glauber_targ = sum_Npart_Optical_Glauber_targ +   &
                                                         Npart_Optical_Glauber_targ
                    sum_Overlap_function_Optical_Glauber =  &
                    sum_Overlap_function_Optical_Glauber +  &
                        Overlap_function_Optical_Glauber
            ! For the case of random impact-parameter (b) sampling method.
            else if( i_b_sampling_method == 2 ) then
                !28 #! psno, ave. b, N_part and N_bin =
                read(99,"(A200)") comment_line(28)
                !29
                read(99,*) b_sampling_method, ave_b_param,  &
                           Npart_Optical_Glauber_proj, Npart_Optical_Glauber_targ,  &
                           Ncoll_Optical_Glauber_2
                    sum_ave_b_param = sum_ave_b_param + ave_b_param
                    sum_Npart_Optical_Glauber_proj = sum_Npart_Optical_Glauber_proj +   &
                                                         Npart_Optical_Glauber_proj
                    sum_Npart_Optical_Glauber_targ = sum_Npart_Optical_Glauber_targ +   &
                                                         Npart_Optical_Glauber_targ
                    sum_Ncoll_Optical_Glauber_2 = sum_Ncoll_Optical_Glauber_2 + &
                                                      Ncoll_Optical_Glauber_2
            end if
            !30 #!-----------------------------------------------------------------------------
            !31 #! event averaged energy of gamma after partonic initiation, partonic cascade,
            !32 #!  hadronization and end of event =
            do j=30,32,1
                read(99,"(A200)") comment_line(j)
            end do
            !33
            read(99,*) E_gamma_1, E_gamma_2, E_gamma_3, E_gamma_4
                sum_E_gamma_1 = sum_E_gamma_1 + E_gamma_1
                sum_E_gamma_2 = sum_E_gamma_2 + E_gamma_2
                sum_E_gamma_3 = sum_E_gamma_3 + E_gamma_3
                sum_E_gamma_4 = sum_E_gamma_4 + E_gamma_4
            !34 #!-----------------------------------------------------------------------------
            !35 #! # of successful, blocked and all collision in parton cascade =
            do j=34,35,1
                read(99,"(A200)") comment_line(j)
            end do
            !36
            read(99,*) n_coll_parton_rescattering_success,  &
                       n_coll_parton_rescattering_block,    &
                       n_coll_parton_rescattering_total
                sum_n_coll_parton_rescattering_success =    &
                sum_n_coll_parton_rescattering_success +    &
                    n_coll_parton_rescattering_success
                sum_n_coll_parton_rescattering_block =  &
                sum_n_coll_parton_rescattering_block +  &
                    n_coll_parton_rescattering_block
                sum_n_coll_parton_rescattering_total =  &
                sum_n_coll_parton_rescattering_total +  &
                    n_coll_parton_rescattering_total
            !37 #! average collision # in parton cascade =
            read(99,"(A200)") comment_line(37)
            !38
            read(99,*) n_coll_parton_rescattering
                sum_n_coll_parton_rescattering = sum_n_coll_parton_rescattering +   &
                                                     n_coll_parton_rescattering
            !39 #! # of scaterring processes in parton cascade
            read(99,"(A200)") comment_line(39)
            !40-42
            read(99,*) ( n_process_parton_rescattering(k), k=1,3,1 )
            read(99,*) ( n_process_parton_rescattering(k), k=4,6,1 )
            read(99,*) ( n_process_parton_rescattering(k), k=7,9,1 )
                sum_n_process_parton_rescattering = sum_n_process_parton_rescattering +    &
                                                        n_process_parton_rescattering
            !43 #! average frequency of the occurring of each inela. in hadron cascade (at the end of the file)
            !44 #! el. and inel. coll. # and sum in hadron cascade=
            do j=43,44,1
                read(99,"(A200)") comment_line(j)
            end do
            !45
            read(99,*) n_coll_hadron_rescattering_elastic,      &
                       n_coll_hadron_rescattering_inelastic,    &
                       n_coll_hadron_rescattering_total
                sum_n_coll_hadron_rescattering_elastic =    &
                sum_n_coll_hadron_rescattering_elastic +    &
                    n_coll_hadron_rescattering_elastic
                sum_n_coll_hadron_rescattering_inelastic =  &
                sum_n_coll_hadron_rescattering_inelastic +  &
                    n_coll_hadron_rescattering_inelastic
                sum_n_coll_hadron_rescattering_total =  &
                sum_n_coll_hadron_rescattering_total +  &
                    n_coll_hadron_rescattering_total
            !46 #!-----------------------------------------------------------------------------
            !47 #! default parj1, parj2, parj3, parj4, parj21 =
            do j=46,47,1
                read(99,"(A200)") comment_line(j)
            end do
            !48
            read(99,*) parj1_default, parj2_default, parj3_default, parj4_default, parj21_default
            !49 #! Eff-parj1, parj2, parj3, parj4, parj21, keff =
            read(99,"(A200)") comment_line(49)
            !50
            read(99,*) parj1_effective, parj2_effective, parj3_effective,   &
                       parj4_effective, parj21_effective, kapa_effective
                sum_parj1_effective = sum_parj1_effective + parj1_effective
                sum_parj2_effective = sum_parj2_effective + parj2_effective
                sum_parj3_effective = sum_parj3_effective + parj3_effective
                sum_parj4_effective = sum_parj4_effective + parj4_effective
                sum_parj21_effective = sum_parj21_effective + parj21_effective
                sum_kapa_effective = sum_kapa_effective + kapa_effective
            !51 #! averaged # of gluon in a string when kjp22=1,3
            read(99,"(A200)") comment_line(51)
            !52
            read(99,*) n_gluon_in_a_string_single_multi_string
                sum_n_gluon_in_a_string_single_multi_string =   &
                sum_n_gluon_in_a_string_single_multi_string +   &
                    n_gluon_in_a_string_single_multi_string
            !53 #! event averaged value of the factor related to # 
            !54 #!  of gluons and hardest gluon in a string, event 
            !55 #!  averaged transverse momentum of hardest gluon,
            !56 #!  event averaged # strings when kjp22=1,3 =
            do j=53,56,1
                read(99,"(A200)") comment_line(j)
            end do
            !57
            read(99,*) xi_factor_single_string, pT_hardest_gluon_single_string, n_string_single_string
                sum_xi_factor_single_string = sum_xi_factor_single_string + &
                                                  xi_factor_single_string
                sum_pT_hardest_gluon_single_string = sum_pT_hardest_gluon_single_string +   &
                                                         pT_hardest_gluon_single_string
                sum_n_string_single_string = sum_n_string_single_string +   &
                                                 n_string_single_string
            !58 #!-----------------------------------------------------------------------------
            !59 #! times & sum=
            do j=58,59,1
                read(99,"(A200)") comment_line(j)
            end do
            !60
            read(99,*) time_NN_collision, time_parton_rescattering, time_hadron_rescattering
                sum_time_NN_collision = sum_time_NN_collision + time_NN_collision
                sum_time_parton_rescattering = sum_time_parton_rescattering + time_parton_rescattering
                sum_time_hadron_rescattering = sum_time_hadron_rescattering + time_hadron_rescattering
            !61 #!-----------------------------------------------------------------------------
            !62 #! q, qbar, charge thrown away =
            do j=61,62,1
                read(99,"(A200)") comment_line(j)
            end do
            !63
            read(99,*) n_q_thrown, n_qbar_thrown, charge_thrown
                sum_n_q_thrown    = sum_n_q_thrown    + n_q_thrown
                sum_n_qbar_thrown = sum_n_qbar_thrown + n_qbar_thrown
                sum_charge_thrown = sum_charge_thrown + charge_thrown
            !64 #! 3-momentum and energy thrown away =
            read(99,"(A200)") comment_line(64)
            !65
            read(99,*) px_thrown, py_thrown, pz_thrown, E_thrown
                sum_px_thrown = sum_px_thrown + px_thrown
                sum_py_thrown = sum_py_thrown + py_thrown
                sum_pz_thrown = sum_pz_thrown + pz_thrown
                sum_E_thrown  = sum_E_thrown  + E_thrown
            !66 Null line.
            !67 #!-----------------------------------------------------------------------------
            !68 #! multiplicity of negative, positive particles and sums, partial & full =
            do j=66,68,1
                read(99,"(A200)") comment_line(j)
            end do
            !69-70
            read(99,*) multiplicity_minus_partial, multiplicity_plus_partial, multiplicity_total_partial
            read(99,*) multiplicity_minus_full, multiplicity_plus_full, multiplicity_total_full
                sum_multiplicity_minus_partial = sum_multiplicity_minus_partial + multiplicity_minus_partial
                sum_multiplicity_plus_partial  = sum_multiplicity_plus_partial  + multiplicity_plus_partial
                sum_multiplicity_minus_full    = sum_multiplicity_minus_full    + multiplicity_minus_full
                sum_multiplicity_plus_full     = sum_multiplicity_plus_full     + multiplicity_plus_full
            !71 #!-----------------------------------------------------------------------------
            !72 #! particle multiplicity, partial =
            do j=71,72,1
                read(99,"(A200)") comment_line(j)
            end do
            !73   pi+K+p name_1 ... name_20
            read(99,*) shebang_sign, ( name_specie(j), j=1,21,1 )
            !74
            read(99,*) ( multiplicity_specie_partial(j), j=1,21,1 )
                sum_multiplicity_specie_partial = sum_multiplicity_specie_partial + multiplicity_specie_partial
            !75 #! particle multiplicity, full    =
            read(99,"(A200)") comment_line(75)
            !76
            read(99,*) ( multiplicity_specie_full(j), j=1,21,1 )
                sum_multiplicity_specie_full = sum_multiplicity_specie_full + multiplicity_specie_full
            !77-79 Null line.
            do j=1,3,1
                read(99,*)
            end do

            ! Read out the data of 6 distributions, both hadrons and partons in 
            !   the partial and full phase-spaces.
            
            ! For hadrons:
            ! #!*******************|    Hadron  Distribution  Output    |******************!#
            read(99,"(A200)") comment_line(200)
            i_line = 200
            distr_hadron_partial = 0D0
            distr_hadron_full = 0D0
            loop_hadron: do i_distr=1,6,1
                ! Partial phase-space.
                ! #!-----------------------------------------------------------------------------
                ! #! ID of distribution m2=           x    xxxx              
                ! #! partial phase-space,  x.xx < pT <  x.xx  (nominal cuts of 1-st setting in usu.dat)
                do j=1,3,1
                    i_line = i_line + 1
                    read(99,"(A200)") comment_line(i_line)
                end do
                ! #! x    pi+K+p  name_1 ... name_20
                read(99,*)   ! Jumps out this line.
                ! Reads out data of partial phase-space.
                do i_bin=1,40,1
                    read(99,*) x_coor(i_bin,i_distr),   &
                        ( distr_hadron_partial(i_specie,i_bin,i_distr), i_specie=1,21,1 )
                end do
                ! Full phase-space.
                ! #! ID of distribution m2=           x    xxxx              
                ! #! full phase-space
                ! #! x    pi+K+p  name_1 ... name_20
                do j=1,2,1
                    i_line = i_line + 1
                    read(99,"(A200)") comment_line(i_line)
                end do
                ! #! x    pi+K+p  name_1 ... name_20
                read(99,*)   ! Jumps out this line.
                do i_bin=1,40,1
                    read(99,*) x_coor(i_bin,i_distr),   &
                        ( distr_hadron_full(i_specie,i_bin,i_distr), i_specie=1,21,1 )
                end do
            end do loop_hadron
            sum_distr_hadron_partial = sum_distr_hadron_partial + distr_hadron_partial
            sum_distr_hadron_full    = sum_distr_hadron_full    + distr_hadron_full

            ! For partons:
            ! Jumps out 3 null lines and 2 comment lines.
            ! #!-----------------------------------------------------------------------------
            ! #! multiplicity of negative, positive quark, sums and gluon, partial & full =
            do j=1,5,1
                read(99,*)
            end do
            read(99,*) multiplicity_quark_partial, multiplicity_qbar_partial,   &
                       multiplicity_q_qbar_total_partial, multiplicity_gluon_partial
            read(99,*) multiplicity_quark_full, multiplicity_qbar_full,     &
                       multiplicity_q_qbar_total_full, multiplicity_gluon_full
                sum_multiplicity_quark_partial = sum_multiplicity_quark_partial + multiplicity_quark_partial
                sum_multiplicity_qbar_partial  = sum_multiplicity_qbar_partial  + multiplicity_qbar_partial
                sum_multiplicity_gluon_partial = sum_multiplicity_gluon_partial + multiplicity_gluon_partial
                sum_multiplicity_quark_full    = sum_multiplicity_quark_full    + multiplicity_quark_full
                sum_multiplicity_qbar_full     = sum_multiplicity_qbar_full     + multiplicity_qbar_full
                sum_multiplicity_gluon_full    = sum_multiplicity_gluon_full    + multiplicity_gluon_full
            ! #!-----------------------------------------------------------------------------
            ! #! particle multiplicity, partial =
            ! #! g      u+d+s + anti- ... u, d, s, ...
            read(99,*)  ! Jumps out this line.
            read(99,*)
            read(99,*)
            read(99,*) ( multiplicity_parton_partial(j), j=1,14,1 )
                sum_multiplicity_parton_partial = sum_multiplicity_parton_partial + multiplicity_parton_partial
            ! #! particle multiplicity, full    =
            read(99,*)  ! Jumps out this line.
            read(99,*) ( multiplicity_parton_full(j), j=1,14,1 )
                sum_multiplicity_parton_full = sum_multiplicity_parton_full + multiplicity_parton_full
            ! Null lines.
            read(99,*)  ! Jumps out this line.
            read(99,*)
            read(99,*)
            ! #!*******************|    Parton  Distribution  Output    |******************!#
            read(99,"(A200)") comment_line(300)
            i_line = 300
            distr_parton_partial = 0D0
            distr_parton_full = 0D0
            loop_parton: do i_distr=1,6,1
                ! Partial phase-space.
                ! #!-----------------------------------------------------------------------------
                ! #! ID of distribution m2=           x    xxxx              
                ! #! partial phase-space,  x.xx < pT <  x.xx  (nominal cuts of 1-st setting in usu.dat)
                do j=1,3,1
                    i_line = i_line + 1
                    read(99,"(A200)") comment_line(i_line)
                end do
                ! #! x    g     u+d+s + anti-   u, d, s, ...
                read(99,*)   ! Jumps out this line.
                ! Reads out data of partial phase-space.
                do i_bin=1,40,1
                    read(99,*) x_coor(i_bin,i_distr),   &
                        ( distr_parton_partial(i_specie,i_bin,i_distr), i_specie=1,14,1 )
                end do
                ! Full phase-space.
                ! #! ID of distribution m2=           x    xxxx              
                ! #! full phase-space
                do j=1,2,1
                    i_line = i_line + 1
                    read(99,"(A200)") comment_line(i_line)
                end do
                ! #! x    g     u+d+s + anti-   u, d, s, ...
                read(99,*)   ! Jumps out this line.
                do i_bin=1,40,1
                    read(99,*) x_coor(i_bin,i_distr),   &
                        ( distr_parton_full(i_specie,i_bin,i_distr), i_specie=1,14,1 )
                end do
            end do loop_parton
            sum_distr_parton_partial = sum_distr_parton_partial + distr_parton_partial
            sum_distr_parton_full    = sum_distr_parton_full    + distr_parton_full

            ! 3 null lines.
            ! #! average frequency of the occurring of each inela. in hadron cascade =
            do j=1,4,1
                read(99,*)
            end do
            ! Reads out data.
            n_process_hadron_rescattering = 0D0
            do j=1,60,1
                read(99,*) ( n_process_hadron_rescattering(k,j), k=1,6,1 )
            end do
            sum_n_process_hadron_rescattering = sum_n_process_hadron_rescattering +     &
                                                    n_process_hadron_rescattering
!Lei20230819E----------

            !Lei20230820    3 null lines and 1 comment line.
            ! #! Inv. dN/dpT, dN/dpT, dN/dy, dN/deta of h+-, p and f
            IF(.FALSE.)THEN
            read(99,"(3/)")
            distr_h = 0D0
            do j=1,40,1
                read(99,*) distr_h(j,1,1), distr_h(j,1,2), &
                           distr_h(j,2,1), distr_h(j,2,2), &
                           distr_h(j,3,1), distr_h(j,3,2), &
                           distr_h(j,4,1), distr_h(j,4,2)
            end do
                sum_distr_h = sum_distr_h + distr_h
            !Lei20230820B-
            ! #! partial/full multiplicity of h+-:
            read(99,"(1/)")
            read(99,*) mult_h(1), mult_h(2)
                sum_mult_h = sum_mult_h + mult_h
            !Lei20230820E-
            ENDIF

            close(99)
        end do loop_file


!Lei20230820 Event-(file-) averaging.
        sum_Ncoll_in_coll_list = sum_Ncoll_in_coll_list / n_file
        sum_Npart_in_coll_list = sum_Npart_in_coll_list / n_file
        sum_Ncoll_max_in_coll_list = sum_Ncoll_max_in_coll_list / n_file
        sum_Ncoll_call_PYTHIA = sum_Ncoll_call_PYTHIA / n_file
        sum_Ncoll_not_call_PYTHIA = sum_Ncoll_not_call_PYTHIA / n_file
        sum_Npart_real = sum_Npart_real / n_file
        sum_Ncoll_over_Npart_Optical_Glauber_multi_string =     &
            sum_Ncoll_over_Npart_Optical_Glauber_multi_string / n_file
        sum_Ncoll_Optical_Glauber = sum_Ncoll_Optical_Glauber / n_file
        sum_Npart_call_PYTHIA = sum_Npart_call_PYTHIA / n_file
        sum_Ncoll_nn_call_PYTHIA = sum_Ncoll_nn_call_PYTHIA / n_file
        sum_Ncoll_pp_call_PYTHIA = sum_Ncoll_pp_call_PYTHIA / n_file
        sum_Ncoll_np_call_PYTHIA = sum_Ncoll_np_call_PYTHIA / n_file
        sum_Ncoll_call_PYTHIA_2 = sum_Ncoll_call_PYTHIA_2 / n_file
        sum_Ncoll_lp_call_PYTHIA = sum_Ncoll_lp_call_PYTHIA / n_file
        sum_ave_b_param = sum_ave_b_param / n_file
        sum_avneu = sum_avneu / n_file
        sum_Npart_Optical_Glauber_proj = sum_Npart_Optical_Glauber_proj / n_file
        sum_Npart_Optical_Glauber_targ = sum_Npart_Optical_Glauber_targ / n_file
        sum_Overlap_function_Optical_Glauber =  &
            sum_Overlap_function_Optical_Glauber / n_file
        sum_Ncoll_Optical_Glauber_2 = sum_Ncoll_Optical_Glauber_2 / n_file
        sum_E_gamma_1 = sum_E_gamma_1 / n_file
        sum_E_gamma_2 = sum_E_gamma_2 / n_file
        sum_E_gamma_3 = sum_E_gamma_3 / n_file
        sum_E_gamma_4 = sum_E_gamma_4 / n_file
        sum_n_coll_parton_rescattering_success = sum_n_coll_parton_rescattering_success / n_file
        sum_n_coll_parton_rescattering_block = sum_n_coll_parton_rescattering_block / n_file
        sum_n_coll_parton_rescattering_total = sum_n_coll_parton_rescattering_total / n_file
        sum_n_coll_parton_rescattering = sum_n_coll_parton_rescattering / n_file
        sum_n_process_parton_rescattering = sum_n_process_parton_rescattering / n_file
        sum_n_coll_hadron_rescattering_elastic = sum_n_coll_hadron_rescattering_elastic / n_file
        sum_n_coll_hadron_rescattering_inelastic = sum_n_coll_hadron_rescattering_inelastic / n_file
        sum_n_coll_hadron_rescattering_total = sum_n_coll_hadron_rescattering_total / n_file
        sum_parj1_effective = sum_parj1_effective / n_file
        sum_parj2_effective = sum_parj2_effective / n_file
        sum_parj3_effective = sum_parj3_effective / n_file
        sum_parj4_effective = sum_parj4_effective / n_file
        sum_parj21_effective = sum_parj21_effective / n_file
        sum_kapa_effective = sum_kapa_effective / n_file
        sum_n_gluon_in_a_string_single_multi_string = sum_n_gluon_in_a_string_single_multi_string / n_file
        sum_xi_factor_single_string = sum_xi_factor_single_string / n_file
        sum_pT_hardest_gluon_single_string = sum_pT_hardest_gluon_single_string / n_file
        sum_n_string_single_string = sum_n_string_single_string / n_file
        sum_time_NN_collision = sum_time_NN_collision / n_file
        sum_time_parton_rescattering = sum_time_parton_rescattering / n_file
        sum_time_hadron_rescattering = sum_time_hadron_rescattering / n_file
        sum_n_q_thrown = sum_n_q_thrown / n_file
        sum_n_qbar_thrown = sum_n_qbar_thrown / n_file
        sum_charge_thrown = sum_charge_thrown / n_file
        sum_px_thrown = sum_px_thrown / n_file
        sum_py_thrown = sum_py_thrown / n_file
        sum_pz_thrown = sum_pz_thrown / n_file
        sum_E_thrown = sum_E_thrown / n_file
        sum_multiplicity_minus_partial = sum_multiplicity_minus_partial / n_file
        sum_multiplicity_plus_partial = sum_multiplicity_plus_partial / n_file
        sum_multiplicity_minus_full = sum_multiplicity_minus_full / n_file
        sum_multiplicity_plus_full = sum_multiplicity_plus_full / n_file
        sum_multiplicity_specie_partial = sum_multiplicity_specie_partial / n_file
        sum_multiplicity_specie_full = sum_multiplicity_specie_full / n_file

        sum_distr_hadron_partial = sum_distr_hadron_partial / n_file
        sum_distr_hadron_full = sum_distr_hadron_full / n_file

        sum_multiplicity_quark_partial = sum_multiplicity_quark_partial / n_file
        sum_multiplicity_qbar_partial = sum_multiplicity_qbar_partial / n_file
        sum_multiplicity_gluon_partial = sum_multiplicity_gluon_partial / n_file
        sum_multiplicity_quark_full = sum_multiplicity_quark_full / n_file
        sum_multiplicity_qbar_full = sum_multiplicity_qbar_full / n_file
        sum_multiplicity_gluon_full = sum_multiplicity_gluon_full / n_file
        sum_multiplicity_parton_partial = sum_multiplicity_parton_partial / n_file
        sum_multiplicity_parton_full = sum_multiplicity_parton_full / n_file

        sum_distr_parton_partial = sum_distr_parton_partial / n_file
        sum_distr_parton_full = sum_distr_parton_full / n_file

        sum_n_process_hadron_rescattering = sum_n_process_hadron_rescattering / n_file

        sum_distr_h = sum_distr_h / n_file
        sum_mult_h = sum_mult_h / n_file


!Lei20230820 output.
        open(100,file="ave_rms.out",status="unknown")
! Header.
!1 #!***************************************************************************!#
!2 #!*********************|    PACIAE  Analysis  Output    |********************!#
!3 #!***************************************************************************!#
!4
        do j=1,4,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
! Records the current real date and time.
!5 #! Now is   HH:MM:SS   DD/MM/YYYY
        call DATE_AND_TIME(VALUES=n_current_date_and_time)
        do i=1,8,1
            write(c_date_and_time(i),"(I4)") n_current_date_and_time(i)
        enddo
        write(100,*) "#! Now is   " //TRIM(ADJUSTL( c_date_and_time(5) )) //    &
                                ":" //TRIM(ADJUSTL( c_date_and_time(6) )) //    &
                                ":" //TRIM(ADJUSTL( c_date_and_time(7) )) //    &
                              "   " //TRIM(ADJUSTL( c_date_and_time(3) )) //    &
                                "/" //TRIM(ADJUSTL( c_date_and_time(2) )) //    &
                                "/" //TRIM(ADJUSTL( c_date_and_time(1) ))
!6 #! Seed (PYTHIA default=19780503) =   xxxxxxxxx
!7
!8 #!-----------------------------------------------------------------------------
!9 #! parp81, parp82, bp, mstp82 =
        do j=6,9,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!10
        write(100,*) PARP81, PARP82, b_param_current, MSTP82
!11 #!-----------------------------------------------------------------------------
!12 #! MC Glauber-like <N_coll>, <N_part> =
        do j=11,12,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!13
        write(100,*) sum_Ncoll_in_coll_list, sum_Npart_in_coll_list
!14 #! largest ave. # of NN collision pairs =
        write(100,*) TRIM(ADJUSTL( comment_line(14) ))
!15
        write(100,*) sum_Ncoll_max_in_coll_list
!16 #! ave. # of NN collision pairs calling pythia, not calling pythia =
        write(100,*) TRIM(ADJUSTL( comment_line(16) ))
!17
        write(100,*) sum_Ncoll_call_PYTHIA, sum_Ncoll_not_call_PYTHIA
!18 #! ave. # of wounded nucleons in parini =
        write(100,*) TRIM(ADJUSTL( comment_line(18) ))
!19
        write(100,*) sum_Npart_real
!20 #! colli. # suffered by projectile nucleon in target nucleus
        write(100,*) TRIM(ADJUSTL( comment_line(20) ))
        !21
        write(100,*) sum_Ncoll_over_Npart_Optical_Glauber_multi_string
!22 #! event averaged N_bin
        write(100,*) TRIM(ADJUSTL( comment_line(22) ))
!23
        write(100,*) sum_Ncoll_Optical_Glauber
!24 #! (Npart)mini-jet, Nnn, Npp=
        write(100,*) TRIM(ADJUSTL( comment_line(24) ))
!25
        write(100,*) sum_Npart_call_PYTHIA, sum_Ncoll_nn_call_PYTHIA, sum_Ncoll_pp_call_PYTHIA
!26 #! Nnp, Ntot, Nep=
        write(100,*) TRIM(ADJUSTL( comment_line(26) ))
!27
        write(100,*) sum_Ncoll_np_call_PYTHIA, sum_Ncoll_call_PYTHIA_2, sum_Ncoll_lp_call_PYTHIA
!28-29
! For the case of fixed impact-parameter (b).
        if( i_b_sampling_method == 0 ) then
            write(100,*)
            write(100,*)
! For the case of random impact-parameter (b) sampling method.
        else if( i_b_sampling_method == 1 ) then
        !28 #! event averaged b, avneu, Npart_p, Npart_t, T_pt=
            write(100,*) TRIM(ADJUSTL( comment_line(28) ))
        !29
            write(100,*) sum_ave_b_param, sum_avneu,  &
                        sum_Npart_Optical_Glauber_proj, sum_Npart_Optical_Glauber_targ,  &
                        sum_Overlap_function_Optical_Glauber
! For the case of random impact-parameter (b) sampling method.
        else if( i_b_sampling_method == 2 ) then
        !28 #! psno, ave. b, N_part and N_bin =
            write(100,*) TRIM(ADJUSTL( comment_line(28) ))
        !29
            write(100,*) INT(i_b_sampling_method), sum_ave_b_param,  &
                        sum_Npart_Optical_Glauber_proj, sum_Npart_Optical_Glauber_targ,  &
                        sum_Ncoll_Optical_Glauber_2
        end if
!30 #!-----------------------------------------------------------------------------
!31 #! event averaged energy of gamma after partonic initiation, partonic cascade,
!32 #!  hadronization and end of event =
        do j=30,32,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!33
        write(100,*) sum_E_gamma_1, sum_E_gamma_2, sum_E_gamma_3, sum_E_gamma_4
!34 #!-----------------------------------------------------------------------------
!35 #! # of successful, blocked and all collision in parton cascade =
        do j=34,35,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!36
        write(100,*) sum_n_coll_parton_rescattering_success,  &
                    sum_n_coll_parton_rescattering_block,    &
                    sum_n_coll_parton_rescattering_total
!37 #! average collision # in parton cascade =
        write(100,*) TRIM(ADJUSTL( comment_line(37) ))
!38
        write(100,*) sum_n_coll_parton_rescattering
!39 #! # of scaterring processes in parton cascade
        write(100,*) TRIM(ADJUSTL( comment_line(39) ))
!40-42
        write(100,*) ( sum_n_process_parton_rescattering(k), k=1,3,1 )
        write(100,*) ( sum_n_process_parton_rescattering(k), k=4,6,1 )
        write(100,*) ( sum_n_process_parton_rescattering(k), k=7,9,1 )
!43 #! average frequency of the occurring of each inela. in hadron cascade (at the end of the file)
!44 #! el. and inel. coll. # and sum in hadron cascade=
        do j=43,44,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!45
        write(100,*) sum_n_coll_hadron_rescattering_elastic,      &
                     sum_n_coll_hadron_rescattering_inelastic,    &
                     sum_n_coll_hadron_rescattering_total
!46 #!-----------------------------------------------------------------------------
!47 #! default parj1, parj2, parj3, parj4, parj21 =
        do j=46,47,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!48
        write(100,*) parj1_default, parj2_default, parj3_default, parj4_default, parj21_default
!49 #! Eff-parj1, parj2, parj3, parj4, parj21, keff =
        write(100,*) TRIM(ADJUSTL( comment_line(49) ))
!50
        write(100,*) sum_parj1_effective, sum_parj2_effective, sum_parj3_effective,   &
                     sum_parj4_effective, sum_parj21_effective, sum_kapa_effective
!51 #! averaged # of gluon in a string when kjp22=1,3
        write(100,*) TRIM(ADJUSTL( comment_line(51) ))
!52
        write(100,*) sum_n_gluon_in_a_string_single_multi_string
!53 #! event averaged value of the factor related to # 
!54 #!  of gluons and hardest gluon in a string, event 
!55 #!  averaged transverse momentum of hardest gluon,
!56 #!  event averaged # strings when kjp22=1,3 =
        do j=53,56,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!57
        write(100,*) sum_xi_factor_single_string, sum_pT_hardest_gluon_single_string, sum_n_string_single_string
!58 #!-----------------------------------------------------------------------------
!59 #! times & sum=
        do j=58,59,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
        !60
        write(100,*) sum_time_NN_collision, sum_time_parton_rescattering, sum_time_hadron_rescattering,     &
                     sum_time_NN_collision + sum_time_parton_rescattering + sum_time_hadron_rescattering
!61 #!-----------------------------------------------------------------------------
!62 #! q, qbar, charge thrown away =
        do j=61,62,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
        !63
        write(100,*) sum_n_q_thrown, sum_n_qbar_thrown, sum_charge_thrown
!64 #! 3-momentum and energy thrown away =
        write(100,*) TRIM(ADJUSTL( comment_line(64) ))
!65
        write(100,*) sum_px_thrown, sum_py_thrown, sum_pz_thrown, sum_E_thrown
!66 Null line.
!67 #!-----------------------------------------------------------------------------
!68 #! multiplicity of negative, positive particles and sums, partial & full =
        do j=66,68,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
        !69-70
        write(100,*) sum_multiplicity_minus_partial, sum_multiplicity_plus_partial,     &
                     sum_multiplicity_minus_partial + sum_multiplicity_plus_partial
        write(100,*) sum_multiplicity_minus_full, sum_multiplicity_plus_full,           &
                     sum_multiplicity_minus_full + sum_multiplicity_plus_full
!71 #!-----------------------------------------------------------------------------
!72 #! particle multiplicity, partial =
        do j=71,72,1
            write(100,*) TRIM(ADJUSTL( comment_line(j) ))
        end do
!73   pi+K+p name_1 ... name_20
        write(100,*) "#! ", ( name_specie(j), j=1,21,1 )
!74
        write(100,*) ( sum_multiplicity_specie_partial(j), j=1,21,1 )
!75 #! particle multiplicity, full    =
        write(100,*) TRIM(ADJUSTL( comment_line(75) ))
!76
        write(100,*) ( sum_multiplicity_specie_full(j), j=1,21,1 )
!77-79 Null line.
        do j=1,3,1
            write(100,*)
        end do

! Outputs the data of 6 distributions, both hadrons and partons in 
!   the partial and full phase-spaces.

! For hadrons:
        ! #!*******************|    Hadron  Distribution  Output    |******************!#
        i_line = 200
        write(100,*) TRIM(ADJUSTL( comment_line(i_line) ))
        loop_hadron_output: do i_distr=1,6,1
            ! Partial phase-space.
            ! #!-----------------------------------------------------------------------------
            ! #! ID of distribution m2=           x    xxxx              
            ! #! partial phase-space,  x.xx < pT <  x.xx  (nominal cuts of 1-st setting in usu.dat)
            do j=1,3,1
                i_line = i_line + 1
                write(100,*) TRIM(ADJUSTL( comment_line(i_line) ))
            end do
            ! #! x    pi+K+p  name_1 ... name_20
            write(100,*) "#!       " // id_abscissa(i_distr) // "               ", (name_specie(j),j=1,21,1)
            ! Outputs data of partial phase-space.
            do i_bin=1,40,1
                write(100,*) x_coor(i_bin,i_distr),   &
                    ( sum_distr_hadron_partial(i_specie,i_bin,i_distr), i_specie=1,21,1 )
            end do
            ! Full phase-space.
            ! #! ID of distribution m2=           x    xxxx              
            ! #! full phase-space
            ! #! x    pi+K+p  name_1 ... name_20
            do j=1,2,1
                i_line = i_line + 1
                write(100,*) TRIM(ADJUSTL( comment_line(i_line) ))
            end do
            ! #! x    pi+K+p  name_1 ... name_20
            write(100,*) "#!       " // id_abscissa(i_distr) // "               ", (name_specie(j),j=1,21,1)
            do i_bin=1,40,1
                write(100,*) x_coor(i_bin,i_distr),   &
                    ( sum_distr_hadron_full(i_specie,i_bin,i_distr), i_specie=1,21,1 )
            end do
        end do loop_hadron_output

! For partons:
        ! 3 null lines.
        do j=1,3,1
            write(100,*)
        end do
        write(100,*) "#!-----------------------------------------------------------------------------"
        write(100,*) "#! multiplicity of negative, positive quark, sums and gluon, partial & full ="
        write(100,*) sum_multiplicity_quark_partial, sum_multiplicity_qbar_partial,     &
                     sum_multiplicity_quark_partial + sum_multiplicity_qbar_partial,    &
                     sum_multiplicity_gluon_partial
        write(100,*) sum_multiplicity_quark_full, sum_multiplicity_qbar_full,       &
                     sum_multiplicity_quark_full + sum_multiplicity_qbar_full,      &
                     sum_multiplicity_gluon_full
        write(100,*) "#!-----------------------------------------------------------------------------"
        write(100,*) "#! particle multiplicity, partial ="
        write(100,*) "#! g                        u+d+s + anti-             " //      &
                     "u                         ubar                      " //      &
                     "d                         dbar                      " //      &
                     "s                         sbar                      " //      &
                     "c                         cbar                      " //      &
                     "b                         bbar                      " //      &
                     "t                         tbar"
        write(100,*) ( sum_multiplicity_parton_partial(j), j=1,14,1 )
        write(100,*) "#! particle multiplicity, full    ="
        write(100,*) ( sum_multiplicity_parton_full(j), j=1,14,1 )
        ! 3 null lines.
        do j=1,3,1
            write(100,*)
        end do
        ! #!*******************|    Parton  Distribution  Output    |******************!#
        i_line = 300
        write(100,*) TRIM(ADJUSTL( comment_line(i_line) ))
        loop_parton_output: do i_distr=1,6,1
            ! Partial phase-space.
            ! #!-----------------------------------------------------------------------------
            ! #! ID of distribution m2=           x    xxxx              
            ! #! partial phase-space,  x.xx < pT <  x.xx  (nominal cuts of 1-st setting in usu.dat)
            do j=1,3,1
                i_line = i_line + 1
                write(100,*) TRIM(ADJUSTL( comment_line(i_line) ))
            end do
            ! #! x    g     u+d+s + anti-   u, d, s, ...
            write(100,*) "#!       " // id_abscissa(i_distr) // "               " //    &
                         "g                         u+d+s + anti-             " //      &
                         "u                         ubar                      " //      &
                         "d                         dbar                      " //      &
                         "s                         sbar                      " //      &
                         "c                         cbar                      " //      &
                         "b                         bbar                      " //      &
                         "t                         tbar"
            ! Outputs data of partial phase-space.
            do i_bin=1,40,1
                write(100,*) x_coor(i_bin,i_distr),   &
                    ( sum_distr_parton_partial(i_specie,i_bin,i_distr), i_specie=1,14,1 )
            end do
            ! Full phase-space.
            ! #! ID of distribution m2=           x    xxxx              
            ! #! full phase-space
            do j=1,2,1
                i_line = i_line + 1
                write(100,*) TRIM(ADJUSTL( comment_line(i_line) ))
            end do
            ! #! x    g     u+d+s + anti-   u, d, s, ...
            write(100,*) "#!       " // id_abscissa(i_distr) // "               " //    &
                         "g                         u+d+s + anti-             " //      &
                         "u                         ubar                      " //      &
                         "d                         dbar                      " //      &
                         "s                         sbar                      " //      &
                         "c                         cbar                      " //      &
                         "b                         bbar                      " //      &
                         "t                         tbar"
            do i_bin=1,40,1
                write(100,*) x_coor(i_bin,i_distr),   &
                    ( sum_distr_parton_full(i_specie,i_bin,i_distr), i_specie=1,14,1 )
            end do
        end do loop_parton_output

        ! 3 null lines.
        do j=1,3,1
            write(100,*)
        end do
        write(100,*) "#! average frequency of the occurring of each inela. in hadron cascade ="
        ! Outputs data.
        do j=1,60,1
            write(100,*) ( sum_n_process_hadron_rescattering(k,j), k=1,6,1 )
        end do
!Lei20230819E----------

!Lei20230821B---
        stop   !Lei20230912
        ! 3 null lines.
        do j=1,3,1
            write(100,*)
        end do
        write(100,*) "#! Inv. dN/dpT, dN/dpT, dN/dy, dN/deta of h+-, p and f"
        do j=1,40,1
            write(100,*) sum_distr_h(j,1,1), sum_distr_h(j,1,2), &
                         sum_distr_h(j,2,1), sum_distr_h(j,2,2), &
                         sum_distr_h(j,3,1), sum_distr_h(j,3,2), &
                         sum_distr_h(j,4,1), sum_distr_h(j,4,2)
        end do
        ! #! partial/full multiplicity of h+-:   
        write(100,*)
        write(100,*) "#! partial/full multiplicity of h+-:"
        write(100,*) sum_mult_h(1), sum_mult_h(2)
!Lei20230821B---

        close(100)



        stop
    end program