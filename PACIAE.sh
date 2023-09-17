#!/bin/bash
#SBATCH -J PACIAE_SLURM       # Task running name. Arbitrary.
#SBATCH -N 1                  # Number of nodes. Modify as needed.
#SBATCH --ntasks-per-node=1   # Number of cores used per node. Modify as needed.
#SBATCH -p middle             # Name of sequence. Modify as needed.

# Above statements are for normal LINUX system in PC and the SLURM scheduling 
#   system on the computing cluster or super-computer.
# Do not delete them but modify them as needed.

# For LSF scheduling system on the computing clusters or super-computers.
APP_NAME=Gsx_normal  # Name of sequence. Modify as needed.
NP=18                # Total number of cores used. Modify as needed.
NP_PER_NODE=18       # Number of cores used per node. Modify as needed.
RUN="RAW"            # Additional option. Not required to be modified usually.



################################################################################
################################################################################
################################################################################
###                                                                          ###
###   PPPPPPP       AAAAA       CCCCCCC    IIIIIII     AAAAA     EEEEEEEEE   ###
###   P      p     A     A     C       C      I       A     A    E           ###
###   P       p   A       A   C               I      A       A   E           ###
###   P      p    A       A   C               I      A       A   E           ###
###   PPPPPPP     AAAAAAAAA   C               I      AAAAAAAAA   EEEEEEEEE   ###
###   P           A       A   C               I      A       A   E           ###
###   P           A       A   C               I      A       A   E           ###
###   P           A       A    C       C      I      A       A   E           ###
###   P           A       A     CCCCCCC    IIIIIII   A       A   EEEEEEEEE   ###
###                                                                          ###
################################################################################
################################################################################
#                                                                              #
#                                                                              #
################################################################################
#                                                                              #
# This is a universal toy SHELL-script for PACIAE running on normal LINUX      #
#   system in the personal-computer and task submitting on SLURM and LSF       #
#   scheduling systems in the computing cluster and the sunper-computer.       #
#                                                                              #
# Mainly for NN, NA and AA collisions.                                         #
#                                                                              #
# This shell script will generate usu.dat and Makefile automatically. Then     #
#   run / submit tasks automatically, too.                                     #
# A series of folders will be generated to store different type of files,      #
#   including src (source folder, stores the source code files), bin (binary   #
#   folder, stores binary executable files), log (log records folder, stores   #
#   running informations) , and etc (editable text configuration folder,       #
#   stores configuration and other files). The PACIAE running folder with the  #
#   name of collisions system is also generated for easier management.         #
#                                                                              #
# It makes the (pseudo-)SPMD (single program, multiple data) possible to       #
#   achieve the (pseudo-)parallel running.                                     #
# Total N events are split up to multiple separate (N/number_of_CPU) events    #
#   and run simultaneously on multiple processors with different input random  #
#   seed based on the real-time clock of the machine.                          #
# This random number generator seed is set in PACIAE internal code (main.f).   #
#                                                                              #
# One needs to note that it's just a pseudo-parallism. In other words, this    #
#   script just performs PACIAE running one-by-one automatically instead of    #
#   manual runs. One also needs another program to aggregate and average all   #
#   the results (in rms.out) after the whole PACIAE running finished.          #
#                                                                              #
# This script has been tested on UBUNTU 20.04 in PC, on SLURM in CCNU          #
#   FARM-cluster (middle partition), and on LSF in National SuperComputing     #
#   Center in Shenzhen (NSCCSZ, Gsx_normal partition).                         #
#                                                                              #
################################################################################
#                                                                              #
# How to use:                                                                  #
#                                                                              #
#   1. First and formost, give this script file "executable permission" by     #
#        typing command "chmod +x PACIAE.sh". In addition, one needs "make"    #
#        tool. Install "make" using the command, for example on UBUNTU,        #
#        "sudo apt install make".                                              #
#                                                                              #
#   2. Modify the variabls needed.                                             #
#        2.1 Modify the names of source files for Makefile. They shoud be the  #
#              same as the real ones (except for name_x, it's user-defined).   #
#        2.2 Modify the variables: the number of program-running (n_run,       #
#              preferably not larger than the total number of CPU cores of the #
#              computer), the number of events per program-running (n_eve),    #
#              and other main setups of incident particles and collisions etc. #
#              "(D=)" means default value.                                     #
#            Here, the variables are aliases for the ones in usu.dat.          #
#            More detailed setups can be modified in the closely following     #
#              usu.dat directly. Search keywords "usu.dat" to find them.       #
#        2.3 On SLURM/LSF scheduling systems in the computing cluster and the  #
#              sunper-computer. One needs modify the additional statements of  #
#              the settings at the beginning of file. (SBATCH, APP_NAME...)    #
#                                                                              #
#   3. Execution.                                                              #
#        3.1 On normal LINUX, to run this script by typing command             #
#              "./PACIAE.sh".                                                  #
#            The more recommended command is                                   #
#              "time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log",           #
#              which stores the screen information to a log file.              #
#            Also, one can execute them in the backgroud using                 #
#              "nohup ./PACIAE.sh &" or                                        #
#              "nohup time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log &"    #
#       3.2 On SLURM: "sbatch PACIAE.sh". On LSF: "bsub PACIAE.sh"             #
#                                                                              #
#                                           By An-Ke Lei at CCNU on 17/10/2022 #
#                                                                              #
#                                      Last updated by An-Ke Lei on 11/09/2023 #
################################################################################
################################################################################
################################################################################



# Gets current date and time.
echo
current_date_and_time=$(date "+%Y%m%d%H%M%S")
echo "Now is $(date +"%Y-%m-%d %T")"
# mv ./usu.dat ./default_usu.dat
echo



################################################################################
#   Modify the following statements as needed.                                 #
################################################################################



################################################################################
################################################################################
####################             Files name             ########################
# For Makefile.
name_x="xPaciae.x"   # Name of the executable. User-defined.
name_f1="main_30"    # Names of the source files.
name_f2="parini_30"
name_f3="parcas_30"
name_f4="sfm_30"
name_f5="coales_30"
name_f6="hadcas_30"
name_f7="analy"
name_f8="p_30"
# name_f9="stahad_30"   # Not used now.
# name_f10="eps09"   # Not used now.
name_xRms="xRms_analysis.x"   #Lei20230608
name_rms="rms_analysis"   #Lei20230608
####################             Files name             ########################
################################################################################
################################################################################



################################################################################
################################################################################
####################      Variable's alias for usu.dat      ####################
#*******************************************************************************
# Setups of the program-running.
# Total number of events in all program-runnings = n_run * n_eve.
n_run=10     # Number of program-running, i.e. total number of CPU cores used.
# n_run=${NP}   # Uncommented this line if submit jods on LSF system.
# ((n_run=SLURM_JOB_NUM_NODES*SLURM_NTASKS_PER_NODE))   # Uncommented for SLURM.
n_eve=10    # neve, number of events per program-running (per core).
n_out=1     # (D=neve/10) nout, outputs per n_out events.
n_osc=0     # (D=0) nosc, OSCAR-1997A/1999A/2013A output or not (1/2/3/0)
i_random=1  # (D=1) adj(26), random generator seed in PYTHIA:
            #         =0, default PYTHIA seed (19780503), can be used for debug
            #         =1, seed from the real-time system clock (h-min-s-ms)

#*******************************************************************************
# Setups of the incident particles.
na_proj=197 # nap, nucleon number of projectile
nz_proj=79  # nzp, proton number of projectile
na_targ=197 # nat, nucleon number of target
nz_targ=79  # nzt, proton number of target
            # For NN, NA (AN), and AB collisions
            #   p+p   : 1,1,1,1;
            #   p+pbar: 1,1,1,-1; pbar+p: 1,-1,1,1
            #   p+n   : 1,1,1, 0; n+p   : 1, 0,1,1
            #   n+n   : 1,0,1, 0;
            #   p+Pb  : 1,1,208,82; Pb+p: 208,82,1,1;
            #   Au+Au : 197,79,197,79; Pb+Pb: 208,82,208,82;
            #   Xe+Xe : 129,54,129,54; U + U: 238,92,238,92;
            #   Ag+Ag : 108,47,108,47; Cu+Cu: 63,29,63,29;
            #   Ru+Ru : 96,44,96,44; Zr+Zr: 96,40,96,40.

#*******************************************************************************
# Setups of the simulation.
i_sim_mode=3    # (D=3) iMode, simulation mode: 
                #              =1, PACIAE hadronic simulation
                #              =2, PYTHIA-like hadronic simulation
                #              =3, PACIAE parton-hadron cascade simulation

#*******************************************************************************
# Setups of the collision.
b_min=0 #8.99 #0 #8.76  # bmin, min b parameter
b_max=3.72 #10.19 #3.72 #9.93  # bmax, max b parameter
i_frame=1   # (D=1) ifram, collision frame, 1=collider, 0=fixed target.
energy=200  # (GeV) win, colliding CMS energy SQRT(s_NN) / incident momentum
i_channel=8 # (D=0) nchan, =0: inelastic (INEL)
            #               1: Non Single Difractive (NSD)
            #               2: Drell-Yan process
            #               3: J/psi production
            #               4: heavy-flavor production
            #               5: direct photon
            #               6: soft only
            #               7: W+/- production
            #               8: pythia default (msel=1)
            #               9: Z0 production
i_stage=4   # (D=4) adj(40), =1, stops event after parini
            #                 2, ... parcas
            #                 3, ... hadronization with coal from parini
            #                 4, ... whole simulation
i_tune=0    # (D=0) i_tune. i.e. MSTP(5) in PYTHIA, tune number
##------------------------------------------------------------------------------
# Setups for A-loop.
prob_Delta_decay=0.9 # (D=0.9) decpro, Delta decay probability in A-loop

#*******************************************************************************
# Setups of the parton initialization (parini).
k_pythia=2.5    # (D=1.5) adj(10), K factor multiplying the differential 
                #                    cross sections for hard parton-parton 
                #                    process in PYTHIA.
sig_kPerp=2     # (D=2.) adj(39), width of parton promordial transverse 
                #                 momentum k_perp inside hadron.

#*******************************************************************************
# Setups of the parton cascade/rescattering (parcas).
k_parcas=1      # (D=1.) adj(1), K factor multiplying on the differential 
                #                 cross-section in parton cascade.
                #                =0, without parton cascade.

#*******************************************************************************
# Setups of the hadronization (sfm/coal).
i_hadronization=0   # (D=0) adj(12), model for hadronization
                    #                =0, string fragmentation (sfm)
                    #                =1, Monte-Carlo coalescence (coal)
                    #                =2, coal with gluon splitting and 
                    #                    deexcitation before parcas
                    #                =3, coal with gluon splitting but 
                    #                    without deexcitation
i_fragment_function=0   # (D=0) adj(29), fragmentation function:
                        #         i.e. MSTJ(11). Choice of longitudinal
                        #         fragmentation function, i.e. how large a 
                        #         fraction of the energy available a 
                        #         newly-created hadron/qqbar takes:
                        #         =1: Lund symmetric fragmentation function
                        #         =2: Field–Feynman(FF) + Peterson/SLAC(P/S)
                        #         =3: LUND + Peterson/SLAC
                        #         =4: LUND + Bowler, default in PYTHIA
                        #         =5: as = 4, but interpolate for c and b
                        #         =0: as = 4 for sfm / random z for coal
p_perp_min=1.9  # (D=1.9) adj1(9), PARP(81) in PYTHIA, effective minimum 
                #         transverse momentum p_erp_min of multiple 
                #         interactions if MSTP(82)=1.
p_perp0=2.5  # (D=2.) parp82, PARP(82) in PYTHIA, regularization scale p_erp_0 
             #        of the transverse-momentum spectrum for multiple 
             #        interactions with MSTP(82) >= 2.
pT_width=0.45  # (D=0.36) adj(34), i.e. PARJ(21) in PYTHIA, width of pT sampling
##------------------------------------------------------------------------------
## String fragmentation parameter (sfm).
# i_fragment_function=4 # (D=4)
a_lund=0.3  # (D=0.3)  adj(6), alpha in the LUND string fragmentation function.
b_lund=0.1 # (D=0.58) adj(7), beta.
##------------------------------------------------------------------------------
## Coalescence parameter (Coal).
# i_fragment_function=0 # (D=0)
i_pT_samp=1 # (D=3) i_pT, pT sampling method of the daughter qqbar pair in coal:
            #         = 1, Gaussian px and py with width PARJ(21)
            #         = 2, Exponential px and py with width PARJ(21)
            #         = 3, Exponential pT with width PARJ(21)
            #         = 4, random pT from mother
            #         = 5, random px and random py from mother, different factor
            #         = 6, random (px and py) from mother, the same factor
            #         = 7, random (px and py) from mother, the same random 
            #              factor as z which related to adj1(29).
            #         = 8, random pT from mother, the same random 
            #              factor as z which related to adj1(29).
i_pT_max=0  # (D=0) i_pT_max, with/without max value constraint in pT sampling.
eDeex=2 # (D=2.) adj1(17), deexcitation threshold energy.
i_phase_constraint=0    # (D=0) adj(21), phase-space constraint in coalescence.
                        #                =0, without
                        #                =1, with complete phase-space
                        #                =2, with position-space only
                        #                =3, with momentum-space only

#*******************************************************************************
# Setups of the hadron cascade/rescattering (hadcas).
i_hadcas=1  # (D=1) kjp21, = 0, without hadron cascade; 1, with
ratio_xSec_inel_over_tot=0.85   # (D=0.85) x_ratio, ratio of inela. cross 
                                #              section to total cross section 
                                #              of hadron-hadron collisions in 
                                #              A-loop and hadcas.
                                #              = 0.85 for B- & C-loop
                                #              = automatically calculated at 
                                #                E_CMS < 3 GeV in A-loop


################################################################################
#   The above setups are enough to run program normally.                       #
#   The following statements are optional. Usually, they are good enough.      #
################################################################################


#*******************************************************************************
#  The initial geometry and effects.
b_samp=2    # (D=2) psno, b parameter sampling method for [b_min,b_max].
            #             =0, fixed b
            #              1, systematic sampling
            #              2, random sampling
i_overlap=1 # (D=1) adj(30), initial nuclon distribution.
            #                =0, without more requirements
            #                =1, distributes the participant nucleons in 
            #                    the overlapping region forcely
i_shadow=1  # (D=1) adj(5), choice of nuclear shadowing for partons.
            #               =0, without nuclear shadowing,
            #               =1, Wang's nuclear shadowing (PLB 527(2002)85).
i_cme=1     # (D=1) adj(23), choice of chiral magnetic effect(CME) for partons.
            #                =0, without CME-induced charge separation mechanism
            #                =1, with
parecc=0.       # (D=0.) parecc, parameter converting initial spatial space 
                #               eccentricity to final momentum space ellipticity
pT_perturb=0.   # (D=0.) smadel, small perturbation of pT ellipse from circle 
                #                   when generates transverse momentum 
                #                   according to a Gaussian distribution.
##------------------------------------------------------------------------------
# Max simulated volum.
volum_factor=200    # (D=200) para10, factor of effective rescattering region size 
                #                 in the partonic and hadronic rescattering. 
                #                 ( r_max = para10 * MAX(r_A,r_B) )
ex_volum_factor=1   # (D=1) adj1(28), extra factor of effective rescattering region 
                #                 size in the hadronic rescattering.
                #                 ( r_max = para10 * MAX(r_A,r_B) * adj1(28) )

#*******************************************************************************
# Processes in parcas.
i_inelastic=0   # (D=0) iparres, consider inelastic processes or not
                #                =0, elastic processes 1, 2, 3, 5, 8, and 9 only
                #                =1, elastic + inelastic
i_inel_proc=7   # (D=7) i_inel_proc, when i_inelastic=1 (iparres=1),
                #                  =6, with all the inelastic processes 4, 6, 7
                #                  =7, only the inelastic process 7
i_t_shower=0    # (D=0) i_time_shower, when i_inelastic=1 (iparres=1),
                #                =0, without final state time-like parton shower
                #                =1, with

#*******************************************************************************
# Options in string fragmentation.
i_string_tension=4  # (D=4) kjp22, selects the form of string rension in sfm.
                    #              =1, variable single string tension
                    #               2, variable multiple string tension
                    #               3, variable (single+multiple) string tension
                    #               4, default constant string tension
cp0=1.              # (D=1.) cp0, cr0, parameters in parameterization 
cr0=0.2             # (D=0.2)            of multiple string effect
seco=0.05           # (D-0.05) seco, correction of PARJ(1) in popcorn mechanism
parj4=0.05          # (D=0.05) dparj4, PARJ(4) in PYTHIA

#*******************************************************************************
# Options in coalescence.
i_deex=3    # (D=3) i_deex, the deexcitation mode used in coal
            #               =1, light-cone variable mode
            #               =2, energy mode
            #               =3, light-cone variable mode, local pT compensation,
            #                   and sample z for qqbar.
            #               =4, light-cone variable mode, local pT compensation,
            #                   and sample z for q0-q or q0-qbar.
i_deex_gen=0    # (D=0) i_deex_gen, the deexcitation generation of newly qqbar
                #         =0, means no deexcitation for any newly produced qqbar
                #         =1, means just do deexcitation for the directly 
                #             proudced qqbar pairs (1-st daughters) from 
                #             original mother quarks (Orig mothers)
                #         =2, means do deexcitation for "1-st daughters" from 
                #             "Orig mothers" and the subsequent qqbar pairs 
                #             produced from "1-st daughters". (2-nd daughters)
                #         ...
                #         =999, always do deexcitation for newly produced qqbar
n_deexc_time_per_q=999 # (D=999) adj(16), the number of deexcitation per q/qbar.
i_mass=3    # (D=3) i_mass, mass definetion of quarks used in "break_f"
            #               =1, kinematical mass
            #               =2, current algebra mass
            #               =3, constituent mass (no top mass defined)
prob_ratio=(1 1 0.3 0 0 0)   # (D= 1 1 0.3 0 0) Probability ratio of 
                             #  u:d:s:c:b:t for qqbar sampling.

#*******************************************************************************
# Time accuracy.
dt_parini=0.00001   # (D=0.00001) ddt, min distinguishble collision time used in 
                    #                    the partonic initiation
dt_parcas=0.03      # (D=0.03) adj(19), time accuracy used in the parton cascade
dt_hadcas=0.1       # (D=0.1) adj(11), time accuracy used in the hadron cascade

#*******************************************************************************
# Statistical histogram bins in analy. The "asd(i)"" in usu.dat.
bin_1=0.35   # (D=0.35) y bin
bin_2=0.25   # (D=0.5)  pT bin (inv. pT)
bin_3=0.35   # (D=0.35) eta bin
bin_4=0.25   # (D=0.3)  mT bin
bin_5=25     # (D=25)   multiplicity bin
bin_6=0.25   # (D=0.5)  pT bin (dN/dpT)
i_y_or_eta=1    # (D=1) yOrEta, selects y or eta in partial-space statistics.
                #       = 0 , y
                #       = 1 , eta
# For 20 particles, use the same cuts.
pT_low=0.2  # pT lower cut
pT_upp=5.   # pT upper cut
rap_low=0.2    # y/eta lower cut
rap_upp=1.4     # y/eta upper cut

################################################################################
# The following lines do not need to change. Normally automatical calculation. #
################################################################################


#*******************************************************************************
# Total number of events in all program-runnings.
((tot_eve=n_run*n_eve))

# The "ipden" and "itden" in PACIAE for different systems.
i_proj=1
i_targ=1
if [[ "${na_proj}" = "1" ]]; then      # For N.
    i_proj=0
elif [[ "${na_proj}" > "1" ]]; then    # For A.
    i_proj=1
fi
if [[ "${na_targ}" = "1" ]]; then      # For N.
    i_targ=0
elif [[ "${na_targ}" > "1" ]]; then    # For A.
    i_targ=1
fi
if [[ "${i_proj}" = "0" && "${i_targ}" = "0" ]]; then
    b_samp=0
fi
# If others, please select them manully.
# i_proj=1  # ipden =0,  if projectile is nucleon (anti-nucleon)
            #       =1,  if nucleus
            #       =2,  for e+e-
            #       =11, if projectile is e- (e+)
            #       =12, if projectile is nu_e (nu_ebar)
            #       =13, if projectile is mu- (mu+)
            #       =14, if projectile is nu_mu (nu_mubar)
            #       =15, if projectile is tau- (tau+)
            #       =16, if projectile is nu_tau (nu_taubar)
# i_targ=1  # itden =0,  if target is nucleon (anti-nucleon)
            #       =1,  if nucleus
            #       =2,  for e+e-
            #...
            # for eA, nu_eA, etc.:
            # e^-A:     nap=1, nzp=-1, ipden=11, itden=1,
            # e^+A:     nap=1, nzp=1,  ipden=11, itden=1,
            # nu_eA:    nap=1, nzp=-1, ipden=12, itden=1,
            # nu_ebarA: nap=1, nzp=1,  ipden=12, itden=1,
####################      Variable's alias for usu.dat      ####################
################################################################################
################################################################################


# More options could be modified on the following "usu.dat" directly.
################################################################################
################################################################################
####################               usu.dat              ########################
echo "${n_eve},${n_out},${n_osc}                  ! neve,nout,nosc"     > usu.dat
echo "${na_proj},${nz_proj},${na_targ},${nz_targ} ! nap,nzp,nat,nzt"    >> usu.dat
echo "${dt_parini},${ratio_xSec_inel_over_tot},${b_min},${b_max},10      ! ddt,x_ratio,bmin,bmax,nmax"      >> usu.dat
echo "${i_hadcas},${i_frame},1.2,${volum_factor},1     ! kjp21,ifram,para7,para10,kjp20"    >> usu.dat
echo "3.1416,${i_proj},${i_targ}                  ! pio,ipden,itden"    >> usu.dat
echo "20,6,2                                 ! ispmax,isdmax,iflmax"    >> usu.dat
echo "0,0,0,0,0,0,0,0,0,0   ! KF code of particles not to decay"        >> usu.dat
echo "0,0,0,0,0,0,0,0,0,0   ! KF code of particles not to decay"        >> usu.dat
echo "211,-211,321,-321,3122,-3122,3312,-3312,3334,-3334   ! KF code of particles to be analyzed online"    >> usu.dat
echo "2212,-2212,2112,-2112,3212,-3212,3112,3222,333,443   ! KF code of particles to be analyzed online"    >> usu.dat
echo "${bin_1},${bin_2},${bin_3},${bin_4},${bin_5},${bin_6}  ! asd(i), i=1,6"           >> usu.dat
echo "${rap_low},${rap_upp}                        ! afl(j=1,i=1,1), afl(j=1,i=1,2) "   >> usu.dat
echo "${pT_low},${pT_upp}                      ! afl(j=1,i=2,1), afl(j=1,i=2,2) "       >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "${rap_low},${rap_upp}"                                                            >> usu.dat
echo "${pT_low},${pT_upp}"                                                              >> usu.dat
echo "2.7,${i_y_or_eta},${energy}               ! parp21,yOrEta,win"    >> usu.dat
echo "0.,0.5,1,1,${i_channel}       ! ttaup,taujp,iabsb,iabsm,nchan"    >> usu.dat
echo "7.2,4.,${b_samp},40.,20.,0.,0.1 ! para13,para14,psno,para15,para16,ajpsi,vneum"   >> usu.dat
echo "40,40,15,10                 ! para1_1,para1_2,para2,para4" >> usu.dat
echo "${i_deex},${i_deex_gen},${i_pT_samp},${i_pT_max},0.77,0.05,0.005,${p_perp0},${i_tune}   ! i_deex,i_deex_gen,i_pT,i_pT_max,a_FF,a_PS_c,a_PS_b,parp82,i_tune" >> usu.dat
echo "1,${i_inel_proc},${i_t_shower},${i_sim_mode},${prob_Delta_decay},2        ! mstu21,i_inel_proc,i_time_shower,,iMode,decpro,itorw" >> usu.dat
echo "${k_parcas},0.47,0.4,1000,${i_shadow},${a_lund},${b_lund},4,${p_perp_min},${k_pythia}    ! adj1(1)- adj1(10)  " >> usu.dat
echo "${dt_hadcas},${i_hadronization},30,45,1.,${n_deexc_time_per_q},${eDeex},0,${dt_parcas},1          ! adj1(11)- adj1(20) " >> usu.dat
echo "${i_phase_constraint},4.,${i_cme},0.15,0.4,${i_random},800000.,${ex_volum_factor},${i_fragment_function},${i_overlap} ! adj1(21)- adj1(30) " >> usu.dat
echo "0.1,0.3,0.4,${pT_width},1,0,100.,3.,${sig_kPerp},${i_stage}    ! adj1(31)- adj1(40) " >> usu.dat
echo "${i_string_tension},2,2,0.025,0   ! kjp22,kjp23,kjp24,parp78,mstptj" >> usu.dat
echo "${parecc},${i_inelastic},${pT_perturb},${parj4},${cp0},${cr0},${seco}   ! parecc,iparres,smadel,dparj4,cp0,cr0,seco" >> usu.dat
echo "${i_mass},${prob_ratio[0]},${prob_ratio[1]},${prob_ratio[2]},${prob_ratio[3]},${prob_ratio[4]},${prob_ratio[5]} ! i_mass, prob_ratio_q"   >> usu.dat
####################               usu.dat              ########################
################################################################################
####################                                    ########################
####################        Annotation of usu.dat       ########################
####################                                    ########################
################################################################################
echo ""                                                                                     >> usu.dat
echo "######################        Annotation of usu.dat         ####################"     >> usu.dat
echo "# neve,nout,nosc (D=xxx, xxx/10, 0)"                                                  >> usu.dat
echo "# neve: events number to be generated"                                                >> usu.dat
echo "# nout: output the event per nout events"                                             >> usu.dat
echo "# nosc: OSCAR standard output (oscar.out)"                                            >> usu.dat
echo "#       = 0 : no OSCAR output, see subroutine oscar"                                  >> usu.dat
echo "#       = 1 : OSCAR1997A (final_id_p_x, just PACIAE final output)"                    >> usu.dat
echo "#       = 2 : OSCAR1999A (full_event_history)"                                        >> usu.dat
echo "#       = 3 : OSCAR2013A (full_event_history, dummy now)"                             >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# nap,nzp,nat,nzt (D=197, 79, 197, 79 or 208, 82, 208, 82)"                           >> usu.dat
echo "# for NN, NA(AN), AA, etc."                                                           >> usu.dat
echo "#  nap(nzp): nucleons (protons) number of projectile"                                 >> usu.dat
echo "#  nat(nzt): nucleons (protons) number of target"                                     >> usu.dat
echo "#            for NN along with ipden=itden=0"                                         >> usu.dat
echo "#             p+p   : 1, 1, 1, 1 ;"                                                   >> usu.dat
echo "#             p+pbar: 1, 1, 1,-1 ;"                                                   >> usu.dat
echo "#             pbar+p: 1,-1, 1, 1 ;"                                                   >> usu.dat
echo "#             p+n   : 1, 1, 1, 0 ;"                                                   >> usu.dat
echo "#             n+p   : 1, 0, 1, 1 ;"                                                   >> usu.dat
echo "#             n+n   : 1, 0, 1, 0 ;"                                                   >> usu.dat
echo "#            for pA(Ap) along with ipden=0, itden=1 (ipden=1, itden=0)"               >> usu.dat
echo "#             p+Pb  : 1, 1, 208, 82;"                                                 >> usu.dat
echo "#             Pb+p  : 208, 82, 1, 1;"                                                 >> usu.dat
echo "#            for A+B along with ipden=itden=1"                                        >> usu.dat
echo "#             Au+Au: 197, 79, 197, 79; Pb+Pb: 208, 82, 208, 82;"                      >> usu.dat
echo "#             Xe+Xe: 129, 54, 129, 54; U + U: 238, 92, 238, 92;"                      >> usu.dat
echo "#             Ag+Ag: 108, 47, 108, 47; Cu+Cu:  63, 29,  63, 29;"                      >> usu.dat
echo "#             Ru+Ru:  96, 44,  96, 44; Zr+Zr:  96, 40,  96, 40."                      >> usu.dat
echo "# for eA, nu_eA, etc."                                                                >> usu.dat
echo "#  e^-A:     nap=1, nzp=-1, ipden=11, itden=1,"                                       >> usu.dat
echo "#  e^+A:     nap=1, nzp= 1, ipden=11, itden=1,"                                       >> usu.dat
echo "#  nu_eA:    nap=1, nzp=-1, ipden=12, itden=1,"                                       >> usu.dat
echo "#  nu_ebarA: nap=1, nzp= 1, ipden=12, itden=1."                                       >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# ddt,x_ratio,bmin,bmax,nmax (D=0.00001, 0.85, xxx, xxx, 10)"                         >> usu.dat
echo "#  ddt: minimum distinguishble collision time interval used in "                      >> usu.dat
echo "#       partonic initiation in parini.f"                                              >> usu.dat
echo "#  x_ratio: param(6), ratio of inel. cross section to total cross "                   >> usu.dat
echo "#           section of hadron-hadron scattering, automatically "                      >> usu.dat
echo "#           calculated at E_CMS < 3 GeV in A-loop"                                    >> usu.dat
echo "#  bmin: minimum impact parameters, "                                                 >> usu.dat
echo "#  bmax: maximum impact parameters,"                                                  >> usu.dat
echo "#  nmax: the number of intervals segmented in [bmin,bmax] when psno=1"                >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# kjp21,ifram,para7,para10,kjp20 (D=1, 1, 1.2, 200, 1)"                               >> usu.dat
echo "#  kjp21: =0, without hadron rescattering"                                            >> usu.dat
echo "#         =1, with hadron rescattering"                                               >> usu.dat
echo "#  ifram: choice collision system type"                                               >> usu.dat
echo "#         =0, fixed target"                                                           >> usu.dat
echo "#         =1, collider"                                                               >> usu.dat
echo "#  para7: proper formation time in rest-frame of particle"                            >> usu.dat
echo "#  para10: largest allowed size of partonic (hadronic) rescattering"                  >> usu.dat
echo "#          region which is product of para10 and target radius"                       >> usu.dat
echo "#  kjp20: choice the cross sections in hadron rescattering (hadcas.f)"                >> usu.dat
echo "#         =1, constant cross sections "                                               >> usu.dat
echo "#         =0, energy dependent cross sections"                                        >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# pio,ipden,itden (D=3.1416, 1, 1 for A+B)"                                           >> usu.dat
echo "#  pio: pi=3.1416"                                                                    >> usu.dat
echo "#  ipden: =0, if projectile is nucleon (anti-nucleon)"                                >> usu.dat
echo "#         =1, if projectile is nucleus"                                               >> usu.dat
echo "#         =2, for e+e-"                                                               >> usu.dat
echo "#         =11, if projectile is e- (e+)"                                              >> usu.dat
echo "#         =12, if projectile is nu_e (nu_ebar)"                                       >> usu.dat
echo "#         =13, if projectile is mu- (mu+)"                                            >> usu.dat
echo "#         =14, if projectile is nu_mu (nu_mubar)"                                     >> usu.dat
echo "#         =15, if projectile is tau- (tau+)"                                          >> usu.dat
echo "#         =16, if projectile is nu_tau (nu_taubar)"                                   >> usu.dat
echo "#  itden: =0, if target is nucleon (anti-nucleon)"                                    >> usu.dat
echo "#         =1, if projectile is nucleus"                                               >> usu.dat
echo "#         =2, for e+e-"                                                               >> usu.dat
echo "#         ..."                                                                        >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# ispmax,isdmax,iflmax   (D=20, 5, 2)"                                                >> usu.dat
echo "#  ispmax: maximum # of different particle pieces to be considered"                   >> usu.dat
echo "#  isdmax: maximum # of different distributions to be considered"                     >> usu.dat
echo "#  iflmax: maximum # of windows to be set, =0 means no window at all"                 >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# ispkf(i,i=1,ispmax):"                                                               >> usu.dat
echo "# KF code: particle code used in PYTHIA and PACIAE, "                                 >> usu.dat
echo "#          (list at the end and see detail in reference: arXiv:hep-ph/0603175)"       >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# asd(i=1,isdmax): interval of the i-th distribution"                                 >> usu.dat
echo "#  for pp, pbarp, pA(Ap), AB etc."                                                    >> usu.dat
echo "#      (D=0.35, 0.5, 0.35, 0.3, 25)"                                                  >> usu.dat
echo "#      i=1: rapidity distribution (dN/dy v.s. y)"                                     >> usu.dat
echo "#       =2: invariant transverse monmentum distribution (1/pT*dN/dpT v.s. pT)"        >> usu.dat
echo "#       =3: pesudorapidity distribution (dN/deta v.s. eta)"                           >> usu.dat
echo "#       =4: transeverse mass distribution (1/mT*dN/dmT v.s. mT)"                      >> usu.dat
echo "#       =5: event-wise multiplicity distribution"                                     >> usu.dat
echo "#  for ep, nu_ep, etc."                                                               >> usu.dat
echo "#      i=1: Q^2=-q^2 (fq2 in code) distribution"                                      >> usu.dat
echo "#       =2: W^2 (w21) distribution"                                                   >> usu.dat
echo "#       =3: y (yyl) distribution"                                                     >> usu.dat
echo "#       =4: p_h (pph) distribution"                                                   >> usu.dat
echo "#       =5: z (zl) distribution"                                                      >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# afl(j,i,1): lower-boundary of i-th window for j-th particle"                        >> usu.dat
echo "# afl(j,i,2): upper-boundary  of i-th window for j-th particle"                       >> usu.dat
echo "#  for pp, pbarp, pA(Ap), AB etc."                                                    >> usu.dat
echo "#      i=1, rapidity/pesudorapidity window (D= -1.,1. )"                              >> usu.dat
echo "#       =2, transverse monmentum           (D= 0.,50. )"                              >> usu.dat
echo "#  for ep, nu_ep, etc."                                                               >> usu.dat
echo "#      i=1, Q^2=-q^2 window"                                                          >> usu.dat
echo "#       =2, W^2"                                                                      >> usu.dat
echo "#       =3, y"                                                                        >> usu.dat
echo "#       =4, p_h (haron momentum)"                                                     >> usu.dat
echo "#       =5: z"                                                                        >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# parp21,yOrEta,win (D=2.7, 1, xxx)"                                                  >> usu.dat
echo "#  parp21: lowest CM energy running 'pythia'"                                         >> usu.dat
echo "#  yOrEta: select y or eta in partial phase-space statistics (analy.f)"               >> usu.dat
echo "#          = 0 , y"                                                                   >> usu.dat
echo "#          = 1 , eta"                                                                 >> usu.dat
echo "#  win = cms energy if ifram=1 (collider)"                                            >> usu.dat
echo "#      = incident momentum if ifram=0 (fixed target)"                                 >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# ttaup,taujp,iabsb,iabsm,nchan (D=0., 0.5, 1, 1, 8)"                                 >> usu.dat
echo "#  ttaup: proper formation time of particles generated in hadronic rescattering"      >> usu.dat
echo "#  taujp: formation time of J/psi"                                                    >> usu.dat
echo "#  iabsb: =0, without J/psi (psi') + baryon"                                          >> usu.dat
echo "#         =1, with J/psi (psi') + baryon"                                             >> usu.dat
echo "#  iabsm: =0, without J/psi (psi') + meson"                                           >> usu.dat
echo "#         =1, with J/psi (psi') + meson"                                              >> usu.dat
echo "#  nchan: to choose which subset of parton-parton subprocesses to include in"         >> usu.dat
echo "#         the generration"                                                            >> usu.dat
echo "#         =0, inelastic (INEL)"                                                       >> usu.dat
echo "#         =1, Non-Single Difractive (NSD)"                                            >> usu.dat
echo "#         =2, Drell-Yan process"                                                      >> usu.dat
echo "#         =3, J/psi production"                                                       >> usu.dat
echo "#         =4, heavy-flavor production"                                                >> usu.dat
echo "#         =5, direct photon"                                                          >> usu.dat
echo "#         =6, soft only"                                                              >> usu.dat
echo "#         =7, W+/- production"                                                        >> usu.dat
echo "#         =8: PYTHIA default (msel=1)"                                                >> usu.dat
echo "#         =9: Z0 production"                                                          >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# para13,para14,psno,para15,para16,ajpsi,vneum (D=7.2, 4., 2, 40., 20., 0, 0.1)"      >> usu.dat
echo "#  para13: total cross-section of J/psi + n"                                          >> usu.dat
echo "#  para14: total cross-section of J/psi + meson"                                      >> usu.dat
echo "#  psno: =0 fixed impact parameter"                                                   >> usu.dat
echo "#        =1 impact parameter is sampled by systematic sampling method"                >> usu.dat
echo "#        =2 randomly sampled impact parameter"                                        >> usu.dat
echo "#  para15: total cross-section of psi' + n"                                           >> usu.dat
echo "#  para16: total cross-section of psi' + meson"                                       >> usu.dat
echo "#  ajpsi: not used now"                                                               >> usu.dat
echo "#  vneum: relevant to average binary collision number, now it is recalculated"        >> usu.dat
echo "#         in program"                                                                 >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# para1_1,para1_2,para2,para4 (D= 40, 40, 25, 10 for Au+Au; 70, 70... for Pb+Pb)"     >> usu.dat
echo "#  para1_1: nn total cross section used in parton initiation"                         >> usu.dat
echo "#  para1_2: nn total cross section used in hadron cascade"                            >> usu.dat
echo "#  para2: total cross-section of pi-nucleon"                                          >> usu.dat
echo "#  para4: total cross-section of pi-pi"                                               >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# i_deex,i_deex_gen,i_pT,i_pT_max,a_FF,aPS_c,aPS_b,parp82,i_tune"                     >> usu.dat
echo "#  (D=3, 0, 1, 0, 0.77, 0.05, 0.005, 2., 0)"                                          >> usu.dat
echo "#  i_deex: the deexcitation mode used in coal"                                        >> usu.dat
echo "#          = 1, light-cone variable mode"                                             >> usu.dat
echo "#          = 2, energy mode"                                                          >> usu.dat
echo "#          = 3, light-cone variable mode, with local pT compensation "                >> usu.dat
echo "#               and sampling z for qqbar"                                             >> usu.dat
echo "#          = 4, energy mode, with local pT compensation "                             >> usu.dat
echo "#               and sampling z for qqbar"                                             >> usu.dat
echo "#  i_deex_gen: the deexcitation generation of newly produced qqbar in coal"           >> usu.dat
echo "#              = 0, means no deexcitation for any newly produced qqbar pairs"         >> usu.dat
echo "#              = 1, means just do deexcitation for the directly proudced qqbar "      >> usu.dat
echo "#                         pairs (1-st daughters) from original mother quarks "        >> usu.dat
echo "#                         (Orig mothers)"                                             >> usu.dat
echo "#              = 2, means do deexcitation for "1-st daughters" from "                 >> usu.dat
echo "#                         "Orig mothers" and the subsequent qqbar pairs "             >> usu.dat
echo "#                         produced from "1-st daughters". (2-nd daughters)"           >> usu.dat
echo "#                ..."                                                                 >> usu.dat
echo "#              = 999, always do deexcitation for newly produced qqbar"                >> usu.dat
echo "#  i_pT: the pT sampling method of the daughter qqbar pair in coal"                   >> usu.dat
echo "#        = 1, Gaussian px and py with width PARJ(21)"                                 >> usu.dat
echo "#        = 2, Exponential px and py with width PARJ(21)"                              >> usu.dat
echo "#        = 3, Exponential pT with width PARJ(21)"                                     >> usu.dat
echo "#        = 4, random pT from mother"                                                  >> usu.dat
echo "#        = 5, random px and random py from mother, different random factors"          >> usu.dat
echo "#        = 6, random (px and py) from mother, the same random factor"                 >> usu.dat
echo "#        = 7, random (px and py) from mother, the same random factor as "             >> usu.dat
echo "#             z which related to adj1(29)"                                            >> usu.dat
echo "#  i_pT_max: whether the sampled pT in coal deexitation is greater than the "         >> usu.dat
echo "#            the mother quark or not."                                                >> usu.dat
echo "#  a_FF: parameter for light hadron in Field-Feynman function, i.e. u, d, and s "     >> usu.dat
echo "#        hadron --PARJ(51), (52), and (53)--, set them equal"                         >> usu.dat
echo "#  aPS_c: -PARJ(54), parameter for charm-hadron in Petersono/SLAC, note the minus"    >> usu.dat
echo "#  aPS_b: -PARJ(55), parameter for bottom-hadron in P/S function"                     >> usu.dat
echo "#  parp82: PARP(82) in PYTHIA, regularization scale p_erp_0 of the "                  >> usu.dat
echo "#          transverse-momentum spectrum for multiple interactions with "              >> usu.dat
echo "#          MSTP(82) >= 2."                                                            >> usu.dat
echo "#  i_tune: MSTP(5), tune number of PYTHIA. = 350, Perugia 2011 tune."                 >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# mstu21,i_inel_proc,i_time_shower,iMode,decpro,itorw (D=1, 7, 0, 3, 0.9, 2)"         >> usu.dat
echo "#  mstu21: parameter mstu(21) in PYTHIA"                                              >> usu.dat
echo "#  i_inel_proc: = 6, with inelastic processes 4, 6, and 7 if iparres=1 (in parcas.f)" >> usu.dat
echo "#               = 7, with inelastic process 7 only if iparres=1 (in parcas.f)"        >> usu.dat
echo "#  i_time_shower: = 0, w/o final state time-like parton shower if iparres=1"          >> usu.dat
echo "#                 = 1, w/ final state time-like parton shower if iparres=1"           >> usu.dat
echo "#  iMode: =1, low energy simulation A-loop"                                           >> usu.dat
echo "#         =2, PYTHIA-like simulation B-loop"                                          >> usu.dat
echo "#         =3, PACIAE simulation C-loop"                                               >> usu.dat
echo "#  decpro: is Delta decay probability in low energy A-loop"                           >> usu.dat
echo "#  itorw: =1, executing pyevnt"                                                       >> usu.dat
echo "#         =2, executing pyevnw"                                                       >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# adj1(i), i=1,40: switches and/or parameters"                                        >> usu.dat
echo "# --------------------------------------------------------------------------------"   >> usu.dat
echo "#     D= 1-10 :  1., 0.47, 0.4, 1000,   1, 0.3,    0.58,  4,  1.9, 1.5"               >> usu.dat
echo "#        11-20: 0.1,    0,  30,   45,  1.,  1.,      2.,  0, 0.03, 1"                 >> usu.dat
echo "#        21-30:   0,   4.,   1, 0.15, 0.4,   1, 800000., 1.,    4, 0"                 >> usu.dat
echo "#        31-40: 0.1,  0.3, 0.4, 0.36,   1,   0,    100., 3.,   2., 4."                >> usu.dat
echo "# --------------------------------------------------------------------------------"   >> usu.dat
echo "#     i=1: K factor in parton rescattering."                                          >> usu.dat
echo "#       2: alpha_s, effective coupling constant in parton rescattering."              >> usu.dat
echo "#       3: mu^2 (tcut in program), the regulation factor introduced in "              >> usu.dat
echo "#          the parton-parton differential cross section (parcas)."                    >> usu.dat
echo "#       4: idw, the number of intervals in the numerical integration."                >> usu.dat
echo "#       5: =1: with Wang's nuclear shadowing (PLB 527(2002)85),"                      >> usu.dat
echo "#          =0: without nuclear shadowing."                                            >> usu.dat
echo "#       6: alpha in the LUND string fragmentation function (PARJ(41) in PYTHIA)."     >> usu.dat
echo "#       7: beta in the LUND string fragmentation function (PARJ(42) in PYTHIA)."      >> usu.dat
echo "#       8: MSTP(82) in PYTHIA."                                                       >> usu.dat
echo "#       9: PARP(81) in PYTHIA."                                                       >> usu.dat
echo "#       10: K factor (PARP(31) in PYTHIA)."                                           >> usu.dat
echo "#       11: time accuracy used in the hadron rescattering."                           >> usu.dat
echo "#       12: model for hadronization:"                                                 >> usu.dat
echo "#           =0, string fragmentation,"                                                >> usu.dat
echo "#           =1, Monte Carlo coalescence model"                                        >> usu.dat
echo "#       13: dimension of meson table considered in coalescence model."                >> usu.dat
echo "#       14: dimension of baryon table considered coalescence model."                  >> usu.dat
echo "#       15: string tension of qqbar simple string."                                   >> usu.dat
echo "#       16: number of loops in the deexcitation of energetic quark in the "           >> usu.dat
echo "#           Monte Carlo coalescence model."                                           >> usu.dat
echo "#       17: the threshold energy in the deexcitation of energetic quark in "          >> usu.dat
echo "#           the Monte Carlo coalescence model."                                       >> usu.dat
echo "#       18: =0, rest partons hadronize by string fragmentation,"                      >> usu.dat
echo "#           =1, rest partons hadronize by coalescence."                               >> usu.dat
echo "#       19: time accuracy used in the parton rescattering."                           >> usu.dat
echo "#       20: the optional parton-parton cross section in the parton rescattering:"     >> usu.dat
echo "#           =0, LO pQCD parton-parton cross section,"                                 >> usu.dat
echo "#           =1, keeping only leading divergent terms in the LO pQCD parton-parton "   >> usu.dat
echo "#               cross section (B. Zhang),"                                            >> usu.dat
echo "#           =2, the same as 0 but flat scattering angle distribution is assumed,"     >> usu.dat
echo "#           =3, the same as 1 but flat scattering angle distribution is assumed."     >> usu.dat
echo "#       21: with or without phase space constraint in the Monte Carlo coalescence model:"     >> usu.dat
echo "#           =0, without phase space constraint,"                                      >> usu.dat
echo "#           =1, with complete phase space constraint,"                                >> usu.dat
echo "#           =2, with spatial phase space constraint only,"                            >> usu.dat
echo "#           =3, with momentum phase space constraint only."                           >> usu.dat
echo "#       22: critical value (D=4) of the product of radii in position and momentum "   >> usu.dat
echo "#           phase spaces."                                                            >> usu.dat
echo "#       23: switch for chiral magnetic effect (CME):"                                 >> usu.dat
echo "#           =0: CME off,"                                                             >> usu.dat
echo "#           =1: CME on."                                                              >> usu.dat
echo "#       24: the virtuality cut ('tl0' in program) in the time-like radiation in parton "      >> usu.dat
echo "#           rescattering."                                                            >> usu.dat
echo "#       25: Lambda_QCD in parton rescattering."                                       >> usu.dat
echo "#       26: selection of random number seed:"                                         >> usu.dat
echo "#           =0, default PYTHIA seed (19780503), can be used for debug,"               >> usu.dat
echo "#           =other, seed from the real-time clock."                                   >> usu.dat
echo "#       27: largest momentum allowed for produced particle."                          >> usu.dat
echo "#       28: concerned to the largest position allowed for produced particle in "      >> usu.dat
echo "#            hadcas, it will be recalculated in program running "                     >> usu.dat
echo "#           ( drmax=para10*dmax1(rnt,rnp) )."                                         >> usu.dat
echo "#       29: For sfm in PYTHIA, it is MSTJ(11). Choice of longitudinal "               >> usu.dat
echo "#            fragmentation function, i.e. how large a fraction of the energy "        >> usu.dat
echo "#            available a newly-created hadron takes:"                                 >> usu.dat
echo "#           =1: Lund symmetric fragmentation function, see PARJ(41) - PARJ(45),"      >> usu.dat
echo "#           =2: Field-Feynman + Peterson/SLAC, see PARJ(51) PARJ(59),"                >> usu.dat
echo "#           =3: Lund + Peterson/SLAC (light flavor + heavier),"                       >> usu.dat
echo "#           =4: default PYTHIA. Lund + Bowler,"                                       >> usu.dat
echo "#           =5: as = 4, but interpolate for c and b; see PARJ(46) and PARJ(47)."      >> usu.dat
echo "#           For coal, sampling daughter parton energy fraction z taking from "        >> usu.dat
echo "#            mother in 'funcz' for coal:"                                             >> usu.dat
echo "#           =1: by Lund string fragmentation function,"                               >> usu.dat
echo "#           =2: by Field-Feynmman fragmentation function,"                            >> usu.dat
echo "#           =3: by Peterson/SLAC fragmentation function,"                             >> usu.dat
echo "#           =4: null,"                                                                >> usu.dat
echo "#           =5: null."                                                                >> usu.dat
echo "#       30: =1, distribute the participant nucleons in overlapping region forcely,"   >> usu.dat
echo "#           =0, without more requirements."                                           >> usu.dat
echo "#       31: PARJ(1) in PYTHIA."                                                       >> usu.dat
echo "#       32: PARJ(2) in PYTHIA."                                                       >> usu.dat
echo "#       33: PARJ(3) in PYTHIA."                                                       >> usu.dat
echo "#       34: PARJ(21) in PYTHIA, width of px/py/pT sampling in PYPTDI/paptdi."         >> usu.dat
echo "#       35: MSTP(91) in PYTHIA, parton transverse momentum (k_perp) distribution "    >> usu.dat
echo "#           inside hadron:"                                                           >> usu.dat
echo "#           =1: Gaussian,"                                                            >> usu.dat
echo "#           =2: exponential."   >> usu.dat
echo "#       36: with or without phenomenological parton energy loss in parton rescattering:"      >> usu.dat
echo "#           =0, without,"                                                             >> usu.dat
echo "#           =1, with."                                                                >> usu.dat
echo "#       37: the coefficient in phenomenological parton energy loss."                  >> usu.dat
echo "#       38: p_T cut in phenomenological parton energy loss."                          >> usu.dat
echo "#       39: PARP(91) (D=2.), width of Gaussian parton k_perp distribution in hadron " >> usu.dat
echo "#           if MSTP(91)=1,"                                                           >> usu.dat
echo "#           PARP(92) (D=0.4), width of Exponential k_perp distribution in hadron "    >> usu.dat
echo "#           if MSTP(91)=2 ( ~ PARP(92)/SQRT(6) )."                                    >> usu.dat
echo "#       40: optional event stopping point"                                            >> usu.dat
echo "#           =1, after parton initiation,"                                             >> usu.dat
echo "#           =2, after parton rescattering,"                                           >> usu.dat
echo "#           =3, after hadronization with coalescence from parton initiation directly,"        >> usu.dat
echo "#           =4, after the whole simulation."                                          >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# kjp22,kjp23,kjp24,parp78,mstptj (D=4, 2, 2, 0.025, 0)"                              >> usu.dat
echo "#  kjp22: =1, variable single string tension and PARJ(1) etc."                        >> usu.dat
echo "#         =2, variable multiple string tension and PARJ(1) etc."                      >> usu.dat
echo "#         =3, variable (single+multiple) string tension and PARJ(1) etc."             >> usu.dat
echo "#         =4, default string tension and PARJ(1) etc."                                >> usu.dat
echo "#  kjp23: optional model for the calculation of participant nucleon echo (npart)"     >> usu.dat
echo "#         =1, geometric model"                                                        >> usu.dat
echo "#         =2, Glauber model"                                                          >> usu.dat
echo "#  kjp24: optional distribution in Glauber model"                                     >> usu.dat
echo "#         =1, sharp sphere"                                                           >> usu.dat
echo "#         =2, Woods-Saxon"                                                            >> usu.dat
echo "#  parp78: parameter controling amount of colour reconnection in final state"         >> usu.dat
echo "#  mstptj: =0, input MSTP(111) (MSTJ(1)) for pp, pA (AP), and AA (for e+e-) in"       >> usu.dat
echo "#              PACIAE simulation developed from partonic initial stage, "             >> usu.dat
echo "#              to partonic rescattering, hadronization, and to hadronic"              >> usu.dat
echo "#              rescttering stage"                                                     >> usu.dat
echo "#          =1, PYTHIA like simulation without partonic & hadronic "                   >> usu.dat
echo "#              rescatterings but with setting of kjp21=0"                             >> usu.dat
echo "#          It will be set automatically in program now."                              >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# parecc,iparres,smadel,dparj4,cp0,cr0,seco (D=0., 0, 0., 0.05, 1., 0.2, 0.05)"       >> usu.dat
echo "#  parecc: proportional factor between initial spatial space eccentricity and final " >> usu.dat
echo "#          momentum space ellipticity"                                                >> usu.dat
echo "#  iparres: =0 consider elastic parton-parton collisions only in parton rescattering" >> usu.dat
echo "#           =1 otherwise"                                                             >> usu.dat
echo "#  smadel: small perpurbation of ellipse from circle"                                 >> usu.dat
echo "#  dparj4: default PARJ(4)"                                                           >> usu.dat
echo "#  cp0,cr0: parameters in parameterization of multiple string effect"                 >> usu.dat
echo "#  seco: parameter in popcorn mechanism for correction of PARJ(1)"                    >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# i_mass, prob_ratio_q (D=3, 1, 1, 0.3, 0, 0, 0)"                                     >> usu.dat
echo "# i_mass: mass definetion of quarks used in \"break_f\""                              >> usu.dat
echo "# prob_ratio_q: probability ratio u-ubar:d-dbar:s-sbar:c-cbar:b-bbar:t-tbar."         >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "# cumulent sum of flavor generation probability, for an example: "                    >> usu.dat
echo "#  p_31=(1-amd/sms)/2., p_32=(1-amu/sms)/2., p_33=(1-ams/sms)/2.;"                    >> usu.dat
echo "#  amd (amu, ams): d (u,s) quark constituent mass, sms=amd+amu+ams;"                  >> usu.dat
echo "#  p_31+p_32+p33=1;"                                                                  >> usu.dat
echo "#  cumulent sum: csp_31=p_31, csp_32=p_31+p_32"                                       >> usu.dat
echo "#  u and d have same constituent mass and probability"                                >> usu.dat
echo "######################        Annotation of usu.dat         ####################"     >> usu.dat
echo "################################################################################"     >> usu.dat
echo ""                                                                                     >> usu.dat
echo ""                                                                                     >> usu.dat
echo ""                                                                                     >> usu.dat
echo "################################################################################"     >> usu.dat
echo "################################################################################"     >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "#                     List of KF codes in program"                                    >> usu.dat
echo "#"                                                                                    >> usu.dat
echo "#             1    d                            -1    dbar            "               >> usu.dat
echo "#             2    u                            -2    ubar            "               >> usu.dat
echo "#             3    s                            -3    sbar            "               >> usu.dat
echo "#             4    c                            -4    cbar            "               >> usu.dat
echo "#             5    b                            -5    bbar            "               >> usu.dat
echo "#             6    t                            -6    tbar            "               >> usu.dat
echo "#             7    b'                           -7    b'bar           "               >> usu.dat
echo "#             8    t'                           -8    t'bar           "               >> usu.dat
echo "#            11    e-                          -11    e+              "               >> usu.dat
echo "#            12    nu_e                        -12    nu_ebar         "               >> usu.dat
echo "#            13    mu-                         -13    mu+             "               >> usu.dat
echo "#            14    nu_mu                       -14    nu_mubar        "               >> usu.dat
echo "#            15    tau-                        -15    tau+            "               >> usu.dat
echo "#            16    nu_tau                      -16    nu_taubar       "               >> usu.dat
echo "#            17    tau'-                       -17    tau'+           "               >> usu.dat
echo "#            18    nu'_tau                     -18    nu'_taubar      "               >> usu.dat
echo "#            21    g               "                                                  >> usu.dat
echo "#            22    gamma           "                                                  >> usu.dat
echo "#            23    Z0              "                                                  >> usu.dat
echo "#            24    W+                          -24    W-              "               >> usu.dat
echo "#            25    h0              "                                                  >> usu.dat
echo "#            32    Z'0             "                                                  >> usu.dat
echo "#            33    Z\"0             "                                                 >> usu.dat
echo "#            34    W'+                         -34    W'-             "               >> usu.dat
echo "#            35    H0              "                                                  >> usu.dat
echo "#            36    A0              "                                                  >> usu.dat
echo "#            37    H+                          -37    H-              "               >> usu.dat
echo "#            39    Graviton        "                                                  >> usu.dat
echo "#            41    R0                          -41    Rbar0           "               >> usu.dat
echo "#            42    LQ_ue                       -42    LQ_uebar        "               >> usu.dat
echo "#          2101    ud_0                      -2101    ud_0bar         "               >> usu.dat
echo "#          3101    sd_0                      -3101    sd_0bar         "               >> usu.dat
echo "#          3201    su_0                      -3201    su_0bar         "               >> usu.dat
echo "#          4101    cd_0                      -4101    cd_0bar         "               >> usu.dat
echo "#          4201    cu_0                      -4201    cu_0bar         "               >> usu.dat
echo "#          4301    cs_0                      -4301    cs_0bar         "               >> usu.dat
echo "#          5101    bd_0                      -5101    bd_0bar         "               >> usu.dat
echo "#          5201    bu_0                      -5201    bu_0bar         "               >> usu.dat
echo "#          5301    bs_0                      -5301    bs_0bar         "               >> usu.dat
echo "#          5401    bc_0                      -5401    bc_0bar         "               >> usu.dat
echo "#          1103    dd_1                      -1103    dd_1bar         "               >> usu.dat
echo "#          2103    ud_1                      -2103    ud_1bar         "               >> usu.dat
echo "#          2203    uu_1                      -2203    uu_1bar         "               >> usu.dat
echo "#          3103    sd_1                      -3103    sd_1bar         "               >> usu.dat
echo "#          3203    su_1                      -3203    su_1bar         "               >> usu.dat
echo "#          3303    ss_1                      -3303    ss_1bar         "               >> usu.dat
echo "#          4103    cd_1                      -4103    cd_1bar         "               >> usu.dat
echo "#          4203    cu_1                      -4203    cu_1bar         "               >> usu.dat
echo "#          4303    cs_1                      -4303    cs_1bar         "               >> usu.dat
echo "#          4403    cc_1                      -4403    cc_1bar         "               >> usu.dat
echo "#          5103    bd_1                      -5103    bd_1bar         "               >> usu.dat
echo "#          5203    bu_1                      -5203    bu_1bar         "               >> usu.dat
echo "#          5303    bs_1                      -5303    bs_1bar         "               >> usu.dat
echo "#          5403    bc_1                      -5403    bc_1bar         "               >> usu.dat
echo "#          5503    bb_1                      -5503    bb_1bar         "               >> usu.dat
echo "#           111    pi0             "                                                  >> usu.dat
echo "#           211    pi+                        -211    pi-             "               >> usu.dat
echo "#           221    eta             "                                                  >> usu.dat
echo "#           311    K0                         -311    Kbar0           "               >> usu.dat
echo "#           130    K_L0            "                                                  >> usu.dat
echo "#           310    K_S0            "                                                  >> usu.dat
echo "#           321    K+                         -321    K-              "               >> usu.dat
echo "#           331    eta'            "                                                  >> usu.dat
echo "#           411    D+                         -411    D-              "               >> usu.dat
echo "#           421    D0                         -421    Dbar0           "               >> usu.dat
echo "#           431    D_s+                       -431    D_s-            "               >> usu.dat
echo "#           441    eta_c           "                                                  >> usu.dat
echo "#           511    B0                         -511    Bbar0           "               >> usu.dat
echo "#           521    B+                         -521    B-              "               >> usu.dat
echo "#           531    B_s0                       -531    B_sbar0         "               >> usu.dat
echo "#           541    B_c+                       -541    B_c-            "               >> usu.dat
echo "#           551    eta_b           "                                                  >> usu.dat
echo "#           113    rho0            "                                                  >> usu.dat
echo "#           213    rho+                       -213    rho-            "               >> usu.dat
echo "#           223    omega           "                                                  >> usu.dat
echo "#           313    K*0                        -313    K*bar0          "               >> usu.dat
echo "#           323    K*+                        -323    K*-             "               >> usu.dat
echo "#           333    phi             "                                                  >> usu.dat
echo "#           413    D*+                        -413    D*-             "               >> usu.dat
echo "#           423    D*0                        -423    D*bar0          "               >> usu.dat
echo "#           433    D*_s+                      -433    D*_s-           "               >> usu.dat
echo "#           443    J/psi           "                                                  >> usu.dat
echo "#           513    B*0                        -513    B*bar0          "               >> usu.dat
echo "#           523    B*+                        -523    B*-             "               >> usu.dat
echo "#           533    B*_s0                      -533    B*_sbar0        "               >> usu.dat
echo "#           543    B*_c+                      -543    B*_c-           "               >> usu.dat
echo "#           553    Upsilon         "                                                  >> usu.dat
echo "#         10113    b_10            "                                                  >> usu.dat
echo "#         10213    b_1+                     -10213    b_1-            "               >> usu.dat
echo "#         10223    h_1             "                                                  >> usu.dat
echo "#         10313    K_10                     -10313    K_1bar0         "               >> usu.dat
echo "#         10323    K_1+                     -10323    K_1-            "               >> usu.dat
echo "#         10333    h'_1            "                                                  >> usu.dat
echo "#         10413    D_1+                     -10413    D_1-            "               >> usu.dat
echo "#         10423    D_10                     -10423    D_1bar0         "               >> usu.dat
echo "#         10433    D_1s+                    -10433    D_1s-           "               >> usu.dat
echo "#         10443    h_1c            "                                                  >> usu.dat
echo "#         10513    B_10                     -10513    B_1bar0         "               >> usu.dat
echo "#         10523    B_1+                     -10523    B_1-            "               >> usu.dat
echo "#         10533    B_1s0                    -10533    B_1sbar0        "               >> usu.dat
echo "#         10543    B_1c+                    -10543    B_1c-           "               >> usu.dat
echo "#         10553    h_1b            "                                                  >> usu.dat
echo "#         10111    a_00            "                                                  >> usu.dat
echo "#         10211    a_0+                     -10211    a_0-            "               >> usu.dat
echo "#         10221    f_0             "                                                  >> usu.dat
echo "#         10311    K*_00                    -10311    K*_0bar0        "               >> usu.dat
echo "#         10321    K*_0+                    -10321    K*_0-           "               >> usu.dat
echo "#         10331    f'_0            "                                                  >> usu.dat
echo "#         10411    D*_0+                    -10411    D*_0-           "               >> usu.dat
echo "#         10421    D*_00                    -10421    D*_0bar0        "               >> usu.dat
echo "#         10431    D*_0s+                   -10431    D*_0s-          "               >> usu.dat
echo "#         10441    chi_0c          "                                                  >> usu.dat
echo "#         10511    B*_00                    -10511    B*_0bar0        "               >> usu.dat
echo "#         10521    B*_0+                    -10521    B*_0-           "               >> usu.dat
echo "#         10531    B*_0s0                   -10531    B*_0sbar0       "               >> usu.dat
echo "#         10541    B*_0c+                   -10541    B*_0c-          "               >> usu.dat
echo "#         10551    chi_0b          "                                                  >> usu.dat
echo "#         20113    a_10            "                                                  >> usu.dat
echo "#         20213    a_1+                     -20213    a_1-            "               >> usu.dat
echo "#         20223    f_1             "                                                  >> usu.dat
echo "#         20313    K*_10                    -20313    K*_1bar0        "               >> usu.dat
echo "#         20323    K*_1+                    -20323    K*_1-           "               >> usu.dat
echo "#         20333    f'_1            "                                                  >> usu.dat
echo "#         20413    D*_1+                    -20413    D*_1-           "               >> usu.dat
echo "#         20423    D*_10                    -20423    D*_1bar0        "               >> usu.dat
echo "#         20433    D*_1s+                   -20433    D*_1s-          "               >> usu.dat
echo "#         20443    chi_1c          "                                                  >> usu.dat
echo "#         20513    B*_10                    -20513    B*_1bar0        "               >> usu.dat
echo "#         20523    B*_1+                    -20523    B*_1-           "               >> usu.dat
echo "#         20533    B*_1s0                   -20533    B*_1sbar0       "               >> usu.dat
echo "#         20543    B*_1c+                   -20543    B*_1c-          "               >> usu.dat
echo "#         20553    chi_1b          "                                                  >> usu.dat
echo "#           115    a_20            "                                                  >> usu.dat
echo "#           215    a_2+                       -215    a_2-            "               >> usu.dat
echo "#           225    f_2             "                                                  >> usu.dat
echo "#           315    K*_20                      -315    K*_2bar0        "               >> usu.dat
echo "#           325    K*_2+                      -325    K*_2-           "               >> usu.dat
echo "#           335    f'_2            "                                                  >> usu.dat
echo "#           415    D*_2+                      -415    D*_2-           "               >> usu.dat
echo "#           425    D*_20                      -425    D*_2bar0        "               >> usu.dat
echo "#           435    D*_2s+                     -435    D*_2s-          "               >> usu.dat
echo "#           445    chi_2c          "                                                  >> usu.dat
echo "#           515    B*_20                      -515    B*_2bar0        "               >> usu.dat
echo "#           525    B*_2+                      -525    B*_2-           "               >> usu.dat
echo "#           535    B*_2s0                     -535    B*_2sbar0       "               >> usu.dat
echo "#           545    B*_2c+                     -545    B*_2c-          "               >> usu.dat
echo "#           555    chi_2b          "                                                  >> usu.dat
echo "#        100443    psi'            "                                                  >> usu.dat
echo "#        100553    Upsilon'        "                                                  >> usu.dat
echo "#          3122    Lambda0                   -3122    Lambdabar0      "               >> usu.dat
echo "#          4122    Lambda_c+                 -4122    Lambda_cbar-    "               >> usu.dat
echo "#          4132    Xi_c0                     -4132    Xi_cbar0        "               >> usu.dat
echo "#          4232    Xi_c+                     -4232    Xi_cbar-        "               >> usu.dat
echo "#          5122    Lambda_b0                 -5122    Lambda_bbar0    "               >> usu.dat
echo "#          5132    Xi_b-                     -5132    Xi_bbar+        "               >> usu.dat
echo "#          5232    Xi_b0                     -5232    Xi_bbar0        "               >> usu.dat
echo "#          5142    Xi_bc0                    -5142    Xi_bcbar0       "               >> usu.dat
echo "#          5242    Xi_bc+                    -5242    Xi_bcbar-       "               >> usu.dat
echo "#          5342    Omega_bc0                 -5342    Omega_bcbar0    "               >> usu.dat
echo "#          2112    n0                        -2112    nbar0           "               >> usu.dat
echo "#          2212    p+                        -2212    pbar-           "               >> usu.dat
echo "#          3112    Sigma-                    -3112    Sigmabar+       "               >> usu.dat
echo "#          3212    Sigma0                    -3212    Sigmabar0       "               >> usu.dat
echo "#          3222    Sigma+                    -3222    Sigmabar-       "               >> usu.dat
echo "#          3312    Xi-                       -3312    Xibar+          "               >> usu.dat
echo "#          3322    Xi0                       -3322    Xibar0          "               >> usu.dat
echo "#          4112    Sigma_c0                  -4112    Sigma_cbar0     "               >> usu.dat
echo "#          4212    Sigma_c+                  -4212    Sigma_cbar-     "               >> usu.dat
echo "#          4222    Sigma_c++                 -4222    Sigma_cbar--    "               >> usu.dat
echo "#          4312    Xi'_c0                    -4312    Xi'_cbar0       "               >> usu.dat
echo "#          4322    Xi'_c+                    -4322    Xi'_cbar-       "               >> usu.dat
echo "#          4332    Omega_c0                  -4332    Omega_cbar0     "               >> usu.dat
echo "#          4412    Xi_cc+                    -4412    Xi_ccbar-       "               >> usu.dat
echo "#          4422    Xi_cc++                   -4422    Xi_ccbar--      "               >> usu.dat
echo "#          4432    Omega_cc+                 -4432    Omega_ccbar-    "               >> usu.dat
echo "#          5112    Sigma_b-                  -5112    Sigma_bbar+     "               >> usu.dat
echo "#          5212    Sigma_b0                  -5212    Sigma_bbar0     "               >> usu.dat
echo "#          5222    Sigma_b+                  -5222    Sigma_bbar-     "               >> usu.dat
echo "#          5312    Xi'_b-                    -5312    Xi'_bbar+       "               >> usu.dat
echo "#          5322    Xi'_b0                    -5322    Xi'_bbar0       "               >> usu.dat
echo "#          5332    Omega_b-                  -5332    Omega_bbar+     "               >> usu.dat
echo "#          5412    Xi'_bc0                   -5412    Xi'_bcbar0      "               >> usu.dat
echo "#          5422    Xi'_bc+                   -5422    Xi'_bcbar-      "               >> usu.dat
echo "#          5432    Omega'_bc0                -5432    Omega'_bcba     "               >> usu.dat
echo "#          5442    Omega_bcc+                -5442    Omega_bccbar-   "               >> usu.dat
echo "#          5512    Xi_bb-                    -5512    Xi_bbbar+       "               >> usu.dat
echo "#          5522    Xi_bb0                    -5522    Xi_bbbar0       "               >> usu.dat
echo "#          5532    Omega_bb-                 -5532    Omega_bbbar+    "               >> usu.dat
echo "#          5542    Omega_bbc0                -5542    Omega_bbcbar0   "               >> usu.dat
echo "#          1114    Delta-                    -1114    Deltabar+       "               >> usu.dat
echo "#          2114    Delta0                    -2114    Deltabar0       "               >> usu.dat
echo "#          2214    Delta+                    -2214    Deltabar-       "               >> usu.dat
echo "#          2224    Delta++                   -2224    Deltabar--      "               >> usu.dat
echo "#          3114    Sigma*-                   -3114    Sigma*bar+      "               >> usu.dat
echo "#          3214    Sigma*0                   -3214    Sigma*bar0      "               >> usu.dat
echo "#          3224    Sigma*+                   -3224    Sigma*bar-      "               >> usu.dat
echo "#          3314    Xi*-                      -3314    Xi*bar+         "               >> usu.dat
echo "#          3324    Xi*0                      -3324    Xi*bar0         "               >> usu.dat
echo "#          3334    Omega-                    -3334    Omegabar+       "               >> usu.dat
echo "#          4114    Sigma*_c0                 -4114    Sigma*_cbar0    "               >> usu.dat
echo "#          4214    Sigma*_c+                 -4214    Sigma*_cbar-    "               >> usu.dat
echo "#          4224    Sigma*_c++                -4224    Sigma*_cbar--   "               >> usu.dat
echo "#          4314    Xi*_c0                    -4314    Xi*_cbar0       "               >> usu.dat
echo "#          4324    Xi*_c+                    -4324    Xi*_cbar-       "               >> usu.dat
echo "#          4334    Omega*_c0                 -4334    Omega*_cbar0    "               >> usu.dat
echo "#          4414    Xi*_cc+                   -4414    Xi*_ccbar-      "               >> usu.dat
echo "#          4424    Xi*_cc++                  -4424    Xi*_ccbar--     "               >> usu.dat
echo "#          4434    Omega*_cc+                -4434    Omega*_ccbar-   "               >> usu.dat
echo "#          4444    Omega*_ccc++              -4444    Omega*_cccbar-  "               >> usu.dat
echo "#          5114    Sigma*_b-                 -5114    Sigma*_bbar+    "               >> usu.dat
echo "#          5214    Sigma*_b0                 -5214    Sigma*_bbar0    "               >> usu.dat
echo "#          5224    Sigma*_b+                 -5224    Sigma*_bbar-    "               >> usu.dat
echo "#          5314    Xi*_b-                    -5314    Xi*_bbar+       "               >> usu.dat
echo "#          5324    Xi*_b0                    -5324    Xi*_bbar0       "               >> usu.dat
echo "#          5334    Omega*_b-                 -5334    Omega*_bbar+    "               >> usu.dat
echo "#          5414    Xi*_bc0                   -5414    Xi*_bcbar0      "               >> usu.dat
echo "#          5424    Xi*_bc+                   -5424    Xi*_bcbar-      "               >> usu.dat
echo "#          5434    Omega*_bc0                -5434    Omega*_bcbar0   "               >> usu.dat
echo "#          5444    Omega*_bcc+               -5444    Omega*_bccbar-  "               >> usu.dat
echo "#          5514    Xi*_bb-                   -5514    Xi*_bbbar+      "               >> usu.dat
echo "#          5524    Xi*_bb0                   -5524    Xi*_bbbar0      "               >> usu.dat
echo "#          5534    Omega*_bb-                -5534    Omega*_bbbar+   "               >> usu.dat
echo "#          5544    Omega*_bbc0               -5544    Omega*_bbcbar0  "               >> usu.dat
echo "#          5554    Omega*_bbb-               -5554    Omega*_bbbbar+  "               >> usu.dat
echo "#       1000001    ~d_L                   -1000001    ~d_Lbar         "               >> usu.dat
echo "#       1000002    ~u_L                   -1000002    ~u_Lbar         "               >> usu.dat
echo "#       1000003    ~s_L                   -1000003    ~s_Lbar         "               >> usu.dat
echo "#       1000004    ~c_L                   -1000004    ~c_Lbar         "               >> usu.dat
echo "#       1000005    ~b_1                   -1000005    ~b_1bar         "               >> usu.dat
echo "#       1000006    ~t_1                   -1000006    ~t_1bar         "               >> usu.dat
echo "#       1000011    ~e_L-                  -1000011    ~e_L+           "               >> usu.dat
echo "#       1000012    ~nu_eL                 -1000012    ~nu_eLbar       "               >> usu.dat
echo "#       1000013    ~mu_L-                 -1000013    ~mu_L+          "               >> usu.dat
echo "#       1000014    ~nu_muL                -1000014    ~nu_muLbar      "               >> usu.dat
echo "#       1000015    ~tau_1-                -1000015    ~tau_1+         "               >> usu.dat
echo "#       1000016    ~nu_tauL               -1000016    ~nu_tauLbar     "               >> usu.dat
echo "#       1000021    ~g              "                                                  >> usu.dat
echo "#       1000022    ~chi_10         "                                                  >> usu.dat
echo "#       1000023    ~chi_20         "                                                  >> usu.dat
echo "#       1000024    ~chi_1+                -1000024    ~chi_1-         "               >> usu.dat
echo "#       1000025    ~chi_30         "                                                  >> usu.dat
echo "#       1000035    ~chi_40         "                                                  >> usu.dat
echo "#       1000037    ~chi_2+                -1000037    ~chi_2-         "               >> usu.dat
echo "#       1000039    ~Gravitino      "                                                  >> usu.dat
echo "#       2000001    ~d_R                   -2000001    ~d_Rbar         "               >> usu.dat
echo "#       2000002    ~u_R                   -2000002    ~u_Rbar         "               >> usu.dat
echo "#       2000003    ~s_R                   -2000003    ~s_Rbar         "               >> usu.dat
echo "#       2000004    ~c_R                   -2000004    ~c_Rbar         "               >> usu.dat
echo "#       2000005    ~b_2                   -2000005    ~b_2bar         "               >> usu.dat
echo "#       2000006    ~t_2                   -2000006    ~t_2bar         "               >> usu.dat
echo "#       2000011    ~e_R-                  -2000011    ~e_R+           "               >> usu.dat
echo "#       2000012    ~nu_eR                 -2000012    ~nu_eRbar       "               >> usu.dat
echo "#       2000013    ~mu_R-                 -2000013    ~mu_R+          "               >> usu.dat
echo "#       2000014    ~nu_muR                -2000014    ~nu_muRbar      "               >> usu.dat
echo "#       2000015    ~tau_2-                -2000015    ~tau_2+         "               >> usu.dat
echo "#       2000016    ~nu_tauR               -2000016    ~nu_tauRbar     "               >> usu.dat
echo "#       3000111    pi_tc0          "                                                  >> usu.dat
echo "#       3000211    pi_tc+                 -3000211    pi_tc-          "               >> usu.dat
echo "#       3000221    pi'_tc0         "                                                  >> usu.dat
echo "#       3000331    eta_tc0         "                                                  >> usu.dat
echo "#       3000113    rho_tc0         "                                                  >> usu.dat
echo "#       3000213    rho_tc+                -3000213    rho_tc-         "               >> usu.dat
echo "#       3000223    omega_tc        "                                                  >> usu.dat
echo "#       3100021    V8_tc           "                                                  >> usu.dat
echo "#       3100111    pi_22_1_tc      "                                                  >> usu.dat
echo "#       3200111    pi_22_8_tc      "                                                  >> usu.dat
echo "#       3100113    rho_11_tc       "                                                  >> usu.dat
echo "#       3200113    rho_12_tc       "                                                  >> usu.dat
echo "#       3300113    rho_21_tc       "                                                  >> usu.dat
echo "#       3400113    rho_22_tc       "                                                  >> usu.dat
echo "#       4000001    d*                     -4000001    d*bar           "               >> usu.dat
echo "#       4000002    u*                     -4000002    u*bar           "               >> usu.dat
echo "#       4000011    e*-                    -4000011    e*bar+          "               >> usu.dat
echo "#       4000012    nu*_e0                 -4000012    nu*_ebar0       "               >> usu.dat
echo "#       5000039    Graviton*       "                                                  >> usu.dat
echo "#       9900012    nu_Re           "                                                  >> usu.dat
echo "#       9900014    nu_Rmu          "                                                  >> usu.dat
echo "#       9900016    nu_Rtau         "                                                  >> usu.dat
echo "#       9900023    Z_R0            "                                                  >> usu.dat
echo "#       9900024    W_R+                   -9900024    W_R-            "               >> usu.dat
echo "#       9900041    H_L++                  -9900041    H_L--           "               >> usu.dat
echo "#       9900042    H_R++                  -9900042    H_R--           "               >> usu.dat
echo "#       9900110    rho_diff0       "                                                  >> usu.dat
echo "#       9900210    pi_diffr+              -9900210    pi_diffr-       "               >> usu.dat
echo "#       9900220    omega_di        "                                                  >> usu.dat
echo "#       9900330    phi_diff        "                                                  >> usu.dat
echo "#       9900440    J/psi_di        "                                                  >> usu.dat
echo "#       9902110    n_diffr0               -9902110    n_diffrbar0     "               >> usu.dat
echo "#       9902210    p_diffr+               -9902210    p_diffrbar-     "               >> usu.dat
echo "#       9900443    cc~[3S18]       "                                                  >> usu.dat
echo "#       9900441    cc~[1S08]       "                                                  >> usu.dat
echo "#       9910441    cc~[3P08]       "                                                  >> usu.dat
echo "#       9900553    bb~[3S18]       "                                                  >> usu.dat
echo "#       9900551    bb~[1S08]       "                                                  >> usu.dat
echo "#       9910551    bb~[3P08]       "                                                  >> usu.dat
####################               usu.dat              ########################
################################################################################



################################################################################
#   The following statements are not required to be modified.                  #
################################################################################



################################################################################
################################################################################
####################              Makefile              ########################
echo "# This is a toy Makefile for PACIAE."     > Makefile
echo "# By An-Ke Lei at CCNU on 17/10/2022"     >> Makefile
echo "# Last updated"                           >> Makefile
echo "#         by An-Ke Lei on 17/10/2022"     >> Makefile
echo "# How to use:"                            >> Makefile
echo "#   1. Type \"make\" command to compile and build PACIAE running file (${name_x})."   >> Makefile
echo "#   2. Type \"make clean\" command to clean the *.o , *.mod and *.x files. "          >> Makefile
echo "#   Usually, second command is not required to be used."      >> Makefile
echo                                                                >> Makefile
echo                                                                >> Makefile
echo "# The name of the executeble file."                           >> Makefile
# echo "target := xPaciae.x"                                        >> Makefile
echo "target := ${name_x}"                                          >> Makefile
echo                                                                >> Makefile
echo "# The names of source files."                                 >> Makefile
echo "#   Modify the following name to the needed one."             >> Makefile
echo "src_1 := ${name_f1}"                                          >> Makefile
echo "src_2 := ${name_f2}"                                          >> Makefile
echo "src_3 := ${name_f3}"                                          >> Makefile
echo "src_4 := ${name_f4}"                                          >> Makefile
echo "src_5 := ${name_f5}"                                          >> Makefile
echo "src_6 := ${name_f6}"                                          >> Makefile
echo "src_7 := ${name_f7}"                                          >> Makefile
echo "src_8 := ${name_f8}"                                          >> Makefile
echo "# src_9 := ${name_f9}"                                        >> Makefile
echo "# scr_10 := ${name_f10}"                                      >> Makefile
echo                                                                >> Makefile
echo "# The intermediate files."                                    >> Makefile
echo "objects := \$(src_1).o \$(src_2).o \$(src_3).o \\"            >> Makefile
echo "           \$(src_4).o \$(src_5).o \$(src_6).o \\"            >> Makefile
echo "           \$(src_7).o \$(src_8).o # \$(src_9).o"             >> Makefile
echo "#          \$(src_10).o"                                      >> Makefile
echo                                                                >> Makefile
#Lei20230608B-------------------------------------------------------------------
echo "targ_rms := ${name_xRms}"                                     >> Makefile
echo "src_rms := ${name_rms}"                                       >> Makefile
echo "# The rms intermediate files."                                >> Makefile
echo "obj_rms := \$(src_rms).o"                                     >> Makefile
#Lei20230608E-------------------------------------------------------------------
echo "# Compiler"                                                   >> Makefile
echo "compiler := gfortran"                                         >> Makefile
echo "# Compiling flags"                                            >> Makefile
echo "comp_flags := -g -fbounds-check"                              >> Makefile
echo "# comp_flags := -g -Wall -fbounds-check"                      >> Makefile
echo "# comp_flags = -g -Wall -fbounds-Check -mcmodel=large"        >> Makefile
echo                                                                >> Makefile
echo "# Generating"                                                 >> Makefile
#Lei20230608B-------------------------------------------------------------------
echo "all : \$(target) \$(targ_rms)"                                >> Makefile
echo "# Generating xRms"                                            >> Makefile
echo "\$(targ_rms) : \$(obj_rms)"                                   >> Makefile
echo "	\$(compiler) \$(comp_flags) -o \$(targ_rms) -O \$(obj_rms)" >> Makefile
#Lei20230608E-------------------------------------------------------------------
echo "\$(target) : \$(objects)"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -o \$(target) -O \$(objects)"   >> Makefile
echo                                                                >> Makefile
echo "\$(src_1).o: \$(src_1).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_1).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_2).o: \$(src_2).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_2).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_3).o: \$(src_3).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_3).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_4).o: \$(src_4).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_4).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_5).o: \$(src_5).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_5).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_6).o: \$(src_6).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_6).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_7).o: \$(src_7).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_7).f"              >> Makefile
echo                                                                >> Makefile
echo "\$(src_8).o: \$(src_8).f"                                     >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_8).f"              >> Makefile
echo                                                                >> Makefile
echo "# \$(src_9).o: \$(src_9).f"                                   >> Makefile
echo "#	\$(compiler) \$(comp_flags) -c -O \$(src_9).f"              >> Makefile
echo                                                                >> Makefile
echo "# \$(src_10).o: \$(src_10).f"                                 >> Makefile
echo "# 	\$(compiler) \$(comp_flags) -c -O \$(src_10).f"         >> Makefile
echo                                                                >> Makefile
#Lei20230608B-------------------------------------------------------------------
echo "\$(src_rms).o: \$(src_rms).f90"                               >> Makefile
echo "	\$(compiler) \$(comp_flags) -c -O \$(src_rms).f90"          >> Makefile
#Lei20230608E-------------------------------------------------------------------
echo                                                                >> Makefile
echo "# Cleans .o, .mod and .x files."                              >> Makefile
echo ".PHONY : clean"                                               >> Makefile
echo "clean : "                                                     >> Makefile
echo "	rm -rf *.o *.mod *.x"                                       >> Makefile
####################              Makefile              ########################
################################################################################
################################################################################



################################################################################
################################################################################
####################           Command lines            ########################
#*******************************************************************************
# Creates src folder.
if [ -d "./src/" ]; then
    echo "Folder src already exists."
else
    mkdir ./src
fi
f_file=$(ls *.f 2> /dev/null | wc -l)
if [ "${f_file}" != "0" ] ; then
    mv -f *.f ./src
fi
mv -f Makefile ./src
#Lei20230608B-------------------------------------------------------------------
if [ -f "./${name_rms}.f90" ]; then
    mv -f ./${name_rms}.f90 ./src
fi
#Lei20230608E-------------------------------------------------------------------
wait
# Enters src folder then compiles and builds executable file.
cd ./src
echo
make
echo
cd ../
##------------------------------------------------------------------------------
# Creates bin folder.
if [ -d "./bin/" ]; then
    echo "Folder bin already exists."
else
    mkdir ./bin
fi
cp -f ./src/${name_x} ./bin
cp -f ./src/${name_xRms} ./bin   #Lei20230608
##------------------------------------------------------------------------------
# Creates doc folder.
if [ -d "./doc/" ]; then
    echo "Folder doc already exists."
else
    mkdir ./doc
fi
if [ -f "./README.md" ]; then
   mv -f ./README.md ./doc
fi
if [ -f "./LICENSE" ]; then
   mv -f ./LICENSE ./doc
fi
##------------------------------------------------------------------------------
# Creates etc folder.
if [ -d "./etc/" ]; then
    echo "Folder etc already exists."
else
    mkdir ./etc
fi
# if [ -f "./EPS09*" ]; then
#    mv -f EPS09* ./etc
# fi
if [ -f "./default_usu.dat*" ]; then
    mv -f ./default_usu.dat ./etc
fi
mv -f usu.dat ./etc
##------------------------------------------------------------------------------
# Creates log folder.
if [ -d "./log/" ]; then
    echo "Folder log already exists."
else
    mkdir ./log
fi
# Good for normal LINUX and SLURM, but not for LSF.
# TODO(Lei20221016): need to improve for LSF.
out_file=$(ls *.out 2> /dev/null | wc -l)
if [ "${out_file}" != "0" ] ; then
    mv -f *.out ./log
fi
log_file=$(ls *.log 2> /dev/null | wc -l)
if [ "${log_file}" != "0" ] ; then
    mv -f *.log ./log
fi
# mv -f *.out *.log ./log
##------------------------------------------------------------------------------
# Creates sim folder.
if [ -d "./sim/" ]; then
    echo "Folder sim already exists."
else
    mkdir ./sim
fi
cd ./sim   # The simulation will be performed in sim folder.

#*******************************************************************************
# File name
if [[ "${na_proj}" = "1" && "${nz_proj}" = "1" ]]; then
    proj_name="p"
elif [[ "${na_proj}" = "1" && "${nz_proj}" = "-1" ]]; then
    proj_name="pbar"
elif [[ "${na_proj}" = "1" && "${nz_proj}" = "0" ]]; then
    proj_name="n"
elif [ "${na_proj}" = "63" ]; then
    proj_name="Cu${na_proj}"
elif [[ "${na_proj}" = "96" && "${nz_proj}" = "44" ]]; then
    proj_name="Ru${na_proj}"
elif [[ "${na_proj}" = "96" && "${nz_proj}" = "40" ]]; then
    proj_name="Zr${na_proj}"
elif [ "${na_proj}" = "108" ]; then
    proj_name="Ag${na_proj}"
elif [ "${na_proj}" = "129" ]; then
    proj_name="Xe${na_proj}"
elif [ "${na_proj}" = "197" ]; then
    proj_name="Au${na_proj}"
elif [ "${na_proj}" = "208" ]; then
    proj_name="Pb${na_proj}"
elif [ "${na_proj}" = "238" ]; then
    proj_name="U${na_proj}"
else
    proj_name="${na_proj}"
fi
##------------------------------------------------------------------------------
# File name
if [[ "${na_targ}" = "1" && "${nz_targ}" = "1" ]]; then
    targ_name="p"
elif [[ "${na_targ}" = "1" && "${nz_targ}" = "-1" ]]; then
    targ_name="pbar"
elif [[ "${na_targ}" = "1" && "${nz_targ}" = "0" ]]; then
    targ_name="n"
elif [ "${na_targ}" = "63" ]; then
    targ_name="Cu${na_targ}"
elif [[ "${na_targ}" = "96" && "${nz_targ}" = "44" ]]; then
    targ_name="Ru${na_targ}"
elif [[ "${na_targ}" = "96" && "${nz_targ}" = "40" ]]; then
    targ_name="Zr${na_targ}"
elif [ "${na_targ}" = "108" ]; then
    targ_name="Ag${na_targ}"
elif [ "${na_targ}" = "129" ]; then
    targ_name="Xe${na_targ}"
elif [ "${na_targ}" = "197" ]; then
    targ_name="Au${na_targ}"
elif [ "${na_targ}" = "208" ]; then
    targ_name="Pb${na_targ}"
elif [ "${na_targ}" = "238" ]; then
    targ_name="U${na_targ}"
else
    targ_name="${na_targ}"
fi
##------------------------------------------------------------------------------
# File name
if [ "${i_frame}" = "1" ]; then
    coll_name="COLL"
elif [ "${i_frame}" = "0" ]; then
    coll_name="FXT"
else
    coll_name="${i_frame}"
fi
##------------------------------------------------------------------------------
# Creates files.
if [ -d "./${proj_name}_${targ_name}_${coll_name}/" ]; then
    echo "Folder ${proj_name}_${targ_name}_${coll_name} already exists."
else
    mkdir ./${proj_name}_${targ_name}_${coll_name}
fi
cd ./${proj_name}_${targ_name}_${coll_name}

#*******************************************************************************
# iStg: i_stage
# iChl: i_channel
# Loop name
if [ "${i_sim_mode}" = "1" ]; then
    loop_name="A"
elif [ "${i_sim_mode}" = "2" ]; then
    loop_name="B"
elif [ "${i_sim_mode}" = "3" ]; then
    loop_name="C"
else
    loop_name="${i_sim_mode}"
fi
dir_energy="${energy}GeV_iChl${i_channel}_iStg${i_stage}_loop${loop_name}"
if [ "${i_sim_mode}" = "1" ]; then
    dir_energy="${energy}GeV_loop${loop_name}"
fi
if [ -d "./${dir_energy}/" ]; then
    echo "Folder ${dir_energy} already exists."
else
    mkdir ./${dir_energy}
fi
cd ./${dir_energy}

#*******************************************************************************
# b: impact parameter
# iO: i_overlap
# aL: a_lund
# bL: b_lund
# aFF: a_FF
# K : k_pythia (in parini)
#     kI: k_pythia (in parini)
#     kC: k_parcas
# iH: i_hadronization
# iHc: i_hadcas
# iC : i_channel
# iStr: i_string_tension
# iPc : i_phase_constraint
# iFF : i_fragment_function
# iM : i_cMe
# pTp : pT_perturb
# vF : volum_factor
# xVf : ex_volum_factor
# iCr : i_coord_recover
# Default
dir="./b${b_min}_${b_max}_aL${a_lund}_bL${b_lund}_eDeex${eDeex}_K${k_pythia}"
# dir="./b${b_min}_${b_max}_aL${a_lund}_bL${b_lund}_kI${k_pythia}_kC${k_parcas}_pTp${pT_perturb}_iO${i_overlap}_iH${i_hadronization}_iHc${i_hadcas}_iStr${i_string_tension}_iP${i_phase_constraint}_iFF${i_fragment_function}_iM${i_cme}"

if [[ "${na_proj}" = "1" && "${na_targ}" = "1" ]]; then
# Elementary NN collision
    if [[ "${i_hadronization}" = "0" ]]; then
        dir="./Sfm_aL${a_lund}_bL${b_lund}_pP0${p_perp0}_pPm${p_perp_min}_sig${pT_width}_K${k_pythia}_Kc${k_parcas}_iHc${i_hadcas}_iTune${i_tune}"
    elif [[ "${i_hadronization}" != "0" ]]; then
        dir="./Coal${i_hadronization}_pP0${p_perp0}_pPm${p_perp_min}_iFF${i_fragment_function}_ipT${i_pT_samp}_iPm${i_pT_max}_iPc${i_phase_constraint}_eDeex${eDeex}_nDq${n_deexc_time_per_q}_iDg${i_deex_gen}_iDx${i_deex}_sig${pT_width}_sigP${sig_kPerp}_K${k_pythia}_Kc${k_parcas}_iHc${i_hadcas}"
    # else
    fi
    if [ "${i_sim_mode}" = "1" ]; then
        dir="./xRatio${ratio_xSec_inel_over_tot}_pDelta${prob_Delta_decay}_iHc${i_hadcas}"
    fi
# elif [[  ]]; then
else
# pA, AA collisions
    if [[ "${i_hadronization}" = "0" ]]; then
        dir="./b${b_min}_${b_max}_sfm_aL${a_lund}_bL${b_lund}_pP0${p_perp0}_pPm${p_perp_min}_sig${pT_width}_K${k_pythia}_Kc${k_parcas}_iO${i_overlap}_iHc${i_hadcas}_iTune${i_tune}"
    elif [[ "${i_hadronization}" != "0" ]]; then
        dir="./b${b_min}_${b_max}_coal${i_hadronization}_pP0${p_perp0}_pPm${p_perp_min}_iFF${i_fragment_function}_ipT${i_pT_samp}_iPm${i_pT_max}_iPc${i_phase_constraint}_eDeex${eDeex}_nDq${n_deexc_time_per_q}_iDg${i_deex_gen}_iDx${i_deex}_sig${pT_width}_sigP${sig_kPerp}_K${k_pythia}_Kc${k_parcas}_iO${i_overlap}_iHc${i_hadcas}"
    # else
    fi
    if [ "${i_sim_mode}" = "1" ]; then
        dir="./b${b_min}_${b_max}_xRatio${ratio_xSec_inel_over_tot}_pDelta${prob_Delta_decay}_iO${i_overlap}_iHc${i_hadcas}"
    fi
fi


if [ -d "${dir}/" ]; then
    echo "Folder ${dir} already exists."
else
    mkdir ${dir}
fi
#ls -a

if [ -d "${dir}/${tot_eve}_events/" ]; then
    echo "Folder ${tot_eve}_events already exists."
else
    mkdir ${dir}/${tot_eve}_events
fi
cd ${dir}/${tot_eve}_events
echo
pwd
echo

# Not required to be modified. Pre-setups for CPU/program-running.
i_cpu=0      # Do not change i_cpu=0. ID od CPU (from 0)
i_run=1      # DO not change i_run=1. Counter for program-running.

while [ ${i_run} -le ${n_run} ]
do
    if [ -d "./PACIAE_${i_run}/" ]; then
        echo
    else
        mkdir ./PACIAE_${i_run}
    fi
    # ls -a
    cp -f ../../../../../bin/${name_x} ./PACIAE_${i_run}/${i_run}_${name_x}
    cp -f ../../../../../etc/usu.dat ./PACIAE_${i_run}/
    # cp -f ../../../../etc/EPS09* ./PACIAE_${i_run}/
    cd ./PACIAE_${i_run}/
    # ls -a
    nohup time ./${i_run}_${name_x} > screen.log &
# Binds one running task to one CPU.
    # nohup time taskset -c ${i_cpu} ./${i_run}_${name_x} > screen.log &
    cd ..
    # ls -a

    echo "PACIAE-${i_run}"

    ((i_cpu=i_cpu+1))
    ((i_run=i_run+1))
done
((i_run=i_run-1))

wait
# sleep 30s   #Lei20230608

#Lei20230608B-----------------------------------------------------------
if [ -d "./rms_file/" ]; then
    echo
else
    mkdir ./rms_file
fi
cp -f ../../../../../bin/${name_xRms} ./rms_file
cd ./rms_file
echo "${n_run},${n_eve}" > input_rms_analysis.dat
echo "${i_sim_mode},${b_samp}" >> input_rms_analysis.dat   #Lei20230819 i_sim_mode, b_samp

# Not required to be modified. Pre-setups for CPU/program-running.
i_run=1      # DO not change i_run=1. Counter for program-running.

while [ ${i_run} -le ${n_run} ]
do
    # ls -a
    cp -f ../PACIAE_${i_run}/rms.out ./rms_${i_run}.out
    echo "rms-${i_run}"
    ((i_run=i_run+1))
done
((i_run=i_run-1))
wait
nohup time ./${name_xRms} > screen.log &
wait
cd ..
#Lei20230608E-----------------------------------------------------------

cd ../../../../..
echo
ls -a
echo
echo "PACIAE-script running ended!"
echo
####################           Command lines            ########################
################################################################################
################################################################################



################################################################################
################################################################################
################################################################################
#                                                                              #
#                                                                              #
################################################################################
################################################################################
#                                                                              #
#         PPPPPPP       AAAAA      CCCCCCC     IIIIIII     AAAAA     EEEEEEEEE #
#        P      p     A     A    C       C       I       A     A    E          #
#       P       p   A       A   C               I      A       A   E           #
#      P      p    A       A   C               I      A       A   E            #
#     PPPPPPP     AAAAAAAAA   C               I      AAAAAAAAA   EEEEEEEEE     #
#    P           A       A   C               I      A       A   E              #
#   P           A       A   C               I      A       A   E               #
#  P           A       A    C       C      I      A       A   E                #
# P           A       A     CCCCCCC    IIIIIII   A       A   EEEEEEEEE         #
#                                                                              #
################################################################################
################################################################################
################################################################################
