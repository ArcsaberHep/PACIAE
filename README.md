<!-- This is a README file for usage of PACIAE.
     Written by Markdown language.
                 By Anke at CCNU on 10/16/2022 
                                               
                    Last updated on 08/13/2023 
 -->

# The parton and hadron cascade model PACIAE

 ***PACIAE*** model (***Parton And-hadron China Institute of Atomic Energy***) is a multipurpose Monte Carlo event generator developed to describe a wide range of collisions, including lepton-lepton, lepton-hadron, lepton-nucleus, hadron-hadron, hadron-nucleus, and nucleus-nucleus collisions. It is built based on PYTHIA-6.428 and incorporates parton and hadron rescattering stages to take care of the nuclear medium effects.

## Installation

Just install the file and decompress it. Then you can get the PACIAE source code directly.

## Usage

We encourage users to run the program on LINUX.
There are two ways to run the program.
 1. Use the normal compiler to complie the soruce code and run the PACIAE program. Take GFortran on LINUX as an example:
    - Compile and link the programs by the command:
        ```
            gfortran -O -C *.f -o paciae.x
        ```
    - Modify the input file of usu.dat according to your wish.
    - Run the program by the command:
        ```
            time ./paciae.x
        ```
      One could use the following command to run the program in the background and record the log information.
        ```
            nohup time ./paciae.x > paciae.log &
        ```
 2. Use the PACIAE.sh shell-script to compile and link the code, generate usu.dat file and run program automatically.
    - Modify the PACIAE.sh file as needed.
    - Give executable permission to the PACIAE.sh file by command:
        ```
            chmod +x PACIAE.sh
        ```
    - Run the PACIAE.sh script by the command:
        ```
            time ./PACIAE.sh
        ```
         The more recommended command is:
         ```
            time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log
         ``` 
         which records the screen information to a log file with date and time.
         Of course, run them in the background:
         ```
            nohup time ./PACIAE.sh &
         ```
         ```
            nohup time ./PACIAE.sh | tee $(date "+%Y%m%d%H%M%S").log &
        ```
    It is also worth mentioning that tasks can be submitted to computer clusters and supercomputers using the PACIAE.sh script. (It is currently only available in SLRUM and LSF scheduling systems) More detailed information and usage please read the PACIAE.sh file.

## Maintainers

[@ArcsaberHep](https://github.com/ArcsaberHep)
[@ArcsaberkxL](https://github.com/ArcsaberkxL)

## Contributing

Feel free to dive in! Any bug reports, comments and suggestions are welcome. Please do not hesitate to contact us.

## Contributors

 - [An-Ke Lei](https://inspirehep.net/authors/1965068), ankeleihep@gmail.com or ankelei@mails.ccnu.edu.cn <!-- Key Laboratory of Quark and Lepton Physics (MOE) and Institute of Particle Physics, Central China Normal University, Wuhan 430079, China. -->
 - [Zhi-Lei She](https://inspirehep.net/authors/1903611), shezhilei@cug.edu.cn <!-- School of Mathematical and Physical Sciences, Wuhan Textile University, Wuhan 430200, China --> 
 - [Dai-Mei Zhou](https://inspirehep.net/authors/1030208), zhoudm@mail.ccnu.edu.cn <!-- Key Laboratory of Quark and Lepton Physics (MOE) and Institute of Particle Physics, Central China Normal University, Wuhan 430079, China. -->
 - [Yu-Liang Yan](https://inspirehep.net/authors/1051028), yanyl@ciae.ac.cn <!-- China Institute of Atomic Energy, P.O. Box 275 (10), Beijing, 102413,China. -->
 - [Ben-Hao Sa](https://inspirehep.net/authors/990834), sabh@ciae.ac.cn (model founder) <!-- China Institute of Atomic Energy, P.O. Box 275 (10), Beijing, 102413,China. -->

## License

[GPL v2.0](LICENSE)

## Released papers

 - **Recent version:** Revisiting the centrality definition and observable centrality dependence of relativistic heavy-ion collisions in PACIAE model, [Comput. Phys. Commun. 284 (2023) 108615](https://doi.org/10.1016/j.cpc.2022.108615) or [arXiv:2212.04087 [nucl-th]](https://doi.org/10.48550/arXiv.2212.04087)
<br/>

 - **Past version:** PACIAE 2.2.1: An updated issue of the parton and hadron cascade model PACIAE 2.2, [Comput. Phys. Commun. 274 (2022) 108289](https://doi.org/10.1016/j.cpc.2022.108289) <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

 - Announcement for the replacement of the PACIAE 2.1 and PACIAE 2.2 series, [Comput. Phys. Commun. 224 (2018) 417-418](https://doi.org/10.1016/j.cpc.2017.10.006) <!-- or [arXiv: [nucl-th]](https://doi.org/) -->

 - An upgraded issue of the parton and hadron cascade model, PACIAE 2.2, [Comput. Phys. Commun. 193 (2015) 89-94](https://doi.org/10.1016/j.cpc.2015.01.022) or [arXiv:1412.7579 [nucl-th]](https://doi.org/10.48550/arXiv.1412.7579)

 - PACIAE 2.1: An updated issue of the parton and hadron cascade model PACIAE 2.0, [Comput. Phys. Commun. 184 (2013) 1476-1479](https://doi.org/10.1016/j.cpc.2012.12.026) or [arXiv:1206.4795 [nucl-th]](https://doi.org/10.48550/arXiv.1206.4795)

 - PACIAE 2.0: An updated parton and hadron cascade model (program) for the relativistic nuclear collisions, [Comput. Phys. Commun. 183 (2012) 333-346](https://doi.org/10.1016/j.cpc.2011.08.021) or [arXiv:1104.1238 [nucl-th]](https://doi.org/10.48550/arXiv.1104.1238)

 - (PACIAE 1.0) Charge particle universal rapidity scaling in e+e-, pbar + p and Au + Au collisions at relativistic energies and its partonic origin, [J. Phys. G 32 (2006) 243-250](https://doi.org/10.1088/0954-3899/32/3/001); Influence of the partonic Pauli blocking on the hadronic final state in relativistic nucleus-nucleus collisions, [Phys. Rev. C 70 (2004) 034904](https://doi.org/10.1103/PhysRevC.70.034904)
<br/>

 - **Ancient version:** (JPCIAE) J/psi dynamical suppression in a hadron and string cascade model (or Formation time effect on J/psi dynamical nuclear suppression), [Phys. Rev. C 59 (1999) 2728-2733](https://doi.org/10.1103/PhysRevC.59.2728) or [arXiv:nucl-th/9803033](https://arxiv.org/abs/nucl-th/9803033); J/psi normal and anomalous suppressions in a hadron and string cascade model, [J.Phys.G 25 (1999) 1123-1133](https://doi.org/10.1088/0954-3899/25/6/302) or [arXiv:nucl-th/9809020](https://arxiv.org/abs/nucl-th/9809020); Inclusive and direct photons in S + Au collisions at 200A GeV/c, [Phys. Rev. C 61 (2000) 064905](https://doi.org/10.1103/PhysRevC.61.064905) or [arXiv:nucl-th/9904035](https://arxiv.org/abs/nucl-th/9904035).
<br/>

 - **Origin version:** LUCIAE 3.0: A New version of a computer program for firecracker model and rescattering in relativistic heavy ion collisions, [Comput. Phys. Commun. 116 (1999) 353](https://doi.org/10.1016/S0010-4655(98)00138-6) or [arXiv:nucl-th/9804001](https://doi.org/10.48550/arXiv.nucl-th/9804001)

 - An Event generator for the firecracker model and the rescattering in high-energy pA and AA collisions: LUCIAE version 2.0, [Comput. Phys. Commun. 90 (1995) 121-140](https://doi.org/10.1016/0010-4655(95)00066-O); Final state interactions in the (nuclear) FRITIOF string interaction scenario, [Z. Phys. C 70 (1996) 499](https://doi.org/10.1007/s002880050127).
<br/>

 - **Special version:** HYDRO-PACIAE, a hydrodynamic and transport hybrid model for ultra-relativistic heavy ion collisions, [J.Phys.G 40 (2013) 025102](https://doi.org/10.1088/0954-3899/40/2/025102) or [arXiv:1110.6704 [nucl-th]](https://doi.org/10.48550/arXiv.1110.6704)

## Relevant papers

- PYTHIA 6.4 Physics and Manual, [JHEP 05 (2006) 026](https://doi.org/10.1088/1126-6708/2006/05/026) or [arXiv:hep-ph/0603175 [hep-ph]](https://doi.org/10.48550/arXiv.hep-ph/0603175)

## Update Notes:

<!----------------------------------------------------------------------------->
...(waiting for updating)

<!----------------------------------------------------------------------------->
#### <font color=red> 08/2023: </font> In version PACIAE 2.3 ###

- In subroutine "oscar" of "main_23.f", a bug was fixed for the case nosc=0.
- In subroutine "tabh" of "coales_23.f", the probabilities of mixing-state measons were corrected.
- PACIAE.sh updated. A "sim" folder (simulation) would be generated to storing the results of simulations.


<!----------------------------------------------------------------------------->
#### <font color=red> 07/2023: </font> In version PACIAE 2.3 ###

- In "main_23.f" and "parini_23.f", the parameters "dtt" and "smtj3" were renamed/replaced as "x_ratio" and "decpro" directly.
- In "parini_23.f", several subroutines about leading-proton reconstruction were introduced. **TODO:** It still need detailed and complete treatment.
- In "PAPTDI", "deexcitation_E", and "deexcitation_EP" of "coales_23.f", minor optimizations were made.
- In "coales_23.f", two new implementations from "deexcitation_EP" were introduced, in which the excited quark-antiquark would be local pT compensation. **TODO:** It still need detailed and complete treatment.
- "usu.dat" updated.
- "PACIAE.sh" updated.
- Redundant statements deleted / modified.
- Bug fixed.

<!----------------------------------------------------------------------------->
#### <font color=red> 06/2023: </font> In version PACIAE 2.3 ###

- In "main_23.f", three input parameter "parp82", "i_coord_recover", "i_tune" were introduced. "parp82" corresponds to PARP(82) in PYTHIA. "i_coord_recover" controls recovering the positions of partons after parton rescattering to those before or not. "i_tune" corresponds to MSTP(5) in PYTHIA which gives easy access to different tunes encoded in PYTHIA.
- In "main_23.f", the 4-position of one of the rest partons would be assigned to the first parton in the string. This treatment would give random 3-coordinates to produced hadrons that surround the first parton after PYTHIA SFM, i.e. more random position distribution for produced hadrons.
- In "main_23.f", the output format of "rms0.out" was optimized.
- The subroutine "oscar" of "main_23.f", the OSCAR output was optimized.
- In "main_23.f", "parini_23.f", and "analy.f", the subroutine "PASTAT" was introduced to output parton-parton level cross-sections to "main.out" file. However, one should note that it is just for the basic check.
- In "parini_23.f", the leading-proton reconstruction was improved. **TODO:** It still need detailed and complete treatment.
- In subroutine "xevent" of "parini_23.f", now when simulate pA and Ap collision, the junction-type NN event would be re-generated.
- In subroutine "bream" of "parini_23.f", the broken quarks would be massless if breaking diquark is massless from PYTHIA treatment.
- In subroutines "updtlp" and "updatl" of "parini_23.f", minor optimizations were made.
- In subroutines "ctlcre_par", "ctlcre_para", "his_p", and "updpli_p" of "parcas_23.f", minor optimizations were made.
- In subroutines "ctlcre_h" of "hadcas_23.f", one minor optimization was made.
- "usu.dat" updated.
- "PACIAE.sh" updated.
- Redundant statements deleted / modified.
- Bug fixed.


<!----------------------------------------------------------------------------->
#### <font color=red> 05/2023: </font> In version PACIAE 2.3 ###

- In "main_23.f" and "analy.f", the internal online analyzing module was improved and optimized. The 5 distributions was extended to 6, i.e. the transverse momentum spectra dN/dpT. The format of user output file "rms.out" was optimized, too.
- In ""main_23.f" and "parini_23.f", the diffractive NN event without parton generation after PYTHIA calling would not be thrown away now.
- In subroutine "scat" of "parini_23.f", a minor bug of "m1 -> mm1" was fixed.
- In "parcas_23.f", the max simulated volume was introduced, controlled by "PARAM(10)" (para10).
- In "parcas_23.f", the inelastic processes and time-like radiation was fixed and improved to consider special pure-gluon strings (glue) and junction-type strings properly. **TODO:** It still need detailed and complete treatment.
- In "analy.f", the usbroutine "analy_parton" and "stati_parton" were introduced to analyze partons. The subroutines "output_hadron_distribution" and "output_parton_distribution" were introduced to output analyzing information of hadrons and partons, respectively. The analyzing module was improved and optimized to 6 distributions.
- "usu.dat" updated.
- "PACIAE.sh" updated.
- Redundant statements deleted / modified.
- Bug fixed.

<!----------------------------------------------------------------------------->
#### <font color=red> 04/2023: </font> In version PACIAE 2.3 ###

- In "main_23.f", a energy-dependent "smtj3" was introduced for low-energy loop-A at $\sqrt{s_{NN}} < 3$ GeV.
- In "parini_23.f", the loop-A was improved.
- In "parini_23.f", the hadron-hadron cross-sections inspired by additive quark model (AQM) wee introduced.
- In all of ".f" code files, many "DO-ENDDO" statements and some statements of variable initialization were optimized to save running time (memory optimization).
- In "main_23.f", different code blocks were separated by dashed lines "-----" with corresponding comments to explain their functions, i.e. "what we are doing here".
- In "main_23.f", three subroutines "rest_hadronization", "rest_sfm" and "rest_coal" were introduced to hadronize the rest partons that failed during the normal hadronization process. The re-hadronization is set after normal hadronization now, instead of calling "pa2evnt" after the whole simulation as before. "pa2evnt" was also fixed and improved again.
- In "main_23.f" and "parini_23.f" (near "PYINIT" calling), the re-generation of the NN binary collision were introduced to deal with case where charge or 4-momentum was not conserved or any errors occurred after calling PYTHIA.
- In "main_23.f", the Coalescence Mode 2 was introduced. In this mode, the gluon splitting and energetic quark deexcitation would be performed before parton rescattering. It means that there will be no gluons into parton rescattering.
- In "main_23.f", a subroutine "prt_final_info" was introduced to print final 4-momentum information to "main.out".
- In subroutines "decmom_sbe" of "main_23.f" and "decmom" of "parini_23.f", minor bugs were fixed.
- The subroutine "oscar" of "main_23.f" was improved to print event particle history correctly.
- The spectator nucleons were moved after parton rescattering and hadron rescattering correctly avoiding the interactions between them and other hadrons, i.e. "sbh moving" blocs in "main_23.f".
- The input parameter adj1(29) could select fragmentation function for both string fragmentation model (SFM) and deexcitation of coalescence model (Coal) now. A series of parameters of fragmentation functions were introduced.
- In subroutine "scat" of "parini_23.f", gamma66 removing block was moved to do the job correctly.
- In subroutine "ctlcre_par" of "parcas_23.f", minor bug were fixed.
- In "coales_23.f", the coalescence model was improved.
- In subroutine "tabh" of "coales_23.f", the specie of hadron is extended to 200 now.
- In subroutine "coal" of "coales_23.f", the reconstruction of parton list after normal "coal" calling was fixed and improved. A process of "final coalescence try" for the quarks that failed in normal coalescence was introduced.
- "usu.dat" updated.
- "PACIAE.sh" updated.
- Redundant statements deleted / modified.
- Bug fixed.

<!----------------------------------------------------------------------------->
#### <font color=red> 03/2023: </font> In version PACIAE 2.3 ###

- The "stahad_23.f" file has beed removed.
- The simulation mode in PACIAE has been classified into three modes, controlled by switch "mstj2":
    - =1, (low-energy) pure hadron simulation ***loop-A***;
    - =2, PYTHIA-like pure hadron simulation ***loop-B***;
    - =3, parton-hadron cascade simulation ***loop-C***.

- In "main_23.f, parini_23.f", the low-energy simulation loop-A has been improved considering endothermic, exothermic processes and Delta particle decay more carefully. A COMMON BLOCK /delt/ has been added for processing Delta particle in loop-A.
- Parameter "mstj3" was renamed as "smtj3" and "smtj3/10" was set as the probability of the Delta particle decay. Parameter "dtt" is the ratio of the inelastic cross-section to the total one for hadron scattering now (D=0.85 for high energy, 0.1 for low-energy loop-A).
- In "main_23.f", two subroutines "share_p_PYJETS" and "share_p_sbe" are introduced to share the 4-momentum for the conservation adjustment. The array "throe_p" of COMMON BLOCK /sa16/ is used to store the lost 4-momentum now.
- In "main_23.f", the subroutine "pa2ent" has been fixed.
- In "remo" of "parini_23.f", the junctions were removed correctly now.
- In "coales_23.f", the quark energy sorting, i.e. the subroutine "eord", was closed.
- In "coales_23.f", a subroutine "PAPTDI" has been introduced for transverse momentum sampling. The subroutine "ffm" for high-energy deexcitation mechanism has been renamed and improved as "deexcitation" (two implementations "deexcitation_E" and "deexcitation_EP" for energy mode and light-cone variable mode).
- The subroutine "break_glu" of "coales_23.f" has been fixed to break-up gluons correctly.
- The subroutine "funcz" of "coales_23.f" was improved.
- The COMMON BLOCK /sa18/ was re-set as switches and parameters about coalescence.
- The "analy.f" has been improved for output the partial phase-space statistics. The forth distribution is transverse mass distribution now.

- ***In "parini_23.f", "parcas_23.f", "coales_23.f", and "hadcas_23.f", all of the "TAB character" were replaced by corresponding "space character" safely.***
- "usu.dat" updated.
- "PACIAE.sh" updated.
- Redundant statements deleted / modified.
- Bug fixed.

<!----------------------------------------------------------------------------->
#### <font color=red> 02/2023: </font> In version PACIAE 2.3 ###

- In "main_23.f, analy.f", the parameter "parp22" was used to select y/eta in partial phase-space statistics.
- In "main_23.f", The mistake-proofing statements were added.
- In subroutine "oscar" of "main_23.f", the OSCAR1997A/1999A/2013A output was re-wrote. 1997A for final particle information, 1999A for full event history (initial nucleon, initial parton state, final partonic state, initial hadron state and final particle state), and 2013A dummy up to now. Some related statements were added in "main_23.f" and "parini_23.f".
- In "parini_23.f", some pre-statements for statistical hadronization was added.
- In "coal_23.f", the separate position/momentum phase-space constraint was introduced.
- In "p_23.f", the size of COMMON BLOCK /HEPEVT/ was extended to 80000, also with local /PYCBLS/. In subroutine "PYTIME", DATE_AND_TIME was activated. (However, a potential bug is MSTU(5). Still 10000 in p_23.f now. A possible way may be using the compiling option "-fdefault-integer-8" when one uses GFortran.) The modification of PYTHIA6 were indicated at the beginning of the file.
- In "usu.dat", comments are updated and default values are given in (D=) (not yet tuned). In addition, the KF code list was attached at the end of the file directly, which was printed from PYTHIA6.
- A file "stahad_23.f" was introduced for the statistical hadronization in the future. Dummy up to now.
- "PACIAE.sh" updated.
- Redundant statements deleted.
- Bug fixed
- ***In "main_23.f", "sfm_23.f", "stahad_23.f", "analy.f", and "p_23.f", all of the "TAB character" were replaced by corresponding "space character" safely. Partly done in "parini_23.f", "parcas_23.f", and "hadcas_23.f"***

<!----------------------------------------------------------------------------->
### <font color=red> 01/2023: </font> In version PACIAE 2.3 ###

- In "main_23.f" and "parini_23.f", the low-energy simulation was introduced.
- "PACIAE.sh" updated.
- Bug fixed.

<!----------------------------------------------------------------------------->
### <font color=red> 12/2022: </font> In version PACIAE 2.3 ###

- In "parcas_23.f", the working COMMON BLOCK /parlist/ was replaced by /PYJETS/.
- Output display optimized for "rms.out".
- Redundant statements deleted.
- The "LICENSE" file was introduced. (GPL v2.0 based on PYTHIA 8)
- Bug fixed.

<!----------------------------------------------------------------------------->
### <font color=red> 11/2022: </font> In version PACIAE 2.3 ###

- In "parini_23.f", the selections of distributing nucleons into the overlapping region forcedly was introduced, controlled by adj(30).
- In "parcas_23.f", all of six quarks and gluons were included now. (u, d, s, g old)
- "eps09.f" has been removed now.
- Redundant statements deleted.
- Bug fixed.

<!----------------------------------------------------------------------------->
### <font color=red> 10/2022: </font> In version PACIAE 2.3 ###

- In subroutine "main" of "main_23.f", the real-time clock random seed was introduced.
- In subroutine "scat" of "parini_23.f", the long-written statement about executing the binary collision by calling PYEVNW / PYEVNT was replaced by a new subroutine "xevent".
- In "parini_23.f", the CME was introduced in PACIAE 2.3 based on Zhi-Lei's improvement in PACIAE 2.2.1b and PACIAE 2.2.1c. The statements about calculation of the eccentricity were corrected.
- In "main_23.f", "eps09.f" and "p_23.f", the statements with old syntax were re-wrote by more modern-style one.
- A shell-script file "PACIAE.sh" was introduced for automatic compilation, building and running (pseudo-parallel running possible).
- In "usu.dat", the gluon-splitting possibility were introduced for "coal_23.f".
- The "README.md" file was introduce.
- Bug fixed.

<!----------------------------------------------------------------------------->
### <font color=red> 04/2021: </font> In version PACIAE 2.2.1 b and c ###

- The CME was introduced by Zhi-Lei She et al.

<!----------------------------------------------------------------------------->
### <font color=red> 01/2020: </font> In version PACIAE 2.3 ###

- The EPS09 nuclear shadowing was introduced by Liang, whose subroutine is called "shanul_eps09". An extra file named "eps09.f" was added.