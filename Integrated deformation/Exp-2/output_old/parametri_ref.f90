 subroutine parametri(t_final,v_sx_x,v_dx_x,v_t_x,v_b_x,Rhoref,Cpref,Diff,Viscref,alpharef,thick_BC_sx,thick_BC_dx,&
 transition_BC_sx,transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_max,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,Hl,e0,&
 e1,e0v,e1v,rw,OH_0,alpham_dl,alpham_df,esp_m,esp_n,OUT_DIR,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,&
 Gy,Rgas,PF,phase_changes,idratazione,power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,OUT_FILE_nb,OUT_FILE_nl,ndof,&
 ndofT,ngauss,ngaussP,it_max,ip_max,marc_min,marc_max,dx,file_marc,file_grid,media,vel_corr,OUT_FILE_nl_vtu,OUT_FILE_nl_it,runge_4,&
 OUT_FILE_final_distribution,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,&
 NEW_FILE_GRID,NEW_FILE_MARKERS,MARKERS_DISTRIBUTION,MARKERS_CONSTANT,init_num_mat,marc_elem,left_v,right_v,top_v,bottom_v,left_T,&
 right_T,top_T,bottom_T,LAYER_X,LAYER_Y,thick_x,thick_y,n_elem_x,n_elem_y,WEAK_SEED,WEAK_TYPE,WEAK_POS_X,WEAK_POS_Y,WEAK_DIM_X,&
 WEAK_DIM_Y,THICK_AIR,file_world,WEAK_NUMBER,ph_maps_M,ph_maps_UC,ph_maps_LC,ph_maps_O,ph_maps_S,ph_maps_Sp,courant_number,&
 time_step,check_markers,eff_average,t_initial,shear,adiabatic,solver_T,layer_T,T_top,T_bot,y_top,y_bot,radio,cond,parabolic_T,&
 MARKERS_OUT,lagrangian_grid,dim_x_lg,dim_y_lg,old_lgrid,lgrid,serpentine,Press_0,H2O_serp,ocean,weak_plastic,xmin_lg,xmax_lg,&
 ymin_lg,ymax_lg,max_time_step,healing,healing_rate,thermal_act,OUT_MODEL,OUT_PHASE,DIR_ph_changes,H2O_ocean,delta_water,melt_ext,&
 penalty,element_type,tolM,iterM,pre_cond,norm_press_volume,air_layer,Temp_top,Temp_bottom,v_sx_y,v_dx_y,v_t_y,v_b_y,min_time_step,&
 random,directory,file,weakening,UCC,LCC,Mantle,Sed,Serp,UOC,LOC,relaxation,radio_exp,y_hr,v_dep,v_in_sx,v_in_dx,theta_stab)

 Implicit none

 logical, intent(out) :: phase_changes,idratazione,power_law,erosione,melting,equilib,OUT_VEL,OUT_TIME,OUT_TYPE,OUT_FILE_nl,&
 OUT_FILE_nl_vtu,OUT_FILE_nl_it,vel_corr,runge_4,OUT_FILE_final_distribution,NEW_FILE_GRID,NEW_FILE_MARKERS,MARKERS_CONSTANT,&
 left_T,right_T,top_T,bottom_T,WEAK_SEED,courant_number,check_markers,eff_average,shear,adiabatic,solver_T,parabolic_T,MARKERS_OUT,&
 lagrangian_grid,old_lgrid,serpentine,ocean,weak_plastic,healing,penalty,pre_cond,norm_press_volume,air_layer,random,radio_exp,v_dep

 character(len=*), intent(out) :: file_marc,file_grid,media,OUT_DIR,file_world,ph_maps_M,ph_maps_UC,ph_maps_LC,&
 ph_maps_O,ph_maps_S,ph_maps_Sp,lgrid,OUT_MODEL,OUT_PHASE,DIR_ph_changes,directory,file,weakening,UCC,LCC,Mantle,Sed,Serp,UOC,LOC

 integer, intent(out) :: ndof,ndofT,ngauss,ngaussP,it_max,ip_max,marc_min,marc_max,dx,OUT_FILE_nb,init_num_mat,marc_elem,&
 left_v,right_v,top_v,bottom_v,MARKERS_DISTRIBUTION,LAYER_X,LAYER_Y,WEAK_NUMBER,WEAK_TYPE,layer_T,element_type,iterM,relaxation

 double precision, intent(out) :: t_final,v_sx_x,v_dx_x,v_t_x,v_b_x,Rhoref,Diff,Viscref,alpharef,Cpref,thick_BC_sx,thick_BC_dx,&
 H2O_serp,transition_BC_sx,transition_BC_dx,OUT_FILE_tm,H2O_0,H2O_max,kappa_melt,gamma,chi,lambda,DH2O,A1,A2,A3,x_t,y_t,tol_melt,&
 Hl,e0,e1,e0v,e1v,rw,OH_0,alpham_dl,alpham_df,esp_m,esp_n,kappa_f,kappa_fsed,kappa_d,kappa_dsed,tol,rtau,pi,teta,eps,Cn,tolV,Gy,&
 Rgas,PF,g_coeff,gsed_coeff,p_coeff,sealevel,poro_s,poro_sh,zporo_s,zporo_sh,ratio,L_coeff,kappa_s,kappa_sh,THICK_AIR,time_step,&
 t_initial,dim_x_lg,dim_y_lg,Press_0,xmin_lg,xmax_lg,ymin_lg,ymax_lg,max_time_step,healing_rate,thermal_act,H2O_ocean,delta_water,&
 melt_ext,tolM,Temp_top,Temp_bottom,v_sx_y,v_dx_y,v_t_y,v_b_y,min_time_step,y_hr,v_in_sx,v_in_dx,theta_stab

 integer, dimension(:), allocatable :: n_elem_x,n_elem_y
 double precision, dimension(:), allocatable :: thick_x,thick_y,WEAK_POS_X,WEAK_POS_Y,WEAK_DIM_X,WEAK_DIM_Y,T_top,T_bot,y_top,&
 y_bot,radio,cond

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PARAMETRI MODIFICABILI DEL SISTEMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc Parametri generali cccccccccccccccccccccccccccccccc
 air_layer=.false.                                              !Strato di aria
 power_law=.false.                                            !Utilizzo power law
 solver_T=.false.                                             !Risoluzione dell'equazione dell'energia
 shear=.false.                                         !Utilizzo shear heating nell'equazione dell'energia
 adiabatic=.false.                                         !Utilizzo adiabatic heating nell'equazione dell'energia
 healing=.false.                                            !Calcolo strain healing
 phase_changes=.false.                                               !Utilizzo phase_changes per transizioni di fase
 idratazione=.false.                                             !Utilizzo idratazione
 serpentine=.false.                                             !Serpentinizzazione del mantello
 ocean=.false.                                             !Oceanizzazione del melting
 melting=.false.                                            !Utlizzo meccanismo di melting
 erosione=.false.                                            !Utilizzo meccanismo erosione/sedimentazione
 equilib=.false.                                       !Fase di riequilibrazione isostatica
 vel_corr=.true.                                       !Conservative Velocity Interpolation
 runge_4=.false.                                          !Utilizzo 4th-order Runge-Kutta
 norm_press_volume=.false.                                 !Normalizzazione pressione su volume totale
 check_markers=.false.                                     !Aggiunta/eliminazione marcatori - se .false. il numero di marcatori non viene modificato
 media='geometrica'                                        !Media usata per calcolo viscosità elementale
 file_world='Initial_stone.wb'                                   !Nome file World Builder

!cccccccccccccccccccccccccccc Parametri per il time step ccccccccccccccccccccccccccccc
 courant_number=.true.                                            !Utilizzo Courant number per time step
 Cn=0.2D0                                                  !Courant number - compreso tra 0 e 1
 time_step=5.D-03                                                 !Time step se Courant number=.false.
 max_time_step=5.D05                                                 !Massimo time step consentito
 min_time_step=0.D0                                                !Minimo time step consentito
 t_initial=0.D06                                               !Tempo iniziale
 t_final=151.D0                                           !Tempo finale - se OUT_TYPE=.true. in Myr - se OUT_TYPE=.false. in iterazioni massime

!ccccccccccccccccccccccccccccc Parametri per file output ccccccccccccccccccccccccccccc
 OUT_DIR='./output_12'              !Main Directory in cui scrivere i file di output
 OUT_MODEL=''                                                   !Nome della directory per il modello in corso
 OUT_PHASE=''                                         !Nome della directory per la fase evolutiva in corso
 OUT_VEL=.false.                                                       !Output velocità - .true. in cm/yr - .false. in m/s
 OUT_TIME=.false.                                                       !Output tempo - .true. in Myr - .false. in s
 OUT_TYPE=.false.                                                       !Intervallo file di output - .true. per tempo - .false. per numero
 OUT_FILE_nb=1                                                       !Numero di iterazioni per il file di output - se OUT_TYPE=.false.
 OUT_FILE_tm=0.1D06                                                       !Tempo di uscita per il file di output - se OUT_TYPE=.true.
 OUT_FILE_nl=.false.                                                     !Output file per ogni iterazione non-lineare
 OUT_FILE_nl_vtu=.false.                                                     !Output file vtu per ogni iterazione non-lineare
 OUT_FILE_nl_it=.false.                                                     !Output file con numero iterazioni non-lineari
 OUT_FILE_final_distribution=.false.                                        !Output file con distribuzione finale nodi e marcatori
 MARKERS_OUT=.false.                                                  !Output file binari per i marcatori

!cccccccccccccccccccccccc Parametri della griglia Lagrangiana cccccccccccccccccccccccc
 lagrangian_grid=.false.                                                  !Visualizzazione griglia Lagrangiana
 xmin_lg=300.D03                                                        !Dimensione lungo X griglia
 xmax_lg=900.D03                                                        !Dimensione lungo Y griglia
 ymin_lg=300.D03                                                        !Dimensione lungo X griglia
 ymax_lg=600.D03                                                        !Dimensione lungo Y griglia
 dim_x_lg=2.D03                                                        !Dimensione lungo X griglia
 dim_y_lg=2.D03                                                        !Dimensione lungo Y griglia
 old_lgrid=.false.                                                        !Lettura vecchio file griglia
 lgrid='../Output_models/Sesia/Isostatic/vtu_files/Grid.00492.vtu'      !Nome vecchio file griglia

!cccccccccccccccccccccccccccc Parametri fisici del sistema ccccccccccccccccccccccccccc
 Rhoref=3200.D0                                                  !Densità di riferimento
 Cpref=1250.D0                                                  !Calore specifico di riferimento
 Diff=1.D-06                                                  !Coefficiente di diffusione termica
 Viscref=1.D21                                                  !Viscosità di riferimento
 alpharef=3.1D-05                                                  !Coefficiente di espansione termica di riferimento
 Gy=-9.81D0                                                       !Gravità
 Press_0=0.D0                                                   !Pressione superficiale

!cccccccccccccccccccccccccccc Boundary condition velocità cccccccccccccccccccccccccccc
 v_dep=.false.                                            !Velocità dipendente dal tempo
 v_in_sx=-0.2d0                                               !Velocità iniziale lato sinistro
 v_in_dx=0.2d0                                            !Velocità iniziale lato destro

 !left border
 left_v=2                                                 !1-open; 2-free-slip; 3-fixed X and Y; 4-fixed X+free Y; 5-free X+fixed Y
 v_sx_x=0.D0                                                  !Velocità lungo x in cm/yr
 v_sx_y=0.D0                                                  !Velocità lungo y in cm/yr
 thick_BC_sx=1.D03                                          !Spessore litosfera
 transition_BC_sx=0.D03                                        !Spessore transizione buondary condition
 !right border
 right_v=2                                                 !1-open; 2-free-slip; 3-fixed X and Y; 4-fixed X+free Y; 5-free X+fixed Y
 v_dx_x=0.D0                                                  !Velocità lungo x in cm/yr
 v_dx_y=0.D0                                                  !Velocità lungo y in cm/yr
 thick_BC_dx=1.D03                                         !Spessore litosfera
 transition_BC_dx=0.D03                                         !Spessore transizione buondary condition
 !bottom
 bottom_v=2                                                !1-open; 2-free-slip; 3-fixed X and Y; 4-free X+fixed Y; 5-fixed X+free Y
 v_b_x=0.D0                                                  !Velocità lungo x in cm/yr
 v_b_y=0.D0                                                  !Velocità lungo y in cm/yr
 !top
 top_v=2                                                 !1-open; 2-free-slip; 3-fixed X and Y; 4-free X+fixed Y; 5-fixed X+free Y
 v_t_x=0.D0                                                  !Velocità lungo x in cm/yr
 v_t_y=0.D0                                                  !Velocità lungo y in cm/yr

!cccccccccccccccccccccccccc Boundary condition temperatura ccccccccccccccccccccccccccc
 left_T=.false.                                                 !BC lato sinistro
 right_T=.false.                                                 !BC lato destro
 top_T=.true.                                                 !BC top
 Temp_top=273.15D0                                                  !Temperatura fissata al top in K
 bottom_T=.true.                                                 !BC bottom
 Temp_bottom=1803.15D0                                                  !Temperatura fissata al bottom in K

!ccccccccccccccccccccccccccccc Produzione calore interno ccccccccccccccccccccccccccccc
 radio_exp=.false.                                         !Utilizzo produzione radiogenica decrescente con profondità
 y_hr=10.D03                                                 !Profondità decrescita esponenziale in m

!ccccccccccccccccccccccccccc Profilo iniziale temperatura cccccccccccccccccccccccccccc
 parabolic_T=.false.                                        !Incremento profilo termico iniziale non lineare
 LAYER_T=1                                                 !Numero livelli termici con profili non lineari
 allocate(T_top(layer_t),T_bot(layer_t),y_top(layer_t),y_bot(layer_t),radio(layer_t),cond(layer_t)) !!!!!!!!!NON MODIFICARE!!!!!!!!!
 T_top(1)=273.15                                             !Temperatura alla base del livello termico con profilo non lineare in K
 T_bot(1)=823.15                                             !Temperatura al top del livello termico con profilo non lineare in K
 y_top(1)=700.D03                                             !Coordinata della base del livello termico con profilo non lineare in m
 y_bot(1)=665.D03                                             !Coordinata del livello termico con profilo non lineare in m
 radio(1)=1.3D-06                                             !Calore radiogenico del livello termico con profilo non lineare
 cond(1)=2.5D0                                             !Conduttività del livello termico con profilo non lineare

!cccccccccccccccccccccccccccccc Parametri della griglia cccccccccccccccccccccccccccccc
 file_grid='Grid_stone_12.nc'               !Nome del file griglia
! file_grid='../../../../Output_models/test_falcon_13/00/Final_grid.nc'
 NEW_FILE_GRID=.true.                                          !File griglia - .true. crea nuovo file - .false. utilizza file esistente
!Modificare i parametri successivi solo se NEW_FILE_GRID=.true.
 LAYER_X=1                                                !Numero settori con diversa larghezza degli elementi lungo X
 LAYER_Y=1                                                !Numero settori con diversa altezza degli elementi lungo Y
 allocate(thick_x(layer_x),thick_y(layer_y),n_elem_x(layer_x),n_elem_y(layer_y)) !!!!!!!!!NON MODIFICARE!!!!!!!!!
 THICK_X(1)=1.D0                                          !Larghezza settore lungo X in metri
 N_ELEM_X(1)=40                                           !Numero elementi nel settore lungo X
 THICK_Y(1)=1.D0                                            !Altezza settore lungo Y in metri
 N_ELEM_Y(1)=40                                              !Numero elementi nel settore lungo Y

!cccccccccccccccccccccccccccccccc Parametri marcatori cccccccccccccccccccccccccccccccc
 marc_min=9                                             !Marcatori minimi per elemento di minore dimensione - solo se check_markers=.true.
 marc_max=20                                             !Marcatori massimi per elemento di minore dimensione - solo se check_markers=.true.
 THICK_AIR=0.D0                                               !Spessore sticky air
 file_marc='Markers_stone_12.nc'                                        !Nome del file marcatori
! file_marc='../../../../Output_models/test_falcon_13/00/Final_markers.nc'
 NEW_FILE_MARKERS=.true.                                  !File marcatori - .true. crea nuovo file - .false. utilizza file esistente
!Modificare i parametri successivi solo se NEW_FILE_MARKERS=.true.
 init_num_mat=1                                            !Numero iniziale dei materiali
 marc_elem=16                                              !Numero di marcatori per elemento di minore dimensione
 MARKERS_CONSTANT=.true.                                          !Numero marcatori costante per ogni elemento
 MARKERS_DISTRIBUTION=1                                  !Distribuzione marcatori - 1-regolare; 2-random; 3-random seguendo una griglia regolare

!cccccccccccccccccccccccccccccccc Parametri weak seed cccccccccccccccccccccccccccccccc
 WEAK_SEED=.false.                                                   !Esistenza di un weak seed
 WEAK_PLASTIC=.true.                                                  !.true. weak plastico - .false. weak viscoso
 WEAK_TYPE=1                                               !Valore strain iniziale - 1-Fixed; 2-Random; 3-Sinusoidale
 WEAK_NUMBER=1                                                    !Numero di weak seed
 allocate(weak_pos_x(weak_number),weak_pos_y(weak_number),weak_dim_x(weak_number),weak_dim_y(weak_number)) !!!!!!!!!NON MODIFICARE!!!!!!!!!
 WEAK_POS_X(1)=600.D03                                                         !Posizione X del weak seed in metri
 WEAK_POS_Y(1)=675.D03                                                         !Posizione Y del weak seed in metri
 WEAK_DIM_X(1)=10.D03                                                         !Dimensione X del weak seed in metri
 WEAK_DIM_Y(1)=10.D03                                                         !Dimensione Y del weak seed in metri

!cccccccccccccccccccccccccccccccc Parametri power-law cccccccccccccccccccccccccccccccc
 eff_average=.false.                                    !Calcolo viscosità effettiva markers - .true. media armonica - .false. taglio visc. min e max
 it_max=100                                                       !Iterazioni massime power-law
 tol=1.D-03                                                  !Tolleranza convergenza power-law
 e0=0.05D0                                                  !Strain iniziale per weakening plastico
 e1=1.05D0                                                  !Strain finale per weakening plastico
 e0v=5.D0                                                  !Strain iniziale per weakening viscoso
 e1v=10.D0                                                  !Strain finale per weakening viscoso

!ccccccccccccccccccccccccccccccccc Parametri healing ccccccccccccccccccccccccccccccccc
 healing_rate=1.D-14                                              !Velocità di strain healing
 thermal_act=100.D0                                               !Thermal activation constant

!ccccccccccccccccccccccccccccccccc Parametri densità ccccccccccccccccccccccccccccccccc
 delta_water=1.5D-02                                         !Coefficiente dipendenza della densità da H2O

!cccccccccccccccccccccccccccccccccc Parametri random ccccccccccccccccccccccccccccccccc
 random=.false.                                             !Utilizzo valori random per le proprietà dei materiali
 directory='Random_data/Random_parameters/'          !Directory con file random
 file='Random.081.txt'                               !Nome del file
 weakening='Weakening'                               !Directory con parametri random del weakening
 LCC='LCC'                               !Directory con parametri random della crosta continentale inferiore
 UCC='UCC'                               !Directory con parametri random del crosta continentale superiore
 Mantle='Mantle'                               !Directory con parametri random del mantello
 Sed='Sediments'                               !Directory con parametri random dei sedimenti
 Serp='Serpentine'                               !Directory con parametri random del serpentino
 LOC='LOC'                               !Directory con parametri random della crosta oceanica inferiore
 UOC='UOC'                               !Directory con parametri random della crosta oceanica superiore

!ccccccccccccccccccccccccccccccccc Parametri phase_changes cccccccccccccccccccccccccccccccc
 DIR_ph_changes='Phase_changes'
 ph_maps_M='pyrolite.txt'                          !Nome file phase_changes per il mantello anidro
 ph_maps_UC='ucc.txt'                          !Nome file phase_changes per il la crosta continentale superiore
 ph_maps_LC='lcc.txt'                          !Nome file phase_changes per la crosta continental inferiore
 ph_maps_O='wet_basalt_fixed.txt'                          !Nome file phase_changes per la crosta oceanica
 ph_maps_S='sediments.txt'                          !Nome file phase_changes per i sedimenti
 ph_maps_Sp='serpentinite.txt'                          !Nome file phase_changes per il mantello idrato

!ccccccccccccccccccccccccccccc Parametri per idratazione ccccccccccccccccccccccccccccc
 H2O_0=0.D0                                                  !Quantità iniziale di acqua nel mantello in wt%
 H2O_max=0.2D0                                                  !Quantità massima di acqua nel mantello in wt%
 H2O_serp=6.8D0                                                  !Quantità massima di acqua nel serpentino in wt%
 H2O_ocean=5.D0                                                  !Quantità massima di acqua nell'oceano in wt%
 rw=1.D0                                               !Esponente per la diminuzione della viscosità in presenza di acqua
 OH_0=0.062D0                                               !Limite quantità di acqua per comportamento anidro

!ccccccccccccccccccccccccccccccc Parametri per melting ccccccccccccccccccccccccccccccc
 kappa_melt=60.D0                                 !Costante per il calcolo di deltaT
 gamma=0.54D0                                  !Esponente per il calcolo di deltaT
 chi=12.D0                                    !Costante per il calcolo di acqua massima nel melt
 lambda=0.6D0                                   !Esponente per il calcolo di acqua massima nel melt
 DH2O=0.01D0                                   !Coefficiente di partizione dell'acqua
 A1=1120.66061D0                                !Coefficiente per calcolo del solidus - Kelley et al., 2010
 A2=132.899012D0                                !Coefficiente per calcolo del solidus - Kelley et al., 2010
 A3=-5.1404654D0                                !Coefficiente per calcolo del solidus - Kelley et al., 2010
 x_t=-221.34D0                                 !Coefficiente per calcolo della differenza solidus-liquidus - Kelley et al., 2010
 y_t=536.86D0                                 !Coefficiente per calcolo della differenza solidus-liquidus - Kelley et al., 2010
 tol_melt=1.D-02                                   !Tolleranza per la convergenza del root finder - Falsi point method
 Hl=3.D05                                   !Calore latente
 alpham_dl=-45.D0                                            !Esponente dislocation per la componente di melt
 alpham_df=-30.D0                                            !Esponente diffusion per la componente di melt
 melt_ext=10.D0                                        !Melt extraction in wt%

!cccccccccccccccccccccccccccccc Parametri per erosione ccccccccccccccccccccccccccccccc
!See https://fastscape.org/fastscapelib-fortran/#_fastscape_set_marine_parameters for a complete explanation of all parameters

 kappa_f=2.D-06                                        !Bedrock river incision (SPL) rate parameter in m/yr
 kappa_fsed=2.D-06                                     !Sediments river incision (SPL) rate parameter in m/yr
 esp_m=0.4D0                                                 !Drainage area exponent in the SPL
 esp_n=1.D0                                                 !Slope exponent in the SPL
 kappa_d=1.D-02                                                 !Bedrock transport coefficient (or diffusivity) for hillslope processes in m2/yr
 kappa_dsed=1.D-02                                                 !Sediments transport coefficient (or diffusivity) for hillslope processes in m2/yr
 g_coeff=0.5D0                                                  !Bedrock dimensionless deposition/transport coefficient for the enriched SPL
 gsed_coeff=0.5D0                                                  !Sediments dimensionless deposition/transport coefficient for the enriched SPL
 p_coeff=-1.D0                                                   !Slope exponent for multi-direction flow

!ccccccccccccccccccccccccccc Parametri per sedimentazione cccccccccccccccccccccccccccc
!See https://fastscape.org/fastscapelib-fortran/#_fastscape_set_marine_parameters for a complete explanation of all parameters

 sealevel=700.D03                                        !Sea level in meters
 poro_s=0.49D0                                           !Reference/surface porosity for sand
 poro_sh=0.63D0                                           !Reference/surface porosity for shale 
 zporo_s=3.7D03                                           !e-folding depth for exponential porosity law for sand
 zporo_sh=1.96D03                                           !e-folding depth for exponential porosity law for shale
 ratio=0.5D0                                              !Sand-shale ratio for material leaving the continent
 L_coeff=1.D02                                              !Averaging depth/thickness needed to solve the sand-shale equation in meters
 kappa_s=3.D02                                             !Marine transport coefficient (diffusivity) for sand in m2/yr
 kappa_sh=kappa_s                                             !Marine transport coefficient (diffusivity) for shale in m2/yr

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     PARAMETRI NON MODIFICABILI DEL SISTEMA
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 ndof=2                                                            !Gradi di libertà velocità
 ndofT=1                                                            !Gradi di libertà temperatura
 ngauss=4                                                            !Punti di integrazione velocità
 ngaussP=1                                                       !Punti di integrazione pressione
 penalty=.true.                                                    !Utilizzo penalty method
 pre_cond=.true.                                                    !Utilizzo preconditioner matrix
 element_type=1                                                   !Tipo elemento utilizzato - 1-Q1xP0; 2-Q1+xQ1; 3-Q2xQ1
 dx=11                                                       !Numero nodi markers-chain per elemento
 rtau=1.D0                                                  !Coefficiente correzione convezione
 pi=3.14159265359D0                                             !Pi greco
 teta=0.5D0                                                  !Coefficiente derivata temporale
 eps=10.D-10                                                  !Valore infinitesimale
 Rgas=8.314462618D0                                             !Costante gas
 ip_max=10                                                       !Iterazioni massime divergenza velocità
 tolV=1.D-02                                                  !Tolleranza convergenza divergenza velocità
 PF=1.D06                                                  !Penalty factor
 tolM=1.D-05                                                 !Tolleranza Uzawa mixed formulation
 iterM=500                                                      !Iterazioni massime mixed formulation
 relaxation=50                                             !Percentuale di memory relaxation per MUMPS
 theta_stab=0.5d0                                           !Valore stabilizzazione free surface

 end subroutine parametri
