!
!
!
!                            Human Health Risk Calculation Code   
!
!	 Exposure / Individual Risk/ Aggregate Risk/ Cumulative Risk / Radiological risk 
!
!    Non Cancer & Cancer Effects of Chemicals 
!
!	 Local-Time Matrix of Environmental Concentration, Exposure & Risk
!
!
      PROGRAM HHRISK
!
!
!    *****************************************************************************************************************
!    *****************************************************************************************************************
!
!	                              HEALTH HUMAN RISK MODEL 
!
!    
!
!         	                     ----- MODELO GERAL ------ 
!
!                                      Version 7.0
!        
!
!    *****************************************************************************************************************				 
!    *****************************************************************************************************************
!    *****************************************************************************************************************
!
!
!
!
!
!							INPUT-INPUT-INPUT
!			
!	  - INPUT files: 
!                     Datachemical.inp
!                     Dataexp.inp
!                     SCENARY.inp
!                     concentration.inp
!							
!    	
!                     ________________________________________________________________			   
!
!					----------------------------------------------------------------
!
!                             THE USE OF THIS PROGRAM MAY BE CONSULTED WITH:
!
!
!					        Jordan Brizi Neris(1), Lívia Correia &  Fermin G. Velasco(2)  
!
!
!      					   (1) FEDERAL UNIVERSITY OF SÃO CARLOS, SÃO CARLOS, SÃO PAULO
!                          (2) STATE UNIVERSITY OF SANTA CRUZ, ILHEUS, BAHIA
!								 EMAIL: fermin@uesc.br, jordanbrizi@gmail.com
!
!                     ------------------------------------------------------------------      
!                   	_________________________________________________________________
!            
!
!
!
!*********************************************************************************************************************                                                 
!*********************************************************************************************************************
!
!
! 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
	LOGICAL W,KEY_SD, MUTAGENIC,KEY_ANALYZE
    INTEGER SCENAR
	REAL*8 Kd,IR
    PARAMETER (NTYPECONC=15)   ! O 15 SÃO OS 15 TIPOS DE CONCENTRAÇOES !!!! ATENÇAO!!!!!!!!!
	PARAMETER (NIDADE=9)	 ! O 9 TIPO DE IDADES AVALIADAS !!!! ATENÇAO!!!!!!!!!
	PARAMETER (NVIAS=14)   	 ! O NÚMERO DE VIAS DE EXPOSIÇÃO AVALIADAS !!!! ATENÇAO!!!!!!!!!
    PARAMETER (NVP=0)   	 ! NVP = 0 	o file Variable_Parameters.out não será exibido, NVP = 1 	o file Variable_Parameters.out será exibido
    CHARACTER(LEN=50)  :: CHEMICAL(500),NOMES(500),POLLUTANT(500),RAD_POL(4),nomevias(NVIAS)
	CHARACTER(LEN=10)  :: KKINFORMATION(6),NRISKTYPE(NVIAS)
	CHARACTER(LEN=20)  :: TYPE_POLLUTANT(500)
	CHARACTER(LEN=14)  :: SCENARIES(4),TYPE_CHEMICAL(2)
!
	CHARACTER aspas*1
!
	REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: HQ
		REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: SD_HQ
	REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: CR
		REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: SD_CR
!
	REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: HIag
		REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: SD_HIag
	REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: CRag
		REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: SD_CRag
!
	REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: HIag_tot	
		REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SD_HIag_tot			  
	REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: CRag_tot
		REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SD_CRag_tot					   
!
	REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: HIcum
		REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SD_HIcum
	REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: CRcum
		REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SD_CRcum
	REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: HIag_ac
		REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: SD_HIag_ac
	REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: CRag_ac
		REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: SD_CRag_ac
	REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: HIcum_ac
		REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SD_HIcum_ac
	REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: CRcum_ac
		REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SD_CRcum_ac
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: HIcum_tot
		REAL*8, DIMENSION (:,:), ALLOCATABLE :: SD_HIcum_tot
	REAL*8, DIMENSION (:,:), ALLOCATABLE :: CRcum_tot
		REAL*8, DIMENSION (:,:), ALLOCATABLE :: SD_CRcum_tot
!
	REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: PORCENTAGEM
!
!
	REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: HQ_tot
		REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: CR_tot
!
	REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: PORC_VIA_NC
	    REAL*8, DIMENSION (:,:,:,:), ALLOCATABLE :: PORC_VIA_C
!
!	
!
    DIMENSION RfD(500,6,3),SF(500,6,3), W(NVIAS),BAF(500,14), SD_BAF(500,14) 
    DIMENSION SD_RfD(500,6,3),SD_SF(500,6,3) 
	DIMENSION BTF(500,15),SD_BTF(500,15)
	  DIMENSION ED(3,NIDADE), BW(NIDADE)
	  DIMENSION SD_ED(3,NIDADE), SD_BW(NIDADE)
	DIMENSION PC(500),SD_PC(500)
    DIMENSION ABS_(500), Kd(500),fw(500)
    DIMENSION SD_ABS(500), SD_Kd(500),SD_fw(500)
	DIMENSION CF(3),SD_CF(3),FI(9),SD_FI(9),Fa(4),SD_Fa(4),Fp(4),SD_Fp(4)
	DIMENSION IR(NIDADE,20),SD_IR(NIDADE,20),ET(NIDADE,3),SD_ET(NIDADE,3),SA(NIDADE,2),SD_SA(NIDADE,2)
	DIMENSION EV(NIDADE,2),SD_EV(NIDADE,2),AF(NIDADE),SD_AF(NIDADE)
	DIMENSION EF_INI(NVIAS,NIDADE),SD_EF_INI(NVIAS,NIDADE),AT_INI(2,NVIAS,NIDADE),SD_AT_INI(2,NVIAS,NIDADE)
	DIMENSION J_INI_AGE(NIDADE)
!
      DIMENSION MUTAGENIC(500,3),KEY_ANALYZE(3),AMOLAR(4),VIDAMEIA(4)
!
	DIMENSION TIMESP(0:1000) ! ACOMODA NUMERO DE PERIODOS DE TEMPO CONSIDERADOS  !   
!
!
!
!
!
!
!
!	  
!                                                 'NUMTIME' (NUMERO DE DIAS, MESES, ANOS NAS QUAIS FORAM FEITAS
!                                                 AS ANÁLISES, DARÁ O NÚMERO DE LINHAS DA MATRIZ "CONCEN" E 
!                                                 'NULOCAL' (NÚMERO DE LOCAIS ONDE FORAM COLETADAS AS AMOSTRAS.
!                                                 DARÁ O NÚMERO DE COLUNAS DA MATRIZ "CONCEN".	
!
!
!
!*********************************************************************************************************************	
!
!
!                                  APRESENTAÇÃO 
!
!
      write(*,*)
      write(*,'( //)')
      write(*,*)'    **************************************************'
      write(*,*)'    *                                                *'
      write(*,*)'    *                    HHRISK                      *'   
      write(*,*)'    *                                                *'
      write(*,*)'    ** *  ** * ** ** ** ** ** ** ** * ** * ** ** ** **'
      write(*,*)'    *                 Version  2.0                   *'
      write(*,*)'    *                                                *'
      write(*,*)'    *   Compiled with FORTRAN PowerStation v4.0      *'
      write(*,*)'    *             UFSCar- Sao Carlos/SP              *'
      write(*,*)'    *                 - May 2020-                    *'
      write(*,*)'    **************************************************'
      write(*,*)'    *                                                *' 
      write(*,*)'    *                 GENERAL MODEL                  *'
      write(*,*)'    *                                                *'
      write(*,*)'    *           --HEALTH HUMAN RISK MODEL--          *'
      write(*,*)'    *                                                *'
      write(*,*)'    *                    Authors:                    *'
      write(*,*)'    *                                                *'
      write(*,*)'    *          J. B. Neris & F. G.Velasco            *'
      write(*,*)'    *                                                *'
	  write(*,*)'    *                                                *'
      write(*,*)'    *                                                *'
      write(*,*)'    **************************************************'      
	  WRITE(*,*)
!
!
!
!
!
!
!*********************************************************************************************************************
!                                           BLOCO DE INICIALIZAÇÃO
!*********************************************************************************************************************
!    
!
!
!_____________________________________________________________________________________________________________________
!
!	                                  INICIALIZANDO MATRIZES DE CÁLCULO
!	
!
!		                                  INPUT -INPUT-INPUT-INPUT
!_____________________________________________________________________________________________________________________
!
	  open (UNIT=99, file='Results\WARNINGS.txt')
!
	  WRITE(99,'("********************************************************************************************************************************************")')
	  WRITE(99,'("                                                         !!!!!!!!!WARNINGS!!!!!!!!!")')
	  WRITE(99,'("********************************************************************************************************************************************")')
	  WRITE(99,*)
!
!
      CALL READSCENARY(NIDADE,NVIAS,SCENAR,NCHEM,NLOCAL,NDURATION,KEY_SD,KEY_ANALYZE,W,CF,SD_CF,FI,SD_FI,IR,SD_IR,FA,SD_FA,FP,SD_FP,&
	  ET,SD_ET,SA,SD_SA,AF,SD_AF,EV,SD_EV,NRISKTYPE,NINFORMATION,EF_INI,SD_EF_INI,AT_INI,SD_AT_INI)
!
!
      IF((SCENAR.EQ.4).AND.(KEY_ANALYZE(1).EQV..TRUE.))THEN
!
	  WRITE(99,'('' Scenario 4 - In Natura - does not allow calculation of human health non-radiological risks. '')')
	  WRITE(99,*)
	  WRITE(99,'('' The calculation of non-radiological risk will not be performed '')')
	  WRITE(99,*)
	  WRITE(99,'('' For more information, enable key --- Information = 1 --- in the Scenary input file '')')
	  WRITE(99,*)
      ENDIF 
!
!
      IF(NINFORMATION.EQ.1)THEN
      KKINFORMATION(6)='y'
!
      IF((KKINFORMATION(6).EQ.'y').OR.(KKINFORMATION(6).EQ.'Y').OR.(KKINFORMATION(6).EQ.'yes').OR.(KKINFORMATION(6).EQ.'Yes').OR.(KKINFORMATION(6).EQ.'YES'))THEN
!
      DO KUK=1,5
      KKINFORMATION(KUK)='YES'
	  ENDDO
!
      ENDIF
!
      CALL INFORMATION(KKINFORMATION)
	  ENDIF
!
!
      CALL READDATABASE (NVIAS,NIDADE,POLLUTANT,TYPE_POLLUTANT,RAD_POL,RfD,SF,SD_BTF,BTF,ED,AT_INI,PC,BW,ABS_,Kd,fw,NPOL,&
	  SD_RfD,SD_SF,SD_ED,SD_AT_INI,SD_PC,SD_ABS,SD_Kd,SD_fw,SD_BW,BAF,SD_BAF,MUTAGENIC,AMOLAR,VIDAMEIA)	   
!
!
      IF(SCENAR.EQ.1)THEN
	  NTIME_TEMPORANEO=MAX(ED(1,1),ED(1,2),ED(1,3),ED(1,4),ED(1,5),ED(1,6),ED(1,7),ED(1,8),ED(1,9))
!	  NAVARAGE=MAX(AT(1,1,1),AT(1,1,2),AT(1,1,3),AT(1,1,4),AT(1,1,5),AT(1,1,6),AT(1,1,7),AT(1,1,8),AT(1,1,9))
	  ELSEIF(SCENAR.EQ.2)THEN
	  NTIME_TEMPORANEO=MAX(ED(2,7),ED(2,8),ED(2,9))
!	  NAVARAGE=MAX(AT(1,2,7),AT(1,2,8),AT(1,2,9))
	  ELSEIF(SCENAR.EQ.3)THEN
	  NTIME_TEMPORANEO=MAX(ED(3,1),ED(3,2),ED(3,3),ED(3,4),ED(3,5),ED(3,6),ED(3,7),ED(3,8),ED(3,9))
!	  NAVARAGE=MAX(AT(1,3,1),AT(1,3,2),AT(1,3,3),AT(1,3,4),AT(1,3,5),AT(1,3,6),AT(1,3,7),AT(1,3,8),AT(1,3,9))
	  ELSEIF(SCENAR.EQ.4)THEN
	  NTIME_TEMPORANEO=NDURATION
!	  NAVARAGE=1.0
	  ENDIF
!
!	  
      N_VARIASAO=1
	  NTIME=NTIME_TEMPORANEO/N_VARIASAO
!
!
	  KOP=1
!
!
!      WRITE(*,*)NTIME,NDURATION
!
      IF(NTIME.LT.NDURATION)THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' The number of time concentrations given in the Concentration.inp file is larger than necessary. '')')
	  WRITE(*,*)
	  WRITE(*,'('' decrease the number of concentrations given in '')')
	  WRITE(*,'('' the concentration sheet '')')
	  WRITE(*,*)
	  WRITE(*,'('' ---------------------------------------------------------------------------------------------- '')')
	  WRITE(*,*)
	  WRITE(*,'('' For example: Maximum exposure duration (ED) value among the 9 provided (for each age group) '')')
	  WRITE(*,'(''              in dataexp sheet = 30 years '')')
	  WRITE(*,*)
	  WRITE(*,'(''              provided 31 concentrations in each matrix in concentration sheet '')')
	  WRITE(*,*)
	  WRITE(*,'(''              ERROR! '')')
	  WRITE(*,*)
	  WRITE(*,'(''              A maximum of 30 concentration should be provided in each concentration matrix '')')
!
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
!
	ALLOCATE (HQ(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL))
		ALLOCATE (SD_HQ(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL))
	ALLOCATE (CR(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL))
	    ALLOCATE (SD_CR(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL))		
	ALLOCATE (HIag(NIDADE,NCHEM,0:NTIME,NLOCAL))
	    ALLOCATE (SD_HIag(NIDADE,NCHEM,NTIME,NLOCAL))
	ALLOCATE (CRag(NIDADE,NCHEM,NTIME,NLOCAL))
		ALLOCATE (SD_CRag(NIDADE,NCHEM,NTIME,NLOCAL))
	ALLOCATE (HIcum(NIDADE,0:NTIME,NLOCAL))
		ALLOCATE (SD_HIcum(NIDADE,NTIME,NLOCAL))
	ALLOCATE (CRcum(NIDADE,NTIME,NLOCAL))
    	ALLOCATE (SD_CRcum(NIDADE,NTIME,NLOCAL))
!
	ALLOCATE (HIag_tot(NIDADE,NCHEM,NLOCAL))			             
    	ALLOCATE (SD_HIag_tot(NIDADE,NCHEM,NLOCAL))			               
	ALLOCATE (CRag_tot(NIDADE,NCHEM,NLOCAL))				           
		ALLOCATE (SD_CRag_tot(NIDADE,NCHEM,NLOCAL))				               
!
	ALLOCATE (HIag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL))
		ALLOCATE (SD_HIag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL))
	ALLOCATE (CRag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL))
	    ALLOCATE (SD_CRag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL))
	ALLOCATE (HIcum_ac(NIDADE,0:NTIME,NLOCAL))
		ALLOCATE (SD_HIcum_ac(NIDADE,0:NTIME,NLOCAL))
	ALLOCATE (CRcum_ac(NIDADE,0:NTIME,NLOCAL))
		ALLOCATE (SD_CRcum_ac(NIDADE,0:NTIME,NLOCAL))
	ALLOCATE (HIcum_tot(NIDADE,NLOCAL))
		ALLOCATE (SD_HIcum_tot(NIDADE,NLOCAL))
	ALLOCATE (CRcum_tot(NIDADE,NLOCAL))
		ALLOCATE (SD_CRcum_tot(NIDADE,NLOCAL))
!
	ALLOCATE (PORCENTAGEM(2,NIDADE,NCHEM,NTIME,NLOCAL))	
!
!
	ALLOCATE (HQ_tot(NIDADE,NVIAS,NTIME,NLOCAL))
	    ALLOCATE (CR_tot(NIDADE,NVIAS,NTIME,NLOCAL))
!
    ALLOCATE (PORC_VIA_NC(NIDADE,NVIAS,NTIME,NLOCAL))
	    ALLOCATE (PORC_VIA_C(NIDADE,NVIAS,NTIME,NLOCAL))
!
!
!
	  SCENARIES(1)='AGRICULTURAL'
	  SCENARIES(2)='INDUSTRIAL'
	  SCENARIES(3)='RESIDENTIAL'
	  SCENARIES(4)='IN NATURE'
!
	  TYPE_CHEMICAL(1)='NON-CUMULATIVE'
	  TYPE_CHEMICAL(2)='CUMULATIVE'
!
!
       PRINT*,'_____________________________________________________________'
       PRINT*,'   SCENARY           N. CHEMICALS    N. TIMES     N. LOCALS  '
	   PRINT*,'_____________________________________________________________'
       WRITE(*,'(3X,A13,9X,I4,10X,I4,9X,I4,10X,A11)') SCENARIES(SCENAR), NCHEM, NTIME,NLOCAL
	   PRINT*,'_____________________________________________________________'
	   WRITE(*,*)
	   WRITE(*,*)
!
!
      IF((KEY_ANALYZE(1).EQV..TRUE.).AND.(SCENAR.LE.3))THEN
!
      IF((SCENAR.EQ.1).OR.(SCENAR.EQ.3))THEN
	  INICIO=1
	  ELSEIF(SCENAR.EQ.2)THEN
	  INICIO=7
	  ENDIF
!
!
	  CALL EXPOSURE (W,CF,SD_CF,FI,SD_FI,IR,SD_IR,FA,SD_FA,FP,SD_FP,ET,SD_ET,SA,SD_SA,AF,SD_AF,EV,SD_EV,&
	  NVP,NIDADE,NVIAS,SCENAR,NCHEM,NTIME,NDURATION,VARIASAO,NLOCAL,NTYPECONC,HQ,SD_HQ,CR,SD_CR,CHEMICAL,TIMESP,NRISKTYPE,EF_INI,SD_EF_INI,AT_INI,SD_AT_INI,INICIO) 
!
!
!
      CALL RISCOag(NIDADE,NVIAS,CHEMICAL,POLLUTANT,TYPE_POLLUTANT,TYPE_CHEMICAL,NPOL,NCHEM,NTIME,NLOCAL,HQ,SD_HQ,CR,SD_CR,HIag,SD_HIag,CRag,SD_CRag,HIag_ac,SD_HIag_ac,&
	  CRag_ac,SD_CRag_ac,HIag_tot,SD_HIag_tot,CRag_tot,SD_CRag_tot,HQ_tot,CR_tot,ED,SCENAR,INICIO)
!
!
      CALL RISCOcum(NIDADE,NCHEM,NTIME,NLOCAL,HIag_ac,SD_HIag_ac,CRag,SD_CRag,HIcum,SD_HIcum,CRcum,SD_CRcum,&
	  HIcum_ac,SD_HIcum_ac,CRcum_ac,SD_CRcum_ac,HIcum_tot,SD_HIcum_tot,CRcum_tot,SD_CRcum_tot,ED,SCENAR,INICIO)
!
      aspas = '"'
!	  
      J_INI_AGE(1)=1
	  J_INI_AGE(2)=2
	  J_INI_AGE(3)=3
	  J_INI_AGE(4)=6
	  J_INI_AGE(5)=11
	  J_INI_AGE(6)=16
	  J_INI_AGE(7)=18
	  J_INI_AGE(8)=21
	  J_INI_AGE(9)=65
!
!
	  lstop0=0
	  istop0=0
	  jstop0=0
	  kstop0=0
	  lstop1=0
	  istop1=0
	  kstop1=0
!
      DO i=1,NCHEM
	  DO l=INICIO,NIDADE
	  IF(ED(SCENAR,l).NE.0.0)THEN
      NTIMEexp=ED(SCENAR,l)
	  DO K=1,NLOCAL
      DO j=1,NTIMEexp
!
	  IF(HIag(l,i,j,k).NE.HIag(l,i,j-1,k))THEN
	  lstop0=l
	  istop0=i
	  jstop0=j
	  kstop0=k
	  ENDIF
!
      ENDDO	 ! fim ciclo "j"
!
	  IF(CRag_tot(l,i,k).NE.0.0)THEN
	  lstop1=l
	  istop1=i
	  kstop1=k
	  ENDIF
!
	  ENDDO	  ! fim ciclo "k"
	  ENDIF
	  ENDDO	  ! fim ciclo "l"
	  ENDDO	  ! fim ciclo "i"
!
!
	  lstop2=0
	  jstop2=0
	  kstop2=0
	  lstop3=0
	  kstop3=0
!
	  DO l=INICIO,NIDADE
	  IF(ED(SCENAR,l).NE.0.0)THEN
      NTIMEexp=ED(SCENAR,l)
	  DO K=1,NLOCAL
      DO j=1,NTIMEexp
!
	  IF(HIcum(l,j,k).NE.HIcum(l,j-1,k))THEN
	  lstop2=l
	  jstop2=j
	  kstop2=k
	  ENDIF
!
	  ENDDO	  ! fim ciclo "j"
!
      IF(CRcum_tot(l,k).NE.0.0) THEN
	  lstop3=l
	  kstop3=k
	  ENDIF
!
	  ENDDO	  ! fim ciclo "k"
	  ENDIF
	  ENDDO	  ! fim ciclo "l"
!
!
      CALL NOMINATION(CHEMICAL,NCHEM,NOMES)
!
!*******************************************************************************************************************************************
!
!		                                            	PRINT -	PROTOCOL   
!
!********************************************************************************************************************************************
!
!
! 
	  open (UNIT=44, file='Results\Summary Aggregate Risk.json')
!
!
      WRITE(44,'("{")')
!
	  WRITE(44,'(A1,"Summarized aggregate non-carcinogenic risks",A1,": [")')aspas,aspas
!
!
	  LKCHEM=NCHEM
      NEW_POL=NPOL
!  
	  DO i=1,LKCHEM
!
!
!
      DO ii=1,NEW_POL
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
!
!
      DO l=INICIO,NIDADE
!
      IF(ED(SCENAR,l).NE.0.0)THEN
!
      NTIMEexp=ED(SCENAR,l)
!
!
	  DO K=1,NLOCAL
!
	  IF((SCENAR.EQ.2).AND.(l.LE.6))THEN
	  HIag_tot(l,i,k)=0.0
	  SD_HIag_tot(l,i,k)=0.0
	  CRag_tot(l,i,k)=0.0
	  SD_CRag_tot(l,i,k)=0.0
	  ENDIF
!
!
      DO j=1,NTIMEexp
!
!
      IF(l.EQ.1)THEN
	  N_AGE=j
      ELSEIF(l.EQ.2)THEN
	  N_AGE=j+1
      ELSEIF(l.EQ.3)THEN
	  N_AGE=j+2
      ELSEIF(l.EQ.4)THEN
	  N_AGE=j+5
      ELSEIF(l.EQ.5)THEN
	  N_AGE=j+10
      ELSEIF(l.EQ.6)THEN
	  N_AGE=j+15
      ELSEIF(l.EQ.7)THEN
	  N_AGE=j+17
      ELSEIF(l.EQ.8)THEN
	  N_AGE=j+20
      ELSEIF(l.EQ.9)THEN
	  N_AGE=j+64
	  ENDIF
!
!
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
	  IF(HIag(l,i,j,k).NE.HIag(l,i,j-1,k))THEN
!
      WRITE(44,'("{")')
!
      write(44,'(A1,"Chemical species",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NOMES(i),aspas
      write(44,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(44,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(44,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(44,'(A1,"HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIag(l,i,j,k)
      write(44,'(A1,"HI error",A1,":",1x,ES12.5) )') aspas,aspas,SD_HIag(l,i,j,k)
!
!      write(*,*) j,k,l
!
      IF((j.EQ.jstop0).and.(k.EQ.kstop0).and.(l.EQ.lstop0).and.(i.EQ.istop0))THEN
	  WRITE(44,'("}")')
	  ELSE
	  WRITE(44,'("},")')
	  ENDIF
!
	  ENDIF		
!
!
	  ELSE
!
	  IF(HIag(l,i,j,k).NE.HIag(l,i,j-1,k))THEN
!
      WRITE(44,'("{")')
!
      write(44,'(A1,"Chemical species",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NOMES(i),aspas
      write(44,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l) 
      write(44,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(44,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(44,'(A1,"HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIag(l,i,j,k)
      write(44,'(A1,"HI error",A1,":",1X,"null") )') aspas,aspas
!
!
      IF((j.EQ.jstop0).and.(k.EQ.kstop0).and.(l.EQ.lstop0).and.(i.EQ.istop0))THEN
	  WRITE(44,'("}")')
	  ELSE
	  WRITE(44,'("},")')
	  ENDIF
!
	  ENDIF 
!
      ENDIF
!
!
	  ENDDO	   !FIM DO CICLO "j"
!	             
!
	  ENDDO	   !FIM DO CICLO "K"
!
      ENDIF ! FIM do IF ED(SCENAR,l)....	             
!
      ENDDO  ! FECHA O "DO" DE l
!
	  ENDDO	  ! FIM DO CICLO "i"
!	       
	  WRITE(44,'("],")')      
!

!     tabela carcinogenicos
!
	  WRITE(44,'(A1,"Summarized aggregate carcinogenic risks",A1,": [")')aspas,aspas
!
	  LKCHEM=NCHEM
      NEW_POL=NPOL
!  
	  DO i=1,LKCHEM
!
      DO ii=1,NEW_POL
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
!
!
      DO l=INICIO,NIDADE
!
      IF(l.EQ.1)THEN
	  M_AGE=NTIMEexp
      ELSEIF(l.EQ.2)THEN
	  M_AGE=NTIMEexp+1
      ELSEIF(l.EQ.3)THEN
	  M_AGE=NTIMEexp+2
      ELSEIF(l.EQ.4)THEN
	  M_AGE=NTIMEexp+5
      ELSEIF(l.EQ.5)THEN
	  M_AGE=NTIMEexp+10
      ELSEIF(l.EQ.6)THEN
	  M_AGE=NTIMEexp+15
      ELSEIF(l.EQ.7)THEN
	  M_AGE=NTIMEexp+17
      ELSEIF(l.EQ.8)THEN
	  M_AGE=NTIMEexp+20
      ELSEIF(l.EQ.9)THEN
	  M_AGE=NTIMEexp+64
	  ENDIF
!
      IF(ED(SCENAR,l).NE.0.0)THEN
!
      NTIMEexp=ED(SCENAR,l)
!
!
	  DO K=1,NLOCAL
!
!
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
	  IF(CRag_tot(l,i,k).NE.0.0)THEN
!
	  WRITE(44,'("{")')
!
      write(44,'(A1,"Chemical species",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NOMES(i),aspas
      write(44,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l) 
      write(44,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(44,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,M_AGE
      write(44,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRag_tot(l,i,k)
      write(44,'(A1,"CR error",A1,":",1x,ES12.5) )') aspas,aspas,SD_CRag_tot(l,i,k)
!
      IF((k.EQ.kstop1).and.(l.EQ.lstop1).and.(i.EQ.istop1))THEN
	  WRITE(44,'("}")')	             
	  ELSE
	  WRITE(44,'("},")')	             
	  ENDIF
!
!
      ENDIF		
!
	  ELSE
!
	  IF(CRag_tot(l,i,k).NE.0.0)THEN
!
	  WRITE(44,'("{")')
!
      write(44,'(A1,"Chemical species",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NOMES(i),aspas
      write(44,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l) 
      write(44,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(44,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,M_AGE
      write(44,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRag_tot(l,i,k)
      write(44,'(A1,"CR error",A1,":",1X,"null") )') aspas,aspas
!
      IF((k.EQ.kstop1).and.(l.EQ.lstop1).and.(i.EQ.istop1))THEN
	  WRITE(44,'("}")')	             
	  ELSE
	  WRITE(44,'("},")')	             
	  ENDIF  
!
	  ENDIF 
!
      ENDIF
!
!
	  ENDDO	   !FIM DO CICLO "K"
!
      ENDIF ! FIM do IF ED(SCENAR,l)....	             
!
      ENDDO  ! FECHA O "DO" DE l
!
	  ENDDO	  ! FIM DO CICLO "i"
!
	  WRITE(44,'("]")')
!
	  WRITE(44,'("}")')
!
!
      CLOSE (44)
!
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!
!
	  open (UNIT=37, file='Results\Extensive Aggregate Risk.json')


!
!
!		 *******************************************************************************************
!!
!			PRINT -	PROTOCOL    AGREGATE SP AND ACUMULATE RISK
!
!		*******************************************************************************************
!
      WRITE(37,'("{")')
!
!
      IICHEM=NCHEM
	  NEW_POLI=NPOL
!  
	  DO i=1,IICHEM
!
	  WRITE(37,'(A1,"HI and CR values for",1x,A10,A1,": [")')aspas,NOMES(i),aspas
!
!
      DO ii=1,NEW_POLI
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
!
!
      DO l=INICIO,NIDADE
!
      IF(ED(SCENAR,l).NE.0.0)THEN
!
      NTIMEexp=ED(SCENAR,l)
!																						
!
	  DO k=1,NLOCAL
	  DO j=1,NTIMEexp
!
	  WRITE(37,'("{")')
!
      IF(l.EQ.1)THEN
	  N_AGE=j
      ELSEIF(l.EQ.2)THEN
	  N_AGE=j+1
      ELSEIF(l.EQ.3)THEN
	  N_AGE=j+2
      ELSEIF(l.EQ.4)THEN
	  N_AGE=j+5
      ELSEIF(l.EQ.5)THEN
	  N_AGE=j+10
      ELSEIF(l.EQ.6)THEN
	  N_AGE=j+15
      ELSEIF(l.EQ.7)THEN
	  N_AGE=j+17
      ELSEIF(l.EQ.8)THEN
	  N_AGE=j+20
      ELSEIF(l.EQ.9)THEN
	  N_AGE=j+64
	  ENDIF
!
	  IF((SCENAR.EQ.2).AND.(l.LE.6))THEN
	  HIag(l,i,j,k)=0.0
	  SD_HIag(l,i,j,k)=0.0
	  HIag_ac(l,i,j,k)=0.0
	  SD_HIag_ac(l,i,j,k)=0.0
	  CRag(l,i,j,k)=0.0
	  SD_CRag(l,i,j,k)=0.0
	  CRag_ac(l,i,j,k)=0.0
	  SD_CRag_ac(l,i,j,k)=0.0
	  ENDIF
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
	  IF(CRag(l,i,j,k).NE.0.0)THEN
!
      write(37,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(37,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(37,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(37,'(A1,"HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIag(l,i,j,k)
      write(37,'(A1,"HI error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_HIag(l,i,j,k)
      write(37,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRag(l,i,j,k)
      write(37,'(A1,"CR error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_CRag(l,i,j,k)
      write(37,'(A1,"Accumulated CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRag_ac(l,i,j,k)
      write(37,'(A1,"Accumulated CR error",A1,":",1x,ES12.5) )') aspas,aspas,SD_CRag_ac(l,i,j,k)
!
	  ELSE
!
      write(37,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(37,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(37,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(37,'(A1,"HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIag(l,i,j,k)
      write(37,'(A1,"HI error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_HIag(l,i,j,k)
      write(37,'(A1,"CR value",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"CR error",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"Accumulated CR value",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"Accumulated CR error",A1,":",1x,"null") )') aspas,aspas
!
	  ENDIF
!
	  ELSE
!
!
      IF(CRag(l,i,j,k).NE.0.0)THEN
!
      write(37,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(37,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(37,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(37,'(A1,"HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIag(l,i,j,k)
      write(37,'(A1,"HI error",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRag(l,i,j,k)
      write(37,'(A1,"CR error",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"Accumulated CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRag_ac(l,i,j,k)
      write(37,'(A1,"Accumulated CR error",A1,":",1x,"null") )') aspas,aspas
!
	  ELSE
!					 
      write(37,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(37,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(37,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(37,'(A1,"HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIag(l,i,j,k)
      write(37,'(A1,"HI error",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"CR value",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"CR error",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"Accumulated CR value",A1,":",1x,"null",",") )') aspas,aspas
      write(37,'(A1,"Accumulated CR error",A1,":",1x,"null") )') aspas,aspas
!
      ENDIF
!
	  ENDIF
!
!
!
      IF((j.EQ.NTIMEexp).and.(K.EQ.NLOCAL).and.(l.EQ.lstop0))THEN
	  WRITE(37,'("}")')	             
	  ELSE
	  WRITE(37,'("},")')	             
	  ENDIF
!
	  ENDDO	   !FIM DO CICLO "j"
!
	  ENDDO	   !FIM DO CICLO "K"
!
      ENDIF ! FIM do IF ED(SCENAR,l)....	             
!
      ENDDO  ! FECHA O "DO" DE l
!
      IF(i.EQ.NCHEM)THEN
	  WRITE(37,'("]")')
	  ELSE
	  WRITE(37,'("],")')
	  ENDIF
!
	  ENDDO	  ! FIM DO CICLO "i"
!
!
	  WRITE(37,'("}")')

!
      CLOSE (37)
!
!
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!
!
!
	  open (UNIT=38, file='Results\Extensive Cumulative Risk.json')


!
!
!		 *******************************************************************************************
!!
!			PRINT -	PROTOCOL    AGREGATE SP AND ACUMULATE RISK
!
!		*******************************************************************************************
!
      WRITE(38,'("{")')
!
	  WRITE(38,'(A1,"Carcinogenic and non-carcinogenic risks in extensive form",A1,": [")')aspas,aspas

!
!
      DO l=INICIO,NIDADE
!
      IF(ED(SCENAR,l).NE.0.0)THEN
!
      NTIMEexp=ED(SCENAR,l)
!																						
!
	  DO k=1,NLOCAL
	  DO j=1,NTIMEexp
!
	  WRITE(38,'("{")')
!
!
      IF(l.EQ.1)THEN
	  N_AGE=j
      ELSEIF(l.EQ.2)THEN
	  N_AGE=j+1
      ELSEIF(l.EQ.3)THEN
	  N_AGE=j+2
      ELSEIF(l.EQ.4)THEN
	  N_AGE=j+5
      ELSEIF(l.EQ.5)THEN
	  N_AGE=j+10
      ELSEIF(l.EQ.6)THEN
	  N_AGE=j+15
      ELSEIF(l.EQ.7)THEN
	  N_AGE=j+17
      ELSEIF(l.EQ.8)THEN
	  N_AGE=j+20
      ELSEIF(l.EQ.9)THEN
	  N_AGE=j+64
	  ENDIF
!
	  IF((SCENAR.EQ.2).AND.(l.LE.6))THEN
	  HIcum(l,j,k)=0.0
	  SD_HIcum(l,j,k)=0.0
	  HIcum_ac(l,j,k)=0.0
	  SD_HIcum_ac(l,j,k)=0.0
	  CRcum(l,j,k)=0.0
	  SD_CRcum(l,j,k)=0.0
	  CRcum_ac(l,j,k)=0.0
	  SD_CRcum_ac(l,j,k)=0.0
	  ENDIF
!
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
      IF(CRcum(l,j,k).NE.0.0)THEN
!
      write(38,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(38,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(38,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(38,'(A1,"HItot value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIcum(l,j,k)
      write(38,'(A1,"HItot error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_HIcum(l,j,k)
      write(38,'(A1,"CRcum value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRcum(l,j,k)
      write(38,'(A1,"CRcum error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_CRcum(l,j,k)
      write(38,'(A1,"Accumulated CRcum value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRcum_ac(l,j,k)
      write(38,'(A1,"Accumulated CRcum error",A1,":",1x,ES12.5))') aspas,aspas,SD_CRcum_ac(l,j,k)
!
	  ELSE
!
      write(38,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(38,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(38,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(38,'(A1,"HItot value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIcum(l,j,k)
      write(38,'(A1,"HItot error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_HIcum(l,j,k)
      write(38,'(A1,"CRcum value",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"CRcum error",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"Accumulated CRcum value",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"Accumulated CRcum error",A1,":",1x,"null"))') aspas,aspas
!
	  ENDIF
!
	  ELSE
!
      IF(CRcum(l,j,k).NE.0.0)THEN
!
      write(38,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(38,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(38,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(38,'(A1,"HItot value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIcum(l,j,k)
      write(38,'(A1,"HItot error",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"CRcum value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRcum(l,j,k)
      write(38,'(A1,"CRcum error",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"Accumulated CRcum value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRcum_ac(l,j,k)
      write(38,'(A1,"Accumulated CRcum error",A1,":",1x,"null"))') aspas,aspas
!
	  ELSE
!
      write(38,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(38,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(38,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(38,'(A1,"HItot value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIcum(l,j,k)
      write(38,'(A1,"HItot error",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"CRcum value",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"CRcum error",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"Accumulated CRcum value",A1,":",1x,"null",",") )') aspas,aspas
      write(38,'(A1,"Accumulated CRcum error",A1,":",1x,"null"))') aspas,aspas
!
	  ENDIF
!
	  ENDIF
!
!
      IF((j.EQ.NTIMEexp).and.(K.EQ.NLOCAL).and.(l.EQ.lstop2))THEN
	  WRITE(38,'("}")')	             
	  ELSE
	  WRITE(38,'("},")')	             
	  ENDIF
!
	  ENDDO	   !FIM DO CICLO "j"
!
	  ENDDO	   !FIM DO CICLO "K"
!
      ENDIF ! FIM do IF ED(SCENAR,l)....	             
!
      ENDDO  ! FECHA O "DO" DE l
!
	  WRITE(38,'("]")')
!
	  WRITE(38,'("}")')

!
      close(38)
!
!
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!*******************************************************************************************************************************************
!
!
!		*******************************************************************************************
!!
!			PRINT -	PROTOCOL    Summary Cumulative RISK
!
!		*******************************************************************************************
!
!
	  open (UNIT=66, file='Results\Summary Cumulative Risk.json')																				 																				 
!
!
      WRITE(66,'("{")')
!
	  WRITE(66,'(A1,"Summarized cumulative non-carcinogenic risks",A1,": [")')aspas,aspas
!
!
      DO l=INICIO,NIDADE
!
      IF(ED(SCENAR,l).NE.0.0)THEN
!
	  NTIMEexp=ED(SCENAR,l)
!
!
	  DO k=1,NLOCAL
	  DO j=1,NTIMEexp
!
      IF(l.EQ.1)THEN
	  N_AGE=j
      ELSEIF(l.EQ.2)THEN
	  N_AGE=j+1
      ELSEIF(l.EQ.3)THEN
	  N_AGE=j+2
      ELSEIF(l.EQ.4)THEN
	  N_AGE=j+5
      ELSEIF(l.EQ.5)THEN
	  N_AGE=j+10
      ELSEIF(l.EQ.6)THEN
	  N_AGE=j+15
      ELSEIF(l.EQ.7)THEN
	  N_AGE=j+17
      ELSEIF(l.EQ.8)THEN
	  N_AGE=j+20
      ELSEIF(l.EQ.9)THEN
	  N_AGE=j+64
	  ENDIF
!
	  IF((SCENAR.EQ.2).AND.(l.LE.6))THEN
	  HIcum(l,j,k)=0.0
	  SD_HIcum(l,j,k)=0.0
	  HIcum_ac(l,j,k)=0.0
	  SD_HIcum_ac(l,j,k)=0.0
	  CRcum(l,j,k)=0.0
	  SD_CRcum(l,j,k)=0.0
	  CRcum_ac(l,j,k)=0.0
	  SD_CRcum_ac(l,j,k)=0.0
	  ENDIF
!
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
!
	  IF(HIcum(l,j,k).NE.HIcum(l,j-1,k))THEN
!
      WRITE(66,'("{")')
!
      write(66,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(66,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(66,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(66,'(A1,"HItot value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIcum(l,j,k)
      write(66,'(A1,"HItot error",A1,":",1x,ES12.5) )') aspas,aspas,SD_HIcum(l,j,k)
!
      IF((j.EQ.jstop2).and.(k.EQ.kstop2).and.(l.EQ.lstop2))THEN
	  WRITE(66,'("}")')
	  ELSE
	  WRITE(66,'("},")')
	  ENDIF
!
	  ENDIF
!
	  ELSE
!
	  IF(HIcum(l,j,k).NE.HIcum(l,j-1,k))THEN
!
      WRITE(66,'("{")')
!
      write(66,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(66,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(66,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(66,'(A1,"HItot value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIcum(l,j,k)
      write(66,'(A1,"HItot error",A1,":",1x,"null") )') aspas,aspas
!
      IF((j.EQ.jstop2).and.(k.EQ.kstop2).and.(l.EQ.lstop2))THEN
	  WRITE(66,'("}")')
	  ELSE
	  WRITE(66,'("},")')
	  ENDIF
!
	  ENDIF
!
	  ENDIF
!
	  ENDDO	   !FIM DO CICLO "j"
!	             
!
	  ENDDO	   !FIM DO CICLO "K"
!
      ENDIF ! FIM do IF ED(SCENAR,l)....	             
!
      ENDDO  ! FECHA O "DO" DE l
!
!	       
	  WRITE(66,'("],")')      
!
!     tabela carcinogenicos
!
	  WRITE(66,'(A1,"Summarized cumulative carcinogenic risks",A1,": [")')aspas,aspas
!
!
      DO l=INICIO,NIDADE
!
!
      IF(l.EQ.1)THEN
	  M_AGE=NTIMEexp
      ELSEIF(l.EQ.2)THEN
	  M_AGE=NTIMEexp+1
      ELSEIF(l.EQ.3)THEN
	  M_AGE=NTIMEexp+2
      ELSEIF(l.EQ.4)THEN
	  M_AGE=NTIMEexp+5
      ELSEIF(l.EQ.5)THEN
	  M_AGE=NTIMEexp+10
      ELSEIF(l.EQ.6)THEN
	  M_AGE=NTIMEexp+15
      ELSEIF(l.EQ.7)THEN
	  M_AGE=NTIMEexp+17
      ELSEIF(l.EQ.8)THEN
	  M_AGE=NTIMEexp+20
      ELSEIF(l.EQ.9)THEN
	  M_AGE=NTIMEexp+64
	  ENDIF
!
      IF(ED(SCENAR,l).NE.0.0)THEN
!
      NTIMEexp=ED(SCENAR,l)
!
!
	  DO K=1,NLOCAL
!
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
      IF(CRcum_tot(l,k).NE.0.0) THEN
!
      WRITE(66,'("{")')
!
      write(66,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l) 
      write(66,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(66,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,M_AGE
      write(66,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRcum_tot(l,k)
      write(66,'(A1,"CR error",A1,":",1x,ES12.5) )') aspas,aspas,SD_CRcum_tot(l,k)
!
      IF((k.EQ.kstop3).and.(l.EQ.lstop3))THEN
	  WRITE(66,'("}")')	             
	  ELSE
	  WRITE(66,'("},")')	             
	  ENDIF
!
	  ENDIF
!
	  ELSE
!
      IF(CRcum_tot(l,k).NE.0.0)THEN
!
      WRITE(66,'("{")')
!
      write(66,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l) 
      write(66,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(66,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,M_AGE
      write(66,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CRcum_tot(l,k)
      write(66,'(A1,"CR error",A1,":",1x,"null") )') aspas,aspas
!
      IF((k.EQ.kstop3).and.(l.EQ.lstop3))THEN
	  WRITE(66,'("}")')	             
	  ELSE
	  WRITE(66,'("},")')	             
	  ENDIF
!
	  ENDIF
!
	  ENDIF
!
!
	  ENDDO
!
      ENDIF ! FIM do IF ED(SCENAR,l)....
!
      ENDDO    ! FIM DO CICLO "l"
!
	  WRITE(66,'("]")')	             
!
	  WRITE(66,'("}")')	             
!
      close(66)
!
!
!
!
!
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
!									   NOVO DADOS DE SAÍDA (ANÁLISES COMPLEMENTARES)
!
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
	  open (UNIT=55, file='Results\Complementary Analyzes.json')
!
!
!
	  NTIMEexp=NTIME
!
!
!
	  JIJ=NCHEM
!
!	  
      DO l=1,NIDADE
	  DO i=1,JIJ
	  DO JOP=1,NTIMEexp
	  DO KP=1,NLOCAL
	  IF(CRcum(l,JOP,KP).NE.0.0)THEN
      PORCENTAGEM(1,l,i,JOP,KP)=(CRag(l,i,JOP,KP)/CRcum(l,JOP,KP))*100
	  ELSE
	  PORCENTAGEM(1,l,i,JOP,KP)=0.0
	  ENDIF
	  IF(HIcum(l,JOP,KP).NE.0.0)THEN
      PORCENTAGEM(2,l,i,JOP,KP)=(HIag(l,i,JOP,KP)/HIcum(l,JOP,KP))*100
	  ELSE
	  PORCENTAGEM(2,l,i,JOP,KP)=0.0
	  ENDIF
	  ENDDO
	  ENDDO
	  ENDDO
	  ENDDO
!
	  WRITE(55,'("{")')	 
!
	  WRITE(55,'(A1,"Non-carcinogenic risks contribution (%) of each chemical species",A1,": [")')aspas,aspas            
! 
!
      DO KO=1,NLOCAL
!	
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"1 to <2",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,1,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,1,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"2 to <3",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,2,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,2,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"3 to <6",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,3,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,3,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"6 to <11",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,4,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,4,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"11 to <16",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,5,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,5,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"16 to <18",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,6,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,6,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"18 to <21",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,7,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,7,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"21 to <65",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,8,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,8,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
      ifinal=64+NTIMEexp
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"65 to ",I3,A1,",")') aspas,aspas,aspas,ifinal,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,9,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(2,9,IOO,1,KO)
	  ENDIF																					
	  ENDDO
	  IF(KO.NE.NLOCAL)THEN
	  WRITE(55,'("},")')!
	  ELSE
	  WRITE(55,'("}")')!
	  ENDIF
!
	  ENDDO	  ! FIM DO KO=1,NLOCAL
!
	  WRITE(55,'("],")')
!
!
!
!
	  WRITE(55,'(A1,"Carcinogenic risks contribution (%) of each chemical species",A1,": [")')aspas,aspas            
!
      DO KO=1,NLOCAL
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"1 to <2",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,1,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,1,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"2 to <3",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,2,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,2,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"3 to <6",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,3,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,3,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"6 to <11",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,4,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,4,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"11 to <16",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,5,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,5,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"16 to <18",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,6,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,6,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"18 to <21",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,7,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,7,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"21 to <65",A1,",")') aspas,aspas,aspas,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,8,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,8,IOO,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
      ifinal=64+NTIMEexp
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"65 to ",I3,A1,",")') aspas,aspas,aspas,ifinal,aspas
	  DO IOO=1,JIJ
	  IF(IOO.NE.JIJ)THEN
      write(55,'(A1,A10,A1,":",1X,ES10.3,",")') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,9,IOO,1,KO)
	  ELSE
      write(55,'(A1,A10,A1,":",1X,ES10.3)') aspas,NOMES(IOO),aspas,PORCENTAGEM(1,9,IOO,1,KO)
	  ENDIF																					
	  ENDDO
	  IF(KO.NE.NLOCAL)THEN
	  WRITE(55,'("},")')!
	  ELSE
	  WRITE(55,'("}")')!
	  ENDIF
!
	  ENDDO	  ! FIM DO KO=1,NLOCAL
!
	  WRITE(55,'("],")')
!
!
!
	  WRITE(55,'(A1,"Contribution (%) of each pathway to the non-carcinogenic risks",A1,": [")')aspas,aspas            
!
!
      DO l=1,NIDADE
      DO n=1,NVIAS
	  DO j=1,NTIMEexp
	  DO k=1,NLOCAL
!
	  IF(HIcum(l,j,k).NE.0.0)THEN
      PORC_VIA_NC(l,n,j,k)=(HQ_tot(l,n,j,k)/HIcum(l,j,k))*100
	  ELSE
	  PORC_VIA_NC(l,n,j,k)=0.0
	  ENDIF
	  IF(CRcum(l,j,k).NE.0.0)THEN
	  PORC_VIA_C(l,n,j,k)=(CR_tot(l,n,j,k)/CRcum(l,j,k))*100
	  ELSE
	  PORC_VIA_C(l,n,j,k)=0.0
	  ENDIF
!
      ENDDO
	  ENDDO
	  ENDDO
	  ENDDO
!
      nomevias(1)='Soil ingestion'
	  nomevias(2)='Fruit ingestion'
	  nomevias(3)='Meat ingestion'
	  nomevias(4)='Milk ingestion'
	  nomevias(5)='Water ingestion'
	  nomevias(6)='Vegetables ingestion'
	  nomevias(7)='Fish ingestion'
	  nomevias(8)='Bird ingestion'
	  nomevias(9)='Eggs ingestion'
	  nomevias(10)='Grains ingestion'
	  nomevias(11)='Particulate inhalation'
	  nomevias(12)='Steam inhalation'
	  nomevias(13)='Dermal water'
	  nomevias(14)='Dermal soil'
!
!
!
      DO KO=1,NLOCAL
!
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"1 to <2",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(1,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(1,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"2 to <3",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(2,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(2,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"3 to <6",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(3,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(3,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"6 to <11",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(4,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(4,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"11 to <16",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(5,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(5,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"16 to <18",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(6,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(6,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"18 to <21",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(7,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(7,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"21 to <65",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(8,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(8,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,">65",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_NC(9,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_NC(9,KP,1,KO)
	  ENDIF
	  ENDDO
	  IF(KO.NE.NLOCAL)THEN
	  WRITE(55,'("},")')!
	  ELSE
	  WRITE(55,'("}")')!
	  ENDIF
!
	  ENDDO	  ! FIM DO KO=1,NLOCAL
!
	  WRITE(55,'("],")')
!
!
!
!
	  WRITE(55,'(A1,"Contribution of each pathway to the carcinogenic risks",A1,": [")')aspas,aspas            
!
!
!
      DO KO=1,NLOCAL
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"1 to <2",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(1,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(1,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"2 to <3",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(2,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(2,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"3 to <6",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(3,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(3,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"6 to <11",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(4,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(4,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"11 to <16",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(5,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(5,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"16 to <18",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(6,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(6,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"18 to <21",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(7,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(7,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')	
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,"21 to <65",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(8,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(8,KP,1,KO)
	  ENDIF
	  ENDDO
	  WRITE(55,'("},")')
!
	  WRITE(55,'("{")')
      write(55,'(A1,"Local",A1,":",1x,I3,",")') aspas,aspas,KO
      write(55,'(A1,"Age groups",A1,":",1X,A1,">65",A1,",")') aspas,aspas,aspas,aspas
	  DO KP=1,NVIAS
	  IF(KP.NE.14)THEN
      write(55,'(A1,A22,A1,":",1X,ES10.3,",")') aspas,nomevias(KP),aspas,PORC_VIA_C(9,KP,1,KO)
	  ELSE
      write(55,'(A1,A22,A1,":",1X,ES10.3)') aspas,nomevias(KP),aspas,PORC_VIA_C(9,KP,1,KO)
	  ENDIF
	  ENDDO
	  IF(KO.NE.NLOCAL)THEN
	  WRITE(55,'("},")')!
	  ELSE
	  WRITE(55,'("}")')!
	  ENDIF
!
	  ENDDO	  ! FIM DO KO=1,NLOCAL
!
	  WRITE(55,'("]")')
!
	  WRITE(55,'("}")')
!
!								 
      CLOSE (55)
!
!
!
      ENDIF ! fim do IF da chave KEY_ANALYZE
!
!
!-------------------------------------------------------------------------------------------------------
!	                               CÁLCULO DO RISCO RADIOLOGICO
!-------------------------------------------------------------------------------------------------------
!
      IF(KEY_ANALYZE(2).EQV..TRUE.)THEN
!
!
      CALL RADIO(VARIASAO,NDURATION,NCHEM,NTIME,NLOCAL,NTYPECONC,CHEMICAL,RAD_POL,AMOLAR,VIDAMEIA,KEY_SD)
!
      ENDIF
!
!
      IF(KEY_ANALYZE(3).EQV..TRUE.)THEN
!
      WRITE(*,*)
	  WRITE(*,'("____________________________________________________________________________________")')   
	  WRITE(*,'("                            ECOLOGICAL RISK ASSESMENT ")')   
	  WRITE(*,'("____________________________________________________________________________________")')   
!
	  WRITE(99,'("____________________________________________________________________________________________________________________________________________")')
	  WRITE(99,'("                                                          ECOLOGICAL RISK ASSESMENT")')
	  WRITE(99,'("____________________________________________________________________________________________________________________________________________")')
      WRITE(99,*)
!
      CALL ECOLOGICAL(SCENAR,VARIASAO,NDURATION,NCHEM,NTIME,NLOCAL,NTYPECONC,CHEMICAL,SCENARIES)
!
      ENDIF
!
!
!
       WRITE(*,*)
	   WRITE(*,*)
	   PRINT*,'THE CODE FINISHED ALL CALCULATIONS'
       WRITE(*,*)
	   WRITE(*,*)
!
	STOP
	END
!
!
!
!*********************************************************************************************************************
!=====================================================================================================================
!                     ----------------------------------------------------------------------
!                               VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!                                     VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!                                            VVVVVVVVVVVVVVVVVVVV
!                                                  VVVVVVVV
!                                             
!                                             FIM DO MAIN PROGRAM
!
!
!---------------------------------------------------------------------------------------------------------------------
!
!
!
      SUBROUTINE READSCENARY(NIDADE,NVIAS,SCENAR,NCHEM,NLOCAL,NDURATION,KEY_SD,KEY_ANALYZE,W,CF,SD_CF,FI,SD_FI,IR,SD_IR,FA,SD_FA,FP,SD_FP,&
	  ET,SD_ET,SA,SD_SA,AF,SD_AF,EV,SD_EV,NRISKTYPE,NINFORMATION,EF_INI,SD_EF_INI,AT_INI,SD_AT_INI)
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL W,KEY_SD,KEY_ANALYZE
	  REAL*8 IR
      INTEGER SCENAR
 	  CHARACTER(LEN=10)  :: NRISKTYPE(NVIAS)
	  DIMENSION TIMESP(0:2000) ! ACOMODA NUMERO DE PERIODOS DE TEMPO CONSIDERADOS  ! 
	  DIMENSION CF(3),SD_CF(3),FI(9),SD_FI(9),Fa(4),SD_Fa(4),Fp(4),SD_Fp(4)
	  DIMENSION IR(NIDADE,20),SD_IR(NIDADE,20),ET(NIDADE,3),SD_ET(NIDADE,3),SA(NIDADE,2),SD_SA(NIDADE,2),EV(NIDADE,2),SD_EV(NIDADE,2),AF(NIDADE),SD_AF(NIDADE)
	  DIMENSION W(NVIAS),KEY_ANALYZE(3),EF_INI(NVIAS,NIDADE),SD_EF_INI(NVIAS,NIDADE),AT_INI(2,NVIAS,NIDADE),SD_AT_INI(2,NVIAS,NIDADE)
!
!
!	  SUBROUTINA usada para ler os dados de dados como as concentrações iniciais da água, solo ou ar durante vários
!     tempos, além de ler o cenário analisado, o tempo de exposição, a frequência de exposição da população estudada
!     e os caminhos possíveis que o contaminante poderia percorrer antes de chegar nos humanos (esses caminhos podem
!     ser escolhidos, os não desejados podem ser desligados e o cálculo relacionado a eles não será executado). 
!
      DO KH=1,NVIAS
	  W(KH)=.FALSE.
	  NRISKTYPE(KH)='Chronic'
	  ENDDO
!
      DO KO=1,3
      KEY_ANALYZE(KO)=.FALSE.
	  ENDDO
!
	  DO JJJ=1,NVIAS
	  DO KKK=1,NIDADE
	  EF_INI(JJJ,KKK)=0.0
	  SD_EF_INI(JJJ,KKK)=0.0
	  DO III=1,2
	  AT_INI(III,JJJ,KKK)=0.0
	  SD_AT_INI(III,JJJ,KKK)=0.0
	  ENDDO
	  ENDDO
	  ENDDO
!
      DO LK=1,NIDADE
!
	  AF(LK)=0.0
	  SD_AF(LK)=0.0
!
      DO K=1,3
	  ET(LK,K)=0.0
	  SD_ET(LK,K)=0.0
	  ENDDO
!
      DO I=1,2
	  SA(LK,I)=0.0
	  SD_SA(LK,I)=0.0
      EV(LK,I)=0.0
	  SD_EV(LK,I)=0.0
	  ENDDO
!
	  DO L=1,20
	  IR(LK,L)=0.0
	  SD_IR(LK,L)=0.0
	  ENDDO
!
      ENDDO
!    
      DO I=1,3
	  CF(I)=0.0
	  SD_CF(I)=0.0
	  ENDDO

	  DO K=1,9
	  FI(K)=0.0
	  SD_FI(K)=0.0
	  ENDDO

	  DO J=1,4
	  Fa(J)=0.0
	  SD_Fa(J)=0.0
	  Fp(J)=0.0
	  SD_Fp(J)=0.0
	  ENDDO
!
!
      OPEN(UNIT=1,STATUS='OLD',FILE='SCENARY.prn')
!
!	 
      READ(1,12)
	  READ(1,*) SCENAR,NCHEM,NDURATION,NLOCAL,NINFORMATION
!
      IF((SCENAR.GT.4).OR.(SCENAR.LT.1))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' SCENARY PARAMETER MUST BE Agricultural, Industrial, Residential OR In natura '')')
	  WRITE(*,*)
	  WRITE(*,'('' For more information, enable key information --- Help --- in the Scenary Sheet '')')
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
      IF((NINFORMATION.LT.0).OR.(NINFORMATION.GT.1))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Information PARAMETER MUST BE Help OR Disabled '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
!
      IF((NCHEM.LE.0).OR.(NDURATION.LE.0).OR.(NLOCAL.LE.0))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Num. of Chemicals, Num. of Times AND Num. of Sites PARAMETERS MUST GREATER OR EQUAL TO 1 '')')
	  WRITE(*,*)
	  WRITE(*,'('' For more information, enable key information --- Help --- in the Scenary Sheet '')')
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  READ(1,*)
	  READ(1,*)
	  READ(1,*)KEY_SD 
!      
      DO KOP=1,3
	  READ(1,*)
	  READ(1,*)
	  READ(1,*)KEY_ANALYZE(KOP)
	  ENDDO
!      
      IF((SCENAR.EQ.4).AND.(KEY_ANALYZE(1).EQV..TRUE.))THEN
	  WRITE(99,*)
	  WRITE(99,'('' WARNING!!!!   WARNING!!!!   WARNING!!!!   WARNING!!!!   WARNING!!!!   WARNING!!!!   '')')
	  WRITE(99,'('' Scenario 4 - In Natura - does not allow calculation of human health non-radiological risks. '')')
	  WRITE(99,*)
	  WRITE(99,'('' The calculation of human health risk will not be performed!!! '')')
	  WRITE(99,'('' For more information, enable key information --- Help --- in the Scenary Sheet '')')
	  WRITE(99,*)
      ENDIF  
!                    
      	     
!         
12      FORMAT (7/)
19      FORMAT (5/)
20      FORMAT (9/)
!
35     	FORMAT (2/)
36     	FORMAT (3/)
37     	FORMAT (7/)
38		FORMAT (4/)
39     	FORMAT (6/)
41     	FORMAT (1/)
!
!
	  IF (SCENAR.EQ.1) THEN	   ! ESSA PARTE SERÁ LIDA SOMENTE SE O CENÁRIO ESCOLHIDO FOR O AGRICULTURA
!    
	  READ(1,19) 
! 
!---------------------------------------------------------------------------------------------------------------    
	  READ(1,*) W(1)   ! PARAMETRO CUJO VALOR É 0 OU 1 QUE DEFINE O "WAY 1" (CAMINHO 1) DO CENÁRIO AGRICULTURA 
	  IF(W(1).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(1) 
	  READ(1,41)
	  READ(1,*) CF(1),SD_CF(1),FI(1),SD_FI(1) 	 ! FI CAMINHO 1
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,1),SD_IR(I,1),EF_INI(1,I),SD_EF_INI(1,I),AT_INI(1,1,I),SD_AT_INI(1,1,I)		! IR SOLO
	  ENDDO
!
     IF(CF(1).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF(FI(1).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0  ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR1=0
	 KSTOPEF1=0
	 KSTOPAT1=0
	 DO L=1,NIDADE
     IF((IR(L,1).EQ.0.0).AND.(KSTOPIR1.EQ.0))THEN
	 KSTOPIR1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(1,L).EQ.0.0).AND.(KSTOPEF1.EQ.0))THEN
	 KSTOPEF1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,1,L).EQ.0.0).AND.(KSTOPAT1.EQ.0))THEN
	 KSTOPAT1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
      IF((NRISKTYPE(1).NE.'Chronic').AND.(NRISKTYPE(1).NE.'Subchronic').AND.(NRISKTYPE(1).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 1 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(2)
	  IF(W(2).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(2) 
	  READ(1,41)
	  READ(1,*) FI(2),SD_FI(2)	 ! FI CAMINHO 2
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,2),SD_IR(I,2),EF_INI(2,I),SD_EF_INI(2,I),AT_INI(1,2,I),SD_AT_INI(1,2,I)		! IR SOLO
	  ENDDO
!
     IF(FI(2).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR2=0
	 KSTOPEF2=0
	 KSTOPAT2=0
	 DO L=1,NIDADE
     IF((IR(L,2).EQ.0.0).AND.(KSTOPIR2.EQ.0))THEN
	 KSTOPIR2=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(2,L).EQ.0.0).AND.(KSTOPEF2.EQ.0))THEN
	 KSTOPEF2=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,2,L).EQ.0.0).AND.(KSTOPAT2.EQ.0))THEN
	 KSTOPAT2=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF 
!
	 ENDDO
!
      IF((NRISKTYPE(2).NE.'Chronic').AND.(NRISKTYPE(2).NE.'Subchronic').AND.(NRISKTYPE(2).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 2 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(3)
	  IF(W(3).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(3) 
	  READ(1,41)
	  READ(1,*) FI(3),SD_FI(3)	 ! FI CAMINHO 3
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,3),SD_IR(I,3),EF_INI(3,I),SD_EF_INI(3,I),AT_INI(1,3,I),SD_AT_INI(1,3,I)		! IR SOLO
	  ENDDO
	  READ(1,36)
	  READ(1,*) FA(1),SD_FA(1),FP(1),SD_FP(1)
	  READ(1,35)
	  READ(1,*)IR(1,4),SD_IR(1,4),IR(1,5),SD_IR(1,5)		! IR FOOD BTF
	  READ(1,35)
	  READ(1,*)IR(1,11),SD_IR(1,11)				
!
     IF(FI(3).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR3=0
	 KSTOPEF3=0
	 KSTOPAT3=0
	 DO L=1,NIDADE
     IF(((IR(L,3).EQ.0.0).OR.(IR(1,4).EQ.0.0).OR.(IR(1,5).EQ.0.0).OR.(IR(1,11).EQ.0.0)).AND.(KSTOPIR3.EQ.0))THEN
	 KSTOPIR3=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(3,L).EQ.0.0).AND.(KSTOPEF3.EQ.0))THEN
	 KSTOPEF3=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,3,L).EQ.0.0).AND.(KSTOPAT3.EQ.0))THEN
	 KSTOPAT3=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	 IF(Fa(1).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fa VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
	 IF(Fp(1).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fp VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
      IF((NRISKTYPE(3).NE.'Chronic').AND.(NRISKTYPE(3).NE.'Subchronic').AND.(NRISKTYPE(3).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 3 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(4)
	  IF(W(4).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(4)
	  READ(1,41)
	  READ(1,*) FI(4),SD_FI(4)	 ! FI CAMINHO 4
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,6),SD_IR(I,6),EF_INI(4,I),SD_EF_INI(4,I),AT_INI(1,4,I),SD_AT_INI(1,4,I)					
	  ENDDO
	  READ(1,36)
	  READ(1,*) FA(2),SD_FA(2),FP(2),SD_FP(2)
	  READ(1,35)
	  READ(1,*)IR(1,7),SD_IR(1,7),IR(1,8),SD_IR(1,8)
	  READ(1,35)
	  READ(1,*)IR(1,12),SD_IR(1,12)		
!
     IF(FI(4).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR4=0
	 KSTOPEF4=0
	 KSTOPAT4=0
	 DO L=1,NIDADE
     IF(((IR(L,6).EQ.0.0).OR.(IR(1,7).EQ.0.0).OR.(IR(1,8).EQ.0.0).OR.(IR(1,12).EQ.0.0)).AND.(KSTOPIR4.EQ.0))THEN
	 KSTOPIR4=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(4,L).EQ.0.0).AND.(KSTOPEF4.EQ.0))THEN
	 KSTOPEF4=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,4,L).EQ.0.0).AND.(KSTOPAT4.EQ.0))THEN
	 KSTOPAT4=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	 IF(Fa(2).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fa VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
	 IF(Fp(2).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fp VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
      IF((NRISKTYPE(4).NE.'Chronic').AND.(NRISKTYPE(4).NE.'Subchronic').AND.(NRISKTYPE(4).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 4 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(5)
	  IF(W(5).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(5)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,9),SD_IR(I,9),EF_INI(5,I),SD_EF_INI(5,I),AT_INI(1,5,I),SD_AT_INI(1,5,I)				
	  ENDDO
!
     KSTOPIR5=0
	 KSTOPEF5=0
	 KSTOPAT5=0
	 DO L=1,NIDADE
     IF((IR(L,9).EQ.0.0).AND.(KSTOPIR5.EQ.0))THEN
	 KSTOPIR5=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 5 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(5,L).EQ.0.0).AND.(KSTOPEF5.EQ.0))THEN
	 KSTOPEF5=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 5 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,5,L).EQ.0.0).AND.(KSTOPAT5.EQ.0))THEN
	 KSTOPAT5=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 5 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
      IF((NRISKTYPE(5).NE.'Chronic').AND.(NRISKTYPE(5).NE.'Subchronic').AND.(NRISKTYPE(5).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 5 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(6)
	  IF(W(6).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(6)
	  READ(1,41)
	  READ(1,*) FI(5),SD_FI(5)	 ! FI CAMINHO 8
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,10),SD_IR(I,10),EF_INI(6,I),SD_EF_INI(6,I),AT_INI(1,6,I),SD_AT_INI(1,6,I)					
	  ENDDO
!
     IF(FI(5).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR6=0
	 KSTOPEF6=0
     KSTOPAT6=0
!
	 DO L=1,NIDADE
     IF((IR(L,10).EQ.0.0).AND.(KSTOPIR6.EQ.0))THEN
	 KSTOPIR6=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(6,L).EQ.0.0).AND.(KSTOPEF6.EQ.0))THEN
	 KSTOPEF6=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,6,L).EQ.0.0).AND.(KSTOPAT6.EQ.0))THEN
	 KSTOPAT6=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
      IF((NRISKTYPE(6).NE.'Chronic').AND.(NRISKTYPE(6).NE.'Subchronic').AND.(NRISKTYPE(6).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 6 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(7)
	  IF(W(7).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(7)
	  READ(1,41)
	  READ(1,*) FI(6),SD_FI(6)	 ! FI CAMINHO 13
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,13),SD_IR(I,13),EF_INI(7,I),SD_EF_INI(7,I),AT_INI(1,7,I),SD_AT_INI(1,7,I)				
	  ENDDO
!
     IF(FI(6).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 7 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR7=0
     KSTOPEF7=0
     KSTOPAT7=0
!
	 DO L=1,NIDADE
     IF((IR(L,13).EQ.0.0).AND.(KSTOPIR7.EQ.0))THEN
	 KSTOPIR7=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 7 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(7,L).EQ.0.0).AND.(KSTOPEF7.EQ.0))THEN
	 KSTOPEF7=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 7 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,7,L).EQ.0.0).AND.(KSTOPAT7.EQ.0))THEN
	 KSTOPAT7=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 7 ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
      IF((NRISKTYPE(7).NE.'Chronic').AND.(NRISKTYPE(7).NE.'Subchronic').AND.(NRISKTYPE(7).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 7 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(8)
	  IF(W(8).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(8)
	  READ(1,41)
	  READ(1,*) FI(7),SD_FI(7)	 ! FI CAMINHO 14
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,14),SD_IR(I,14),EF_INI(8,I),SD_EF_INI(8,I),AT_INI(1,8,I),SD_AT_INI(1,8,I)			
	  ENDDO
	  READ(1,36)
	  READ(1,*) FA(3),SD_FA(3),FP(3),SD_FP(3)
	  READ(1,35)
	  READ(1,*)IR(1,15),SD_IR(1,15),IR(1,16),SD_IR(1,16)
!
     IF(FI(7).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 8 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR8=0
	 KSTOPEF8=0
	 KSTOPAT8=0
	 DO L=1,NIDADE
     IF(((IR(L,14).EQ.0.0).OR.(IR(1,15).EQ.0.0).OR.(IR(1,16).EQ.0.0)).AND.(KSTOPIR8.EQ.0))THEN
	 KSTOPIR8=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 8 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(8,L).EQ.0.0).AND.(KSTOPEF8.EQ.0))THEN
	 KSTOPEF8=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 8 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,8,L).EQ.0.0).AND.(KSTOPAT8.EQ.0))THEN
	 KSTOPAT8=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 8 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	 IF(Fa(3).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fa VALUE IS EQUAL TO 0.0   ---> WAY 8 ")')
	 WRITE(99,*)
	 ENDIF
!
	 IF(Fp(3).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fp VALUE IS EQUAL TO 0.0   ---> WAY 8 ")')
	 WRITE(99,*)
	 ENDIF
!
      IF((NRISKTYPE(8).NE.'Chronic').AND.(NRISKTYPE(8).NE.'Subchronic').AND.(NRISKTYPE(8).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 8 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(9)
	  IF(W(9).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(9)
	  READ(1,41)
	  READ(1,*) FI(8),SD_FI(8)	 ! FI CAMINHO 15
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,17),SD_IR(I,17),EF_INI(9,I),SD_EF_INI(9,I),AT_INI(1,9,I),SD_AT_INI(1,9,I)				
	  ENDDO
	  READ(1,36)
	  READ(1,*) FA(4),SD_FA(4),FP(4),SD_FP(4)
	  READ(1,35)
	  READ(1,*)IR(1,18),SD_IR(1,18),IR(1,19),SD_IR(1,19)		! IR WATER BTF WAY 15
!
     IF(FI(8).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 9 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR9=0
	 KSTOPEF9=0
     KSTOPAT9=0
!
	 DO L=1,NIDADE
     IF(((IR(L,17).EQ.0.0).OR.(IR(1,18).EQ.0.0).OR.(IR(1,19).EQ.0.0)).AND.(KSTOPIR9.EQ.0))THEN
	 KSTOPIR9=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 9 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(9,L).EQ.0.0).AND.(KSTOPEF9.EQ.0))THEN
	 KSTOPEF9=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 9 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,9,L).EQ.0.0).AND.(KSTOPAT9.EQ.0))THEN
	 KSTOPAT9=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 9 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	 IF(Fa(4).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fa VALUE IS EQUAL TO 0.0   ---> WAY 9 ")')
	 WRITE(99,*)
	 ENDIF
!
	 IF(Fp(4).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fp VALUE IS EQUAL TO 0.0   ---> WAY 9 ")')
	 WRITE(99,*)
	 ENDIF
!
      IF((NRISKTYPE(9).NE.'Chronic').AND.(NRISKTYPE(9).NE.'Subchronic').AND.(NRISKTYPE(9).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 9 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(10)
	  IF(W(10).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(10)
	  READ(1,41)
	  READ(1,*) FI(9),SD_FI(9)	 ! FI CAMINHO 16
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,20),SD_IR(I,20),EF_INI(10,I),SD_EF_INI(10,I),AT_INI(1,10,I),SD_AT_INI(1,10,I)				
	  ENDDO
!
     IF(FI(9).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 10 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR10=0
	 KSTOPEF10=0
	 KSTOPAT10=0
	 DO L=1,NIDADE
     IF((IR(L,20).EQ.0.0).AND.(KSTOPIR10.EQ.0))THEN
	 KSTOPIR10=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 10 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(10,L).EQ.0.0).AND.(KSTOPEF10.EQ.0))THEN
	 KSTOPEF10=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 10 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,10,L).EQ.0.0).AND.(KSTOPAT10.EQ.0))THEN
	 KSTOPAT10=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 10 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
      IF((NRISKTYPE(10).NE.'Chronic').AND.(NRISKTYPE(10).NE.'Subchronic').AND.(NRISKTYPE(10).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 10 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,35)
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*) W(11)
	  IF(W(11).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(11)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,1),SD_ET(I,1),EF_INI(11,I),SD_EF_INI(11,I),AT_INI(1,11,I),SD_AT_INI(1,11,I)		
	  ENDDO
!
     KSTOPET1=0
	 KSTOPEF11=0
	 KSTOPAT11=0
	 DO L=1,NIDADE
     IF((ET(L,1).EQ.0.0).AND.(KSTOPET1.EQ.0))THEN
	 KSTOPET1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 11  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(11,L).EQ.0.0).AND.(KSTOPEF11.EQ.0))THEN
	 KSTOPEF11=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 11 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,11,L).EQ.0.0).AND.(KSTOPAT11.EQ.0))THEN
	 KSTOPAT11=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 11 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
      IF((NRISKTYPE(11).NE.'Chronic').AND.(NRISKTYPE(11).NE.'Subchronic').AND.(NRISKTYPE(11).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 11 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*)
	  READ(1,*) W(12)
	  IF(W(12).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(12)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,2),SD_ET(I,2),EF_INI(12,I),SD_EF_INI(12,I),AT_INI(1,12,I),SD_AT_INI(1,12,I)			
	  ENDDO
!
     KSTOPET2=0
	 KSTOPEF12=0
	 KSTOPAT12=0
	 DO L=1,NIDADE
     IF((ET(L,2).EQ.0.0).AND.(KSTOPET2.EQ.0))THEN
	 KSTOPET2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 12  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(12,L).EQ.0.0).AND.(KSTOPEF12.EQ.0))THEN
	 KSTOPEF12=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 12 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,12,L).EQ.0.0).AND.(KSTOPAT12.EQ.0))THEN
	 KSTOPAT12=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 12 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
      IF((NRISKTYPE(12).NE.'Chronic').AND.(NRISKTYPE(12).NE.'Subchronic').AND.(NRISKTYPE(12).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 12 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,35)
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*) W(13)
	  IF(W(13).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(13)
	  READ(1,41)
	  READ(1,*)	CF(2),SD_CF(2)
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,3),SD_ET(I,3),EV(I,1),SD_EV(I,1),SA(I,1),SD_SA(I,1),EF_INI(13,I),SD_EF_INI(13,I),AT_INI(1,13,I),SD_AT_INI(1,13,I)				
	  ENDDO
!
     IF(CF(2).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 13 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPET3=0
	 KSTOPEF13=0
	 KSTOPAT13=0
	 DO L=1,NIDADE
     IF((ET(L,3).EQ.0.0).AND.(KSTOPET3.EQ.0))THEN
	 KSTOPET3=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 13  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(13,L).EQ.0.0).AND.(KSTOPEF13.EQ.0))THEN
	 KSTOPEF13=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 13 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,13,L).EQ.0.0).AND.(KSTOPAT13.EQ.0))THEN
	 KSTOPAT13=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 13 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
     KSTOPEV1=0
	 DO L=1,NIDADE
     IF((EV(L,1).EQ.0.0).AND.(KSTOPEV1.EQ.0))THEN
	 KSTOPEV1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! EV VALUE IS EQUAL TO 0.0   ---> WAY 13  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
     KSTOPSA1=0
	 DO L=1,NIDADE
     IF((SA(L,1).EQ.0.0).AND.(KSTOPSA1.EQ.0))THEN
	 KSTOPSA1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! SA VALUE IS EQUAL TO 0.0   ---> WAY 13  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
      IF((NRISKTYPE(13).NE.'Chronic').AND.(NRISKTYPE(13).NE.'Subchronic').AND.(NRISKTYPE(13).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 13 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*)
	  READ(1,*) W(14)
	  IF(W(14).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(14)
	  READ(1,41)
	  READ(1,*)	CF(3),SD_CF(3)
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)EV(I,2),SD_EV(I,2),AF(I),SD_AF(I),SA(I,2),SD_SA(I,2),EF_INI(14,I),SD_EF_INI(14,I),AT_INI(1,14,I),SD_AT_INI(1,14,I)				
	  ENDDO	
!
     IF(CF(3).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPAF=0
	 KSTOPEF14=0
	 KSTOPAT14=0
	 DO L=1,NIDADE
     IF((AF(L).EQ.0.0).AND.(KSTOPAF.EQ.0))THEN
	 KSTOPAF=1
	 WRITE(99, '(" WARNING!!!!!!!!!! AF VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(14,L).EQ.0.0).AND.(KSTOPEF14.EQ.0))THEN
	 KSTOPEF14=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,14,L).EQ.0.0).AND.(KSTOPAT14.EQ.0))THEN
	 KSTOPAT14=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
     KSTOPEV2=0
	 DO L=1,NIDADE
     IF((EV(L,2).EQ.0.0).AND.(KSTOPEV2.EQ.0))THEN
	 KSTOPEV2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! EV VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
     KSTOPSA2=0
	 DO L=1,NIDADE
     IF((SA(L,2).EQ.0.0).AND.(KSTOPSA2.EQ.0))THEN
	 KSTOPSA2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! SA VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO	
!
      IF((NRISKTYPE(14).NE.'Chronic').AND.(NRISKTYPE(14).NE.'Subchronic').AND.(NRISKTYPE(14).NE.'Acute'))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' Pathway 14 exposure duration must be Chronic, Subchronic or Acute   '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
	  ENDIF
!
!---------------------------------------------------------------------------------------------------------------         
!
!      PAUSE '1_5'
!
10      FORMAT (7x,L9)
!

	 ELSEIF (SCENAR.EQ.2) THEN	   ! ESSA PARTE SERÁ LIDA SOMENTE SE O CENÁRIO ESCOLHIDO FOR O industrial
! 
	  READ(1,12) 
!
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*) W(1)   ! PARAMETRO CUJO VALOR É 0 OU 1 QUE DEFINE O "WAY 1" (CAMINHO 1) DO CENÁRIO AGRICULTURA 
	  IF(W(1).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(1) 
	  READ(1,41)
	  READ(1,*) CF(1),SD_CF(1),FI(1),SD_FI(1) 	 ! FI CAMINHO 1
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,1),SD_IR(I,1),EF_INI(1,I),SD_EF_INI(1,I),AT_INI(1,1,I),SD_AT_INI(1,1,I)		! IR SOLO
	  ENDDO
!
     IF(CF(1).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF(FI(1).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0  ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR1=0
	 KSTOPEF1=0
	 KSTOPAT1=0
	 DO L=1,NIDADE
     IF((IR(L,1).EQ.0.0).AND.(KSTOPIR1.EQ.0))THEN
	 KSTOPIR1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(1,L).EQ.0.0).AND.(KSTOPEF1.EQ.0))THEN
	 KSTOPEF1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,1,L).EQ.0.0).AND.(KSTOPAT1.EQ.0))THEN
	 KSTOPAT1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF

!---------------------------------------------------------------------------------------------------------------         
	  READ(1,35)
!---------------------------------------------------------------------------------------------------------------        
	  READ(1,*) W(11)
	  IF(W(11).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(11)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,1),SD_ET(I,1),EF_INI(11,I),SD_EF_INI(11,I),AT_INI(1,11,I),SD_AT_INI(1,11,I)		
	  ENDDO
!
     KSTOPET1=0
	 KSTOPEF11=0
	 KSTOPAT11=0
	 DO L=1,NIDADE
     IF((ET(L,1).EQ.0.0).AND.(KSTOPET1.EQ.0))THEN
	 KSTOPET1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 11  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(11,L).EQ.0.0).AND.(KSTOPEF11.EQ.0))THEN
	 KSTOPEF11=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 11 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,11,L).EQ.0.0).AND.(KSTOPAT11.EQ.0))THEN
	 KSTOPAT11=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 11 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*)
	  READ(1,*) W(12)
	  IF(W(12).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(12)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,2),SD_ET(I,2),EF_INI(12,I),SD_EF_INI(12,I),AT_INI(1,12,I),SD_AT_INI(1,12,I)			
	  ENDDO
!
     KSTOPET2=0
	 KSTOPEF12=0
	 KSTOPAT12=0
	 DO L=1,NIDADE
     IF((ET(L,2).EQ.0.0).AND.(KSTOPET2.EQ.0))THEN
	 KSTOPET2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 12  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(12,L).EQ.0.0).AND.(KSTOPEF12.EQ.0))THEN
	 KSTOPEF12=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 12 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,12,L).EQ.0.0).AND.(KSTOPAT12.EQ.0))THEN
	 KSTOPAT12=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 12 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,35)
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*) W(14)
	  IF(W(14).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(14)
	  READ(1,41)
	  READ(1,*)	CF(3),SD_CF(3)
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)EV(I,2),SD_EV(I,2),AF(I),SD_AF(I),SA(I,2),SD_SA(I,2),EF_INI(14,I),SD_EF_INI(14,I),AT_INI(1,14,I),SD_AT_INI(1,14,I)				
	  ENDDO	
!
     IF(CF(3).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPAF=0
	 KSTOPEF14=0
	 KSTOPAT14=0
	 DO L=1,NIDADE
     IF((AF(L).EQ.0.0).AND.(KSTOPAF.EQ.0))THEN
	 KSTOPAF=1
	 WRITE(99, '(" WARNING!!!!!!!!!! AF VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(14,L).EQ.0.0).AND.(KSTOPEF14.EQ.0))THEN
	 KSTOPEF14=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,14,L).EQ.0.0).AND.(KSTOPAT14.EQ.0))THEN
	 KSTOPAT14=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
     KSTOPEV2=0
	 DO L=1,NIDADE
     IF((EV(L,2).EQ.0.0).AND.(KSTOPEV2.EQ.0))THEN
	 KSTOPEV2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! EV VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
     KSTOPSA2=0
	 DO L=1,NIDADE
     IF((SA(L,2).EQ.0.0).AND.(KSTOPSA2.EQ.0))THEN
	 KSTOPSA2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! SA VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO	
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
!
!
      ELSEIF(SCENAR.EQ.3) THEN
!
!13    FORMAT (444/,7x,L9)
	  READ(1,20) 
! 
!---------------------------------------------------------------------------------------------------------------    
	  READ(1,*) W(1)   ! PARAMETRO CUJO VALOR É 0 OU 1 QUE DEFINE O "WAY 1" (CAMINHO 1) DO CENÁRIO AGRICULTURA 
	  IF(W(1).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(1) 
	  READ(1,41)
	  READ(1,*) CF(1),SD_CF(1),FI(1),SD_FI(1) 	 ! FI CAMINHO 1
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,1),SD_IR(I,1),EF_INI(1,I),SD_EF_INI(1,I),AT_INI(1,1,I),SD_AT_INI(1,1,I)		! IR SOLO
	  ENDDO
!
     IF(CF(1).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF(FI(1).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0  ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR1=0
	 KSTOPEF1=0
	 KSTOPAT1=0
	 DO L=1,NIDADE
     IF((IR(L,1).EQ.0.0).AND.(KSTOPIR1.EQ.0))THEN
	 KSTOPIR1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(1,L).EQ.0.0).AND.(KSTOPEF1.EQ.0))THEN
	 KSTOPEF1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,1,L).EQ.0.0).AND.(KSTOPAT1.EQ.0))THEN
	 KSTOPAT1=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 1 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(2)
	  IF(W(2).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(2) 
	  READ(1,41)
	  READ(1,*) FI(2),SD_FI(2)	 ! FI CAMINHO 2
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,2),SD_IR(I,2),EF_INI(2,I),SD_EF_INI(2,I),AT_INI(1,2,I),SD_AT_INI(1,2,I)		! IR SOLO
	  ENDDO
!
     IF(FI(2).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR2=0
	 KSTOPEF2=0
	 KSTOPAT2=0
	 DO L=1,NIDADE
     IF((IR(L,2).EQ.0.0).AND.(KSTOPIR2.EQ.0))THEN
	 KSTOPIR2=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(2,L).EQ.0.0).AND.(KSTOPEF2.EQ.0))THEN
	 KSTOPEF2=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,2,L).EQ.0.0).AND.(KSTOPAT2.EQ.0))THEN
	 KSTOPAT2=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 2 ")')
	 WRITE(99,*)
	 ENDIF 
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(3)
	  IF(W(3).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(3) 
	  READ(1,41)
	  READ(1,*) FI(3),SD_FI(3)	 ! FI CAMINHO 3
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,3),SD_IR(I,3),EF_INI(3,I),SD_EF_INI(3,I),AT_INI(1,3,I),SD_AT_INI(1,3,I)		! IR SOLO
	  ENDDO
	  READ(1,36)
	  READ(1,*) FA(1),SD_FA(1),FP(1),SD_FP(1)
	  READ(1,35)
	  READ(1,*)IR(1,4),SD_IR(1,4),IR(1,5),SD_IR(1,5)		! IR FOOD BTF
	  READ(1,35)
	  READ(1,*)IR(1,11),SD_IR(1,11)				
!
     IF(FI(3).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR3=0
	 KSTOPEF3=0
	 KSTOPAT3=0
	 DO L=1,NIDADE
     IF(((IR(L,3).EQ.0.0).OR.(IR(1,4).EQ.0.0).OR.(IR(1,5).EQ.0.0).OR.(IR(1,11).EQ.0.0)).AND.(KSTOPIR3.EQ.0))THEN
	 KSTOPIR3=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(3,L).EQ.0.0).AND.(KSTOPEF3.EQ.0))THEN
	 KSTOPEF3=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,3,L).EQ.0.0).AND.(KSTOPAT3.EQ.0))THEN
	 KSTOPAT3=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	 IF(Fa(1).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fa VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
	 IF(Fp(1).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fp VALUE IS EQUAL TO 0.0   ---> WAY 3 ")')
	 WRITE(99,*)
	 ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(4)
	  IF(W(4).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(4)
	  READ(1,41)
	  READ(1,*) FI(4),SD_FI(4)	 ! FI CAMINHO 4
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,6),SD_IR(I,6),EF_INI(4,I),SD_EF_INI(4,I),AT_INI(1,4,I),SD_AT_INI(1,4,I)					
	  ENDDO
	  READ(1,36)
	  READ(1,*) FA(2),SD_FA(2),FP(2),SD_FP(2)
	  READ(1,35)
	  READ(1,*)IR(1,7),SD_IR(1,7),IR(1,8),SD_IR(1,8)
	  READ(1,35)
	  READ(1,*)IR(1,12),SD_IR(1,12)		
!
     IF(FI(4).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR4=0
	 KSTOPEF4=0
	 KSTOPAT4=0
	 DO L=1,NIDADE
     IF(((IR(L,6).EQ.0.0).OR.(IR(1,7).EQ.0.0).OR.(IR(1,8).EQ.0.0).OR.(IR(1,12).EQ.0.0)).AND.(KSTOPIR4.EQ.0))THEN
	 KSTOPIR4=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(4,L).EQ.0.0).AND.(KSTOPEF4.EQ.0))THEN
	 KSTOPEF4=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,4,L).EQ.0.0).AND.(KSTOPAT4.EQ.0))THEN
	 KSTOPAT4=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	 IF(Fa(2).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fa VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
	 IF(Fp(2).EQ.0.0)THEN
	 WRITE(99, '(" WARNING!!!!!!!!!! Fp VALUE IS EQUAL TO 0.0   ---> WAY 4 ")')
	 WRITE(99,*)
	 ENDIF
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(5)
	  IF(W(5).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(5)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,9),SD_IR(I,9),EF_INI(5,I),SD_EF_INI(5,I),AT_INI(1,5,I),SD_AT_INI(1,5,I)				
	  ENDDO
!
     KSTOPIR5=0
	 KSTOPEF5=0
	 KSTOPAT5=0
	 DO L=1,NIDADE
     IF((IR(L,9).EQ.0.0).AND.(KSTOPIR5.EQ.0))THEN
	 KSTOPIR5=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 5 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(5,L).EQ.0.0).AND.(KSTOPEF5.EQ.0))THEN
	 KSTOPEF5=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 5 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,5,L).EQ.0.0).AND.(KSTOPAT5.EQ.0))THEN
	 KSTOPAT5=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 5 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,*)
	  READ(1,*) W(6)
	  IF(W(6).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(6)
	  READ(1,41)
	  READ(1,*) FI(5),SD_FI(5)	 ! FI CAMINHO 8
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)IR(I,10),SD_IR(I,10),EF_INI(6,I),SD_EF_INI(6,I),AT_INI(1,6,I),SD_AT_INI(1,6,I)					
	  ENDDO
!
     IF(FI(5).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! FI VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPIR6=0
	 KSTOPEF6=0
     KSTOPAT6=0
!
	 DO L=1,NIDADE
     IF((IR(L,10).EQ.0.0).AND.(KSTOPIR6.EQ.0))THEN
	 KSTOPIR6=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME IR VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(6,L).EQ.0.0).AND.(KSTOPEF6.EQ.0))THEN
	 KSTOPEF6=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,6,L).EQ.0.0).AND.(KSTOPAT6.EQ.0))THEN
	 KSTOPAT6=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 6 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------
	  READ(1,35)
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*) W(11)
	  IF(W(11).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(11)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,1),SD_ET(I,1),EF_INI(11,I),SD_EF_INI(11,I),AT_INI(1,11,I),SD_AT_INI(1,11,I)		
	  ENDDO
!
     KSTOPET1=0
	 KSTOPEF11=0
	 KSTOPAT11=0
	 DO L=1,NIDADE
     IF((ET(L,1).EQ.0.0).AND.(KSTOPET1.EQ.0))THEN
	 KSTOPET1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 11  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(11,L).EQ.0.0).AND.(KSTOPEF11.EQ.0))THEN
	 KSTOPEF11=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 11 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,11,L).EQ.0.0).AND.(KSTOPAT11.EQ.0))THEN
	 KSTOPAT11=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 11 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*)
	  READ(1,*) W(12)
	  IF(W(12).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(12)
	  READ(1,*)
	  READ(1,*)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,2),SD_ET(I,2),EF_INI(12,I),SD_EF_INI(12,I),AT_INI(1,12,I),SD_AT_INI(1,12,I)			
	  ENDDO
!
     KSTOPET2=0
	 KSTOPEF12=0
	 KSTOPAT12=0
	 DO L=1,NIDADE
     IF((ET(L,2).EQ.0.0).AND.(KSTOPET2.EQ.0))THEN
	 KSTOPET2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 12  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(12,L).EQ.0.0).AND.(KSTOPEF12.EQ.0))THEN
	 KSTOPEF12=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 12 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,12,L).EQ.0.0).AND.(KSTOPAT12.EQ.0))THEN
	 KSTOPAT12=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 12 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,35)
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*) W(13)
	  IF(W(13).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(13)
	  READ(1,41)
	  READ(1,*)	CF(2),SD_CF(2)
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)ET(I,3),SD_ET(I,3),EV(I,1),SD_EV(I,1),SA(I,1),SD_SA(I,1),EF_INI(13,I),SD_EF_INI(13,I),AT_INI(1,13,I),SD_AT_INI(1,13,I)				
	  ENDDO
!
     IF(CF(2).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 13 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPET3=0
	 KSTOPEF13=0
	 KSTOPAT13=0
	 DO L=1,NIDADE
     IF((ET(L,3).EQ.0.0).AND.(KSTOPET3.EQ.0))THEN
	 KSTOPET3=1
	 WRITE(99, '(" WARNING!!!!!!!!!! ET VALUE IS EQUAL TO 0.0   ---> WAY 13  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(13,L).EQ.0.0).AND.(KSTOPEF13.EQ.0))THEN
	 KSTOPEF13=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 13 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,13,L).EQ.0.0).AND.(KSTOPAT13.EQ.0))THEN
	 KSTOPAT13=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 13 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
     KSTOPEV1=0
	 DO L=1,NIDADE
     IF((EV(L,1).EQ.0.0).AND.(KSTOPEV1.EQ.0))THEN
	 KSTOPEV1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! EV VALUE IS EQUAL TO 0.0   ---> WAY 13  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
     KSTOPSA1=0
	 DO L=1,NIDADE
     IF((SA(L,1).EQ.0.0).AND.(KSTOPSA1.EQ.0))THEN
	 KSTOPSA1=1
	 WRITE(99, '(" WARNING!!!!!!!!!! SA VALUE IS EQUAL TO 0.0   ---> WAY 13  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
	  READ(1,*)
	  READ(1,*) W(14)
	  IF(W(14).EQV..TRUE.)THEN
	  READ(1,*)
	  READ(1,*)NRISKTYPE(14)
	  READ(1,41)
	  READ(1,*)	CF(3),SD_CF(3)
	  READ(1,35)
	  DO I=1,NIDADE
	  READ(1,*)EV(I,2),SD_EV(I,2),AF(I),SD_AF(I),SA(I,2),SD_SA(I,2),EF_INI(14,I),SD_EF_INI(14,I),AT_INI(1,14,I),SD_AT_INI(1,14,I)				
	  ENDDO	
!
     IF(CF(3).EQ.0.0)THEN
	 WRITE(99,'(" WARNING!!!!!!!!!! CF VALUE IS EQUAL TO 0.0  ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
     KSTOPAF=0
	 KSTOPEF14=0
	 KSTOPAT14=0
	 DO L=1,NIDADE
     IF((AF(L).EQ.0.0).AND.(KSTOPAF.EQ.0))THEN
	 KSTOPAF=1
	 WRITE(99, '(" WARNING!!!!!!!!!! AF VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((EF_INI(14,L).EQ.0.0).AND.(KSTOPEF14.EQ.0))THEN
	 KSTOPEF14=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME EF VALUE IS EQUAL TO 0.0   ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
     IF((AT_INI(1,14,L).EQ.0.0).AND.(KSTOPAT14.EQ.0))THEN
	 KSTOPAT14=1
	 WRITE(99,'(" WARNING!!!!!!!!!! SOME AT - non-carc. VALUE IS EQUAL TO 0.0   ---> WAY 14 ")')
	 WRITE(99,*)
	 ENDIF
!
	 ENDDO
!
     KSTOPEV2=0
	 DO L=1,NIDADE
     IF((EV(L,2).EQ.0.0).AND.(KSTOPEV2.EQ.0))THEN
	 KSTOPEV2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! EV VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO
!
     KSTOPSA2=0
	 DO L=1,NIDADE
     IF((SA(L,2).EQ.0.0).AND.(KSTOPSA2.EQ.0))THEN
	 KSTOPSA2=1
	 WRITE(99, '(" WARNING!!!!!!!!!! SA VALUE IS EQUAL TO 0.0   ---> WAY 14  ")')
	 WRITE(99,*)
	 ENDIF
	 ENDDO	
!
	  ENDIF
!---------------------------------------------------------------------------------------------------------------         
!
      ENDIF		! TERMINA O IF DE "CASE - SCENARY"
!--------------------------------------------------------------------------------------------------------
!
	 TIMESP(0)=0.0
!
!	  
      DO JJ=1,NIDADE
      IR(JJ,4)=IR(1,4)
	  SD_IR(JJ,4)=SD_IR(1,4)
	  IR(JJ,5)=IR(1,5)
	  SD_IR(JJ,5)=SD_IR(1,5)
      IR(JJ,7)=IR(1,7)
	  SD_IR(JJ,7)=SD_IR(1,7)
      IR(JJ,8)=IR(1,8)
	  SD_IR(JJ,8)=SD_IR(1,8)
      IR(JJ,11)=IR(1,11)
	  SD_IR(JJ,11)=SD_IR(1,11)
      IR(JJ,12)=IR(1,12)
	  SD_IR(JJ,12)=SD_IR(1,12)
      IR(JJ,15)=IR(1,15)
	  SD_IR(JJ,15)=SD_IR(1,15)
      IR(JJ,16)=IR(1,16)
	  SD_IR(JJ,16)=SD_IR(1,16)
      IR(JJ,18)=IR(1,18)
	  SD_IR(JJ,18)=SD_IR(1,18)
      IR(JJ,19)=IR(1,19)
	  SD_IR(JJ,19)=SD_IR(1,19)
	  ENDDO
!
!
      CLOSE(1)
!
!
      RETURN
	  END
!
!
!**********************************************************************************************
!																																																																																										   
      SUBROUTINE CONCENTRATION(VARIASAO,NDURATION,NCHEM,NTYPECONC,CHEMICAL,CSOIL,CWATER, &
	  CPAR,CSTEAM,CFRUIT,CLEAVES,CBEEF,CMILK,CAVE,CEGG,CFISH,CGRAIN,CWATERDER,CWATEROTHER,CSEDIMENT,&
	  KEYCONC,NTIME,NLOCAL,TIMESP,SD_CSOIL,&
      SD_CWATER,SD_CPAR,SD_CSTEAM,SD_CFRUIT,SD_CLEAVES,SD_CBEEF,SD_CMILK,&
	  SD_CAVE,SD_CEGG,SD_CFISH,SD_CGRAIN,SD_CWATERDER,SD_CWATEROTHER,SD_CSEDIMENT)
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL KEYCONC
	CHARACTER(LEN=50)  :: CHEMICAL(500)
    CHARACTER(LEN=30)  :: CONC_TYPE(NTYPECONC)
	DIMENSION CSOIL(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL)
	DIMENSION CWATER(NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_CWATER(NCHEM,NTIME,NLOCAL)
	DIMENSION CWATERDER(NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_CWATERDER(NCHEM,NTIME,NLOCAL)
    DIMENSION CWATEROTHER(NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_CWATEROTHER(NCHEM,NTIME,NLOCAL)
	DIMENSION CPAR(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CPAR(NCHEM,NTIME,NLOCAL)
	DIMENSION CSTEAM(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CSTEAM(NCHEM,NTIME,NLOCAL)
	DIMENSION CFRUIT(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CFRUIT(NCHEM,NTIME,NLOCAL)
	DIMENSION CLEAVES(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CLEAVES(NCHEM,NTIME,NLOCAL)
	DIMENSION CBEEF(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CBEEF(NCHEM,NTIME,NLOCAL)
	DIMENSION CMILK(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CMILK(NCHEM,NTIME,NLOCAL)
	DIMENSION CAVE(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CAVE(NCHEM,NTIME,NLOCAL)
	DIMENSION CEGG(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CEGG(NCHEM,NTIME,NLOCAL)
	DIMENSION CFISH(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CFISH(NCHEM,NTIME,NLOCAL)
	DIMENSION CGRAIN(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CGRAIN(NCHEM,NTIME,NLOCAL)
	DIMENSION CSEDIMENT(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CSEDIMENT(NCHEM,NTIME,NLOCAL)
!
	DIMENSION CONC_MATRIX(NTYPECONC,NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CONC_MATRIX(NTYPECONC,NCHEM,NTIME,NLOCAL)
!
!
	DIMENSION KEYCONC(NCHEM,NTYPECONC)  ! NTYPECONC SÃO OS 15 TIPOS DE CONCENTRAÇOES ACIMA !!!!
	DIMENSION TIMESP(0:NTIME)
!
      OPEN(UNIT=3,STATUS='OLD',FILE='concentration.prn')
!
!
      DO i=1,NCHEM
	  DO j=1,NDURATION	
	  DO k=1,NLOCAL
!
      CSOIL(i,j,k)=0.0
	  SD_CSOIL(i,j,k)=0.0
	  CWATER(i,j,k)=0.0
	  SD_CWATER(i,j,k)=0.0
	  CWATERDER(i,j,k)=0.0
	  SD_CWATERDER(i,j,k)=0.0
	  CWATEROTHER(i,j,k)=0.0
	  SD_CWATEROTHER(i,j,k)=0.0
	  CPAR(i,j,k)=0.0
	  SD_CPAR(i,j,k)=0.0
	  CSTEAM(i,j,k)=0.0
	  SD_CSTEAM(i,j,k)=0.0
	  CFRUIT(i,j,k)=0.0
	  SD_CFRUIT(i,j,k)=0.0
	  CLEAVES(i,j,k)=0.0
	  SD_CLEAVES(i,j,k)=0.0
	  CBEEF(i,j,k)=0.0
	  SD_CBEEF(i,j,k)=0.0
	  CMILK(i,j,k)=0.0
	  SD_CMILK(i,j,k)=0.0
	  CAVE(i,j,k)=0.0
	  SD_CAVE(i,j,k)=0.0
	  CEGG(i,j,k)=0.0
	  SD_CEGG(i,j,k)=0.0
	  CFISH(i,j,k)=0.0
	  SD_CFISH(i,j,k)=0.0
	  CGRAIN(i,j,k)=0.0
	  SD_CGRAIN(i,j,k)=0.0
	  CSEDIMENT(i,j,k)=0.0
	  SD_CSEDIMENT(i,j,k)=0.0
!
	  DO LE=1,NTYPECONC
!
      CONC_MATRIX(LE,i,j,k)=0.0
	  SD_CONC_MATRIX(LE,i,j,k)=0.0
!
      KEYCONC(i,LE) = .FALSE.
!
      ENDDO
	  ENDDO
	  ENDDO
	  ENDDO
!
      READ(3,63)  ! READ USADOS PARA COLOCAR COMENTARIOS E FAZER ABERTURA DO FILE
!
63     	FORMAT (18/)
!
!      pause '2'
!
      DO i=1,NCHEM
!
!
				
      READ(3,*) CHEMICAL(i)
!
!
!
!      IF((N_MATRIX.GT.13).OR.(N_MATRIX.LT.1))THEN
!	  WRITE(*,*)
!	  WRITE(*,*)
!	  WRITE(*,'('' NUMBER OF MATRIX UTILIZED GREATER THAN 1 AND LESS THAN 13 '')')
!	  WRITE(*,*)
!	  WRITE(*,'('' For more information, enable key --- Information = 1 --- in the Scenary input file '')')
!	  WRITE(*,*)
!	  WRITE(*,'('' THE CODE WILL STOP '')')
!	  WRITE(*,*)
!	  WRITE(*,*)
!	  WRITE(*,*)
!	  stop
!     ENDIF
!
!
	  DO LO=1,NTYPECONC
!
      READ(3,*)
      READ(3,*)	CONC_TYPE(LO)	
!
!
!
!							KEYCONC(i,1)= KC1SOIL
!							KEYCONC(i,2)= KC2DRINKING_WATER
!							KEYCONC(i,3)= KC3PARTICULATE
!							KEYCONC(i,4)= KC4FRUIT
!							KEYCONC(i,5)= KC5LEAVES
!							KEYCONC(i,6)= KC6MEAT
!							KEYCONC(i,7)= KC7MILK
!							KEYCONC(i,8)= KC8AVE
!							KEYCONC(i,9)= KC9EGG
!							KEYCONC(i,10)= KC10FISH
!							KEYCONC(i,11)= KC11GRAIN
!                           KEYCONC(i,12)= KC12STEAM
!                           KEYCONC(i,13)= KC13BATH_WATER
!                           KEYCONC(i,14)= KC14OTHER_WATERS
!                           KEYCONC(i,15)= KC15SEDIMENT
!
	IF((CONC_TYPE(LO).EQ.'GRAIN').OR.(CONC_TYPE(LO).EQ.'FISH').OR.(CONC_TYPE(LO).EQ.'EGG').OR.(CONC_TYPE(LO).EQ.'BIRD').OR.&
	(CONC_TYPE(LO).EQ.'MILK').OR.(CONC_TYPE(LO).EQ.'BEEF').OR.(CONC_TYPE(LO).EQ.'LEAVES').OR.(CONC_TYPE(LO).EQ.'FRUIT').OR.&
	(CONC_TYPE(LO).EQ.'STEAM').OR.(CONC_TYPE(LO).EQ.'PARTICULATE').OR.(CONC_TYPE(LO).EQ.'BATH_WATER').OR.(CONC_TYPE(LO).EQ.'DRINKING_WATER').OR.&
	(CONC_TYPE(LO).EQ.'OTHER_WATERS').OR.(CONC_TYPE(LO).EQ.'SOIL').OR.(CONC_TYPE(LO).EQ.'SEDIMENTS'))THEN
!
!
	DO j=1,NDURATION	  
	READ(3,*) (CONC_MATRIX(LO,i,j,k), k=1,NLOCAL)
	READ(3,*) (SD_CONC_MATRIX(LO,i,j,k), k=1,NLOCAL)	
	ENDDO
!
	IF (CONC_TYPE(LO).EQ.'SOIL') THEN
	KEYCONC(i,1)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CSOIL (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CSOIL (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'DRINKING_WATER') THEN
	KEYCONC(i,2)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CWATER (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CWATER (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'BATH_WATER') THEN
	KEYCONC(i,13)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CWATERDER (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CWATERDER (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'OTHER_WATERS') THEN
	KEYCONC(i,14)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CWATEROTHER (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CWATEROTHER (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO

!
	ELSEIF (CONC_TYPE(LO).EQ.'PARTICULATE') THEN
	KEYCONC(i,3)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CPAR (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CPAR (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'STEAM') THEN
	KEYCONC(i,12)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CSTEAM (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CSTEAM (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'FRUIT') THEN
	KEYCONC(i,4)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CFRUIT (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CFRUIT (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'LEAVES') THEN
	KEYCONC(i,5)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CLEAVES (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CLEAVES (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'MEAT') THEN
	KEYCONC(i,6)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CBEEF(i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CBEEF (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'MILK') THEN
	KEYCONC(i,7)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CMILK (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CMILK (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'BIRD') THEN
	KEYCONC(i,8)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CAVE (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CAVE (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'EGG') THEN
	KEYCONC(i,9)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CEGG (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CEGG (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'FISH') THEN
	KEYCONC(i,10)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CFISH (i,j,k)=CONC_MATRIX(LO,i,j,k)
	SD_CFISH (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'GRAIN') THEN
	KEYCONC(i,11)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CGRAIN (i,j,k)=CONC_MATRIX(LO,i,j,k)
    SD_CGRAIN (i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ELSEIF (CONC_TYPE(LO).EQ.'SEDIMENTS') THEN
	KEYCONC(i,15)= .TRUE.
	DO j=1,NDURATION
	DO k=1,NLOCAL
	CSEDIMENT(i,j,k)=CONC_MATRIX(LO,i,j,k)
    SD_CSEDIMENT(i,j,k)=SD_CONC_MATRIX(LO,i,j,k)
	ENDDO
	ENDDO
!
	ENDIF  !FIM DO IF ATRIBUIÇÃO DAS CONCENTRAÇÃO AS MATRIZES CERTAS
!
    ELSE
!
    GOTO 101 
!
    ENDIF	! FIM DO IF CASO NÃO EXISTA A MATRIZ INFORMADA
!
    ENDDO  !ACABA O CICLO "LO"
!
101 CONTINUE
!
    ENDDO  !ACABA O CICLO DOS METAIS
!
!
!
    TIMESP(0)=0.0
!
    DO KKK=1,NTIME
	TIMESP(KKK)=TIMESP(KKK-1)+VARIASAO
	ENDDO
!
    IF(NDURATION.EQ.1)THEN
!
	DO i=1,NCHEM
	DO k=1,NLOCAL
	DO j=2,NTIME
	CSOIL (i,j,k)= CSOIL (i,1,k)
	SD_CSOIL(i,j,k)=SD_CSOIL(i,1,k)
	CWATER (i,j,k)= CWATER (i,1,k)
	SD_CWATER(i,j,k)=SD_CWATER(i,1,k)
	CWATERDER (i,j,k)= CWATERDER (i,1,k)
	SD_CWATERDER(i,j,k)=SD_CWATERDER(i,1,k)
	CWATEROTHER (i,j,k)= CWATEROTHER (i,1,k)
	SD_CWATEROTHER(i,j,k)=SD_CWATEROTHER(i,1,k)
	CPAR (i,j,k)= CPAR (i,1,k)
	SD_CPAR(i,j,k)=SD_CPAR(i,1,k)
	CSTEAM (i,j,k)= CSTEAM (i,1,k)
	SD_CSTEAM(i,j,k)=SD_CSTEAM(i,1,k)
    CFRUIT (i,j,k)= CFRUIT (i,1,k)
	SD_CFRUIT(i,j,k)=SD_CFRUIT(i,1,k)
	CLEAVES (i,j,k)= CLEAVES (i,1,k)
	SD_CLEAVES(i,j,k)=SD_CLEAVES(i,1,k)
    CBEEF (i,j,k)= CBEEF (i,1,k)
	SD_CBEEF(i,j,k)=SD_CBEEF(i,1,k)
    CMILK (i,j,k)= CMILK (i,1,k)
	SD_CMILK(i,j,k)=SD_CMILK(i,1,k)
    CAVE (i,j,k)= CAVE (i,1,k)
	SD_CAVE(i,j,k)=SD_CAVE(i,1,k)
    CEGG(i,j,k)= CEGG(i,1,k)
	SD_CEGG(i,j,k)=SD_CEGG(i,1,k)
	CFISH (i,j,k)= CFISH (i,1,k)
	SD_CFISH(i,j,k)=SD_CFISH(i,1,k)
	CGRAIN (i,j,k)= CGRAIN (i,1,k)
	SD_CGRAIN(i,j,k)=SD_CGRAIN(i,1,k)
	CSEDIMENT(i,j,k)= CSEDIMENT(i,1,k)
	SD_CSEDIMENT(i,j,k)=SD_CSEDIMENT(i,1,k)
	ENDDO
	ENDDO
	ENDDO
!
    ELSEIF ((NDURATION.NE.1).AND.(NDURATION.LT.NTIME))THEN
!
	DO i=1,NCHEM
	DO k=1,NLOCAL
	DO j=NDURATION+1,NTIME
	CSOIL (i,j,k)= CSOIL (i,NDURATION,k)
	SD_CSOIL(i,j,k)=SD_CSOIL(i,NDURATION,k)
	CWATER (i,j,k)= CWATER (i,NDURATION,k)
	SD_CWATER(i,j,k)=SD_CWATER(i,NDURATION,k)
	CWATERDER (i,j,k)= CWATERDER (i,NDURATION,k)
	SD_CWATERDER(i,j,k)=SD_CWATERDER(i,NDURATION,k)
	CWATEROTHER (i,j,k)= CWATEROTHER (i,NDURATION,k)
	SD_CWATEROTHER(i,j,k)=SD_CWATEROTHER(i,NDURATION,k)
	CPAR (i,j,k)= CPAR (i,NDURATION,k)
	SD_CPAR(i,j,k)=SD_CPAR(i,NDURATION,k)
	CSTEAM (i,j,k)= CSTEAM (i,NDURATION,k)
	SD_CSTEAM(i,j,k)=SD_CSTEAM(i,NDURATION,k)
    CFRUIT (i,j,k)= CFRUIT (i,NDURATION,k)
	SD_CFRUIT(i,j,k)=SD_CFRUIT(i,NDURATION,k)
	CLEAVES (i,j,k)= CLEAVES (i,NDURATION,k)
	SD_CLEAVES(i,j,k)=SD_CLEAVES(i,NDURATION,k)
    CBEEF (i,j,k)= CBEEF (i,NDURATION,k)
	SD_CBEEF(i,j,k)=SD_CBEEF(i,NDURATION,k)
    CMILK (i,j,k)= CMILK (i,NDURATION,k)
	SD_CMILK(i,j,k)=SD_CMILK(i,NDURATION,k)
    CAVE (i,j,k)= CAVE (i,NDURATION,k)
	SD_CAVE(i,j,k)=SD_CAVE(i,NDURATION,k)
    CEGG(i,j,k)= CEGG(i,NDURATION,k)
	SD_CEGG(i,j,k)=SD_CEGG(i,NDURATION,k)
	CFISH (i,j,k)= CFISH (i,NDURATION,k)
	SD_CFISH(i,j,k)=SD_CFISH(i,NDURATION,k)
	CGRAIN (i,j,k)= CGRAIN (i,NDURATION,k)
	SD_CGRAIN(i,j,k)=SD_CGRAIN(i,NDURATION,k)
	CSEDIMENT(i,j,k)= CSEDIMENT(i,NDURATION,k)
	SD_CSEDIMENT(i,j,k)=SD_CSEDIMENT(i,NDURATION,k)
	ENDDO
	ENDDO
	ENDDO
!
!
    ENDIF
!
      CLOSE(3)
!
!
      RETURN 
	  END
!
!
!***********************************************************************************************
!***********************************************************************************************
!
      SUBROUTINE READDATABASE (NVIAS,NIDADE,POLLUTANT,TYPE_POLLUTANT,RAD_POL,RfD,SF,SD_BTF,&
      BTF,ED,AT,PC,BW,ABS_,Kd,fw,NPOL,SD_RfD,SD_SF,SD_ED,SD_AT,SD_PC,SD_ABS,SD_Kd,SD_fw,SD_BW,BAF,SD_BAF,MUTAGENIC,AMOLAR,VIDAMEIA)  
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
	  CHARACTER(LEN=50)  :: POLLUTANT(500),RAD_POL(4)
	  CHARACTER(LEN=20)  :: TYPE_POLLUTANT(500)
	  LOGICAL MUTAGENIC
	  REAL*8 Kd
      DIMENSION RfD(500,6,3),SF(500,6,3),BAF(500,14), SD_BAF(500,14)
      DIMENSION SD_RfD(500,6,3),SD_SF(500,6,3),SF_MEDIA(500,6,3),SF_DESVIO(500,6,3),SIGMA_L(500,6,3),SD_NI(500,6,3)  
	  DIMENSION BTF(500,15),SD_BTF(500,15)
	  DIMENSION PC(500),SD_PC(500),ABS_(500), Kd(500), fw(500)
      DIMENSION	SD_ABS(500), SD_Kd(500),SD_fw(500)
!
	  DIMENSION ED(3,NIDADE), AT(2,NVIAS,NIDADE), BW(NIDADE)
	  DIMENSION SD_ED(3,NIDADE), SD_AT(2,NVIAS,NIDADE), SD_BW(NIDADE)
!
      DIMENSION MUTAGENIC(500,3),AMOLAR(4),VIDAMEIA(4)
!
!
!
      OPEN(UNIT=4,STATUS='OLD',FILE='Datachemical.prn')     
!
!	
      READ(4,*)
      READ(4,*)
      READ(4,*)  ! READ USADOS PARA COLOCAR COMENTARIOS E FAZER ABERTURA DO FILE
      READ(4,*)
      READ(4,*)
      READ(4,*)
	  READ(4,*)
      READ(4,*)
	  READ(4,*)
      READ(4,*) NPOL,KEY_SF
!
!     LEITURA DE POLUENTES NAO CANCERIGENOS   
!
      !
	  do i=1,500
	  TYPE_POLLUTANT(i)="NON-CUMULATIVE"
      DO LOK=1,15
	  BTF(i,LOK)=0.0
	  SD_BTF(i,LOK)=0.0
	  ENDDO
      DO IK=1,14
	  BAF(i,IK)=0.0
	  SD_BAF(i,IK)=0.0
	  ENDDO
	  DO ku=1,3
      DO j=1,6
      RfD(i,j,ku)=0.0
	  SD_RfD(i,j,ku)=0.0
	  SF(i,j,ku)=0.0
	  SD_SF(i,j,ku)=0.0
	  ENDDO	
	  ENDDO
	  enddo
!
!76    FORMAT (L7,L9)
!
!
      DO i=1,NPOL
!
!
      READ(4,*)
	  READ(4,*)
	  READ(4,*) POLLUTANT(i)
!
!
      READ(4,*)
      READ(4,*)
	  READ(4,*)	MUTAGENIC(i,1),MUTAGENIC(i,2)
!
      MUTAGENIC(i,3)=.FALSE.
!
!
!    WRITE(*,*) POLLUTANT(i)
!
      READ(4,*)
	  READ(4,*)
      READ(4,*)		  
	  READ(4,*) (BAF(i,lk), lk=1,14)
!
      DO JUJ=1,14
      IF((BAF(i,JUJ).GT.1.0).OR.(BAF(i,JUJ).LT.0.0))THEN
	  WRITE(*,*)
	  WRITE(*,'('' Chemical =  '')')POLLUTANT(i)
	  WRITE(*,*)
	  WRITE(*,'('' BAF parameter must be less than or equal to 1.000 and greater than or equal to zero '')')
	  WRITE(*,*)
	  WRITE(*,'('' For more information, enable key information --- Help --- in the Scenary Sheet '')')
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
	  ENDDO
!
	  READ(4,*)
      READ(4,*)
      READ(4,*) (SD_BAF(i,LP), LP=1,14)
!
      DO IIO=5,6
	  BAF(i,IIO)=1.0		! forçando o BAF numero 5 e 6 que correspondem a contato dérmico de água e solo a ser igual a 1 porque o PC e ABS presentes nas equações das doses ja fazer a função dos BAF
	  SD_BAF(i,IIO)=0.0		! forçando a incerteza dos BAF numero 5 e 6 que correspondem a contato dérmico de água e solo a ser igual a zero porque o PC e ABS presentes nas equações das doses ja fazer a função dos BAF
	  ENDDO
!
      READ(4,*)
      READ(4,*)
	  READ(4,*)
!
!
	  DO k=1,3
	  READ(4,*) (RfD(i,KO,k), KO=1,6)
	  ENDDO
!
!
      READ(4,*)
      READ(4,*)
!
	  DO k=1,3
	  READ(4,*) (SD_RfD(i,j,k), j=1,6)
	  ENDDO
!
      READ(4,*)
      READ(4,*)
      READ(4,*)
!
	DO k=1,3				     ! j=1 oral SOLO, j=2 inalaçao PARTICULAS, j=3 dermico SOLO, j=4 oral água, j=5 inalação VAPORES, j=6, dermico água
	READ(4,*) (SF(i,JJ,k), JJ=1,6)	! k=1 cronica, k=2 sub-cronica, k=3 aguda (exposição)
	ENDDO
!
!
	IF((KEY_SF.LT.0).OR.(KEY_SF.GT.1))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' key SF uncertainties MUST BE Disable OR Active '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
	ENDIF
!
    IF(KEY_SF.EQ.1)THEN
      READ(4,*)
      READ(4,*)
      READ(4,*)
	  READ(4,*)
      READ(4,*)
!
	DO K=1,3
	do j=1,6
	IF(SF(i,j,k).GT.0.0)THEN
	SF_MEDIA(i,j,k)=SF(i,j,k)/3.425
	SF_DESVIO(i,j,k)=SF(i,j,k)/2.060
	SIGMA_L(i,j,k)=SQRT(((SF_DESVIO(i,j,k)/SF_MEDIA(i,j,k))**2)+1)
	SD_NI(i,j,k)=LOG(SF_MEDIA(i,j,k))-(SIGMA_L(i,j,k)**2/2)
	SD_SF(i,j,k)=SIGMA_L(i,j,k)*EXP(SD_NI(i,j,k))
	ELSE
	SD_SF(i,j,k)=0.0
	ENDIF
	enddo
	ENDDO
!
    ELSE IF(KEY_SF.EQ.0)THEN
      READ(4,*)
      READ(4,*)
!
	DO k=1,3				     ! j=1 oral SOLO, j=2 inalaçao PARTICULAS, j=3 dermico SOLO, j=4 oral água, j=5 inalação VAPORES, j=6, dermico água
	READ(4,*) (SD_SF(i,JU,k), JU=1,6)	! k=1 cronica, k=2 sub-cronica, k=3 aguda (exposição)
	ENDDO
!
!
!
!
	ENDIF  ! FIM DO IF(KEY_SF)
!
!
      READ(4,*)
	  READ(4,*)
      READ(4,*) PC(i), ABS_(i), fw(i), Kd(i)
!
!
      READ(4,*)
      READ(4,*) SD_PC(i), SD_ABS(i), SD_fw(i), SD_Kd(i)
!
!
!                                     FATORES DE BIOTRANSFERENCIA - BTF
!
      READ(4,*)
	  READ(4,*)
	  DO m=1,15
      READ(4,*) BTF(i,m),SD_BTF(i,m)
	  ENDDO
!	 
!     WRITE(*,*) POLLUTANT(I)
!     WRITE(*,*) (RFD(1,JT,1), JT=1,6)
!
      ENDDO	 
!
     close(4)
!
!    DATABASE RADIOLÓGICO COM AS MASSAS MOLARES E OS TEMPOS DE MEIA VIDA DOS RADIONÚCLIDEOS (U-238, Th-232, Ra-226 E K-40)
!
	 RAD_POL(1)="K"
	 RAD_POL(2)="Th"
	 RAD_POL(3)="U"
	 RAD_POL(4)="Ra"
	 AMOLAR(1)=3.99630E+4
	 AMOLAR(2)=2.32038E+5
	 AMOLAR(3)=2.38028E+5
	 AMOLAR(4)=2.26030E+5
	 VIDAMEIA(1)=4.0366E+16
	 VIDAMEIA(2)=4.4340E+17
	 VIDAMEIA(3)=1.4160E+17
     VIDAMEIA(4)=5.1151E+10
!
!	
!				
!	 LEITURA DOS PARAMETROS DE EXPOSIÇÃO 
!
    OPEN(UNIT=5, FILE='Dataexp.prn') 
!
!		         
!                     ROTAS
!     i=1  ingestao de agua contaminada 
!     i=2  ingestao de agua durante a natação 
!     i=3  contato dermico com agua contaminada 
!     i=4  ingestao com solo contaminado 
!     i=5  contato dermico com solo contaminado  
!     i=6  inalaçao de vapores presentes no ar
!     i=7  inalaçao de particulas presentes no ar
!
!
    READ(5,67)
!
67	FORMAT(10/)
!
    KSTOPED=0
! 
	READ(5,*) ED(1,1),SD_ED(1,1),ED(1,2),SD_ED(1,2),ED(1,3),SD_ED(1,3),ED(1,4),SD_ED(1,4),ED(1,5),SD_ED(1,5),ED(1,6),SD_ED(1,6),ED(1,7),SD_ED(1,7),ED(1,8),SD_ED(1,8),ED(1,9),SD_ED(1,9)	! j= idade da população 
!
	DO i=2,3
	DO K=1,NIDADE			  
	ED(i,K)=ED(1,K)
	SD_ED(i,K)=SD_ED(1,K)
	ENDDO
	ENDDO
!
    READ(5,*)
	READ(5,*)
    READ(5,*)
	READ(5,*)
!
!
	READ(5,*) AT(2,1,1),SD_AT(2,1,1),AT(2,1,2),SD_AT(2,1,2),AT(2,1,3),SD_AT(2,1,3),AT(2,1,4),SD_AT(2,1,4),AT(2,1,5),SD_AT(2,1,5),AT(2,1,6),SD_AT(2,1,6),AT(2,1,7),SD_AT(2,1,7),AT(2,1,8),SD_AT(2,1,8),AT(2,1,9),SD_AT(2,1,9)
!   
    KSTOPAT=0
	DO j=1,NIDADE
	IF((AT(2,1,j).EQ.0.0).AND.(KSTOPAT.EQ.0))THEN
	KSTOPAT=1
	WRITE(99,'('' WARNING!!! SOME VALUES OF AT carcinogenic ARE EQUAL TO ZERO!!!!! '')')
	WRITE(99,*)
	ENDIF
	ENDDO
!
	DO I=2,NVIAS
	DO K=1,NIDADE			  
	AT(2,I,K)=AT(2,1,K)
	SD_AT(2,I,K)=SD_AT(2,1,K)
	ENDDO
	ENDDO
!
    READ(5,*)
	READ(5,*)
	READ(5,*)
	READ(5,*)
!
	READ(5,*) BW(1),SD_BW(1),BW(2),SD_BW(2),BW(3),SD_BW(3),BW(4),SD_BW(4),BW(5),SD_BW(5),BW(6),SD_BW(6),BW(7),SD_BW(7),BW(8),SD_BW(8),BW(9),SD_BW(9)
!
    KSTOPBW=0
    DO I=1,NIDADE
	IF((BW(I).EQ.0.0).AND.(KSTOPBW.EQ.0))THEN
	KSTOPBW=1
	WRITE(99,'('' WARNING!!! SOME VALUE OF BW ARE EQUAL TO ZERO!!!!! '')')
	WRITE(99,*)
	ENDIF
	ENDDO
!
!
!
!
	close(5)
!
!
    RETURN
	END
!
!
!**********************************************************************
!**********************************************************************
!
!
      SUBROUTINE EXPOSURE(W,CF,SD_CF,FI,SD_FI,IR_INI,SD_IR_INI,FA,SD_FA,FP,SD_FP,ET_INI,SD_ET_INI,SA_INI,SD_SA_INI,AF_INI,SD_AF_INI,EV_INI,SD_EV_INI,&
	  NVP,NIDADE,NVIAS,SCENAR,NCHEM,NTIME,NDURATION,VARIASAO,NLOCAL,NTYPECONC,HQ,SD_HQ,CR,SD_CR,CHEMICAL,TIMESP,NRISKTYPE,EF_INI,SD_EF_INI,AT_INI,SD_AT_INI,INICIO)
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL W,MUTAGENIC
	  LOGICAL KEYCONC,KEYCONC_2,KEYCONC_6,KEYCONC_B3,KEYCONC_S3,KEYCONC_W3
	  LOGICAL KEYCONC_M4,KEYCONC_S4,KEYCONC_W4,KEYCONC_F7,KEYCONC_AVE8,KEYCONC_G10,KEYCONC_EGG9
!
	  INTEGER SCENAR
	  REAL*8 Kd,IR, IR_INI
	  CHARACTER(LEN=50)  :: CHEMICAL(500),NOMES(500),POLLUTANT(500),RAD_POL(4)
	  CHARACTER(LEN=20)  ::	TYPE_POLLUTANT(500)
	  CHARACTER(LEN=10)  :: NRISKTYPE(NVIAS) 
!
	  CHARACTER aspas*1
!
!
	DIMENSION CSOIL(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL)
	DIMENSION CWATER(NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_CWATER(NCHEM,NTIME,NLOCAL)
	DIMENSION CWATERDER(NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_CWATERDER(NCHEM,NTIME,NLOCAL)
    DIMENSION CWATEROTHER(NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_CWATEROTHER(NCHEM,NTIME,NLOCAL)
	DIMENSION CPAR(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CPAR(NCHEM,NTIME,NLOCAL)
	DIMENSION CSTEAM(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CSTEAM(NCHEM,NTIME,NLOCAL)
	DIMENSION CFRUIT(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CFRUIT(NCHEM,NTIME,NLOCAL)
	DIMENSION CLEAVES(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CLEAVES(NCHEM,NTIME,NLOCAL)
	DIMENSION CBEEF(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CBEEF(NCHEM,NTIME,NLOCAL)
	DIMENSION CMILK(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CMILK(NCHEM,NTIME,NLOCAL)
	DIMENSION CAVE(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CAVE(NCHEM,NTIME,NLOCAL)
	DIMENSION CEGG(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CEGG(NCHEM,NTIME,NLOCAL)
	DIMENSION CFISH(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CFISH(NCHEM,NTIME,NLOCAL)
	DIMENSION CGRAIN(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CGRAIN(NCHEM,NTIME,NLOCAL)
	DIMENSION CSEDIMENT(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CSEDIMENT(NCHEM,NTIME,NLOCAL)
!
	DIMENSION ADD(NIDADE,2,NVIAS,NCHEM,NTIME,NLOCAL)
	    DIMENSION SD_ADD(NIDADE,2,NVIAS,NCHEM,NTIME,NLOCAL)
!
      DIMENSION RfD(500,6,3),SF(500,6,3), W(NVIAS),BAF(500,14), SD_BAF(500,14)
      DIMENSION SD_RfD(500,6,3), SD_SF(500,6,3) 
	  DIMENSION BTF(500,15),SD_BTF(500,15)
	  DIMENSION SD_ED_EXP(NVIAS,NCHEM)
	  DIMENSION PC(500),SD_PC(500)
      DIMENSION	ABS_(500), Kd(500),fw(500),Fa(4),SD_Fa(4),Fp(4),SD_Fp(4)  
      DIMENSION	SD_ABS(500), SD_Kd(500),SD_fw(500)
	  DIMENSION CF(3),SD_CF(3),FI(9),SD_FI(9)
	  DIMENSION KSTOPRFD(6),KSTOPSF(6),KSTOPBAF(14)
!
	  DIMENSION IR_INI(NIDADE,20),SD_IR_INI(NIDADE,20),ET_INI(NIDADE,3),SD_ET_INI(NIDADE,3),SA_INI(NIDADE,2),SD_SA_INI(NIDADE,2)
	  DIMENSION EV_INI(NIDADE,2),SD_EV_INI(NIDADE,2),AF_INI(NIDADE),SD_AF_INI(NIDADE)
!
	  DIMENSION IR(20),SD_IR(20),ET(3),SD_ET(3),SA(2),SD_SA(2),EV(2),SD_EV(2)
!
	  DIMENSION ED_INI(3,NIDADE), EF_INI(NVIAS,NIDADE), AT_INI(2,NVIAS,NIDADE), BW_INI(NIDADE)
	  DIMENSION SD_ED_INI(3,NIDADE), SD_EF_INI(NVIAS,NIDADE), SD_AT_INI(2,NVIAS,NIDADE), SD_BW_INI(NIDADE)
!
	  DIMENSION EF(NVIAS), AT(2,NVIAS)
	  DIMENSION SD_EF(NVIAS), SD_AT(2,NVIAS),J_INI_AGE(NIDADE)
!
	DIMENSION KEYCONC(NCHEM,NTYPECONC), KSTOPWAY(NVIAS), NSTOPNVP(NCHEM+1,NIDADE),LTSP(20), MUTAGENIC(500,3),ADAF(3) 
!
!
    DIMENSION HQ(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL),CR(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL)
	     DIMENSION SD_HQ(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL), SD_CR(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL)
	DIMENSION TIMESP(0:NTIME), DELTAT(NTIME),LOCALSP(NLOCAL),AMOLAR(4),VIDAMEIA(4)
!
!
!
      KOP=1
!
!
      CALL CONCENTRATION (VARIASAO,NDURATION,NCHEM,NTYPECONC,CHEMICAL,CSOIL,CWATER, &
	  CPAR,CSTEAM,CFRUIT,CLEAVES,CBEEF,CMILK,CAVE,CEGG,CFISH,CGRAIN,CWATERDER,CWATEROTHER,CSEDIMENT,&
	  KEYCONC,NTIME,NLOCAL,TIMESP,SD_CSOIL,&
      SD_CWATER,SD_CPAR,SD_CSTEAM,SD_CFRUIT,SD_CLEAVES,SD_CBEEF,SD_CMILK,&
	  SD_CAVE,SD_CEGG,SD_CFISH,SD_CGRAIN,SD_CWATERDER,SD_CWATEROTHER,SD_CSEDIMENT)
!
!
      CALL READDATABASE (NVIAS,NIDADE,POLLUTANT,TYPE_POLLUTANT,RAD_POL,RfD,SF,SD_BTF,&
      BTF,ED_INI,AT_INI,PC,BW_INI,ABS_,Kd,fw,NPOL,&
	  SD_RfD,SD_SF,SD_ED_INI,SD_AT_INI,SD_PC,SD_ABS,SD_Kd,SD_fw,SD_BW_INI,BAF,SD_BAF,MUTAGENIC,AMOLAR,VIDAMEIA) 
!
!
      WRITE(*,*)
	  WRITE(*,'("____________________________________________________________________________________")')   
	  WRITE(*,'("                           HUMAN HEALTH RISK ASSESMENT ")')   
	  WRITE(*,'("____________________________________________________________________________________")')   
!
      WRITE(99,*)
	  WRITE(99,'("____________________________________________________________________________________________________________________________________________")')
	  WRITE(99,'("                                                         HUMAN HEALTH RISK ASSESMENT")')
	  WRITE(99,'("____________________________________________________________________________________________________________________________________________")')
!
! 
	   TIMESP(0)=0.0
!
!
!      CRIAÇÃO DE MATRIZES
!
       DO l=1,NIDADE
       DO i=1,NCHEM
	   DO j=1,NTIME
	   DO k=1,NLOCAL
	   DO n=1,NVIAS		! NUMERO DE CAMINHOS
	   DO m=1,2
       ADD(l,m,n,i,j,k)=0.0
	      SD_ADD(l,m,n,i,j,k)=0.0
	   HQ(l,n,i,j,k)=0.0
	      SD_HQ(l,n,i,j,k)=0.0
	   CR(l,n,i,j,k)=0.0
	      SD_CR(l,n,i,j,k)=0.0
       SD_ED_EXP(n,i)=0.0
	   ENDDO
	   ENDDO
	   ENDDO
	   ENDDO
	   ENDDO
	   ENDDO
!
!
      IF(NVP.EQ.1)THEN
	  open (UNIT=17, file='Variable_Parameters.out')
!
	  WRITE(17,'("****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,'("                                                                       VARIABLE PARAMETERS")')
	  WRITE(17,'("****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,*)
	  ENDIF
!
      DO l=1,NIDADE
      NSTOPNVP(0,l)=0
	  ENDDO
!
!*******/////////////////////////////////*******************////////////////////////////////
!*******/////////////////////////////////*******************////////////////////////////////
!*******/////////////////////////////////*******************////////////////////////////////
!
!					                     CASO 1 (AGRICULTURAL)
!
!
     SELECT CASE(SCENAR)
	 CASE(1)
!    CHAMADAS A WAY DO CENARIO=1
!
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!----------------------------------- CICLOS POR METAIS -------------------------------------
!
      DO l=1,NIDADE
!
!------------------
	  NTIMEexp=ED_INI(1,l)
!------------------
!
!
!
!
!
		 LJ=NCHEM
		 NEW_POL=NPOL
!
      DO i=1,LJ
!
      NSTOPNVP(i,l)= NSTOPNVP(i-1,l)+1
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1))THEN
	  WRITE(17,*)
	  WRITE(17,*)
	  WRITE(17,'("****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,'("                                                                             CICLE ",I2)')l
	  WRITE(17,'("****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,*)
	  WRITE(17,*)
!
!
	  npk=20
	  nppp=14
!
	  DO lop=1,20
	  LTSP(lop)=lop
	  ENDDO
!
!	  DO loz=1,14
!	  LTSu(loz)=loz
!	  ENDDO
!

	  WRITE(17,'("________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
	  WRITE(17,'("   TIME      ATc     SD_ATc     BW   SD_BW    ET1    SD_ET1     ET2    SD_ET2    ET3    SD_ET3     SA1     SD_SA1     SA2    SD_SA2    EV1  SD_EV1   EV2  SD_EV2     AF      SD_AF     IR                                                                                                                                                                                 SD_IR    ")')
	  WRITE(17,'(175X,<npk>I9,<npk>I9)')(LTSP(kJy),kJy=1,npk),(LTSP(kJy),kJy=1,npk)
	  WRITE(17,'("________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
	  WRITE(17,*)
!
	  ENDIF
! 
      KSTOPBTF1=0
	  KSTOPBTF_6S=0
	  KSTOPBTF_SB3=0
	  KSTOPBTF_SV3=0
	  KSTOPBTF_VB3=0
	  KSTOPBTF_WB3=0
	  KSTOPBTF_SM4=0
	  KSTOPBTF_SV4=0
	  KSTOPBTF_VM4=0
	  KSTOPBTF_WM4=0
	  KSTOPBTF_WF7=0
	  KSTOPBTF_WAVE8=0
	  KSTOPBTF_SAVE8=0
	  KSTOPBTF_WEGG9=0
	  KSTOPBTF_SEGG9=0
	  KSTOPBTF_SG10=0
!
	  KSTOPfw3=0
      KSTOPfw4=0
	  KSTOPfw8=0
	  KSTOPfw9=0
      KSTOPPC13=0
	  KSTOPABS14=0
!
      DO LOL=1,6
      KSTOPRFD(LOL)=0
	  KSTOPSF(LOL)=0
	  ENDDO
!
      DO KO=1,NVIAS
	  KSTOPWAY(KO)=0
	  ENDDO
!
	  DO KOK=1,14
	  KSTOPBAF(KOK)=0
	  ENDDO
!
!      write(*,*) l,EF_INI(1,l)
!
      DO m=1,2
!
!----------------------------------------------------------------------------------------------------------------------
!
!-----------------
!	  AT(m,1)=AT_INI(m,1,l)
!	  SD_AT(m,1)=SD_AT_INI(m,1,l)
!-----------------
!
      iii=0
!
      DO ii=1,NEW_POL
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
	  IF (iii.eq.0) THEN
!
	  IF (l.eq.1) THEN
      WRITE(*,*)
	  WRITE(*,'(" WARNING!!! CHEMICAL NOT EXIST IN Datachemical DATABASE  --->  ",A30)') CHEMICAL(i) 
	  WRITE(*,'(" The HQ and CR values will all be zero for this Chemical species!  ")')
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! CHEMICAL NOT EXIST IN Datachemical RISK DATABASE  --->  ",A30)') CHEMICAL(i)
	  WRITE(99,'(" The HQ and CR values will all be zero for this Chemical species!  ")')
      WRITE(99,*)
!
      ENDIF
!
	  GOTO 53
	  ENDIF
!
!
!---------------------------------						   
	  DO j=1,NTIMEexp			  ! CICLO POR TEMPO
!
	  DELTAT(J)=1.0
!
!
      CALL REDISTRI(iii,l,j,NVIAS,NIDADE,EF_INI,SD_EF_INI,BW_INI,SD_BW_INI,IR_INI,SD_IR_INI,ET_INI,SD_ET_INI,SA_INI,SD_SA_INI,&
	  EV_INI,SD_EV_INI,AF_INI,SD_AF_INI,EF,SD_EF,BW,SD_BW,IR,SD_IR,ET,SD_ET,SA,SD_SA,EV,SD_EV,AF,SD_AF,MUTAGENIC,ADAF,&
	  AT_INI,SD_AT_INI,AT,SD_AT) 
!
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1).AND.(m.EQ.1))THEN
!
      write(17,'( 3x,I3,3x,F9.1,1x,F9.1,3x,F4.1,2x,F4.1,2x,E8.3,5(1x,E8.3),4(1x,E8.3),3x,F3.1,3x,F3.1,5x,F3.1,3x,F3.1,4x,E8.3,1x,E8.3,20(1x,E8.3),20(1x,E8.3))') j,AT(2,1),SD_AT(2,1),BW,SD_BW,ET(1),SD_ET(1),ET(2),SD_ET(2),ET(3),SD_ET(3),SA(1),SD_SA(1),SA(2),SD_SA(2),EV(1),SD_EV(1),EV(2),SD_EV(2),AF,SD_AF,(IR(iou), iou=1,20),(SD_IR(IQW), IQW=1,20)
!
!
      ENDIF
!
	  DO k=1,NLOCAL
!
!      write(*,'(" k= ",i3)') k
!
!
!*******************************************************************************************
!*******************************************************************************************
      IF(W(1).EQV..TRUE.)THEN
	  IF(KEYCONC(i,1).EQV..TRUE.) THEN
!
      EF1=EF(1)
	  SD_EF1=SD_EF(1)
!
	  ED1=DELTAT(J)
	  SD_ED1=0.0
	  SD_ED_EXP(1,i)=SD_ED1
!
      IF(m.EQ.1)THEN
	  AT1=AT(m,1)*DELTAT(J)
	  ELSE
	  AT1=AT(m,1)
	  ENDIF
	  SD_AT1=SD_AT(m,1)
	  R_IR1=IR(1)
	  SD_IR1=SD_IR(1) 
!
	  FI1=FI(1)
	  SD_FI1=SD_FI(1)
!
      CF1=CF(1)
	  SD_CF1=SD_CF(1)
!
!
      CALL DOSEWAY1(NCHEM,NTIME,NLOCAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF1,SD_EF1,ED1,SD_ED1,AT1,SD_AT1,R_IR1,SD_IR1,&
	  FI1,SD_FI1,CF1,SD_CF1,AD,SD_AD)
!
!
      ADD(l,m,1,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!
      SD_ADD(l,m,1,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
!  
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(1).EQ.'Chronic')THEN
	  RFD1=RfD(iii,1,1)
	  SD_RFD1=SD_RfD(iii,1,1)
	  SF1=SF(iii,1,1)
	  SD_SF1=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(1).EQ.'Subchronic')THEN
	  RFD1=RfD(iii,1,2)
	  SD_RFD1=SD_RfD(iii,1,2)
	  SF1=SF(iii,1,2)
	  SD_SF1=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(1).EQ.'Acute')THEN
	  RFD1=RfD(iii,1,3)
	  SD_RFD1=SD_RfD(iii,1,3)
	  SF1=SF(iii,1,3)
	  SD_SF1=SD_SF(iii,1,3)
	  ENDIF
!
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,1,i,j,k).GT.0.0).AND.(RFD1.GT.0.0).AND.(BAF(iii,1).GT.0.0))THEN
	  HQ(l,1,i,j,k)=ADD(l,m,1,i,j,k)*BAF(iii,1)/RFD1   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,1,i,j,k)=HQ(l,1,i,j,k)* SQRT((SD_ADD(l,m,1,i,j,k)/ADD(l,m,1,i,j,k))**2+(SD_RFD1/RFD1)**2+(SD_BAF(iii,1)/BAF(iii,1))**2)
	  ELSE
	  HQ(l,1,i,j,k)=0.0
	  SD_HQ(l,1,i,j,k)=0.0
	  IF((RFD1.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,1,i,j,k).GT.0.0).AND.(SF1.GT.0.0).AND.(BAF(iii,1).GT.0.0))THEN
	  CR(l,1,i,j,k)=ADD(l,m,1,i,j,k)*BAF(iii,1)*SF1*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,1,i,j,k)=CR(l,1,i,j,k)* SQRT((SD_ADD(l,m,1,i,j,k)/ADD(l,m,1,i,j,k))**2+(SD_SF1/SF1)**2+(SD_BAF(iii,1)/BAF(iii,1))**2)
	  ELSE
	  CR(l,1,i,j,k)=0.0
	  SD_CR(l,1,i,j,k)=0.0
	  IF((SF1.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
	  ENDIF
!
	  IF((BAF(iii,1).LE.0.0).AND.(KSTOPBAF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The soil ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(1)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway1 was not calculed, because soil concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY1 WAS NOT CALCULED, BECAUSE SOIL CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
      IF(W(2).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,4).EQV..TRUE.)) THEN
!
      EF2=EF(2)
	  SD_EF2=SD_EF(2)
!
	  ED2=DELTAT(J)
	  SD_ED2=0.0
	  SD_ED_EXP(2,i)=SD_ED2
!
      IF(m.EQ.1)THEN
	  AT2=AT(m,2)*DELTAT(J)
	  ELSE
	  AT2=AT(m,2)
	  ENDIF
	  SD_AT2=SD_AT(m,2)
      R_IRfruit2=IR(2)
	  SD_IRfruit2=SD_IR(2)
!	  
	  BTF_2=BTF(iii,1)
	  SD_BTF_2=SD_BTF(iii,1)
	  KEYCONC_2=KEYCONC(i,4)
	  FI2=FI(2)
	  SD_FI2=SD_FI(2)
!
!
      CALL DOSEWAY2(CHEMICAL,KSTOPBTF1,NCHEM,NTIME,NLOCAL,KEYCONC_2,i,j,k,CSOIL,SD_CSOIL,CFRUIT,SD_CFRUIT,EF2,SD_EF2,ED2,SD_ED2,AT2,SD_AT2,&
	  BW,SD_BW,BTF_2,SD_BTF_2,R_IRfruit2,SD_IRfruit2,FI2,SD_FI2,AD,SD_AD)
!
!
!
      ADD(l,m,2,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!
      SD_ADD(l,m,2,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
!	  
!
!     CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(2).EQ.'Chronic')THEN
	  RFD2=RfD(iii,1,1)
	  SD_RFD2=SD_RfD(iii,1,1)
	  SF2=SF(iii,1,1)
	  SD_SF2=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(2).EQ.'Subchronic')THEN
	  RFD2=RfD(iii,1,2)
	  SD_RFD2=SD_RfD(iii,1,2)
	  SF2=SF(iii,1,2)
	  SD_SF2=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(2).EQ.'Acute')THEN
	  RFD2=RfD(iii,1,3)
	  SD_RFD2=SD_RfD(iii,1,3)
	  SF2=SF(iii,1,3)
	  SD_SF2=SD_SF(iii,1,3)
	  ENDIF     
! 
      IF(m.EQ.1)THEN
      IF((ADD(l,m,2,i,j,k).GT.0.0).AND.(RFD2.GT.0.0).AND.(BAF(iii,8).GT.0.0))THEN
      HQ(l,2,i,j,k)=ADD(l,m,2,i,j,k)*BAF(iii,8)/RFD2   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,2,i,j,k)=HQ(l,2,i,j,k)* SQRT((SD_ADD(l,m,2,i,j,k)/ADD(l,m,2,i,j,k))**2+(SD_RFD2/RFD2)**2+(SD_BAF(iii,8)/BAF(iii,8))**2)
	  ELSE
	  HQ(l,2,i,j,k)=0.0
	  SD_HQ(l,2,i,j,k)=0.0
	  IF((RFD2.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!						 
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,2,i,j,k).GT.0.0).AND.(SF2.GT.0.0).AND.(BAF(iii,8).GT.0.0))THEN
	  CR(l,2,i,j,k)=ADD(l,m,2,i,j,k)*BAF(iii,8)*SF2*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,2,i,j,k)=CR(l,2,i,j,k)* SQRT((SD_ADD(l,m,2,i,j,k)/ADD(l,m,2,i,j,k))**2+(SD_SF2/SF2)**2+(SD_BAF(iii,8)/BAF(iii,8))**2)
	  ELSE
	  CR(l,2,i,j,k)=0.0
	  SD_CR(l,2,i,j,k)=0.0
	  IF((SF2.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
	  ENDIF
!
	  IF((BAF(iii,8).LE.0.0).AND.(KSTOPBAF(8).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(8)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The fruit ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(2)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway2 was not calculed, because soil and fruits, seeds or tubers concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY2 WAS NOT CALCULED, BECAUSE SOIL AND FRUITS, SEEDS OR TUBERS'
      PRINT*, '  CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!*******************************************************************************************
!
      IF(W(5).EQV..TRUE.)THEN
	  IF(KEYCONC(i,2).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF5=EF(5)
	  SD_EF5=SD_EF(5)
!
	  ED5=DELTAT(J)
	  SD_ED5=0.0
	  SD_ED_EXP(5,i)=SD_ED5
!
      IF(m.EQ.1)THEN
	  AT5=AT(m,5)*DELTAT(J)
	  ELSE
	  AT5=AT(m,5)
	  ENDIF
!
	  SD_AT5=SD_AT(m,5)
	  R_IR5=IR(9)
	  SD_IR5=SD_IR(9)
!
!
!
!
      CALL DOSEWAY5(NCHEM,NTIME,NLOCAL,i,j,k,CWATER,SD_CWATER,BW,SD_BW,EF5,SD_EF5,ED5,SD_ED5,AT5,SD_AT5,R_IR5,SD_IR5,AD,SD_AD)
!
!
      ADD(l,m,5,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,5,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(5).EQ.'Chronic')THEN
	  RFD5=RfD(iii,4,1)
	  SD_RFD5=SD_RfD(iii,4,1)
	  SF5=SF(iii,4,1)
	  SD_SF5=SD_SF(iii,4,1)
      ELSEIF(NRISKTYPE(5).EQ.'Subchronic')THEN
	  RFD5=RfD(iii,4,2)
	  SD_RFD5=SD_RfD(iii,4,2)
	  SF5=SF(iii,4,2)
	  SD_SF5=SD_SF(iii,4,2)
      ELSEIF(NRISKTYPE(5).EQ.'Acute')THEN
	  RFD5=RfD(iii,4,3)
	  SD_RFD5=SD_RfD(iii,4,3)
	  SF5=SF(iii,4,3)
	  SD_SF5=SD_SF(iii,4,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,5,i,j,k).GT.0.0).AND.(RFD5.GT.0.0).AND.(BAF(iii,2).GT.0.0))THEN
      HQ(l,5,i,j,k)=ADD(l,m,5,i,j,k)*BAF(iii,2)/RFD5   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,5,i,j,k)=HQ(l,5,i,j,k)* SQRT((SD_ADD(l,m,5,i,j,k)/ADD(l,m,5,i,j,k))**2+(SD_RFD5/RFD5)**2+(SD_BAF(iii,2)/BAF(iii,2))**2)
	  ELSE
	  HQ(l,5,i,j,k)=0.0
	  SD_HQ(l,5,i,j,k)=0.0
	  IF((RFD5.LE.0.0).AND.(KSTOPRFD(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(2)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral water RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,5,i,j,k).GT.0.0).AND.(SF5.GT.0.0).AND.(BAF(iii,2).GT.0.0))THEN
	  CR(l,5,i,j,k)=ADD(l,m,5,i,j,k)*BAF(iii,2)*SF5*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,5,i,j,k)=CR(l,5,i,j,k)* SQRT((SD_ADD(l,m,5,i,j,k)/ADD(l,m,5,i,j,k))**2+(SD_SF5/SF5)**2+(SD_BAF(iii,2)/BAF(iii,2))**2)
	  ELSE
	  CR(l,5,i,j,k)=0.0
	  SD_CR(l,5,i,j,k)=0.0
	  IF((SF5.LE.0.0).AND.(KSTOPSF(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(2)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral water SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,2).LE.0.0).AND.(KSTOPBAF(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(2)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The ingestion water BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(5)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway5 was not calculed, because water (DRINKING_WATER) concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY5 WAS NOT CALCULED, BECAUSE WATER (DRINKING_WATER) CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!
!*******************************************************************************************
!
!
      IF(W(6).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,5).EQV..TRUE.))THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
	  EF6=EF(6)
	  SD_EF6=SD_EF(6)
!
	  ED6=DELTAT(J)
	  SD_ED6=0.0
	  SD_ED_EXP(6,i)=SD_ED6
!
      IF(m.EQ.1)THEN
	  AT6=AT(m,6)*DELTAT(J)
	  ELSE
	  AT6=AT(m,6)
	  ENDIF
!
	  SD_AT6=SD_AT(m,6)
      R_IRveg6=IR(10)
	  SD_IRveg6=SD_IR(10)
!
	  BTF_6S=BTF(iii,7)
	  SD_BTF_6S=SD_BTF(iii,7)
	  KEYCONC_6=KEYCONC(i,5)
	  FI6=FI(5)
	  SD_FI6=SD_FI(5)
!
!
      CALL DOSEWAY6(CHEMICAL,KSTOPBTF_6S,NCHEM,NTIME,NLOCAL,KEYCONC_6,i,j,k,CSOIL,SD_CSOIL,CLEAVES,SD_CLEAVES,EF6,SD_EF6,ED6,SD_ED6,&
	  AT6,SD_AT6,BW,SD_BW,BTF_6S,SD_BTF_6S,R_IRveg6,SD_IRveg6,FI6,SD_FI6,AD,SD_AD)
!
!
      ADD(l,m,6,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,6,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(6).EQ.'Chronic')THEN
	  RFD6=RfD(iii,1,1)
	  SD_RFD6=SD_RfD(iii,1,1)
	  SF6=SF(iii,1,1)
	  SD_SF6=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(6).EQ.'Subchronic')THEN
	  RFD6=RfD(iii,1,2)
	  SD_RFD6=SD_RfD(iii,1,2)
	  SF6=SF(iii,1,2)
	  SD_SF6=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(6).EQ.'Acute')THEN
	  RFD6=RfD(iii,1,3)
	  SD_RFD6=SD_RfD(iii,1,3)
	  SF6=SF(iii,1,3)
	  SD_SF6=SD_SF(iii,1,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,6,i,j,k).GT.0.0).AND.(RFD6.GT.0.0).AND.(BAF(iii,7).GT.0.0))THEN
	  HQ(l,6,i,j,k)=ADD(l,m,6,i,j,k)*BAF(iii,7)/RFD6   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,6,i,j,k)=HQ(l,6,i,j,k)* SQRT((SD_ADD(l,m,6,i,j,k)/ADD(l,m,6,i,j,k))**2+(SD_RFD6/RFD6)**2+(SD_BAF(iii,7)/BAF(iii,7))**2)	 
	  ELSE
	  HQ(l,6,i,j,k)=0.0
	  SD_HQ(l,6,i,j,k)=0.0
	  IF((RFD6.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF

!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,6,i,j,k).GT.0.0).AND.(SF6.GT.0.0).AND.(BAF(iii,7).GT.0.0))THEN
	  CR(l,6,i,j,k)=ADD(l,m,6,i,j,k)*BAF(iii,7)*SF6*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,6,i,j,k)=CR(l,6,i,j,k)* SQRT((SD_ADD(l,m,6,i,j,k)/ADD(l,m,6,i,j,k))**2+(SD_SF6/SF6)**2+(SD_BAF(iii,7)/BAF(iii,7))**2)
	  ELSE
	  CR(l,6,i,j,k)=0.0
	  SD_CR(l,6,i,j,k)=0.0
	  IF((SF6.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,7).LE.0.0).AND.(KSTOPBAF(7).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(7)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The vegetable ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(6)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway6 was not calculed, because soil and vegetables (leaves) concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY6 WAS NOT CALCULED, BECAUSE SOIL AND VEGETABLES (LEAVES) CONCENTRATIONS NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
!
      IF(W(3).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,6).EQV..TRUE.)) THEN
!
      EF3=EF(3)
	  SD_EF3=SD_EF(3)
!
	  ED3=DELTAT(J)
	  SD_ED3=0.0
	  SD_ED_EXP(3,i)=SD_ED3
!
      IF(m.EQ.1)THEN
	  AT3=AT(m,3)*DELTAT(J)
	  ELSE
	  AT3=AT(m,3)
	  ENDIF
!
	  SD_AT3=SD_AT(m,3)
      R_IRcarne3=IR(3)
	  SD_IRcarne3=SD_IR(3)
	  R_IRali3=IR(4)
	  SD_IRali3=SD_IR(4)
	  R_IRagua3=IR(11)
	  SD_IRagua3=SD_IR(11)
	  R_IRsolo3=IR(5)
	  SD_IRsolo3=SD_IR(5)
!
	  BTF_SV3=BTF(iii,2)
	  SD_BTF_SV3=SD_BTF(iii,2)
	  BTF_VB3=BTF(iii,3)
	  SD_BTF_VB3=SD_BTF(iii,3)
	  BTF_SB3=BTF(iii,4)
	  SD_BTF_SB3=SD_BTF(iii,4)
	  BTF_WB3=BTF(iii,8)
	  SD_BTF_WB3=SD_BTF(iii,8)
	  KEYCONC_B3=KEYCONC(i,6)
	  KEYCONC_S3=KEYCONC(i,1)
	  KEYCONC_W3=KEYCONC(i,14)
	  Fa3=Fa(1)
	  SD_Fa3=SD_Fa(1)
	  Fp3=Fp(1)
	  SD_Fp3=SD_Fp(1)
	  FI3=FI(3)
	  SD_FI3=SD_FI(3)
	  fw3=fw(iii)
	  SD_fw3=SD_fw(iii)
!
!
!
      CALL DOSEWAY3(CHEMICAL,KSTOPBTF_SB3,KSTOPBTF_SV3,KSTOPBTF_VB3,KSTOPBTF_WB3,KSTOPfw3,NCHEM,NTIME,NLOCAL,BTF_SV3,SD_BTF_SV3,BTF_VB3,SD_BTF_VB3,BTF_SB3,SD_BTF_SB3,BTF_WB3,SD_BTF_WB3,i,j,k,&
	  CSOIL,SD_CSOIL,CWATEROTHER,SD_CWATEROTHER,CBEEF,SD_CBEEF,EF3,SD_EF3,ED3,SD_ED3,AT3,SD_AT3,BW,SD_BW,KEYCONC_B3,KEYCONC_S3,KEYCONC_W3,&
	  Fa3,SD_Fa3,Fp3,SD_Fp3,R_IRali3,SD_IRali3,R_IRsolo3,SD_IRsolo3,R_IRagua3,SD_IRagua3,R_IRcarne3,SD_IRcarne3,FI3,SD_FI3,fw3,SD_fw3,AD,SD_AD)
!
!
      ADD(l,m,3,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,3,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(3).EQ.'Chronic')THEN
	  RFD3=RfD(iii,1,1)
	  SD_RFD3=SD_RfD(iii,1,1)
	  SF3=SF(iii,1,1)
	  SD_SF3=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(3).EQ.'Subchronic')THEN
	  RFD3=RfD(iii,1,2)
	  SD_RFD3=SD_RfD(iii,1,2)
	  SF3=SF(iii,1,2)
	  SD_SF3=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(3).EQ.'Acute')THEN
	  RFD3=RfD(iii,1,3)
	  SD_RFD3=SD_RfD(iii,1,3)
	  SF3=SF(iii,1,3)
	  SD_SF3=SD_SF(iii,1,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,3,i,j,k).GT.0.0).AND.(RFD3.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
      HQ(l,3,i,j,k)=ADD(l,m,3,i,j,k)*BAF(iii,3)/RFD3   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,3,i,j,k)=HQ(l,3,i,j,k)* SQRT((SD_ADD(l,m,3,i,j,k)/ADD(l,m,3,i,j,k))**2+(SD_RFD3/RFD3)**2+(SD_BAF(iii,9)/BAF(iii,9))**2)
	  ELSE
	  HQ(l,3,i,j,k)=0.0
	  SD_HQ(l,3,i,j,k)=0.0
	  IF((RFD3.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,3,i,j,k).GT.0.0).AND.(SF3.GT.0.0).AND.(BAF(iii,9).GT.0.0))THEN
	  CR(l,3,i,j,k)=ADD(l,m,3,i,j,k)*BAF(iii,9)*SF3*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,3,i,j,k)=CR(l,3,i,j,k)* SQRT((SD_ADD(l,m,3,i,j,k)/ADD(l,m,3,i,j,k))**2+(SD_SF3/SF3)**2+(SD_BAF(iii,9)/BAF(iii,9))**2)  
	  ELSE
	  CR(l,3,i,j,k)=0.0
	  SD_CR(l,3,i,j,k)=0.0
	  IF((SF3.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF      
	  ENDIF
!
	  IF((BAF(iii,9).LE.0.0).AND.(KSTOPBAF(9).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(9)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The beef ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(3)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway3 was not calculed, because soil, water (OTHER_WATERS) and meat concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY3 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'
      PRINT*, '  AND MEAT CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
!
      IF(W(4).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,7).EQV..TRUE.)) THEN
!
      EF4=EF(4)
	  SD_EF4=SD_EF(4)
!
	  ED4=DELTAT(J)
	  SD_ED4=0.0
	  SD_ED_EXP(4,i)=SD_ED4
!
      IF(m.EQ.1)THEN
	  AT4=AT(m,4)*DELTAT(J)
	  ELSE
	  AT4=AT(m,4)
	  ENDIF
!
	  SD_AT4=SD_AT(m,4)
      R_IRleite4=IR(6)
	  SD_IRleite4=SD_IR(6)
	  R_IRali4=IR(7)
	  SD_IRali4=SD_IR(7)
	  R_IRagua4=IR(12)
	  SD_IRagua4=SD_IR(12)
      R_IRsolo4=IR(8)
	  SD_IRsolo4=SD_IR(8)
!
	  BTF_SV4=BTF(iii,2)
	  SD_BTF_SV4=SD_BTF(iii,2)
	  BTF_VM4=BTF(iii,5)
	  SD_BTF_VM4=SD_BTF(iii,5)
	  BTF_SM4=BTF(iii,6)
	  SD_BTF_SM4=SD_BTF(iii,6)
	  BTF_WM4=BTF(iii,9)
	  SD_BTF_WM4=SD_BTF(iii,9)
	  KEYCONC_M4=KEYCONC(i,7)
	  KEYCONC_S4=KEYCONC(i,1)
	  KEYCONC_W4=KEYCONC(i,14)
	  Fa4=Fa(2)
	  SD_Fa4=SD_Fa(2)
	  Fp4=Fp(2)
	  SD_Fp4=SD_Fp(2)
	  FI4=FI(4)
	  SD_FI4=SD_FI(4)
	  fw4=fw(iii)
	  SD_fw4=SD_fw(iii)
!
!
!
      CALL DOSEWAY4(CHEMICAL,KSTOPBTF_SM4,KSTOPBTF_SV4,KSTOPBTF_VM4,KSTOPBTF_WM4,KSTOPfw4,NCHEM,NTIME,NLOCAL,BTF_SV4,SD_BTF_SV4,BTF_VM4,SD_BTF_VM4,BTF_SM4,SD_BTF_SM4,BTF_WM4,SD_BTF_WM4,i,j,k,&
	  CSOIL,SD_CSOIL,CWATEROTHER,SD_CWATEROTHER,CMILK,SD_CMILK,EF4,SD_EF4,ED4,SD_ED4,AT4,SD_AT4,BW,SD_BW,KEYCONC_M4,KEYCONC_S4,KEYCONC_W4,&
	  Fa4,SD_Fa4,Fp4,SD_Fp4,R_IRali4,SD_IRali4,R_IRsolo4,SD_IRsolo4,R_IRagua4,SD_IRagua4,R_IRleite4,SD_IRleite4,&
	  FI4,SD_FI4,fw4,SD_fw4,AD,SD_AD)
!
!
!
      ADD(l,m,4,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,4,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(4).EQ.'Chronic')THEN
	  RfD4=RfD(iii,1,1)
	  SD_RfD4=SD_RfD(iii,1,1)
	  SF4=SF(iii,1,1)
	  SD_SF4=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(4).EQ.'Subchronic')THEN
	  RfD4=RfD(iii,1,2)
	  SD_RfD4=SD_RfD(iii,1,2)
	  SF4=SF(iii,1,2)
	  SD_SF4=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(4).EQ.'Acute')THEN
	  RfD4=RfD(iii,1,3)
	  SD_RfD4=SD_RfD(iii,1,3)
	  SF4=SF(iii,1,3)
	  SD_SF4=SD_SF(iii,1,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,4,i,j,k).GT.0.0).AND.(RfD4.GT.0.0).AND.(BAF(iii,10).GT.0.0))THEN
      HQ(l,4,i,j,k)=ADD(l,m,4,i,j,k)*BAF(iii,10)/RfD4   ! o RfD 4 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,4,i,j,k)=HQ(l,4,i,j,k)* SQRT((SD_ADD(l,m,4,i,j,k)/ADD(l,m,4,i,j,k))**2+(SD_RfD4/RfD4)**2+(SD_BAF(iii,10)/BAF(iii,10))**2)
	  ELSE
	  HQ(l,4,i,j,k)=0.0
	  SD_HQ(l,4,i,j,k)=0.0
	  IF((RfD4.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,4,i,j,k).GT.0.0).AND.(SF4.GT.0.0).AND.(BAF(iii,10).GT.0.0))THEN
	  CR(l,4,i,j,k)=ADD(l,m,4,i,j,k)*BAF(iii,10)*SF4*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,4,i,j,k)=CR(l,4,i,j,k)* SQRT((SD_ADD(l,m,4,i,j,k)/ADD(l,m,4,i,j,k))**2+(SD_SF4/SF4)**2+(SD_BAF(iii,10)/BAF(iii,10))**2)
	  ELSE
	  CR(l,4,i,j,k)=0.0
	  SD_CR(l,4,i,j,k)=0.0
	  IF((SF4.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,10).LE.0.0).AND.(KSTOPBAF(10).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(10)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The milk ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(4)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway4 was not calculed, because soil, water (OTHER_WATERS) and milk concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY4 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'
      PRINT*, '  AND MILK CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
      IF(W(7).EQV..TRUE.)THEN
	  IF((KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,10).EQV..TRUE.)) THEN
!
      EF7=EF(7)
	  SD_EF7=SD_EF(7)
!
	  ED7=DELTAT(J)
	  SD_ED7=0.0
	  SD_ED_EXP(7,i)=SD_ED7
!
      IF(m.EQ.1)THEN
	  AT7=AT(m,7)*DELTAT(J)
	  ELSE
	  AT7=AT(m,7)
	  ENDIF
!
	  SD_AT7=SD_AT(m,7)
      R_IRfish7=IR(13)		 !valor da MPA (2013) foi de 12 kg por ano, transformando-se isso em kg por dia ou kg por refeição da 0.03333!
      SD_IRfish7=SD_IR(13)
!
	  BTF_WF7=BTF(iii,10)
	  SD_BTF_WF7=SD_BTF(iii,10)
	  KEYCONC_F7=KEYCONC(i,10)
	  FI7=FI(6)
	  SD_FI7=SD_FI(6)
!
!
      CWATER_7=CWATEROTHER(i,j,k)
	  SD_CWATER_7=SD_CWATEROTHER(i,j,k)
!
      CALL DOSEWAY7(CHEMICAL,KSTOPBTF_WF7,NCHEM,NTIME,NLOCAL,BTF_WF7,SD_BTF_WF7,i,j,k,CFISH,SD_CFISH,EF7,SD_EF7,ED7,SD_ED7,&
	  AT7,SD_AT7,BW,SD_BW,KEYCONC_F7,R_IRfish7,SD_IRfish7,FI7,SD_FI7,CWATER_7,SD_CWATER_7,AD,SD_AD)
!
!
!
      ADD(l,m,7,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,7,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(7).EQ.'Chronic')THEN
	  RfD7=RfD(iii,1,1)
	  SD_RfD7=SD_RfD(iii,1,1)
	  SF7=SF(iii,1,1)
	  SD_SF7=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(7).EQ.'Subchronic')THEN
	  RfD7=RfD(iii,1,2)
	  SD_RfD7=SD_RfD(iii,1,2)
	  SF7=SF(iii,1,2)
	  SD_SF7=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(7).EQ.'Acute')THEN
	  RfD7=RfD(iii,1,3)
	  SD_RfD7=SD_RfD(iii,1,3)
	  SF7=SF(iii,1,3)
	  SD_SF7=SD_SF(iii,1,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,7,i,j,k).GT.0.0).AND.(RfD7.GT.0.0).AND.(BAF(iii,13).GT.0.0))THEN
      HQ(l,7,i,j,k)=ADD(l,m,7,i,j,k)*BAF(iii,13)/RfD7   ! o RfD 7 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,7,i,j,k)=HQ(l,7,i,j,k)* SQRT((SD_ADD(l,m,7,i,j,k)/ADD(l,m,7,i,j,k))**2+(SD_RfD7/RfD7)**2+(SD_BAF(iii,13)/BAF(iii,13))**2)
	  ELSE
	  HQ(l,7,i,j,k)=0.0
	  SD_HQ(l,7,i,j,k)=0.0
	  IF((RfD7.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,7,i,j,k).GT.0.0).AND.(SF7.GT.0.0).AND.(BAF(iii,13).GT.0.0))THEN
	  CR(l,7,i,j,k)=ADD(l,m,7,i,j,k)*BAF(iii,13)*SF7*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,7,i,j,k)=CR(l,7,i,j,k)* SQRT((SD_ADD(l,m,7,i,j,k)/ADD(l,m,7,i,j,k))**2+(SD_SF7/SF7)**2+(SD_BAF(iii,13)/BAF(iii,13))**2)
	  ELSE
	  CR(l,7,i,j,k)=0.0
	  SD_CR(l,7,i,j,k)=0.0
	  IF((SF7.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
      IF((BAF(iii,13).LE.0.0).AND.(KSTOPBAF(13).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(13)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The fish ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
!
      
	  IF((KSTOPWAY(7).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(7)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway7 was not calculed, because water (OTHER_WATERS) and fish concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY7 WAS NOT CALCULED, BECAUSE WATER (OTHER_WATERS) AND FISH CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!
!*********************************************************************************************************************
!
!
      IF(W(8).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,8).EQV..TRUE.)) THEN
!
      EF8=EF(8)
	  SD_EF8=SD_EF(8)
!
	  ED8=DELTAT(J)
	  SD_ED8=0.0
	  SD_ED_EXP(8,i)=SD_ED8
!
      IF(m.EQ.1)THEN
	  AT8=AT(m,8)*DELTAT(J)
	  ELSE
	  AT8=AT(m,8)
	  ENDIF
!
	  SD_AT8=SD_AT(m,8)
      R_IRave8=IR(14)
	  SD_IRave8=SD_IR(14)
	  R_IRagua8=IR(15)
	  SD_IRagua8=SD_IR(15)
      R_IRsolo8=IR(16)
	  SD_IRsolo8=SD_IR(16)
!
!
	  BTF_SAVE8=BTF(iii,11)
	  SD_BTF_SAVE8=SD_BTF(iii,11)
	  BTF_WAVE8=BTF(iii,12)
	  SD_BTF_WAVE8=SD_BTF(iii,12)
	  KEYCONC_AVE8=KEYCONC(i,8)
	  Fa8=Fa(3)
	  SD_Fa8=SD_Fa(3)
	  Fp8=Fp(3)
	  SD_Fp8=SD_Fp(3)
	  FI8=FI(7)
	  SD_FI8=SD_FI(7)
	  fw8=fw(iii)
	  SD_fw8=SD_fw(iii)
!
!
!
      CALL DOSEWAY8(CHEMICAL,KSTOPBTF_WAVE8,KSTOPBTF_SAVE8,KSTOPfw8,NCHEM,NTIME,NLOCAL,BTF_SAVE8,SD_BTF_SAVE8,BTF_WAVE8,SD_BTF_WAVE8,i,j,k,CSOIL,SD_CSOIL,CWATEROTHER,SD_CWATEROTHER,&
	  CAVE,SD_CAVE,EF8,SD_EF8,ED8,SD_ED8,AT8,SD_AT8,BW,SD_BW,KEYCONC_AVE8,Fa8,SD_Fa8,Fp8,SD_Fp8,R_IRsolo8,SD_IRsolo8,&
	  R_IRagua8,SD_IRagua8,R_IRave8,SD_IRave8,FI8,SD_FI8,fw8,SD_fw8,AD,SD_AD)
!
!

!
      ADD(l,m,8,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,8,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(8).EQ.'Chronic')THEN
	  RfD8=RfD(iii,1,1)
	  SD_RfD8=SD_RfD(iii,1,1)
	  SF8=SF(iii,1,1)
	  SD_SF8=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(8).EQ.'Subchronic')THEN
	  RfD8=RfD(iii,1,2)
	  SD_RfD8=SD_RfD(iii,1,2)
	  SF8=SF(iii,1,2)
	  SD_SF8=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(8).EQ.'Acute')THEN
	  RfD8=RfD(iii,1,3)
	  SD_RfD8=SD_RfD(iii,1,3)
	  SF8=SF(iii,1,3)
	  SD_SF8=SD_SF(iii,1,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,8,i,j,k).GT.0.0).AND.(RfD8.GT.0.0).AND.(BAF(iii,11).GT.0.0))THEN
      HQ(l,8,i,j,k)=ADD(l,m,8,i,j,k)*BAF(iii,11)/RfD8   ! o RfD 8 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,8,i,j,k)=HQ(l,8,i,j,k)* SQRT((SD_ADD(l,m,8,i,j,k)/ADD(l,m,8,i,j,k))**2+(SD_RfD8/RfD8)**2+(SD_BAF(iii,11)/BAF(iii,11))**2)
	  ELSE
	  HQ(l,8,i,j,k)=0.0
	  SD_HQ(l,8,i,j,k)=0.0
	  IF((RfD8.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,8,i,j,k).GT.0.0).AND.(SF8.GT.0.0).AND.(BAF(iii,11).GT.0.0))THEN
	  CR(l,8,i,j,k)=ADD(l,m,8,i,j,k)*BAF(iii,11)*SF8*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,8,i,j,k)=CR(l,8,i,j,k)* SQRT((SD_ADD(l,m,8,i,j,k)/ADD(l,m,8,i,j,k))**2+(SD_SF8/SF8)**2+(SD_BAF(iii,11)/BAF(iii,11))**2)
	  ELSE
	  CR(l,8,i,j,k)=0.0
	  SD_CR(l,8,i,j,k)=0.0
	  IF((SF8.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,11).LE.0.0).AND.(KSTOPBAF(11).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(11)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The ave ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(8).EQ.0).AND.(l.EQ.1))THEN
      KSTOPWAY(8)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway8 was not calculed, because soil, water (OTHER_WATERS) and bird concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY8 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'		 
      PRINT*, '  AND BIRD CONCENTRATIONS NOT EXIST'		 
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!*********************************************************************************************************************
!
      IF(W(9).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,9).EQV..TRUE.)) THEN
!
      EF9=EF(9)
	  SD_EF9=SD_EF(9)
!
	  ED9=DELTAT(J)
	  SD_ED9=0.0
	  SD_ED_EXP(9,i)=SD_ED9
!
      IF(m.EQ.1)THEN
	  AT9=AT(m,9)*DELTAT(J)
	  ELSE
	  AT9=AT(m,9)
	  ENDIF
!
	  SD_AT9=SD_AT(m,9)
      R_IRegg9=IR(17)
	  SD_IRegg9=SD_IR(17)
	  R_IRsolo9=IR(19)
	  SD_IRsolo9=SD_IR(19)
	  R_IRagua9=IR(18)
	  SD_IRagua9=SD_IR(18)
!
!
	  BTF_SEGG9=BTF(iii,13)
	  SD_BTF_SEGG9=SD_BTF(iii,13)
	  BTF_WEGG9=BTF(iii,14)
	  SD_BTF_WEGG9=SD_BTF(iii,14)
	  KEYCONC_EGG9=KEYCONC(i,9)
	  Fa9=Fa(4)
	  SD_Fa9=SD_Fa(4)
	  Fp9=Fp(4)
	  SD_Fp9=SD_Fp(4)
	  FI9=FI(8)
	  SD_FI9=SD_FI(8)
	  fw9=fw(iii)
	  SD_fw9=SD_fw(iii)
!
!
!
      CALL DOSEWAY9(CHEMICAL,KSTOPBTF_WEGG9,KSTOPBTF_SEGG9,KSTOPfw9,NCHEM,NTIME,NLOCAL,BTF_SEGG9,SD_BTF_SEGG9,BTF_WEGG9,SD_BTF_WEGG9,i,j,k,CSOIL,SD_CSOIL,CWATEROTHER,SD_CWATEROTHER,&
	  CEGG,SD_CEGG,EF9,SD_EF9,ED9,SD_ED9,AT9,SD_AT9,BW,SD_BW,KEYCONC_EGG9,Fa9,SD_Fa9,Fp9,SD_Fp9,R_IRsolo9,SD_IRsolo9,&
	  R_IRagua9,SD_IRagua9,R_IRegg9,SD_IRegg9,FI9,SD_FI9,fw9,SD_fw9,AD,SD_AD)
!
!
!
      ADD(l,m,9,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,9,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(9).EQ.'Chronic')THEN
	  RfD9=RfD(iii,1,1)
	  SD_RfD9=SD_RfD(iii,1,1)
	  SF9=SF(iii,1,1)
	  SD_SF9=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(9).EQ.'Subchronic')THEN
	  RfD9=RfD(iii,1,2)
	  SD_RfD9=SD_RfD(iii,1,2)
	  SF9=SF(iii,1,2)
	  SD_SF9=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(9).EQ.'Acute')THEN
	  RfD9=RfD(iii,1,3)
	  SD_RfD9=SD_RfD(iii,1,3)
	  SF9=SF(iii,1,3)
	  SD_SF9=SD_SF(iii,1,3)
	  ENDIF 
!
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,9,i,j,k).GT.0.0).AND.(RfD9.GT.0.0).AND.(BAF(iii,12).GT.0.0))THEN
      HQ(l,9,i,j,k)=ADD(l,m,9,i,j,k)*BAF(iii,12)/RfD9   ! o RfD 9 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,9,i,j,k)=HQ(l,9,i,j,k)* SQRT((SD_ADD(l,m,9,i,j,k)/ADD(l,m,9,i,j,k))**2+(SD_RfD9/RfD9)**2+(SD_BAF(iii,12)/BAF(iii,12))**2)	
	  ELSE
	  HQ(l,9,i,j,k)=0.0
	  SD_HQ(l,9,i,j,k)=0.0
	  IF((RfD9.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,9,i,j,k).GT.0.0).AND.(SF9.GT.0.0).AND.(BAF(iii,12).GT.0.0))THEN
	  CR(l,9,i,j,k)=ADD(l,m,9,i,j,k)*BAF(iii,12)*SF9*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,9,i,j,k)=CR(l,9,i,j,k)* SQRT((SD_ADD(l,m,9,i,j,k)/ADD(l,m,9,i,j,k))**2+(SD_SF9/SF9)**2+(SD_BAF(iii,12)/BAF(iii,12))**2)
	  ELSE
	  CR(l,9,i,j,k)=0.0
	  SD_CR(l,9,i,j,k)=0.0
	  IF((SF9.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,12).LE.0.0).AND.(KSTOPBAF(12).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(12)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The egg ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(9).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(9)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway9 was not calculed, because soil, water (OTHER_WATERS) and egg concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY9 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'
      PRINT*, '  AND EGG CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!*************************************************************************************************
!
!
      IF(W(10).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,11).EQV..TRUE.)) THEN
!
      EF10=EF(10)
	  SD_EF10=SD_EF(10)
!
	  ED10=DELTAT(J)
	  SD_ED10=0.0
	  SD_ED_EXP(10,i)=SD_ED10
!
      IF(m.EQ.1)THEN
	  AT10=AT(m,10)*DELTAT(J)
	  ELSE
	  AT10=AT(m,10)
	  ENDIF
!
	  SD_AT10=SD_AT(m,10)
      R_IRgrain10=IR(20)
      SD_IRgrain10=SD_IR(20)
!
!
	  BTF_SG10=BTF(iii,15)
	  SD_BTF_SG10=SD_BTF(iii,15)
	  KEYCONC_G10=KEYCONC(i,11)
	  FI10=FI(9)
	  SD_FI10=SD_FI(9)
!				
!
!
      CALL DOSEWAY10(CHEMICAL,KSTOPBTF_SG10,NCHEM,NTIME,NLOCAL,BTF_SG10,SD_BTF_SG10,i,j,k,CSOIL,SD_CSOIL,&
	  CGRAIN,SD_CGRAIN,EF10,SD_EF10,ED10,SD_ED10,AT10,SD_AT10,BW,SD_BW,KEYCONC_G10,R_IRgrain10,SD_IRgrain10,FI10,SD_FI10,AD,SD_AD)
!
!
!
      ADD(l,m,10,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,10,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(10).EQ.'Chronic')THEN
	  RfD10=RfD(iii,1,1)
	  SD_RfD10=SD_RfD(iii,1,1)
	  SF10=SF(iii,1,1)
	  SD_SF10=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(10).EQ.'Subchronic')THEN
	  RfD10=RfD(iii,1,2)
	  SD_RfD10=SD_RfD(iii,1,2)
	  SF10=SF(iii,1,2)
	  SD_SF10=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(10).EQ.'Acute')THEN
	  RfD10=RfD(iii,1,3)
	  SD_RfD10=SD_RfD(iii,1,3)
	  SF10=SF(iii,1,3)
	  SD_SF10=SD_SF(iii,1,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,10,i,j,k).GT.0.0).AND.(RfD10.GT.0.0).AND.(BAF(iii,14).GT.0.0))THEN
      HQ(l,10,i,j,k)=ADD(l,m,10,i,j,k)*BAF(iii,14)/RfD10   ! o RfD 10 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,10,i,j,k)=HQ(l,10,i,j,k)* SQRT((SD_ADD(l,m,10,i,j,k)/ADD(l,m,10,i,j,k))**2+(SD_RfD10/RfD10)**2+(SD_BAF(iii,14)/BAF(iii,14))**2)
	  ELSE
	  HQ(l,10,i,j,k)=0.0
	  SD_HQ(l,10,i,j,k)=0.0
	  IF((RfD10.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,10,i,j,k).GT.0.0).AND.(SF10.GT.0.0).AND.(BAF(iii,14).GT.0.0))THEN
	  CR(l,10,i,j,k)=ADD(l,m,10,i,j,k)*BAF(iii,14)*SF10*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,10,i,j,k)=CR(l,10,i,j,k)* SQRT((SD_ADD(l,m,10,i,j,k)/ADD(l,m,10,i,j,k))**2+(SD_SF10/SF10)**2+(SD_BAF(iii,14)/BAF(iii,14))**2)
	  ELSE 
	  CR(l,10,i,j,k)=0.0
	  SD_CR(l,10,i,j,k)=0.0
	  IF((SF10.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,14).LE.0.0).AND.(KSTOPBAF(14).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(14)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The grain ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(10).EQ.0).AND.(l.EQ.1))THEN
      KSTOPWAY(10)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway10 was not calculed, because soil, water (OTHER_WATERS) and grain concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY10 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'
      PRINT*, '  AND GRAIN CONCENTRATION NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
!
      IF(W(11).EQV..TRUE.)THEN
	  IF(KEYCONC(i,3).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF11=EF(11)
	  SD_EF11=SD_EF(11)
!
	  ED11=DELTAT(J)
	  SD_ED11=0.0
	  SD_ED_EXP(11,i)=SD_ED11
!
      IF(m.EQ.1)THEN
	  AT11=AT(m,11)*DELTAT(J)
	  ELSE
	  AT11=AT(m,11)
	  ENDIF
!
	  SD_AT11=SD_AT(m,11)
!
	  ET_AGR_11=ET(1)
	  SD_ET_AGR_11=SD_ET(1)
!
!
!
!                   
!
      CALL DOSEWAY11(NCHEM,NTIME,NLOCAL,i,j,k,EF11,SD_EF11,ED11,SD_ED11,AT11,SD_AT11,CPAR,SD_CPAR,&
	  ET_AGR_11,SD_ET_AGR_11,AD,SD_AD)
!
!
      ADD(l,m,11,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,11,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(11).EQ.'Chronic')THEN
	  RfD11=RfD(iii,2,1)
	  SD_RfD11=SD_RfD(iii,2,1)
	  SF11=SF(iii,2,1)
	  SD_SF11=SD_SF(iii,2,1)
      ELSEIF(NRISKTYPE(11).EQ.'Subchronic')THEN
	  RfD11=RfD(iii,2,2)
	  SD_RfD11=SD_RfD(iii,2,2)
	  SF11=SF(iii,2,2)
	  SD_SF11=SD_SF(iii,2,2)
      ELSEIF(NRISKTYPE(11).EQ.'Acute')THEN
	  RfD11=RfD(iii,2,3)
	  SD_RfD11=SD_RfD(iii,2,3)
	  SF11=SF(iii,2,3)
	  SD_SF11=SD_SF(iii,2,3)
	  ENDIF
!
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,11,i,j,k).GT.0.00).AND.(RfD11.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
      HQ(l,11,i,j,k)=ADD(l,m,11,i,j,k)*BAF(iii,3)/RfD11   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,11,i,j,k)=HQ(l,11,i,j,k)* SQRT((SD_ADD(l,m,11,i,j,k)/ADD(l,m,11,i,j,k))**2+(SD_RfD11/RfD11)**2+(SD_BAF(iii,3)/BAF(iii,3))**2)
	  ELSE
	  HQ(l,11,i,j,k)=0.0
	  SD_HQ(l,11,i,j,k)=0.0
	  IF((RfD11.LE.0.0).AND.(KSTOPRFD(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particulate matter RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,11,i,j,k).GT.0.0).AND.(SF11.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
	  CR(l,11,i,j,k)=ADD(l,m,11,i,j,k)*BAF(iii,3)*SF11*ADAF(2)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,11,i,j,k)=CR(l,11,i,j,k)* SQRT((SD_ADD(l,m,11,i,j,k)/ADD(l,m,11,i,j,k))**2+(SD_SF11/SF11)**2+(SD_BAF(iii,3)/BAF(iii,3))**2)	
	  ELSE
	  CR(l,11,i,j,k)=0.0
	  SD_CR(l,11,i,j,k)=0.0
	  IF((SF11.LE.0.0).AND.(KSTOPSF(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particulate matter SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,3).LE.0.0).AND.(KSTOPBAF(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particulate matter BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(11).EQ.0).AND.(l.EQ.1))THEN
      KSTOPWAY(11)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway11 was not calculed, because particulate concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY11 WAS NOT CALCULED, BECAUSE PARTICULATE CONCENTRATIONS NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
!
      IF(W(12).EQV..TRUE.)THEN
	  IF(KEYCONC(i,12).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF12=EF(12)
	  SD_EF12=SD_EF(12)
!					 
	  ED12=DELTAT(J)
	  SD_ED12=0.0
	  SD_ED_EXP(12,i)=SD_ED12
!
      IF(m.EQ.1)THEN
	  AT12=AT(m,12)*DELTAT(J)
	  ELSE
	  AT12=AT(m,12)
	  ENDIF
!
	  SD_AT12=SD_AT(m,12)
!
	  ET_AGR_12=ET(2)
	  SD_ET_AGR_12=SD_ET(2)
!
!                   
!
      CALL DOSEWAY12(NCHEM,NTIME,NLOCAL,i,j,k,EF12,SD_EF12,ED12,SD_ED12,AT12,SD_AT12,CSTEAM,SD_CSTEAM,&
	  ET_AGR_12,SD_ET_AGR_12,AD,SD_AD)
!
!
      ADD(l,m,12,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,12,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(12).EQ.'Chronic')THEN
	  RfD12=RfD(iii,5,1)
	  SD_RfD12=SD_RfD(iii,5,1)
	  SF12=SF(iii,5,1)
	  SD_SF12=SD_SF(iii,5,1)
      ELSEIF(NRISKTYPE(12).EQ.'Subchronic')THEN
	  RfD12=RfD(iii,5,2)
	  SD_RfD12=SD_RfD(iii,5,2)
	  SF12=SF(iii,5,2)
	  SD_SF12=SD_SF(iii,5,2)
      ELSEIF(NRISKTYPE(12).EQ.'Acute')THEN
	  RfD12=RfD(iii,5,3)
	  SD_RfD12=SD_RfD(iii,5,3)
	  SF12=SF(iii,5,3)
	  SD_SF12=SD_SF(iii,5,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,12,i,j,k).GT.0.00).AND.(RfD12.GT.0.0).AND.(BAF(iii,4).GT.0.0))THEN
      HQ(l,12,i,j,k)=ADD(l,m,12,i,j,k)*BAF(iii,4)/RfD12   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,12,i,j,k)=HQ(l,12,i,j,k)* SQRT((SD_ADD(l,m,12,i,j,k)/ADD(l,m,12,i,j,k))**2+(SD_RfD12/RfD12)**2+(SD_BAF(iii,4)/BAF(iii,4))**2)
	  ELSE
	  HQ(l,12,i,j,k)=0.0
	  SD_HQ(l,12,i,j,k)=0.0
	  IF((RfD12.LE.0.0).AND.(KSTOPRFD(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,12,i,j,k).GT.0.0).AND.(SF12.GT.0.0).AND.(BAF(iii,4).GT.0.0))THEN
	  CR(l,12,i,j,k)=ADD(l,m,12,i,j,k)*BAF(iii,4)*SF12*ADAF(2)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,12,i,j,k)=CR(l,12,i,j,k)* SQRT((SD_ADD(l,m,12,i,j,k)/ADD(l,m,12,i,j,k))**2+(SD_SF12/SF12)**2+(SD_BAF(iii,4)/BAF(iii,4))**2)	
	  ELSE
	  CR(l,12,i,j,k)=0.0
	  SD_CR(l,12,i,j,k)=0.0
	  IF((SF12.LE.0.0).AND.(KSTOPSF(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,4).LE.0.0).AND.(KSTOPBAF(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(12).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(12)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway12 was not calculed, because steam concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY12 WAS NOT CALCULED, BECAUSE THE STEAM CONCENTRATIONS NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!
!*******************************************************************************************
!*******************************************************************************************
!
!
      IF(W(13).EQV..TRUE.)THEN
	  IF(KEYCONC(i,13).EQV..TRUE.)THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
	  EF13=EF(13)
	  SD_EF13=SD_EF(13)
!
	  ED13=DELTAT(J)
	  SD_ED13=0.0
	  SD_ED_EXP(13,i)=SD_ED13
!
      IF(m.EQ.1)THEN
	  AT13=AT(m,13)*DELTAT(J)
	  ELSE
	  AT13=AT(m,13)
	  ENDIF
!
	  SD_AT13=SD_AT(m,13)
	  ET13=ET(3)
	  SD_ET13=SD_ET(3)
	  SA13=SA(1)
	  SD_SA13=SD_SA(1)
!
	  PC13=PC(iii)
	  SD_PC13=SD_PC(iii)
	  CF13=CF(2)
	  SD_CF13=SD_CF(2)
!
	  EV13=EV(1)
	  SD_EV13=SD_EV(1)
! 
!
!
      CALL DOSEWAY13(NCHEM,NTIME,NLOCAL,KSTOPPC13,CHEMICAL,i,j,k,CWATERDER,SD_CWATERDER,BW,SD_BW,EV13,SD_EV13,EF13,SD_EF13,ED13,SD_ED13,ET13,SD_ET13,AT13,SD_AT13,&
	  PC13,SD_PC13,SA13,SD_SA13,CF13,SD_CF13,AD,SD_AD)
!
!
      ADD(l,m,13,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,13,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(13).EQ.'Chronic')THEN
	  RfD13=RfD(iii,6,1)
	  SD_RfD13=SD_RfD(iii,6,1)
	  SF13=SF(iii,6,1)
	  SD_SF13=SD_SF(iii,6,1)
      ELSEIF(NRISKTYPE(13).EQ.'Subchronic')THEN
	  RfD13=RfD(iii,6,2)
	  SD_RfD13=SD_RfD(iii,6,2)
	  SF13=SF(iii,6,2)
	  SD_SF13=SD_SF(iii,6,2)
      ELSEIF(NRISKTYPE(13).EQ.'Acute')THEN
	  RfD13=RfD(iii,6,3)
	  SD_RfD13=SD_RfD(iii,6,3)
	  SF13=SF(iii,6,3)
	  SD_SF13=SD_SF(iii,6,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,13,i,j,k).GT.0.0).AND.(RfD13.GT.0.0).AND.(BAF(iii,6).GT.0.0))THEN
      HQ(l,13,i,j,k)=ADD(l,m,13,i,j,k)*BAF(iii,6)/RfD13   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,13,i,j,k)=HQ(l,13,i,j,k)* SQRT((SD_ADD(l,m,13,i,j,k)/ADD(l,m,13,i,j,k))**2+(SD_RfD13/RfD13)**2+(SD_BAF(iii,6)/BAF(iii,6))**2)
	  ELSE
	  HQ(l,13,i,j,k)=0.0
	  SD_HQ(l,13,i,j,k)=0.0
	  IF((RfD13.LE.0.0).AND.(KSTOPRFD(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal water RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,13,i,j,k).GT.0.0).AND.(SF13.GT.0.0).AND.(BAF(iii,6).GT.0.0))THEN
	  CR(l,13,i,j,k)=ADD(l,m,13,i,j,k)*BAF(iii,6)*SF13*ADAF(3)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,13,i,j,k)=CR(l,13,i,j,k)* SQRT((SD_ADD(l,m,13,i,j,k)/ADD(l,m,13,i,j,k))**2+(SD_SF13/SF13)**2+(SD_BAF(iii,6)/BAF(iii,6))**2)
	  ELSE
	  CR(l,13,i,j,k)=0.0
	  SD_CR(l,13,i,j,k)=0.0
	  IF((SF13.LE.0.0).AND.(KSTOPSF(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal water SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,6).LE.0.0).AND.(KSTOPBAF(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal water BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(13).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(13)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway13 was not calculed, because water (BATH_WATER) concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY13 WAS NOT CALCULED, BECAUSE WATER (BATH_WATER) CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF

!
!*******************************************************************************************
!
      IF(W(14).EQV..TRUE.)THEN
	  IF(KEYCONC(i,1).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF14=EF(14)
	  SD_EF14=SD_EF(14)
!
	  ED14=DELTAT(J)
	  SD_ED14=0.0
	  SD_ED_EXP(14,i)=SD_ED14
!
      IF(m.EQ.1)THEN
	  AT14=AT(m,14)*DELTAT(J)
	  ELSE
	  AT14=AT(m,14)
	  ENDIF
!
	  SD_AT14=SD_AT(m,14)
	  SA_AGR14=SA(2)
	  SD_SA_AGR14=SD_SA(2)
	  AFsoil_AGR14=AF
	  SD_AFsoil_AGR14=SD_AF
!
!
	  ABS14=ABS_(iii)
	  SD_ABS14=SD_ABS(iii)
	  CF14=CF(3)
	  SD_CF14=SD_CF(3)
	  EV14=EV(2)
	  SD_EV14=SD_EV(2)
!
!
      CALL DOSEWAY14(NCHEM,NTIME,NLOCAL,KSTOPABS14,CHEMICAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF14,SD_EF14,ED14,SD_ED14,AT14,SD_AT14,ABS14,SD_ABS14,&
	  AFsoil_AGR14,SD_AFsoil_AGR14,SA_AGR14,SD_SA_AGR14,CF14,SD_CF14,EV14,SD_EV14,AD,SD_AD)
!
!
      ADD(l,m,14,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,14,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(14).EQ.'Chronic')THEN
	  RfD14=RfD(iii,3,1)
	  SD_RfD14=SD_RfD(iii,3,1)
	  SF14=SF(iii,3,1)
	  SD_SF14=SD_SF(iii,3,1)
      ELSEIF(NRISKTYPE(14).EQ.'Subchronic')THEN
	  RfD14=RfD(iii,3,2)
	  SD_RfD14=SD_RfD(iii,3,2)
	  SF14=SF(iii,3,2)
	  SD_SF14=SD_SF(iii,3,2)
      ELSEIF(NRISKTYPE(14).EQ.'Acute')THEN
	  RfD14=RfD(iii,3,3)
	  SD_RfD14=SD_RfD(iii,3,3)
	  SF14=SF(iii,3,3)
	  SD_SF14=SD_SF(iii,3,3)
	  ENDIF
!
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,14,i,j,k).GT.0.0).AND.(RfD14.GT.0.0).AND.(BAF(iii,5).GT.0.0))THEN
      HQ(l,14,i,j,k)=ADD(l,m,14,i,j,k)*BAF(iii,5)/RfD14   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,14,i,j,k)=HQ(l,14,i,j,k)* SQRT((SD_ADD(l,m,14,i,j,k)/ADD(l,m,14,i,j,k))**2+(SD_RfD14/RfD14)**2+(SD_BAF(iii,5)/BAF(iii,5))**2)
	  ELSE
	  HQ(l,14,i,j,k)=0.0
	  SD_HQ(l,14,i,j,k)=0.0
	  IF((RfD14.LE.0.0).AND.(KSTOPRFD(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,14,i,j,k).GT.0.0).AND.(SF14.GT.0.0).AND.(BAF(iii,5).GT.0.0))THEN
	  CR(l,14,i,j,k)=ADD(l,m,14,i,j,k)*BAF(iii,5)*SF14*ADAF(3)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,14,i,j,k)=CR(l,14,i,j,k)* SQRT((SD_ADD(l,m,14,i,j,k)/ADD(l,m,14,i,j,k))**2+(SD_SF14/SF14)**2+(SD_BAF(iii,5)/BAF(iii,5))**2)
	  ELSE
	  CR(l,14,i,j,k)=0.0
	  SD_CR(l,14,i,j,k)=0.0
	  IF((SF14.LE.0.0).AND.(KSTOPSF(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF	 ! IF "M"
!
	  IF((BAF(iii,5).LE.0.0).AND.(KSTOPBAF(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(14).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(14)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway14 was not calculed, because soil concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY14 WAS NOT CALCULED, BECAUSE SOIL CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!
	  ENDDO	  ! ACABA O CICLO de "k"
	  ENDDO	  ! ACABA O CICLO de "j"
	  ENDDO	  ! ACABA O CICLO de "m"
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1).AND.(m.EQ.1))THEN
	  WRITE(17,'("________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
      ENDIF
!
53    CONTINUE
!
      ENDDO   ! ACABA O CICLO de "i" --> espécie química
!
      ENDDO   ! ACABA O CICLO DE "l"
!
!
!
!*******/////////////////////////////////*******************////////////////////////////////
!*******/////////////////////////////////*******************////////////////////////////////
!*******/////////////////////////////////*******************////////////////////////////////
!
!					                     CASO 2 (INDUSTRIAL)
!
	 CASE(2)
!    CHAMADAS A WAY DO CENARIO=2
!	   
!
!*******************************************************************************************
!-------------------------------------------------------------------------------------------
!---------------------------------- CICLOS POR METAIS  -------------------------------------
!
!
      DO l=7,NIDADE
!
!------------------
	  NTIMEexp=ED_INI(2,l)
!------------------
!
!
!
		 LJ=NCHEM
		 MEW_NPOL=NPOL
!
      DO i=1,LJ
!
      NSTOPNVP(i,l)= NSTOPNVP(i-1,l)+1
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1))THEN
	  WRITE(17,*)
	  WRITE(17,*)
	  WRITE(17,'("**********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,'("                                                                             CICLE ",I2)')l
	  WRITE(17,'("**********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,*)
	  WRITE(17,*)
!
	  npk=20
!
	  DO lop=1,20
	  LTSP(lop)=lop
	  ENDDO
!
	  WRITE(17,'("__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
	  WRITE(17,'("   TIME      ATc     SD_ATc     BW   SD_BW    ET1    SD_ET1     ET2    SD_ET2     SA2    SD_SA2     EV    SD_EV      AF      SD_AF     IR                                                                                                                                                                                 SD_IR ")')
	  WRITE(17,'(127X,<npk>I9,<npk>I9)')(LTSP(kJy),kJy=1,npk),(LTSP(kJy),kJy=1,npk)
	  WRITE(17,'("__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
	  WRITE(17,*)
!
	  ENDIF
!
      KSTOPBTF1=0
	  KSTOPABS14=0
!
      DO LOL=1,6
      KSTOPRFD(LOL)=0
	  KSTOPSF(LOL)=0
	  ENDDO
!
      DO KO=1,NVIAS
	  KSTOPWAY(KO)=0
	  ENDDO
!
	  DO jk=1,14
      KSTOPBAF(jk)=0
	  ENDDO
!
      DO m=1,2
!
!
!---------------------------------------------------------------------------------------------------------------------
!
!-----------------
!	  AT(m,2)=AT_INI(m,2,l)
!	  SD_AT(m,2)=SD_AT_INI(m,2,l)
!-----------------
!
!
      iii=0
!
      DO ii=1,MEW_NPOL
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
	  IF (iii.eq.0) THEN
!
	  IF (l.eq.7) THEN
!
      WRITE(*,*)
	  WRITE(*,'(" WARNING!!! CHEMICAL NOT EXIST IN Datachemical DATABASE  --->  ",A30)') CHEMICAL(i)  
	  WRITE(*,'(" The HQ and CR values will all be zero for this Chemical species!  ")')
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! CHEMICAL NOT EXIST IN Datachemical DATABASE  --->  ",A30)') CHEMICAL(i)   
	  WRITE(99,'(" The HQ and CR values will all be zero for this Chemical species!  ")')
      WRITE(99,*)
!
      ENDIF
!
	  GOTO 54
	  ENDIF

!
	  DO j=1,NTIMEexp
!
	  DELTAT(J)=1.0
!
      CALL REDISTRI(iii,l,j,NVIAS,NIDADE,EF_INI,SD_EF_INI,BW_INI,SD_BW_INI,IR_INI,SD_IR_INI,ET_INI,SD_ET_INI,SA_INI,SD_SA_INI,&
	  EV_INI,SD_EV_INI,AF_INI,SD_AF_INI,EF,SD_EF,BW,SD_BW,IR,SD_IR,ET,SD_ET,SA,SD_SA,EV,SD_EV,AF,SD_AF,MUTAGENIC,ADAF,&
	  AT_INI,SD_AT_INI,AT,SD_AT)
!
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1).AND.(m.EQ.1))THEN
!
      write(17,'( 3x,I3,3x,F9.1,1x,F9.1,3x,F4.1,2x,F4.1,2x,E8.3,3(1x,E8.3),2(1x,E8.3),4x,F3.1,5x,F3.1,3x,E8.3,1x,E8.3,20(1x,E8.3),20(1x,E8.3))') j,AT(2,2),SD_AT(2,2),BW,SD_BW,ET(1),SD_ET(1),ET(2),SD_ET(2),SA(2),SD_SA(2),EV(2),SD_EV(2),AF,SD_AF,(IR(iou), iou=1,20),(SD_IR(IQW), IQW=1,20)
!
      ENDIF
!
!
	  DO k=1,NLOCAL
!
!
!*******************************************************************************************
      IF(W(1).EQV..TRUE.)THEN
	  IF(KEYCONC(i,1).EQV..TRUE.) THEN
!
      EF1=EF(1)
	  SD_EF1=SD_EF(1)
!
	  ED1=DELTAT(J)
	  SD_ED1=0.0
	  SD_ED_EXP(1,i)=SD_ED1
!
      IF(m.EQ.1)THEN
	  AT1=AT(m,1)*DELTAT(J)
	  ELSE
	  AT1=AT(m,1)
	  ENDIF
!
!
	  SD_AT1=SD_AT(m,1)
	  R_IR1_5=IR(1)
	  SD_IR1_5=SD_IR(1)
!
	  FI1=FI(1)
	  SD_FI1=SD_FI(1)
	  CF1=CF(1)
	  SD_CF1=SD_CF(1)
!
!
      CALL DOSEWAY1(NCHEM,NTIME,NLOCAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF1,SD_EF1,ED1,SD_ED1,AT1,SD_AT1,R_IR1_5,SD_IR1_5,&
	  FI1,SD_FI1,CF1,SD_CF1,AD,SD_AD)
!
!
      ADD(l,m,1,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!
      SD_ADD(l,m,1,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
!  
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(1).EQ.'Chronic')THEN
	  RfD1=RfD(iii,1,1)
	  SD_RfD1=SD_RfD(iii,1,1)
	  SF1=SF(iii,1,1)
	  SD_SF1=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(1).EQ.'Subchronic')THEN
	  RfD1=RfD(iii,1,2)
	  SD_RfD1=SD_RfD(iii,1,2)
	  SF1=SF(iii,1,2)
	  SD_SF1=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(1).EQ.'Acute')THEN
	  RfD1=RfD(iii,1,3)
	  SD_RfD1=SD_RfD(iii,1,3)
	  SF1=SF(iii,1,3)
	  SD_SF1=SD_SF(iii,1,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,1,i,j,k).GT.0.0).AND.(RfD1.GT.0.0).AND.(BAF(iii,1).GT.0.0))THEN
	  HQ(l,1,i,j,k)=ADD(l,m,1,i,j,k)*BAF(iii,1)/RfD1   ! o RfD 1 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,1,i,j,k)=HQ(l,1,i,j,k)* SQRT((SD_ADD(l,m,1,i,j,k)/ADD(l,m,1,i,j,k))**2+(SD_RfD1/RfD1)**2+(SD_BAF(iii,1)/BAF(iii,1))**2)
	  ELSE
	  HQ(l,1,i,j,k)=0.0
	  SD_HQ(l,1,i,j,k)=0.0
	  IF((RfD1.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,1,i,j,k).GT.0.0).AND.(SF1.GT.0.0).AND.(BAF(iii,1).GT.0.0))THEN
	  CR(l,1,i,j,k)=ADD(l,m,1,i,j,k)*BAF(iii,1)*SF1*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,1,i,j,k)=CR(l,1,i,j,k)* SQRT((SD_ADD(l,m,1,i,j,k)/ADD(l,m,1,i,j,k))**2+(SD_SF1/SF1)**2+(SD_BAF(iii,1)/BAF(iii,1))**2)
	  ELSE
	  CR(l,1,i,j,k)=0.0
	  SD_CR(l,1,i,j,k)=0.0
	  IF((SF1.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
	  ENDIF
!
	  IF((BAF(iii,1).LE.0.0).AND.(KSTOPBAF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The soil ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(1).EQ.0).AND.(l.EQ.7))THEN
	  KSTOPWAY(1)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway1 was not calculed, because soil concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY1 WAS NOT CALCULED, BECAUSE SOIL CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
      IF(W(11).EQV..TRUE.)THEN
	  IF(KEYCONC(i,3).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
!
      EF11=EF(11)
	  SD_EF11=SD_EF(11)
!
	  ED11=DELTAT(J)
	  SD_ED11=0.0
	  SD_ED_EXP(11,i)=SD_ED11
!
      IF(m.EQ.1)THEN
	  AT11=AT(m,11)*DELTAT(J)
	  ELSE
	  AT11=AT(m,11)
	  ENDIF
!
	  SD_AT11=SD_AT(m,11)
!
	  ET_IND_11=ET(1)
	  SD_ET_IND_11=SD_ET(1)
!
!
!
      CALL DOSEWAY11(NCHEM,NTIME,NLOCAL,i,j,k,EF11,SD_EF11,ED11,SD_ED11,AT11,SD_AT11,CPAR,SD_CPAR,&
	  ET_IND_11,SD_ET_IND_11,AD,SD_AD)
!
!
      ADD(l,m,11,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,11,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(11).EQ.'Chronic')THEN
	  RfD11=RfD(iii,2,1)
	  SD_RfD11=SD_RfD(iii,2,1)
	  SF11=SF(iii,2,1)
	  SD_SF11=SD_SF(iii,2,1)
      ELSEIF(NRISKTYPE(11).EQ.'Subchronic')THEN
	  RfD11=RfD(iii,2,2)
	  SD_RfD11=SD_RfD(iii,2,2)
	  SF11=SF(iii,2,2)
	  SD_SF11=SD_SF(iii,2,2)
      ELSEIF(NRISKTYPE(11).EQ.'Acute')THEN
	  RfD11=RfD(iii,2,3)
	  SD_RfD11=SD_RfD(iii,2,3)
	  SF11=SF(iii,2,3)
	  SD_SF11=SD_SF(iii,2,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,11,i,j,k).GT.0.0).AND.(RfD11.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
      HQ(l,11,i,j,k)=ADD(l,m,11,i,j,k)*BAF(iii,3)/RfD11   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,11,i,j,k)=HQ(l,11,i,j,k)* SQRT((SD_ADD(l,m,11,i,j,k)/ADD(l,m,11,i,j,k))**2+(SD_RfD11/RfD11)**2+(SD_BAF(iii,3)/BAF(iii,3))**2)
	  ELSE
	  HQ(l,11,i,j,k)=0.0
	  SD_HQ(l,11,i,j,k)=0.0
	  IF((RfD11.LE.0.0).AND.(KSTOPRFD(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particulate matter RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,11,i,j,k).GT.0.0).AND.(SF11.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
	  CR(l,11,i,j,k)=ADD(l,m,11,i,j,k)*BAF(iii,3)*SF11*ADAF(2)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,11,i,j,k)=CR(l,11,i,j,k)* SQRT((SD_ADD(l,m,11,i,j,k)/ADD(l,m,11,i,j,k))**2+(SD_SF11/SF11)**2+(SD_BAF(iii,3)/BAF(iii,3))**2)	
	  ELSE
	  CR(l,11,i,j,k)=0.0
	  SD_CR(l,11,i,j,k)=0.0
	  IF((SF11.LE.0.0).AND.(KSTOPSF(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particulate matter SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,3).LE.0.0).AND.(KSTOPBAF(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particulate matter BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(11).EQ.0).AND.(l.EQ.7))THEN
	  KSTOPWAY(11)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway11 was not calculed, because particulate concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY11 WAS NOT CALCULED, BECAUSE PARTICULATE CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
!*******************************************************************************************
!
      IF(W(12).EQV..TRUE.)THEN
	  IF(KEYCONC(i,12).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
!
      EF12=EF(12)
	  SD_EF12=SD_EF(12)
!
	  ED12=DELTAT(J)
	  SD_ED12=0.0
	  SD_ED_EXP(12,i)=SD_ED12
!
      IF(m.EQ.1)THEN
	  AT12=AT(m,12)*DELTAT(J)
	  ELSE
	  AT12=AT(m,12)
	  ENDIF
!
	  SD_AT12=SD_AT(m,12)
!
	  ET_IND_12=ET(2)
	  SD_ET_IND_12=SD_ET(2)
!
!
!
      CALL DOSEWAY12(NCHEM,NTIME,NLOCAL,i,j,k,EF12,SD_EF12,ED12,SD_ED12,AT12,SD_AT12,CSTEAM,SD_CSTEAM,&
	  ET_IND_12,SD_ET_IND_12,AD,SD_AD)
!
!
      ADD(l,m,12,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,12,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(12).EQ.'Chronic')THEN
	  RfD12=RfD(iii,5,1)
	  SD_RfD12=SD_RfD(iii,5,1)
	  SF12=SF(iii,5,1)
	  SD_SF12=SD_SF(iii,5,1)
      ELSEIF(NRISKTYPE(12).EQ.'Subchronic')THEN
	  RfD12=RfD(iii,5,2)
	  SD_RfD12=SD_RfD(iii,5,2)
	  SF12=SF(iii,5,2)
	  SD_SF12=SD_SF(iii,5,2)
      ELSEIF(NRISKTYPE(12).EQ.'Acute')THEN
	  RfD12=RfD(iii,5,3)
	  SD_RfD12=SD_RfD(iii,5,3)
	  SF12=SF(iii,5,3)
	  SD_SF12=SD_SF(iii,5,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,12,i,j,k).GT.0.00).AND.(RfD12.GT.0.0).AND.(BAF(iii,4).GT.0.0))THEN
      HQ(l,12,i,j,k)=ADD(l,m,12,i,j,k)*BAF(iii,4)/RfD12   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,12,i,j,k)=HQ(l,12,i,j,k)* SQRT((SD_ADD(l,m,12,i,j,k)/ADD(l,m,12,i,j,k))**2+(SD_RfD12/RfD12)**2+(SD_BAF(iii,4)/BAF(iii,4))**2)
	  ELSE
	  HQ(l,12,i,j,k)=0.0
	  SD_HQ(l,12,i,j,k)=0.0
	  IF((RfD12.LE.0.0).AND.(KSTOPRFD(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,12,i,j,k).GT.0.0).AND.(SF12.GT.0.0).AND.(BAF(iii,4).GT.0.0))THEN
	  CR(l,12,i,j,k)=ADD(l,m,12,i,j,k)*BAF(iii,4)*SF12*ADAF(2)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,12,i,j,k)=CR(l,12,i,j,k)* SQRT((SD_ADD(l,m,12,i,j,k)/ADD(l,m,12,i,j,k))**2+(SD_SF12/SF12)**2+(SD_BAF(iii,4)/BAF(iii,4))**2)	
	  ELSE
	  CR(l,12,i,j,k)=0.0
	  SD_CR(l,12,i,j,k)=0.0
	  IF((SF12.LE.0.0).AND.(KSTOPSF(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,4).LE.0.0).AND.(KSTOPBAF(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
	  IF((KSTOPWAY(12).EQ.0).AND.(l.EQ.7))THEN
	  KSTOPWAY(12)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway12 was not calculed, because steam concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY12 WAS NOT CALCULED, BECAUSE STEAM CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
      IF(W(14).EQV..TRUE.)THEN
	  IF(KEYCONC(i,1).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
!
      EF14=EF(14)
	  SD_EF14=SD_EF(14)
!
	  ED14=DELTAT(J)
	  SD_ED14=0.0
	  SD_ED_EXP(14,i)=SD_ED14
!
      IF(m.EQ.1)THEN
	  AT14=AT(m,14)*DELTAT(J)
	  ELSE
	  AT14=AT(m,14)
	  ENDIF
!
	  SD_AT14=SD_AT(m,14)
!
	  ABS14=ABS_(iii)
	  SD_ABS14=SD_ABS(iii)
!
	  SA14=SA(2)
	  SD_SA14=SD_SA(2)
!
	  AFsoil14=AF
	  SD_AFsoil14=SD_AF
!
	  CF14=CF(3)
	  SD_CF14=SD_CF(3)
	  EV14=EV(2)		
	  SD_EV14=SD_EV(2)
!
!
      CALL DOSEWAY14(NCHEM,NTIME,NLOCAL,KSTOPABS14,CHEMICAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF14,SD_EF14,ED14,SD_ED14,AT14,SD_AT14,ABS14,SD_ABS14,&
	  AFsoil14,SD_AFsoil14,SA14,SD_SA14,CF14,SD_CF14,EV14,SD_EV14,AD,SD_AD)
!
!
      ADD(l,m,14,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,14,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(14).EQ.'Chronic')THEN
	  RfD14=RfD(iii,3,1)
	  SD_RfD14=SD_RfD(iii,3,1)
	  SF14=SF(iii,3,1)
	  SD_SF14=SD_SF(iii,3,1)
      ELSEIF(NRISKTYPE(14).EQ.'Subchronic')THEN
	  RfD14=RfD(iii,3,2)
	  SD_RfD14=SD_RfD(iii,3,2)
	  SF14=SF(iii,3,2)
	  SD_SF14=SD_SF(iii,3,2)
      ELSEIF(NRISKTYPE(14).EQ.'Acute')THEN
	  RfD14=RfD(iii,3,3)
	  SD_RfD14=SD_RfD(iii,3,3)
	  SF14=SF(iii,3,3)
	  SD_SF14=SD_SF(iii,3,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,14,i,j,k).GT.0.0).AND.(RfD14.GT.0.0).AND.(BAF(iii,5).GT.0.0))THEN
      HQ(l,14,i,j,k)=ADD(l,m,14,i,j,k)*BAF(iii,5)/RfD14   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,14,i,j,k)=HQ(l,14,i,j,k)* SQRT((SD_ADD(l,m,14,i,j,k)/ADD(l,m,14,i,j,k))**2+(SD_RfD14/RfD14)**2+(SD_BAF(iii,5)/BAF(iii,5))**2)
	  ELSE
	  HQ(l,14,i,j,k)=0.0
	  SD_HQ(l,14,i,j,k)=0.0
	  IF((RfD14.LE.0.0).AND.(KSTOPRFD(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,14,i,j,k).GT.0.0).AND.(SF14.GT.0.0).AND.(BAF(iii,5).GT.0.0))THEN
	  CR(l,14,i,j,k)=ADD(l,m,14,i,j,k)*BAF(iii,5)*SF14*ADAF(3)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,14,i,j,k)=CR(l,14,i,j,k)* SQRT((SD_ADD(l,m,14,i,j,k)/ADD(l,m,14,i,j,k))**2+(SD_SF14/SF14)**2+(SD_BAF(iii,5)/BAF(iii,5))**2)
	  ELSE
	  CR(l,14,i,j,k)=0.0
	  SD_CR(l,14,i,j,k)=0.0
	  IF((SF14.LE.0.0).AND.(KSTOPSF(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF	 ! IF "M"
!
	  IF((BAF(iii,5).LE.0.0).AND.(KSTOPBAF(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(14).EQ.0).AND.(l.EQ.7))THEN
	  KSTOPWAY(14)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway14 was not calculed, because soil concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY14 WAS NOT CALCULED, BECAUSE THE SOIL CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
      ENDDO
	  ENDDO		 
	  ENDDO		 ! FIM DO CICLO DOS METAIS CANCERIGENOS
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1).AND.(m.EQ.1))THEN
	  WRITE(17,'("__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
      ENDIF
!
54    CONTINUE
!
      ENDDO   ! ACABA O CICLO de "i"
!
      ENDDO   ! ACABA O CICLO de "l"
!
!
!*******/////////////////////////////////*******************////////////////////////////////
!*******/////////////////////////////////*******************////////////////////////////////
!*******/////////////////////////////////*******************////////////////////////////////
!
!					                     CASO 3 (RESIDENTIAL)
!
	 CASE(3)
!    CHAMADAS A WAY DO CENARIO=3
!	  
! 
!*******************************************************************************************
!-------------------------------------------------------------------------------------------
!---------------------------------- CICLOS POR METAIS  -------------------------------------
!
!
      DO l=1,NIDADE
!
!------------------
	  NTIMEexp=ED_INI(3,l)
!------------------
!
!
!
		 LJ=NCHEM
		 JEP_NPOL=NPOL
!
      DO i=1,LJ
!
      NSTOPNVP(i,l)= NSTOPNVP(i-1,l)+1
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1))THEN
	  WRITE(17,*)
	  WRITE(17,*)
	  WRITE(17,'("****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,'("                                                                             CICLE ",I2)')l
	  WRITE(17,'("****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************")')
	  WRITE(17,*)
	  WRITE(17,*)
!
	  npk=20
!
	  DO lop=1,20
	  LTSP(lop)=lop
	  ENDDO
!
	  WRITE(17,'("________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
	  WRITE(17,'("   TIME      ATc     SD_ATc     BW   SD_BW    ET1    SD_ET1     ET2    SD_ET2    ET3    SD_ET3     SA1     SD_SA1     SA2    SD_SA2    EV1  SD_EV1   EV2  SD_EV2     AF      SD_AF     IR                                                                                                                                                                                 SD_IR ")')
	  WRITE(17,'(175X,<npk>I9,<npk>I9)')(LTSP(kJy),kJy=1,npk),(LTSP(kJy),kJy=1,npk)
	  WRITE(17,'("________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
	  WRITE(17,*)
!
!
	  ENDIF
!
      KSTOPBTF1=0
	  KSTOPBTF_6S=0
	  KSTOPBTF_SB3=0
	  KSTOPBTF_SV3=0
	  KSTOPBTF_VB3=0
	  KSTOPBTF_WB3=0
	  KSTOPBTF_SM4=0
	  KSTOPBTF_SV4=0
	  KSTOPBTF_VM4=0
	  KSTOPBTF_WM4=0
!
	  KSTOPfw3=0
	  KSTOPfw4=0
	  KSTOPPC13=0
	  KSTOPABS14=0
!
      DO LOL=1,6
      KSTOPRFD(LOL)=0
	  KSTOPSF(LOL)=0
	  ENDDO
!
      DO KO=1,NVIAS
	  KSTOPWAY(KO)=0
	  ENDDO
!
      DO JK=1,14
      KSTOPBAF(JK)=0
	  ENDDO
!
      DO m=1,2
!
!
!----------------------------------------------------------------------------------------------------------------------!
!-----------------
!	  AT(m,3)=AT_INI(m,3,l)
!	  SD_AT(m,3)=SD_AT_INI(m,3,l)
!-----------------
!
!
      iii=0
!
      DO ii=1,JEP_NPOL
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
	  IF (iii.eq.0) THEN
!
	  IF (l.eq.1) THEN
!
      WRITE(*,*)
	  WRITE(*,'(" WARNING!!! CHEMICAL NOT EXIST IN Datachemical DATABASE  --->  ",A30)') CHEMICAL(i)  
	  WRITE(*,'(" The HQ and CR values will all be zero for this Chemical species!  ")')
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! CHEMICAL NOT EXIST IN Datachemical DATABASE  --->  ",A30)') CHEMICAL(i)  
	  WRITE(99,'(" The HQ and CR values will all be zero for this Chemical species!  ")')
      WRITE(99,*)
!
      ENDIF
!
	  GOTO 56
	  ENDIF
!
!
	  DO j=1,NTIMEexp 
!  
	  DELTAT(J)=1.0
!
!
      CALL REDISTRI(iii,l,j,NVIAS,NIDADE,EF_INI,SD_EF_INI,BW_INI,SD_BW_INI,IR_INI,SD_IR_INI,ET_INI,SD_ET_INI,SA_INI,SD_SA_INI,&
	  EV_INI,SD_EV_INI,AF_INI,SD_AF_INI,EF,SD_EF,BW,SD_BW,IR,SD_IR,ET,SD_ET,SA,SD_SA,EV,SD_EV,AF,SD_AF,MUTAGENIC,ADAF,&
	  AT_INI,SD_AT_INI,AT,SD_AT)  
!
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1).AND.(m.EQ.1))THEN
!
      write(17,'( 3x,I3,3x,F9.1,1x,F9.1,3x,F4.1,2x,F4.1,2x,E8.3,5(1x,E8.3),4(1x,E8.3),3x,F3.1,3x,F3.1,5x,F3.1,3x,F3.1,4x,E8.3,1x,E8.3,20(1x,E8.3),20(1x,E8.3))') j,AT(2,3),SD_AT(2,3),BW,SD_BW,ET(1),SD_ET(1),ET(2),SD_ET(2),ET(3),SD_ET(3),SA(1),SD_SA(1),SA(2),SD_SA(2),EV(1),SD_EV(1),EV(2),SD_EV(2),AF,SD_AF,(IR(iou), iou=1,20),(SD_IR(IQW), IQW=1,20)
!
      ENDIF
!
!
	  DO k=1,NLOCAL

!
!
!*******************************************************************************************
      IF(W(2).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,4).EQV..TRUE.)) THEN
!
      EF2=EF(2)
	  SD_EF2=SD_EF(2)
!
	  ED2=DELTAT(J)
	  SD_ED2=0.0
	  SD_ED_EXP(2,i)=SD_ED2
!
      IF(m.EQ.1)THEN
	  AT2=AT(m,2)*DELTAT(J)
	  ELSE
	  AT2=AT(m,2)
	  ENDIF
!
	  SD_AT2=SD_AT(m,2)
      R_IRfruit2=IR(2)
	  SD_IRfruit2=SD_IR(2)
!
!
	  BTF_2=BTF(iii,1)
	  SD_BTF_2=SD_BTF(iii,1)
	  KEYCONC_2=KEYCONC(i,4)
	  FI2=FI(2)
	  SD_FI2=SD_FI(2)
!
!
!
      CALL DOSEWAY2(CHEMICAL,KSTOPBTF1,NCHEM,NTIME,NLOCAL,KEYCONC_2,i,j,k,CSOIL,SD_CSOIL,CFRUIT,SD_CFRUIT,EF2,SD_EF2,ED2,SD_ED2,AT2,SD_AT2,&
	  BW,SD_BW,BTF_2,SD_BTF_2,R_IRfruit2,SD_IRfruit2,FI2,SD_FI2,AD,SD_AD)
!
!
!
      ADD(l,m,2,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!
      SD_ADD(l,m,2,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
!	  
!
!     CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(2).EQ.'Chronic')THEN
	  RfD2=RfD(iii,1,1)
	  SD_RfD2=SD_RfD(iii,1,1)
	  SF2=SF(iii,1,1)
	  SD_SF2=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(2).EQ.'Subchronic')THEN
	  RfD2=RfD(iii,1,2)
	  SD_RfD2=SD_RfD(iii,1,2)
	  SF2=SF(iii,1,2)
	  SD_SF2=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(2).EQ.'Acute')THEN
	  RfD2=RfD(iii,1,3)
	  SD_RfD2=SD_RfD(iii,1,3)
	  SF2=SF(iii,1,3)
	  SD_SF2=SD_SF(iii,1,3)
	  ENDIF      
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,2,i,j,k).GT.0.0).AND.(RfD2.GT.0.0).AND.(BAF(iii,8).GT.0.0))THEN
      HQ(l,2,i,j,k)=ADD(l,m,2,i,j,k)*BAF(iii,8)/RfD2   ! o RfD 2 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,2,i,j,k)=HQ(l,2,i,j,k)* SQRT((SD_ADD(l,m,2,i,j,k)/ADD(l,m,2,i,j,k))**2+(SD_RfD2/RfD2)**2+(SD_BAF(iii,8)/BAF(iii,8))**2)
	  ELSE
	  HQ(l,2,i,j,k)=0.0
	  SD_HQ(l,2,i,j,k)=0.0
	  IF((RfD2.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!						 
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,2,i,j,k).GT.0.0).AND.(SF2.GT.0.0).AND.(BAF(iii,8).GT.0.0))THEN
	  CR(l,2,i,j,k)=ADD(l,m,2,i,j,k)*BAF(iii,8)*SF2*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,2,i,j,k)=CR(l,2,i,j,k)* SQRT((SD_ADD(l,m,2,i,j,k)/ADD(l,m,2,i,j,k))**2+(SD_SF2/SF2)**2+(SD_BAF(iii,8)/BAF(iii,8))**2)
	  ELSE
	  CR(l,2,i,j,k)=0.0
	  SD_CR(l,2,i,j,k)=0.0
	  IF((SF2.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
	  ENDIF
!
	  IF((BAF(iii,8).LE.0.0).AND.(KSTOPBAF(8).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(8)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The fruit ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(2)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway2 was not calculed, because soil, fruits, seeds or tubers concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY2 WAS NOT CALCULED, BECAUSE SOIL AND FRUITS, SEEDS OR TUBERS'
      PRINT*, '  CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!

      IF(W(5).EQV..TRUE.)THEN
	  IF(KEYCONC(i,2).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
!
      EF5=EF(5)
	  SD_EF5=SD_EF(5)
!
	  ED5=DELTAT(J)
	  SD_ED5=0.0
	  SD_ED_EXP(5,i)=SD_ED5
!
      IF(m.EQ.1)THEN
	  AT5=AT(m,5)*DELTAT(J)
	  ELSE
	  AT5=AT(m,5)
	  ENDIF
!
	  SD_AT5=SD_AT(m,5)
	  R_IR5=IR(9)
	  SD_IR5=SD_IR(9)
!
!
      CALL DOSEWAY5(NCHEM,NTIME,NLOCAL,i,j,k,CWATER,SD_CWATER,BW,SD_BW,EF5,SD_EF5,ED5,SD_ED5,AT5,SD_AT5,R_IR5,SD_IR5,AD,SD_AD)
!
!
      ADD(l,m,5,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,5,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(5).EQ.'Chronic')THEN
	  RfD5=RfD(iii,4,1)
	  SD_RfD5=SD_RfD(iii,4,1)
	  SF5=SF(iii,4,1)
	  SD_SF5=SD_SF(iii,4,1)
      ELSEIF(NRISKTYPE(5).EQ.'Subchronic')THEN
	  RfD5=RfD(iii,4,2)
	  SD_RfD5=SD_RfD(iii,4,2)
	  SF5=SF(iii,4,2)
	  SD_SF5=SD_SF(iii,4,2)
      ELSEIF(NRISKTYPE(5).EQ.'Acute')THEN
	  RfD5=RfD(iii,4,3)
	  SD_RfD5=SD_RfD(iii,4,3)
	  SF5=SF(iii,4,3)
	  SD_SF5=SD_SF(iii,4,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,5,i,j,k).GT.0.0).AND.(RfD5.GT.0.0).AND.(BAF(iii,2).GT.0.0))THEN
      HQ(l,5,i,j,k)=ADD(l,m,5,i,j,k)*BAF(iii,2)/RfD5   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,5,i,j,k)=HQ(l,5,i,j,k)* SQRT((SD_ADD(l,m,5,i,j,k)/ADD(l,m,5,i,j,k))**2+(SD_RfD5/RfD5)**2+(SD_BAF(iii,2)/BAF(iii,2))**2)
	  ELSE
	  HQ(l,5,i,j,k)=0.0
	  SD_HQ(l,5,i,j,k)=0.0
	  IF((RfD5.LE.0.0).AND.(KSTOPRFD(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(2)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral water RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,5,i,j,k).GT.0.0).AND.(SF5.GT.0.0).AND.(BAF(iii,2).GT.0.0))THEN
	  CR(l,5,i,j,k)=ADD(l,m,5,i,j,k)*BAF(iii,2)*SF5*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,5,i,j,k)=CR(l,5,i,j,k)* SQRT((SD_ADD(l,m,5,i,j,k)/ADD(l,m,5,i,j,k))**2+(SD_SF5/SF5)**2+(SD_BAF(iii,2)/BAF(iii,2))**2)
	  ELSE
	  CR(l,5,i,j,k)=0.0
	  SD_CR(l,5,i,j,k)=0.0
	  IF((SF5.LE.0.0).AND.(KSTOPSF(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(2)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral water SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,2).LE.0.0).AND.(KSTOPBAF(2).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(2)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The water ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(5)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway5 was not calculed, because water (DRINKING_WATER) concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY5 WAS NOT CALCULED, BECAUSE WATER (DRINKING_WATER) CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
!
      IF(W(6).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,5).EQV..TRUE.))THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
! 	  
	  EF6=EF(6)
	  SD_EF6=SD_EF(6)
!
	  ED6=DELTAT(J)
	  SD_ED6=0.0
	  SD_ED_EXP(6,i)=SD_ED6
!
      IF(m.EQ.1)THEN
	  AT6=AT(m,6)*DELTAT(J)
	  ELSE
	  AT6=AT(m,6)
	  ENDIF
!
	  SD_AT6=SD_AT(m,6)
	  R_IRveg6=IR(10)
	  SD_IRveg6=SD_IR(10)
!
!
	  BTF_6S=BTF(iii,7)
	  SD_BTF_6S=SD_BTF(iii,7)
	  KEYCONC_6=KEYCONC(i,5)
	  FI6=FI(5)
	  SD_FI6=SD_FI(5)
!
!
      CALL DOSEWAY6(CHEMICAL,KSTOPBTF_6S,NCHEM,NTIME,NLOCAL,KEYCONC_6,i,j,k,CSOIL,SD_CSOIL,CLEAVES,SD_CLEAVES,EF6,SD_EF6,ED6,SD_ED6,&
	  AT6,SD_AT6,BW,SD_BW,BTF_6S,SD_BTF_6S,R_IRveg6,SD_IRveg6,FI6,SD_FI6,AD,SD_AD)
!
!
      ADD(l,m,6,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,6,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(6).EQ.'Chronic')THEN
	  RfD6=RfD(iii,1,1)
	  SD_RfD6=SD_RfD(iii,1,1)
	  SF6=SF(iii,1,1)
	  SD_SF6=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(6).EQ.'Subchronic')THEN
	  RfD6=RfD(iii,1,2)
	  SD_RfD6=SD_RfD(iii,1,2)
	  SF6=SF(iii,1,2)
	  SD_SF6=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(6).EQ.'Acute')THEN
	  RfD6=RfD(iii,1,3)
	  SD_RfD6=SD_RfD(iii,1,3)
	  SF6=SF(iii,1,3)
	  SD_SF6=SD_SF(iii,1,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,6,i,j,k).GT.0.0).AND.(RfD6.GT.0.0).AND.(BAF(iii,7).GT.0.0))THEN
	  HQ(l,6,i,j,k)=ADD(l,m,6,i,j,k)*BAF(iii,7)/RfD6   ! o RfD 6 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,6,i,j,k)=HQ(l,6,i,j,k)* SQRT((SD_ADD(l,m,6,i,j,k)/ADD(l,m,6,i,j,k))**2+(SD_RfD6/RfD6)**2+(SD_BAF(iii,7)/BAF(iii,7))**2)	 
	  ELSE
	  HQ(l,6,i,j,k)=0.0
	  SD_HQ(l,6,i,j,k)=0.0
	  IF((RfD6.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF

!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,6,i,j,k).GT.0.0).AND.(SF6.GT.0.0).AND.(BAF(iii,7).GT.0.0))THEN
	  CR(l,6,i,j,k)=ADD(l,m,6,i,j,k)*BAF(iii,7)*SF6*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,6,i,j,k)=CR(l,6,i,j,k)* SQRT((SD_ADD(l,m,6,i,j,k)/ADD(l,m,6,i,j,k))**2+(SD_SF6/SF6)**2+(SD_BAF(iii,7)/BAF(iii,7))**2)
	  ELSE
	  CR(l,6,i,j,k)=0.0
	  SD_CR(l,6,i,j,k)=0.0
	  IF((SF6.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,7).LE.0.0).AND.(KSTOPBAF(7).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(7)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The vegetable ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(6)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway6 was not calculed, because soil and vegetables (leaves) concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY6 WAS NOT CALCULED, BECAUSE SOIL AND VEGETABLES (LEAVES) CONCENTRATIONS NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!*******************************************************************************************
!
!
      IF(W(3).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,6).EQV..TRUE.)) THEN
!
      EF3=EF(3)
	  SD_EF3=SD_EF(3)
!
	  ED3=DELTAT(J)
	  SD_ED3=0.0
	  SD_ED_EXP(3,i)=SD_ED3
!
      IF(m.EQ.1)THEN
	  AT3=AT(m,3)*DELTAT(J)
	  ELSE
	  AT3=AT(m,3)
	  ENDIF
!
	  SD_AT3=SD_AT(m,3)
      R_IRcarne3=IR(3)
	  SD_IRcarne3=SD_IR(3)
	  R_IRali3=IR(4)
	  SD_IRali3=SD_IR(4)
	  R_IRagua3=IR(11)
	  SD_IRagua3=SD_IR(11)
      R_IRsolo3=IR(5)
	  SD_IRsolo3=SD_IR(5)
!
!
	  BTF_SV3=BTF(iii,2)
	  SD_BTF_SV3=SD_BTF(iii,2)
	  BTF_VB3=BTF(iii,3)
	  SD_BTF_VB3=SD_BTF(iii,3)
	  BTF_SB3=BTF(iii,4)
	  SD_BTF_SB3=SD_BTF(iii,4)
	  BTF_WB3=BTF(iii,8)
	  SD_BTF_WB3=SD_BTF(iii,8)
	  KEYCONC_B3=KEYCONC(i,6)
	  KEYCONC_S3=KEYCONC(i,1)
	  KEYCONC_W3=KEYCONC(i,14)
	  Fa3=Fa(1)
	  SD_Fa3=SD_Fa(1)
	  Fp3=Fp(1)
	  SD_Fp3=SD_Fp(1)
	  FI3=FI(3)
	  SD_FI3=SD_FI(3)
	  fw3=fw(iii)
	  SD_fw3=SD_fw(iii)
!
!
!
      CALL DOSEWAY3(CHEMICAL,KSTOPBTF_SB3,KSTOPBTF_SV3,KSTOPBTF_VB3,KSTOPBTF_WB3,KSTOPfw3,NCHEM,NTIME,NLOCAL,BTF_SV3,SD_BTF_SV3,BTF_VB3,SD_BTF_VB3,BTF_SB3,SD_BTF_SB3,BTF_WB3,SD_BTF_WB3,i,j,k,&
	  CSOIL,SD_CSOIL,CWATEROTHER,SD_CWATEROTHER,CBEEF,SD_CBEEF,EF3,SD_EF3,ED3,SD_ED3,AT3,SD_AT3,BW,SD_BW,KEYCONC_B3,KEYCONC_S3,KEYCONC_W3,&
	  Fa3,SD_Fa3,Fp3,SD_Fp3,R_IRali3,SD_IRali3,R_IRsolo3,SD_IRsolo3,R_IRagua3,SD_IRagua3,R_IRcarne3,SD_IRcarne3,FI3,SD_FI3,fw3,SD_fw3,AD,SD_AD)
!
!
!
      ADD(l,m,3,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,3,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(3).EQ.'Chronic')THEN
	  RfD3=RfD(iii,1,1)
	  SD_RfD3=SD_RfD(iii,1,1)
	  SF3=SF(iii,1,1)
	  SD_SF3=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(3).EQ.'Subchronic')THEN
	  RfD3=RfD(iii,1,2)
	  SD_RfD3=SD_RfD(iii,1,2)
	  SF3=SF(iii,1,2)
	  SD_SF3=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(3).EQ.'Acute')THEN
	  RfD3=RfD(iii,1,3)
	  SD_RfD3=SD_RfD(iii,1,3)
	  SF3=SF(iii,1,3)
	  SD_SF3=SD_SF(iii,1,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,3,i,j,k).GT.0.0).AND.(RfD3.GT.0.0).AND.(BAF(iii,9).GT.0.0))THEN
      HQ(l,3,i,j,k)=ADD(l,m,3,i,j,k)*BAF(iii,9)/RfD3   ! o RfD 3 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,3,i,j,k)=HQ(l,3,i,j,k)* SQRT((SD_ADD(l,m,3,i,j,k)/ADD(l,m,3,i,j,k))**2+(SD_RfD3/RfD3)**2+(SD_BAF(iii,9)/BAF(iii,9))**2)
	  ELSE
	  HQ(l,3,i,j,k)=0.0
	  SD_HQ(l,3,i,j,k)=0.0
	  IF((RfD3.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,3,i,j,k).GT.0.0).AND.(SF3.GT.0.0).AND.(BAF(iii,9).GT.0.0))THEN
	  CR(l,3,i,j,k)=ADD(l,m,3,i,j,k)*BAF(iii,9)*SF3*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,3,i,j,k)=CR(l,3,i,j,k)* SQRT((SD_ADD(l,m,3,i,j,k)/ADD(l,m,3,i,j,k))**2+(SD_SF3/SF3)**2+(SD_BAF(iii,9)/BAF(iii,9))**2)  
	  ELSE
	  CR(l,3,i,j,k)=0.0
	  SD_CR(l,3,i,j,k)=0.0
	  IF((SF3.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF      
	  ENDIF
!
	  IF((BAF(iii,9).LE.0.0).AND.(KSTOPBAF(9).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(9)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The beef ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(3)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway3 was not calculed, because soil, water (OTHER_WATERS) and meat concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY3 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'
      PRINT*, '  MEAT CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
!
      IF(W(4).EQV..TRUE.)THEN
	  IF((KEYCONC(i,1).EQV..TRUE.).OR.(KEYCONC(i,14).EQV..TRUE.).OR.(KEYCONC(i,7).EQV..TRUE.)) THEN
! 	  
	  EF4=EF(4)
	  SD_EF4=SD_EF(4)
!
	  ED4=DELTAT(J)
	  SD_ED4=0.0
	  SD_ED_EXP(4,i)=SD_ED4
!
      IF(m.EQ.1)THEN
	  AT4=AT(m,4)*DELTAT(J)
	  ELSE
	  AT4=AT(m,4)
	  ENDIF
!
	  SD_AT4=SD_AT(m,4)
      R_IRleite4=IR(6)
	  SD_IRleite4=SD_IR(6)
	  R_IRali4=IR(7)
	  SD_IRali4=SD_IR(7)
	  R_IRagua4=IR(12)
	  SD_IRagua4=SD_IR(12)
      R_IRsolo4=IR(8)
	  SD_IRsolo4=SD_IR(8)
!
!
	  BTF_SV4=BTF(iii,2)
	  SD_BTF_SV4=SD_BTF(iii,2)
	  BTF_VM4=BTF(iii,5)
	  SD_BTF_VM4=SD_BTF(iii,5)
	  BTF_SM4=BTF(iii,6)
	  SD_BTF_SM4=SD_BTF(iii,6)
	  BTF_WM4=BTF(iii,9)
	  SD_BTF_WM4=SD_BTF(iii,9)
	  KEYCONC_M4=KEYCONC(i,7)
	  KEYCONC_S4=KEYCONC(i,1)
	  KEYCONC_W4=KEYCONC(i,14)
	  Fa4=Fa(2)
	  SD_Fa4=SD_Fa(2)
	  Fp4=Fp(2)
	  SD_Fp4=SD_Fp(2)
	  FI4=FI(4)
	  SD_FI4=SD_FI(4)
	  fw4=fw(iii)
	  SD_fw4=SD_fw(iii)
!
!
!
      CALL DOSEWAY4(CHEMICAL,KSTOPBTF_SM4,KSTOPBTF_SV4,KSTOPBTF_VM4,KSTOPBTF_WM4,KSTOPfw4,NCHEM,NTIME,NLOCAL,BTF_SV4,SD_BTF_SV4,BTF_VM4,SD_BTF_VM4,BTF_SM4,SD_BTF_SM4,BTF_WM4,SD_BTF_WM4,i,j,k,&
	  CSOIL,SD_CSOIL,CWATEROTHER,SD_CWATEROTHER,CMILK,SD_CMILK,EF4,SD_EF4,ED4,SD_ED4,AT4,SD_AT4,BW,SD_BW,KEYCONC_M4,KEYCONC_S4,KEYCONC_W4,&
	  Fa4,SD_Fa4,Fp4,SD_Fp4,R_IRali4,SD_IRali4,R_IRsolo4,SD_IRsolo4,R_IRagua4,SD_IRagua4,R_IRleite4,SD_IRleite4,&
	  FI4,SD_FI4,fw4,SD_fw4,AD,SD_AD)
!
!
!
      ADD(l,m,4,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,4,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 	
!
      IF(NRISKTYPE(4).EQ.'Chronic')THEN
	  RfD4=RfD(iii,1,1)
	  SD_RfD4=SD_RfD(iii,1,1)
	  SF4=SF(iii,1,1)
	  SD_SF4=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(4).EQ.'Subchronic')THEN
	  RfD4=RfD(iii,1,2)
	  SD_RfD4=SD_RfD(iii,1,2)
	  SF4=SF(iii,1,2)
	  SD_SF4=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(4).EQ.'Acute')THEN
	  RfD4=RfD(iii,1,3)
	  SD_RfD4=SD_RfD(iii,1,3)
	  SF4=SF(iii,1,3)
	  SD_SF4=SD_SF(iii,1,3)
	  ENDIF	
!								  
      IF(m.EQ.1)THEN
	  IF((ADD(l,m,4,i,j,k).GT.0.0).AND.(RfD4.GT.0.0).AND.(BAF(iii,10).GT.0.0))THEN
      HQ(l,4,i,j,k)=ADD(l,m,4,i,j,k)*BAF(iii,10)/RfD4   ! o RfD 4 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,4,i,j,k)=HQ(l,4,i,j,k)* SQRT((SD_ADD(l,m,4,i,j,k)/ADD(l,m,4,i,j,k))**2+(SD_RfD4/RfD4)**2+(SD_BAF(iii,10)/BAF(iii,10))**2)
	  ELSE
	  HQ(l,4,i,j,k)=0.0
	  SD_HQ(l,4,i,j,k)=0.0
	  IF((RfD4.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,4,i,j,k).GT.0.0).AND.(SF4.GT.0.0).AND.(BAF(iii,10).GT.0.0))THEN
	  CR(l,4,i,j,k)=ADD(l,m,4,i,j,k)*BAF(iii,10)*SF4*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,4,i,j,k)=CR(l,4,i,j,k)* SQRT((SD_ADD(l,m,4,i,j,k)/ADD(l,m,4,i,j,k))**2+(SD_SF4/SF4)**2+(SD_BAF(iii,10)/BAF(iii,10))**2)
	  ELSE
	  CR(l,4,i,j,k)=0.0
	  SD_CR(l,4,i,j,k)=0.0
	  IF((SF4.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,10).LE.0.0).AND.(KSTOPBAF(10).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(10)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The milk ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(4)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway4 was not calculed, because soil, water (OTHER_WATERS) and milk concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY4 WAS NOT CALCULED, BECAUSE SOIL, WATER (OTHER_WATERS)'
      PRINT*, '  AND MILK CONCENTRATIONS NOT EXIST'
	  ENDIF
!
!
	  ENDIF
	  ENDIF
!
!
!
!*******************************************************************************************
!
!
!
      IF(W(13).EQV..TRUE.)THEN
	  IF(KEYCONC(i,13).EQV..TRUE.)THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
	  EF13=EF(13)
	  SD_EF13=SD_EF(13)
!
	  ED13=DELTAT(J)
	  SD_ED13=0.0
	  SD_ED_EXP(13,i)=SD_ED13
!
      IF(m.EQ.1)THEN
	  AT13=AT(m,13)*DELTAT(J)
	  ELSE
	  AT13=AT(m,13)
	  ENDIF
!
	  SD_AT13=SD_AT(m,13)
	  ET13=ET(3)
	  SD_ET13=SD_ET(3)
	  SA13=SA(1)
	  SD_SA13=SD_SA(1)
!
	  PC13=PC(iii)
	  SD_PC13=SD_PC(iii)
!
      CF13=CF(2)
	  SD_CF13=SD_CF(2)
!
	  EV13=EV(1)
	  SD_EV13=SD_EV(1)
!
      CALL DOSEWAY13(NCHEM,NTIME,NLOCAL,KSTOPPC13,CHEMICAL,i,j,k,CWATERDER,SD_CWATERDER,BW,SD_BW,EV13,SD_EV13,EF13,SD_EF13,ED13,SD_ED13,ET13,SD_ET13,AT13,SD_AT13,&
	  PC13,SD_PC13,SA13,SD_SA13,CF13,SD_CF13,AD,SD_AD)
!
!
      ADD(l,m,13,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,13,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(13).EQ.'Chronic')THEN
	  RfD13=RfD(iii,6,1)
	  SD_RfD13=SD_RfD(iii,6,1)
	  SF13=SF(iii,6,1)
	  SD_SF13=SD_SF(iii,6,1)
      ELSEIF(NRISKTYPE(13).EQ.'Subchronic')THEN
	  RfD13=RfD(iii,6,2)
	  SD_RfD13=SD_RfD(iii,6,2)
	  SF13=SF(iii,6,2)
	  SD_SF13=SD_SF(iii,6,2)
      ELSEIF(NRISKTYPE(13).EQ.'Acute')THEN
	  RfD13=RfD(iii,6,3)
	  SD_RfD13=SD_RfD(iii,6,3)
	  SF13=SF(iii,6,3)
	  SD_SF13=SD_SF(iii,6,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,13,i,j,k).GT.0.0).AND.(RfD13.GT.0.0).AND.(BAF(iii,6).GT.0.0))THEN
      HQ(l,13,i,j,k)=ADD(l,m,13,i,j,k)*BAF(iii,6)/RfD13   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,13,i,j,k)=HQ(l,13,i,j,k)* SQRT((SD_ADD(l,m,13,i,j,k)/ADD(l,m,13,i,j,k))**2+(SD_RfD13/RfD13)**2+(SD_BAF(iii,6)/BAF(iii,6))**2)
	  ELSE
	  HQ(l,13,i,j,k)=0.0
	  SD_HQ(l,13,i,j,k)=0.0
	  IF((RfD13.LE.0.0).AND.(KSTOPRFD(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal water RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,13,i,j,k).GT.0.0).AND.(SF13.GT.0.0).AND.(BAF(iii,6).GT.0.0))THEN
	  CR(l,13,i,j,k)=ADD(l,m,13,i,j,k)*BAF(iii,6)*SF13*ADAF(3)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,13,i,j,k)=CR(l,13,i,j,k)* SQRT((SD_ADD(l,m,13,i,j,k)/ADD(l,m,13,i,j,k))**2+(SD_SF13/SF13)**2+(SD_BAF(iii,6)/BAF(iii,6))**2)
	  ELSE
	  CR(l,13,i,j,k)=0.0
	  SD_CR(l,13,i,j,k)=0.0
	  IF((SF13.LE.0.0).AND.(KSTOPSF(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal water SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,6).LE.0.0).AND.(KSTOPBAF(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal water BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
	  IF((KSTOPWAY(13).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(13)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway13 was not calculed, because water (BATH_WATER) concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY13 WAS NOT CALCULED, BECAUSE WATER (BATH_WATER) CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!*******************************************************************************************
!
!
      IF(W(1).EQV..TRUE.)THEN
	  IF(KEYCONC(i,1).EQV..TRUE.) THEN
!
!
      EF1=EF(1)
	  SD_EF1=SD_EF(1)
!
	  ED1=DELTAT(J)
	  SD_ED1=0.0
	  SD_ED_EXP(1,i)=SD_ED1
!
      IF(m.EQ.1)THEN
	  AT1=AT(m,1)*DELTAT(J)
	  ELSE
	  AT1=AT(m,1)
	  ENDIF
!
	  SD_AT1=SD_AT(m,1)
	  R_IR1=IR(1)
	  SD_IR1=SD_IR(1)
!
!
	  FI1=FI(1)
	  SD_FI1=SD_FI(1)
!
      CF1=CF(1)
	  SD_CF1=SD_CF(1)
!
!
      CALL DOSEWAY1(NCHEM,NTIME,NLOCAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF1,SD_EF1,ED1,SD_ED1,AT1,SD_AT1,R_IR1,SD_IR1,&
	  FI1,SD_FI1,CF1,SD_CF1,AD,SD_AD)
!
!      WRITE(*,*) AD
!
      ADD(l,m,1,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!
!      WRITE(*,*) ADD(2,m,1,i,j,k)
!
      SD_ADD(l,m,1,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
!  
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(1).EQ.'Chronic')THEN
	  RfD1=RfD(iii,1,1)
	  SD_RfD1=SD_RfD(iii,1,1)
	  SF1=SF(iii,1,1)
	  SD_SF1=SD_SF(iii,1,1)
      ELSEIF(NRISKTYPE(1).EQ.'Subchronic')THEN
	  RfD1=RfD(iii,1,2)
	  SD_RfD1=SD_RfD(iii,1,2)
	  SF1=SF(iii,1,2)
	  SD_SF1=SD_SF(iii,1,2)
      ELSEIF(NRISKTYPE(1).EQ.'Acute')THEN
	  RfD1=RfD(iii,1,3)
	  SD_RfD1=SD_RfD(iii,1,3)
	  SF1=SF(iii,1,3)
	  SD_SF1=SD_SF(iii,1,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,1,i,j,k).GT.0.0).AND.(RfD1.GT.0.0).AND.(BAF(iii,1).GT.0.0))THEN
	  HQ(l,1,i,j,k)=ADD(l,m,1,i,j,k)*BAF(iii,1)/RfD1   ! o RfD 1 representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,1,i,j,k)=HQ(l,1,i,j,k)* SQRT((SD_ADD(l,m,1,i,j,k)/ADD(l,m,1,i,j,k))**2+(SD_RfD1/RfD1)**2+(SD_BAF(iii,1)/BAF(iii,1))**2)
	  ELSE
	  HQ(l,1,i,j,k)=0.0
	  SD_HQ(l,1,i,j,k)=0.0
	  IF((RfD1.LE.0.0).AND.(KSTOPRFD(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,1,i,j,k).GT.0.0).AND.(SF1.GT.0.0).AND.(BAF(iii,1).GT.0.0))THEN
	  CR(l,1,i,j,k)=ADD(l,m,1,i,j,k)*BAF(iii,1)*SF1*ADAF(1)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,1,i,j,k)=CR(l,1,i,j,k)* SQRT((SD_ADD(l,m,1,i,j,k)/ADD(l,m,1,i,j,k))**2+(SD_SF1/SF1)**2+(SD_BAF(iii,1)/BAF(iii,1))**2)
	  ELSE
	  CR(l,1,i,j,k)=0.0
	  SD_CR(l,1,i,j,k)=0.0
	  IF((SF1.LE.0.0).AND.(KSTOPSF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The oral SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
	  ENDIF
!
	  IF((BAF(iii,1).LE.0.0).AND.(KSTOPBAF(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(1)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The soil ingestion BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(1).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(1)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway1 was not calculed, because soil concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY1 WAS NOT CALCULED, BECAUSE SOIL CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
      IF(W(11).EQV..TRUE.)THEN
	  IF(KEYCONC(i,3).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF11=EF(11)
	  SD_EF11=SD_EF(11)
!
	  ED11=DELTAT(J)
	  SD_ED11=0.0
	  SD_ED_EXP(11,i)=SD_ED11
!
      IF(m.EQ.1)THEN
	  AT11=AT(m,11)*DELTAT(J)
	  ELSE
	  AT11=AT(m,11)
	  ENDIF
!
	  SD_AT11=SD_AT(m,11)
!
	  ET_RES_11=ET(1)
	  SD_ET_RES_11=SD_ET(1)
!
!
!
      CALL DOSEWAY11(NCHEM,NTIME,NLOCAL,i,j,k,EF11,SD_EF11,ED11,SD_ED11,AT11,SD_AT11,CPAR,SD_CPAR,&
	  ET_RES_11,SD_ET_RES_11,AD,SD_AD)
!
!
      ADD(l,m,11,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,11,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(11).EQ.'Chronic')THEN
	  RfD11=RfD(iii,2,1)
	  SD_RfD11=SD_RfD(iii,2,1)
	  SF11=SF(iii,2,1)
	  SD_SF11=SD_SF(iii,2,1)
      ELSEIF(NRISKTYPE(11).EQ.'Subchronic')THEN
	  RfD11=RfD(iii,2,2)
	  SD_RfD11=SD_RfD(iii,2,2)
	  SF11=SF(iii,2,2)
	  SD_SF11=SD_SF(iii,2,2)
      ELSEIF(NRISKTYPE(11).EQ.'Acute')THEN
	  RfD11=RfD(iii,2,3)
	  SD_RfD11=SD_RfD(iii,2,3)
	  SF11=SF(iii,2,3)
	  SD_SF11=SD_SF(iii,2,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,11,i,j,k).GT.0.00).AND.(RfD11.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
      HQ(l,11,i,j,k)=ADD(l,m,11,i,j,k)*BAF(iii,3)/RfD11   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,11,i,j,k)=HQ(l,11,i,j,k)* SQRT((SD_ADD(l,m,11,i,j,k)/ADD(l,m,11,i,j,k))**2+(SD_RfD11/RfD11)**2+(SD_BAF(iii,3)/BAF(iii,3))**2)
	  ELSE
	  HQ(l,11,i,j,k)=0.0
	  SD_HQ(l,11,i,j,k)=0.0
	  IF((RfD11.LE.0.0).AND.(KSTOPRFD(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particles RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,11,i,j,k).GT.0.0).AND.(SF11.GT.0.0).AND.(BAF(iii,3).GT.0.0))THEN
	  CR(l,11,i,j,k)=ADD(l,m,11,i,j,k)*BAF(iii,3)*SF11*ADAF(2)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,11,i,j,k)=CR(l,11,i,j,k)* SQRT((SD_ADD(l,m,11,i,j,k)/ADD(l,m,11,i,j,k))**2+(SD_SF11/SF11)**2+(SD_BAF(iii,3)/BAF(iii,3))**2)	
	  ELSE
	  CR(l,11,i,j,k)=0.0
	  SD_CR(l,11,i,j,k)=0.0
	  IF((SF11.LE.0.0).AND.(KSTOPSF(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particles SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,3).LE.0.0).AND.(KSTOPBAF(3).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(3)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation particles RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(11).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(11)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway11 was not calculed, because particulate concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY11 WAS NOT CALCULED, BECAUSE PARTICULATE CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
!
      IF(W(12).EQV..TRUE.)THEN
	  IF(KEYCONC(i,12).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF12=EF(12)
	  SD_EF12=SD_EF(12)
!
	  ED12=DELTAT(J)
	  SD_ED12=0.0
	  SD_ED_EXP(12,i)=SD_ED12
!
      IF(m.EQ.1)THEN
	  AT12=AT(m,12)*DELTAT(J)
	  ELSE
	  AT12=AT(m,12)
	  ENDIF
!
	  SD_AT12=SD_AT(m,12)
!
	  ET_RES_12=ET(2)
	  SD_ET_RES_12=SD_ET(2)
!
!
!
      CALL DOSEWAY12(NCHEM,NTIME,NLOCAL,i,j,k,EF12,SD_EF12,ED12,SD_ED12,AT12,SD_AT12,CSTEAM,SD_CSTEAM,&
	  ET_RES_12,SD_ET_RES_12,AD,SD_AD)
!
!
      ADD(l,m,12,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,12,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(12).EQ.'Chronic')THEN
	  RfD12=RfD(iii,5,1)
	  SD_RfD12=SD_RfD(iii,5,1)
	  SF12=SF(iii,5,1)
	  SD_SF12=SD_SF(iii,5,1)
      ELSEIF(NRISKTYPE(12).EQ.'Subchronic')THEN
	  RfD12=RfD(iii,5,2)
	  SD_RfD12=SD_RfD(iii,5,2)
	  SF12=SF(iii,5,2)
	  SD_SF12=SD_SF(iii,5,2)
      ELSEIF(NRISKTYPE(12).EQ.'Acute')THEN
	  RfD12=RfD(iii,5,3)
	  SD_RfD12=SD_RfD(iii,5,3)
	  SF12=SF(iii,5,3)
	  SD_SF12=SD_SF(iii,5,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,12,i,j,k).GT.0.00).AND.(RfD12.GT.0.0).AND.(BAF(iii,4).GT.0.0))THEN
      HQ(l,12,i,j,k)=ADD(l,m,12,i,j,k)*BAF(iii,4)/RfD12   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,12,i,j,k)=HQ(l,12,i,j,k)* SQRT((SD_ADD(l,m,12,i,j,k)/ADD(l,m,12,i,j,k))**2+(SD_RfD12/RfD12)**2+(SD_BAF(iii,4)/BAF(iii,4))**2)
	  ELSE
	  HQ(l,12,i,j,k)=0.0
	  SD_HQ(l,12,i,j,k)=0.0
	  IF((RfD12.LE.0.0).AND.(KSTOPRFD(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
	  IF((ADD(l,m,12,i,j,k).GT.0.0).AND.(SF12.GT.0.0).AND.(BAF(iii,4).GT.0.0))THEN
	  CR(l,12,i,j,k)=ADD(l,m,12,i,j,k)*BAF(iii,4)*SF12*ADAF(2)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,12,i,j,k)=CR(l,12,i,j,k)* SQRT((SD_ADD(l,m,12,i,j,k)/ADD(l,m,12,i,j,k))**2+(SD_SF12/SF12)**2+(SD_BAF(iii,4)/BAF(iii,4))**2)	
	  ELSE
	  CR(l,12,i,j,k)=0.0
	  SD_CR(l,12,i,j,k)=0.0
	  IF((SF12.LE.0.0).AND.(KSTOPSF(6).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(6)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF
!
	  IF((BAF(iii,4).LE.0.0).AND.(KSTOPBAF(4).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(4)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The inhalation steam BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(12).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(12)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway12 was not calculed, because steam concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY12 WAS NOT CALCULED, BECAUSE STEAM CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!
!
!
      IF(W(14).EQV..TRUE.)THEN
	  IF(KEYCONC(i,1).EQV..TRUE.) THEN	 ! i= numero metal	j= tipo de concentração (exemplo: solo=1, agua=2, ar=3...)
!
      EF14=EF(14)
	  SD_EF14=SD_EF(14)
!
	  ED14=DELTAT(J)
	  SD_ED14=0.0
	  SD_ED_EXP(14,i)=SD_ED14
!
      IF(m.EQ.1)THEN
	  AT14=AT(m,14)*DELTAT(J)
	  ELSE
	  AT14=AT(m,14)
	  ENDIF
!
	  SD_AT14=SD_AT(m,14)
	  SA14=SA(2)
	  SD_SA14=SD_SA(2)
	  AFsoil14=AF
	  SD_AFsoil14=SD_AF
	  EV14=EV(2)
	  SD_EV14=SD_EV(2)
!
!
	  ABS14=ABS_(iii)
	  SD_ABS14=SD_ABS(iii)
!
      CF14=CF(3)
	  SD_CF14=SD_CF(3)
!
      CALL DOSEWAY14(NCHEM,NTIME,NLOCAL,KSTOPABS14,CHEMICAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF14,SD_EF14,ED14,SD_ED14,AT14,SD_AT14,ABS14,SD_ABS14,&
	  AFsoil14,SD_AFsoil14,SA14,SD_SA14,CF14,SD_CF14,EV14,SD_EV14,AD,SD_AD)
!
!
      ADD(l,m,14,i,j,k)=AD     !MATRIZ QUE GURDA OS VALORES DAS DOSES POR METAL, TEMPO E LOCAL
!	  
      SD_ADD(l,m,14,i,j,k)=SD_AD     !MATRIZ QUE GURDA OS VALORES DAS INCERTEZAS DAS DOSES POR METAL, TEMPO E LOCAL
! 
!	  														  
!	  CALCULO DO RISCO 										  
!
      IF(NRISKTYPE(14).EQ.'Chronic')THEN
	  RfD14=RfD(iii,3,1)
	  SD_RfD14=SD_RfD(iii,3,1)
	  SF14=SF(iii,3,1)
	  SD_SF14=SD_SF(iii,3,1)
      ELSEIF(NRISKTYPE(14).EQ.'Subchronic')THEN
	  RfD14=RfD(iii,3,2)
	  SD_RfD14=SD_RfD(iii,3,2)
	  SF14=SF(iii,3,2)
	  SD_SF14=SD_SF(iii,3,2)
      ELSEIF(NRISKTYPE(14).EQ.'Acute')THEN
	  RfD14=RfD(iii,3,3)
	  SD_RfD14=SD_RfD(iii,3,3)
	  SF14=SF(iii,3,3)
	  SD_SF14=SD_SF(iii,3,3)
	  ENDIF
!
      IF(m.EQ.1)THEN
      IF((ADD(l,m,14,i,j,k).GT.0.0).AND.(RfD14.GT.0.0).AND.(BAF(iii,5).GT.0.0))THEN
      HQ(l,14,i,j,k)=ADD(l,m,14,i,j,k)*BAF(iii,5)/RfD14   ! o RfD (iii,1,1) representa: iii=METAL, 1=oral, 1= cronico
	  SD_HQ(l,14,i,j,k)=HQ(l,14,i,j,k)* SQRT((SD_ADD(l,m,14,i,j,k)/ADD(l,m,14,i,j,k))**2+(SD_RfD14/RfD14)**2+(SD_BAF(iii,5)/BAF(iii,5))**2)
	  ELSE
	  HQ(l,14,i,j,k)=0.0
	  SD_HQ(l,14,i,j,k)=0.0
	  IF((RfD14.LE.0.0).AND.(KSTOPRFD(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPRFD(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil RfD is equal to zero!!! HQ value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
!
      ELSEIF(m.EQ.2)THEN
      IF((ADD(l,m,14,i,j,k).GT.0.0).AND.(SF14.GT.0.0).AND.(BAF(iii,5).GT.0.0))THEN
	  CR(l,14,i,j,k)=ADD(l,m,14,i,j,k)*BAF(iii,5)*SF14*ADAF(3)   ! o SF (i,1,1) representa: i=METAL, 1=oral, 1= cronico
	  SD_CR(l,14,i,j,k)=CR(l,14,i,j,k)* SQRT((SD_ADD(l,m,14,i,j,k)/ADD(l,m,14,i,j,k))**2+(SD_SF14/SF14)**2+(SD_BAF(iii,5)/BAF(iii,5))**2)
	  ELSE
	  CR(l,14,i,j,k)=0.0
	  SD_CR(l,14,i,j,k)=0.0
	  IF((SF14.LE.0.0).AND.(KSTOPSF(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPSF(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil SF is equal to zero!!! CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
	  ENDIF
      ENDIF	 ! IF "M"
!
	  IF((BAF(iii,5).LE.0.0).AND.(KSTOPBAF(5).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPBAF(5)=1
	  WRITE(99,*)
	  WRITE(99,'(" Chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" The dermal soil BAF is equal to zero!!! HQ and CR value will be equal to zero too!! ")')
      WRITE(99,*)
	  ENDIF
!
	  ELSE
!
      
	  IF((KSTOPWAY(14).EQ.0).AND.(l.EQ.1))THEN
	  KSTOPWAY(14)=1
!
	  WRITE(99,*)
	  WRITE(99,'(" For chemical species = ",A30)')CHEMICAL(i)
	  WRITE(99,'(" Doseway14 was not calculed, because soil concentrations not exist ")')
      WRITE(99,*)
!
      PRINT*
	  PRINT*, '  FOR CHEMICAL SPECIES ', CHEMICAL(i)
      PRINT*, '  DOSEWAY14 WAS NOT CALCULED, BECAUSE SOIL CONCENTRATION NOT EXIST'
	  ENDIF
!
	  ENDIF
	  ENDIF
!
!*******************************************************************************************
!	  
!
      ENDDO	 
!
!      IF((l.eq.1).and.(m.eq.1))then
!	  write(*,*)adaf(1),adaf(2),adaf(3)
!	  endif
!
	  ENDDO	
!
!      IF((l.eq.1).and.(m.eq.1))then
!	  write(*,*)mutagenic(iii,1),mutagenic(iii,2),mutagenic(iii,3)
!	  endif
!	 
	  ENDDO		 ! FIM DO CICLO DOS METAIS CANCERIGENOS
!
      IF((NVP.EQ.1).AND.(NSTOPNVP(i,l).EQ.1).AND.(m.EQ.1))THEN
	  WRITE(17,'("________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________")')
      ENDIF
!
!
56    CONTINUE
!
      ENDDO   ! ACABA O CICLO de "i"
!
      ENDDO   ! ACABA O CICLO de "l"
!
     END SELECT
!
!
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
!
      aspas = '"'
!	  
      J_INI_AGE(1)=1
	  J_INI_AGE(2)=2
	  J_INI_AGE(3)=3
	  J_INI_AGE(4)=6
	  J_INI_AGE(5)=11
	  J_INI_AGE(6)=16
	  J_INI_AGE(7)=18
	  J_INI_AGE(8)=21
	  J_INI_AGE(9)=65
!
      lstop10=0
	  nstop10=0
	  istop10=0
	  kstop10=0
	  jstop10=0
!
      DO l=INICIO,NIDADE
      IF(ED_INI(SCENAR,l).NE.0.0)THEN
	  NTIMEexp=ED_INI(SCENAR,l)
      DO n=1,NVIAS
      IF(W(n).EQV..TRUE.)THEN
      DO i=1,NCHEM
      DO K=1,NLOCAL
	  DO j=1,NTIMEexp 
!
      lstop10=l
	  nstop10=n
	  istop10=i
	  kstop10=k
	  jstop10=j
!
	  ENDDO	  ! FIM DO j=1,NTIMEexp 
	  ENDDO	  ! FIM DO K=1,NLOCAL 
	  ENDDO	  ! FIM DO CICLO DE 'i' (POR METAL)     
      ENDIF   ! FIM DO IF PARA CADA WAY	 
      ENDDO   ! FIM DO CICLO "n"
      ENDIF   ! FIM do IF ED(SCENAR,l)....
      ENDDO   ! FIM DO CICLO "l"
!
!           PROTOCOLO DE IMPRESSÃO DA SAIDA DE DADOS
!
!      IF(KEY_ANALYZE.EQV..TRUE.)THEN
!
      CALL NOMINATION(CHEMICAL,NCHEM,NOMES)
!
!
!	  ABRE-SE O ARQUIVO DE SAIDA E 'Exposure.out'
	  OPEN(UNIT=33,FILE='Results\Doses, HQ and CR.json')
!
      nnc=NLOCAL
!
	  DO KK=1,NLOCAL
	  LOCALSP(KK)=KK
	  ENDDO
!
!
	  WRITE(33,'("{")')
!
	  WRITE(33,'(A1,"Non-carcinogenic doses and HQ values",A1,": [")')aspas,aspas
!
      DO l=INICIO,NIDADE
!
      IF(ED_INI(SCENAR,l).NE.0.0)THEN
!
!------------------
	  NTIMEexp=ED_INI(SCENAR,l)
!------------------ 
!
      DO n=1,NVIAS
!
      IF(W(n).EQV..TRUE.)THEN
!
!
	  KCHEM=NCHEM
      NPOL=NPOL
!
      DO i=1,KCHEM
!
!
      DO K=1,NLOCAL
!
	  DO j=1,NTIMEexp 
!
      IF(l.EQ.1)THEN
	  N_AGE=j
      ELSEIF(l.EQ.2)THEN
	  N_AGE=j+1
      ELSEIF(l.EQ.3)THEN
	  N_AGE=j+2
      ELSEIF(l.EQ.4)THEN
	  N_AGE=j+5
      ELSEIF(l.EQ.5)THEN
	  N_AGE=j+10
      ELSEIF(l.EQ.6)THEN
	  N_AGE=j+15
      ELSEIF(l.EQ.7)THEN
	  N_AGE=j+17
      ELSEIF(l.EQ.8)THEN
	  N_AGE=j+20
      ELSEIF(l.EQ.9)THEN
	  N_AGE=j+64
	  ENDIF
!
!
      WRITE(33,'("{")')
!
      write(33,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(33,'(A1,"Pathway",A1,":",1x,I2,",") )') aspas,aspas,n
      write(33,'(A1,"Chemical species",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NOMES(i),aspas
	  write(33,'(A1,"Exposure type",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NRISKTYPE(n),aspas
      write(33,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(33,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(33,'(A1,"Dose value",A1,":",1x,ES12.5,",") )') aspas,aspas,ADD(l,1,n,i,j,k)
      write(33,'(A1,"Dose error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ADD(l,1,n,i,j,k)
      write(33,'(A1,"HQ value",A1,":",1x,ES12.5,",") )') aspas,aspas,HQ(l,n,i,j,k)
      write(33,'(A1,"HQ error",A1,":",1x,ES12.5) )') aspas,aspas,SD_HQ(l,n,i,j,k)
!
!
      IF((j.EQ.jstop10).and.(K.EQ.kstop10).and.(i.EQ.istop10).and.(l.EQ.lstop10).and.(n.EQ.nstop10))THEN
	  WRITE(33,'("}")')
	  ELSE
	  WRITE(33,'("},")')
	  ENDIF
!
	  ENDDO	  ! FIM DO j=1,NTIMEexp 
!
	  ENDDO	  ! FIM DO K=1,NLOCAL 
!
	  ENDDO	  ! FIM DO CICLO DE 'i' (POR METAL)     
!
      ENDIF   ! FIM DO IF PARA CADA WAY	 
!
      ENDDO   ! FIM DO CICLO "n"
!
      ENDIF   ! FIM do IF ED(SCENAR,l)....
!
      ENDDO   ! FIM DO CICLO "l"	
!
	  WRITE(33,'("],")')
!
!
	  WRITE(33,'(A1,"Carcinogenic doses and CR values",A1,": [")')aspas,aspas
!
      DO l=INICIO,NIDADE
!
      IF(ED_INI(SCENAR,l).NE.0.0)THEN
!
!------------------
	  NTIMEexp=ED_INI(SCENAR,l)
!------------------ 
!
      DO n=1,NVIAS
!
      IF(W(n).EQV..TRUE.)THEN
!
!
	  KCHEM=NCHEM
      NPOL=NPOL
!
      DO i=1,KCHEM
!
!
      DO K=1,NLOCAL
!
	  DO j=1,NTIMEexp 
!
      IF(l.EQ.1)THEN
	  N_AGE=j
      ELSEIF(l.EQ.2)THEN
	  N_AGE=j+1
      ELSEIF(l.EQ.3)THEN
	  N_AGE=j+2
      ELSEIF(l.EQ.4)THEN
	  N_AGE=j+5
      ELSEIF(l.EQ.5)THEN
	  N_AGE=j+10
      ELSEIF(l.EQ.6)THEN
	  N_AGE=j+15
      ELSEIF(l.EQ.7)THEN
	  N_AGE=j+17
      ELSEIF(l.EQ.8)THEN
	  N_AGE=j+20
      ELSEIF(l.EQ.9)THEN
	  N_AGE=j+64
	  ENDIF
!
!
      WRITE(33,'("{")')
!
      write(33,'(A1,"Initial age",A1,":",1x,I2,",") )') aspas,aspas,J_INI_AGE(l)
      write(33,'(A1,"Pathway",A1,":",1x,I2,",") )') aspas,aspas,n
      write(33,'(A1,"Chemical species",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NOMES(i),aspas
	  write(33,'(A1,"Exposure type",A1,":",1x,A1,A10,A1,",") )') aspas,aspas,aspas,NRISKTYPE(n),aspas
      write(33,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(33,'(A1,"Age",A1,":",1x,I3,",") )') aspas,aspas,N_AGE
      write(33,'(A1,"Dose value",A1,":",1x,ES12.5,",") )') aspas,aspas,ADD(l,2,n,i,j,k)
      write(33,'(A1,"Dose error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ADD(l,2,n,i,j,k)
      write(33,'(A1,"CR value",A1,":",1x,ES12.5,",") )') aspas,aspas,CR(l,n,i,j,k)
      write(33,'(A1,"CR error",A1,":",1x,ES12.5) )') aspas,aspas,SD_CR(l,n,i,j,k)
!
!
      IF((j.EQ.jstop10).and.(K.EQ.kstop10).and.(i.EQ.istop10).and.(l.EQ.lstop10).and.(n.EQ.nstop10))THEN
	  WRITE(33,'("}")')
	  ELSE
	  WRITE(33,'("},")')
	  ENDIF
!
	  ENDDO	  ! FIM DO j=1,NTIMEexp 
!
	  ENDDO	  ! FIM DO K=1,NLOCAL 
!
	  ENDDO	  ! FIM DO CICLO DE 'i' (POR METAL)     
!
      ENDIF   ! FIM DO IF PARA CADA WAY	 
!
      ENDDO   ! FIM DO CICLO "n"
!
      ENDIF   ! FIM do IF ED(SCENAR,l)....
!
      ENDDO   ! FIM DO CICLO "l"	
!
	  WRITE(33,'("]")')
!
!
	  WRITE(33,'("}")')
!
!
      RETURN	! FIM DA SUB-ROTINA EXPOSURE
	  END		
!
!
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
!
!
!
!*******************************************************************************
!*******************************************************************************
!
      SUBROUTINE DOSEWAY14(NCHEM,NTIME,NLOCAL,KSTOPABS14,CHEMICAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF,SD_EF,ED,SD_ED,AT,SD_AT,ABS,SD_ABS,AFsoil,SD_AFsoil,&
	  SA,SD_SA,CF,SD_CF,EV,SD_EV,AD,SD_AD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),SD_CSOIL(NCHEM,NTIME,NLOCAL)
!
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
!
!     PROTEÇÃO
!
    IF(ABS.EQ.0.0)THEN
	ABS_SUB=1.0
	IF(KSTOPABS14.EQ.0)THEN
	KSTOPABS14=1
	WRITE(99,*)
    WRITE(99,'(" WAY 14!!!! ABS VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	ABS_SUB=ABS
	ENDIF
!
	  IF(CSOIL(i,j,k).EQ.0.0)THEN
	  SOIL=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  SOIL=CSOIL(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÃO
!*****************************************************************************************************************************************************************     
!
!
     AD=(CSOIL(i,j,k)*CF*SA*AFsoil*EF*ED*EV*ABS)/(BW*AT)   !DOSE USEPA			   OK!!!!!!!!!!!!!!!!!!
! 
!
     AD2=AD
!
!
     SD_AD=AD2*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_CF/CF)**2+(SD_SA/SA)**2+(SD_ABS/ABS_SUB)**2+(SD_AFsoil/AFsoil)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_EV/EV)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
      RETURN
      END
!!
!*****************************************************************************
!****************************************************************************
!
!
      SUBROUTINE DOSEWAY11(NCHEM,NTIME,NLOCAL,i,j,k,EF,SD_EF,ED,SD_ED,AT,SD_AT,CPAR,SD_CPAR,&
	  ET,SD_ET,AD,SD_AD)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION CPAR(NCHEM,NTIME,NLOCAL),SD_CPAR(NCHEM,NTIME,NLOCAL)
!
!     PROTEÇÃO
!
	  IF(CPAR(i,j,k).EQ.0.0)THEN
	  PAR=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  PAR=CPAR(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÃO
!********************************************************************************************************************************
!
!    FORMULA USADA POR USEPA(1989): vol I PART F
!
!
    AD=(CPAR(i,j,k)*EF*ED*ET)/(AT*24)	! Antes de entrar o AT é multiplicado por ED e o 24 representa o as 24h do dia
!											! Nessa equação AT representa as horas de exposição ao ar durante o anos de exposição
    AD2=AD
!
    SD_AD=AD2*SQRT((SD_CPAR(i,j,k)/PAR)**2+(SD_ET/ET)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_AT/AT)**2)
!
      RETURN
	END
!
!
!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY12(NCHEM,NTIME,NLOCAL,i,j,k,EF,SD_EF,ED,SD_ED,AT,SD_AT,CSTEAM,SD_CSTEAM,&
	  ET,SD_ET,AD,SD_AD)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION CSTEAM(NCHEM,NTIME,NLOCAL),SD_CSTEAM(NCHEM,NTIME,NLOCAL)
!
!     PROTEÇÃO
!
	  IF(CSTEAM(i,j,k).EQ.0.0)THEN
	  STEAM=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  STEAM=CSTEAM(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÃO
!********************************************************************************************************************************
!
!    FORMULA USADA POR USEPA(1989): vol I PART F
!
    AD=(CSTEAM(i,j,k)*EF*ED*ET)/(AT*24)	! Antes de entrar o AT é multiplicado por ED e o 24 representa o as 24h do dia
!											! Nessa equação AT representa as horas de exposição ao ar durante o anos de exposição
    AD2=AD
!
    SD_AD=AD2*SQRT((SD_CSTEAM(i,j,k)/STEAM)**2+(SD_ET/ET)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_AT/AT)**2)
!
      RETURN
	END
!
!
!*********************************************************************************
!*********************************************************************************
!
!
!
      SUBROUTINE DOSEWAY1(NCHEM,NTIME,NLOCAL,i,j,k,CSOIL,SD_CSOIL,BW,SD_BW,EF,SD_EF,ED,SD_ED,AT,SD_AT,R_IR,SD_R_IR,&
	  FI,SD_FI,CF,SD_CF,AD,SD_AD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),SD_CSOIL(NCHEM,NTIME,NLOCAL)
!
!     PROTEÇÃO
!
      IF(CSOIL(i,j,k).EQ.0.0)THEN
	  SOIL=1.0
	  ELSEIF(CSOIL(i,j,k).NE.0.0)THEN
	  SOIL=CSOIL(i,j,k)
	  ENDIF
!
!
!     FIM PROTEÇÃO
!*****************************************************************************************************************************************************
      AD=(CSOIL(i,j,k)*R_IR*FI*EF*ED*CF)/(BW*AT)  !DOSE USEPA 1989
!
      AD2=AD
!
      SD_AD=AD2*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_R_IR/R_IR)**2+(SD_CF/CF)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_BW/BW)**2+(SD_AT/AT)**2+(SD_ED/ED)**2)
!
      RETURN
	  END 
!
!
!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY2(CHEMICAL,KSTOPBTF1,NCHEM,NTIME,NLOCAL,KEYCONC,i,j,k,CSOIL,SD_CSOIL,CFRUIT,SD_CFRUIT,EF,SD_EF,ED,&
	  SD_ED,AT,SD_AT,BW,SD_BW,BTF1,SD_BTF1,R_IRfruit,SD_R_IRfruit,FI,SD_FI,AD,SD_AD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL KEYCONC
!
	  CHARACTER(LEN=50)  :: CHEMICAL(500)	
!
      DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),CFRUIT(NCHEM,NTIME,NLOCAL)
	        DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL),SD_CFRUIT(NCHEM,NTIME,NLOCAL) 
!
! 	 PROTEÇÃO
!
	  IF(CSOIL(i,j,k).EQ.0.0)THEN
	  SOIL=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  SOIL=CSOIL(i,j,k)
	  ENDIF
!
	  IF(CFRUIT(i,j,k).EQ.0.0)THEN
	  FRUIT=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  FRUIT=CFRUIT(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÃO
!*****************************************************************************************************************************************************
!	CALCULANDO O Cvege USANDO-SE O BTF
!
    IF(KEYCONC.EQV..FALSE.) THEN
	CFRUIT(i,j,k)=CSOIL(i,j,k)*BTF1*(1.00-0.85)
!
	  IF(BTF1.EQ.0.0)THEN
	  BTF1_SUB=1.0			! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  IF(KSTOPBTF1.EQ.0)THEN
	  KSTOPBTF1=1
	  WRITE(99,*)
	  WRITE(99,'(" WAY 2!!!! BTF VALUE (Soil-Fruit) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE
	  BTF1_SUB=BTF1
	  ENDIF
!
    IF((SD_CSOIL(i,j,k).NE.0.0).OR.(SD_BTF1.NE.0.0))THEN
	SD_CFRUIT(i,j,k)=CFRUIT(i,j,k)* SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF1/BTF1_SUB)**2)
	ELSE
	SD_CFRUIT(i,j,k)=0.0
	ENDIF
!
	  IF(CFRUIT(i,j,k).EQ.0.0)THEN
	  FRUIT=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  FRUIT=CFRUIT(i,j,k)
	  ENDIF
!
	ENDIF
!*****************************************************************************************************************************************************
!
!
      AD=(CFRUIT(i,j,k)*R_IRfruit*FI*EF*ED)/(BW*AT)  	
!
      AD2=AD	 
!
      SD_AD=AD2*SQRT((SD_CFRUIT(i,j,k)/FRUIT)**2+(SD_R_IRfruit/R_IRfruit)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
!
      RETURN
	  END
!
!
!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY5(NCHEM,NTIME,NLOCAL,i,j,k,CWATER,SD_CWATER,BW,SD_BW,EF,SD_EF,ED,SD_ED,AT,SD_AT,R_IR,SD_IR,AD,SD_AD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION CWATER(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL)
!
!	  PROTEÇÃO
!
	  IF(CWATER(i,j,k).EQ.0.0)THEN
	  WATER=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  WATER=CWATER(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÃO
!************************************************************************************************************************************
!
      AD=(CWATER(i,j,k)*R_IR*EF*ED)/(BW*AT) 
!
      AD2=AD
!
      SD_AD=AD2*SQRT((SD_CWATER(i,j,k)/WATER)**2+(SD_IR/R_IR)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2) 		 			 	
!
      RETURN
	  END
!
!******************************************************************************
!******************************************************************************
!
    SUBROUTINE DOSEWAY6(CHEMICAL,KSTOPBTF_6S,NCHEM,NTIME,NLOCAL,KEYCONC,i,j,k,CSOIL,SD_CSOIL,CLEAVES,SD_CLEAVES,EF,SD_EF,ED,SD_ED,&
	AT,SD_AT,BW,SD_BW,BTF_6S,SD_BTF_6S,R_IRveg,SD_IRveg,FI,SD_FI,AD,SD_AD)

    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	LOGICAL KEYCONC
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
    DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),CLEAVES(NCHEM,NTIME,NLOCAL)
    DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL),SD_CLEAVES(NCHEM,NTIME,NLOCAL) 
!
!   PROTEÇÃO
!
    IF(CSOIL(i,j,k).EQ.0.0)THEN
	SOIL=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SOIL=CSOIL(i,j,k)
	ENDIF
    IF(BTF_6S.EQ.0.0)THEN
	BTF_6S_SUB=1.0
	IF((KSTOPBTF_6S.EQ.0).AND.(KEYCONC.EQV..FALSE.))THEN
	KSTOPBTF_6S=1
	WRITE(99,*)
    WRITE(99,'(" WAY 6!!!! BTF VALUE (Soil-Leaves) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_6S_SUB=BTF_6S
	ENDIF
!
    IF(CLEAVES(i,j,k).EQ.0.0)THEN
	SLEAVES=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SLEAVES=CLEAVES(i,j,k)
	ENDIF
! 
!
!   FIM PROTEÇÃO
!*************************************************************************************************************************************************************
!	CALCULANDO O Cvege USANDO-SE O BTF
!
    IF(KEYCONC.EQV..FALSE.) THEN
!
	CLEAVES(i,j,k)=CSOIL(i,j,k)*BTF_6S*(1.00-0.85)
	IF((SD_CSOIL(i,j,k).NE.0.0).OR.(SD_BTF_6S.NE.0.0))THEN
	SD_CLEAVES(i,j,k)=CLEAVES(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF_6S/BTF_6S_SUB)**2)
	ELSE
	SD_CLEAVES(i,j,k)=0.0
	ENDIF
!
    IF(CLEAVES(i,j,k).EQ.0.0)THEN
	SLEAVES=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SLEAVES=CLEAVES(i,j,k)
	ENDIF
!
	ENDIF
!*************************************************************************************************************************************************************
!
    AD=(CLEAVES(i,j,k)*R_IRveg*FI*EF*ED)/(BW*AT)  !DOSE LIMA		 	
!
    AD2=AD
!
    SD_AD=AD2*SQRT((SD_CLEAVES(i,j,k)/SLEAVES)**2+(SD_IRveg/R_IRveg)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
    RETURN
	END
!
!*******************************************************************************
!******************************************************************************* 
!
!
      SUBROUTINE DOSEWAY13(NCHEM,NTIME,NLOCAL,KSTOPPC13,CHEMICAL,i,j,k,CWATER,SD_CWATER,BW,SD_BW,EV,SD_EV,EF,SD_EF,ED,SD_ED,ET,SD_ET,AT,SD_AT,PC,SD_PC,SA,SD_SA,CF,SD_CF,AD,SD_AD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION CWATER(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL)
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
!
!   PROTEÇÃO
!
    IF(PC.EQ.0.0)THEN
	PC_SUB=1.0
	IF(KSTOPPC13.EQ.0)THEN
	KSTOPPC13=1
	WRITE(99,*)
    WRITE(99,'(" WAY 13!!!! PC VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	PC_SUB=PC
	ENDIF
!
	  IF(CWATER(i,j,k).EQ.0.0)THEN
	  WATER=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  WATER=CWATER(i,j,k)
	  ENDIF
!
!   FIM PROTEÇÃO
!*****************************************************************************************************************************************************************
    AD=(CWATER(i,j,k)*SA*PC*ET*EV*EF*ED*CF)/(BW*AT)  !DOSE CETESB-LIMA			 		
!
    AD2=AD
!
!
    SD_AD=AD2*SQRT((SD_CWATER(i,j,k)/WATER)**2+(SD_SA/SA)**2+(SD_PC/PC_SUB)**2+(SD_CF/CF)**2+(SD_ET/ET)**2+(SD_EV/EV)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
    RETURN
	END
!

!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY3(CHEMICAL,KSTOPBTF_SB3,KSTOPBTF_SV3,KSTOPBTF_VB3,KSTOPBTF_WB3,KSTOPfw3,NCHEM,NTIME,NLOCAL,BTF_SV3,SD_BTF_SV3,BTF_VB3,SD_BTF_VB3,BTF_SB3,SD_BTF_SB3,BTF_WB3,SD_BTF_WB3,i,j,k,&
	  CSOIL,SD_CSOIL,CWATER,SD_CWATER,CBEEF,SD_CBEEF,EF,SD_EF,ED,SD_ED,AT,SD_AT,BW,SD_BW,&
	  KEYCONC_B3,KEYCONC_S3, KEYCONC_W3,Fa,SD_Fa,Fp,SD_Fp,R_IRali,SD_IRali,R_IRsolo,SD_IRsolo,R_IRagua,SD_IRagua,R_IRcarne,SD_IRcarne,&
	  FI,SD_FI,fw,SD_fw,AD,SD_AD)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
	  LOGICAL KEYCONC_B3,KEYCONC_S3,KEYCONC_W3	
      DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),CWATER(NCHEM,NTIME,NLOCAL)
	  DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL)
	  DIMENSION CBEEF1(NCHEM,NTIME,NLOCAL),CBEEF2(NCHEM,NTIME,NLOCAL),CBEEF3(NCHEM,NTIME,NLOCAL),CBEEF(NCHEM,NTIME,NLOCAL)
	  DIMENSION SD_CBEEF1(NCHEM,NTIME,NLOCAL),SD_CBEEF2(NCHEM,NTIME,NLOCAL),SD_CBEEF3(NCHEM,NTIME,NLOCAL),SD_CBEEF(NCHEM,NTIME,NLOCAL)	
!
!	CALCULANDO O Cbeef USANDO-SE O BTF
!
!
!   PROTEÇÕES
	IF(CSOIL(i,j,k).EQ.0.0)THEN
	SOIL=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SOIL=CSOIL(i,j,k)
	ENDIF
    IF(BTF_SB3.EQ.0.0)THEN
	BTF_SB3_SUB=1.0
	IF((KSTOPBTF_SB3.EQ.0).AND.(KEYCONC_B3.EQV..FALSE.).AND.(KEYCONC_S3.EQV..TRUE.).AND.(KEYCONC_W3.EQV..TRUE.))THEN
	KSTOPBTF_SB3=1
	WRITE(99,*)
    WRITE(99,'(" WAY 3!!!! BTF VALUE (Soil-Meat cattle) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_SB3_SUB=BTF_SB3
	ENDIF
!
    IF(BTF_SV3.EQ.0.0)THEN
	BTF_SV3_SUB=1.0
	IF((KSTOPBTF_SV3.EQ.0).AND.(KEYCONC_B3.EQV..FALSE.).AND.(KEYCONC_S3.EQV..TRUE.).AND.(KEYCONC_W3.EQV..TRUE.))THEN
	KSTOPBTF_SV3=1
	WRITE(99,*)
    WRITE(99,'(" WAY 3!!!! BTF VALUE (Soil-Feed plants) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_SV3_SUB=BTF_SV3
	ENDIF
	IF(BTF_VB3.EQ.0.0)THEN
	BTF_VB3_SUB=1.0
	IF((KSTOPBTF_VB3.EQ.0).AND.(KEYCONC_B3.EQV..FALSE.).AND.(KEYCONC_S3.EQV..TRUE.).AND.(KEYCONC_W3.EQV..TRUE.))THEN
	KSTOPBTF_VB3=1
	WRITE(99,*)
    WRITE(99,'(" WAY 3!!!! BTF VALUE (Feed plants-Meat cattle) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_VB3_SUB=BTF_VB3
	ENDIF
!
	IF(CWATER(i,j,k).EQ.0.0)THEN
	WATER=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	WATER=CWATER(i,j,k)
	ENDIF
    IF(BTF_WB3.EQ.0.0)THEN
	BTF_WB3_SUB=1.0
	IF((KSTOPBTF_WB3.EQ.0).AND.(KEYCONC_B3.EQV..FALSE.).AND.(KEYCONC_S3.EQV..TRUE.).AND.(KEYCONC_W3.EQV..TRUE.))THEN
	KSTOPBTF_WB3=1
	WRITE(99,*)
	WRITE(99,'(" WAY 3!!!! BTF VALUE (Water-Beef cattle) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_WB3_SUB=BTF_WB3
	ENDIF
!
    IF(fw.EQ.0.0)THEN
	fw_SUB=1.0
	IF(KSTOPfw3.EQ.0)THEN
	KSTOPfw3=1
	WRITE(99,*)
	WRITE(99,'(" WAY 3!!!! fw VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	fw_SUB=fw
	ENDIF
!
	IF(CBEEF(i,j,k).EQ.0.0)THEN
	BEEF=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BEEF=CBEEF(i,j,k)
	ENDIF
!
!   FIM DAS PROTEÇÕES
!
!    WRITE(*,*)CSOIL(i,j,k),BTF_SB9,R_IRsolo
!	WRITE(*,*)Fa,Fp
!
!    PAUSE
!
!***************************************************************************************************************************************************************************************
    IF(KEYCONC_B3.EQV..FALSE.) THEN
	CBEEF1(i,j,k)=CSOIL(i,j,k)*BTF_SB3*R_IRsolo*Fa*Fp 
	SD_CBEEF1(i,j,k)=	CBEEF1(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF_SB3/BTF_SB3_SUB)**2+(SD_IRsolo/R_IRsolo)**2+(SD_Fa/Fa)**2+(SD_Fp/Fp)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	CBEEF2(i,j,k)=CSOIL(i,j,k)*BTF_SV3*(1.00-0.85)*BTF_VB3*R_IRali*Fa*Fp
	SD_CBEEF2(i,j,k)=	CBEEF2(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF_SV3/BTF_SV3_SUB)**2+(SD_BTF_VB3/BTF_VB3_SUB)**2+(SD_IRali/R_IRali)**2+(SD_Fa/Fa)**2+(SD_Fp/Fp)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	CBEEF3(i,j,k)=CWATER(i,j,k)*R_IRagua*fw*BTF_WB3
	SD_CBEEF3(i,j,k)=	CBEEF3(i,j,k)*SQRT((SD_CWATER(i,j,k)/WATER)**2+(SD_BTF_WB3/BTF_WB3_SUB)**2+(SD_IRagua/R_IRagua)**2+(SD_fw/fw_SUB)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    CBEEF(i,j,k)=CBEEF1(i,j,k)+CBEEF2(i,j,k)+CBEEF3(i,j,k)
	IF((SD_CBEEF1(i,j,k).NE.0.0).OR.(SD_CBEEF2(i,j,k).NE.0.0).OR.(SD_CBEEF3(i,j,k).NE.0.0))THEN
	SD_CBEEF(i,j,k)=SQRT(SD_CBEEF1(i,j,k)**2+SD_CBEEF2(i,j,k)**2+SD_CBEEF3(i,j,k)**2)
	ELSE
	SD_CBEEF(i,j,k)=0.0
	ENDIF
!
	IF(CBEEF(i,j,k).EQ.0.0)THEN
	BEEF=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BEEF=CBEEF(i,j,k)
	ENDIF
!
!***************************************************************************************************************************************************************************************
	ELSE IF(KEYCONC_B3.EQV..TRUE.)THEN
!
	ENDIF
!***************************************************************************************************************************************************************************************
!
!    WRITE(*,*)CBEEF1(i,j,k)
!	PAUSE
!	WRITE(*,*)CBEEF2(i,j,k),CBEEF3(i,j,k)
!
!
    AD=(CBEEF(i,j,k)*R_IRcarne*FI*EF*ED)/(BW*AT) 		 	
!
    AD2=AD
!
    SD_AD=AD2*SQRT((SD_CBEEF(i,j,k)/BEEF)**2+(SD_IRcarne/R_IRcarne)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)		 	
!
    RETURN
	END
!
!
!!
!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY4(CHEMICAL,KSTOPBTF_SM4,KSTOPBTF_SV4,KSTOPBTF_VM4,KSTOPBTF_WM4,KSTOPfw4,NCHEM,NTIME,NLOCAL,BTF_SV4,SD_BTF_SV4,BTF_VM4,SD_BTF_VM4,BTF_SM4,SD_BTF_SM4,BTF_WM4,SD_BTF_WM4,i,j,k,&
	  CSOIL,SD_CSOIL,CWATER,SD_CWATER,CMILK,SD_CMILK,EF,SD_EF,ED,SD_ED,AT,SD_AT,BW,SD_BW,KEYCONC_M4,KEYCONC_S4,KEYCONC_W4,&
	  Fa,SD_Fa,Fp,SD_Fp,R_IRali,SD_IRali,R_IRsolo,SD_IRsolo,R_IRagua,SD_IRagua,R_IRleite,SD_IRleite,FI,SD_FI,fw,SD_fw,AD,SD_AD)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL KEYCONC_M4,KEYCONC_S4, KEYCONC_W4
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
      DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),CWATER(NCHEM,NTIME,NLOCAL)
	  DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL)
	  DIMENSION CMILK1(NCHEM,NTIME,NLOCAL),CMILK2(NCHEM,NTIME,NLOCAL),CMILK3(NCHEM,NTIME,NLOCAL),CMILK(NCHEM,NTIME,NLOCAL)
	  DIMENSION SD_CMILK1(NCHEM,NTIME,NLOCAL),SD_CMILK2(NCHEM,NTIME,NLOCAL),SD_CMILK3(NCHEM,NTIME,NLOCAL),SD_CMILK(NCHEM,NTIME,NLOCAL)
!
!	CALCULANDO O Cbeef USANDO-SE O BTF
!
!   PROTEÇÕES
    IF(CSOIL(i,j,k).EQ.0.0)THEN
	SOIL=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SOIL=CSOIL(i,j,k)
	ENDIF
    IF(BTF_SM4.EQ.0.0)THEN
	BTF_SM4_SUB=1.0
	IF((KSTOPBTF_SM4.EQ.0).AND.(KEYCONC_M4.EQV..FALSE.).AND.(KEYCONC_S4.EQV..TRUE.).AND.(KEYCONC_W4.EQV..TRUE.))THEN
	KSTOPBTF_SM4=1
	WRITE(99,*)
    WRITE(99,'(" WAY 4!!!! BTF VALUE (Soil-Milk bovine) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_SM4_SUB=BTF_SM4
	ENDIF
!
 	IF(BTF_SV4.EQ.0.0)THEN
	BTF_SV4_SUB=1.0
	IF((KSTOPBTF_SV4.EQ.0).AND.(KEYCONC_M4.EQV..FALSE.).AND.(KEYCONC_S4.EQV..TRUE.).AND.(KEYCONC_W4.EQV..TRUE.))THEN
	KSTOPBTF_SV4=1
	WRITE(99,*)
    WRITE(99,'(" WAY 4!!!! BTF VALUE (Soil-Feed plants) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_SV4_SUB=BTF_SV4
	ENDIF
    IF(BTF_VM4.EQ.0.0)THEN
	BTF_VM4_SUB=1.0
	IF((KSTOPBTF_VM4.EQ.0).AND.(KEYCONC_M4.EQV..FALSE.).AND.(KEYCONC_S4.EQV..TRUE.).AND.(KEYCONC_W4.EQV..TRUE.))THEN
	KSTOPBTF_VM4=1
	WRITE(99,*)
    WRITE(99,'(" WAY 4!!!! BTF VALUE (Feed plants-Milk bovine) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_VM4_SUB=BTF_VM4
	ENDIF
!
	IF(CWATER(i,j,k).EQ.0.0)THEN
	WATER=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	WATER=CWATER(i,j,k)
	ENDIF
    IF(BTF_WM4.EQ.0.0)THEN
	BTF_WM4_SUB=1.0
	IF((KSTOPBTF_WM4.EQ.0).AND.(KEYCONC_M4.EQV..FALSE.).AND.(KEYCONC_S4.EQV..TRUE.).AND.(KEYCONC_W4.EQV..TRUE.))THEN
	KSTOPBTF_WM4=1
	WRITE(99,*)
    WRITE(99,'(" WAY 4!!!! BTF VALUE (Water-Milk bovine) USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	BTF_WM4_SUB=BTF_WM4
	ENDIF
!
    IF(fw.EQ.0.0)THEN
	fw_SUB=1.0
	IF(KSTOPfw4.EQ.0)THEN
	KSTOPfw4=1
	WRITE(99,*)
    WRITE(99,'(" WAY 4!!!! fw VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	fw_SUB=fw
	ENDIF
!
!
	IF(CMILK(i,j,k).EQ.0.0)THEN
	SMILK=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SMILK=CMILK(i,j,k)
	ENDIF
!
!   FIM DAS PROTEÇÕES
!
!***************************************************************************************************************************************************************************************
    IF(KEYCONC_M4.EQV..FALSE.) THEN
!
	CMILK1(i,j,k)=CSOIL(i,j,k)*BTF_SM4*R_IRsolo*Fa*Fp
	SD_CMILK1(i,j,k)=	CMILK1(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF_SM4/BTF_SM4_SUB)**2+(SD_IRsolo/R_IRsolo)**2+(SD_Fa/Fa)**2+(SD_Fp/Fp)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	CMILK2(i,j,k)=CSOIL(i,j,k)*BTF_SV4*(1.00-0.85)*BTF_VM4*R_IRali*Fa*Fp
	SD_CMILK2(i,j,k)=	CMILK2(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF_SV4/BTF_SV4_SUB)**2+(SD_BTF_VM4/BTF_VM4_SUB)**2+(SD_IRali/R_IRali)**2+(SD_Fa/Fa)**2+(SD_Fp/Fp)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	CMILK3(i,j,k)=CWATER(i,j,k)*R_IRagua*fw*BTF_WM4
	SD_CMILK3(i,j,k)=	CMILK3(i,j,k)*SQRT((SD_CWATER(i,j,k)/WATER)**2+(SD_BTF_WM4/BTF_WM4_SUB)**2+(SD_IRagua/R_IRagua)**2+(SD_fw/fw_SUB)**2)
!
	CMILK(i,j,k)=CMILK1(i,j,k)+CMILK2(i,j,k)+CMILK3(i,j,k)
!
    IF((SD_CMILK1(i,j,k).NE.0.0).OR.(SD_CMILK2(i,j,k).NE.0.0).OR.(SD_CMILK3(i,j,k).NE.0.0))THEN
	SD_CMILK(i,j,k)=SQRT(SD_CMILK1(i,j,k)**2+SD_CMILK2(i,j,k)**2+SD_CMILK3(i,j,k)**2)
	ELSE
	SD_CMILK(i,j,k)=0.0
	ENDIF
!
	IF(CMILK(i,j,k).EQ.0.0)THEN
	SMILK=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	SMILK=CMILK(i,j,k)
	ENDIF
!
!***************************************************************************************************************************************************************************************
	ELSE IF(KEYCONC_M4.EQV..TRUE.)THEN
!
	ENDIF
!***************************************************************************************************************************************************************************************
!
      AD=(CMILK(i,j,k)*R_IRleite*FI*EF*ED)/(BW*AT)  !DOSE LIMA	
!
      AD2=AD 	
!
      SD_AD=AD2*SQRT((SD_CMILK(i,j,k)/SMILK)**2+(SD_IRleite/R_IRleite)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)		 	
!
      RETURN
	END
!
!
!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY7(CHEMICAL,KSTOPBTF_WF7,NCHEM,NTIME,NLOCAL,BTF_WF7,SD_BTF_WF7,i,j,k,CFISH,SD_CFISH,EF,SD_EF,ED,SD_ED,&
	  AT,SD_AT,BW,SD_BW,KEYCONC_F7,R_IRfish,SD_IRfish,FI,SD_FI,CWATER,SD_CWATER,AD,SD_AD)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
	  LOGICAL KEYCONC_F7	
      DIMENSION CFISH(NCHEM,NTIME,NLOCAL),SD_CFISH(NCHEM,NTIME,NLOCAL)
!
!	CALCULANDO O Cbeef USANDO-SE O BTF
!
!   PROTEÇÕES
!
	IF(CWATER.EQ.0.0)THEN
	WATER=1.0
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	WATER=CWATER
	ENDIF
!
	  IF(BTF_WF7.EQ.0.0)THEN
	  BTF_WF7_SUB=1.0
	  IF((KSTOPBTF_WF7.EQ.0).AND.(KEYCONC_F7.EQV..FALSE.))THEN
	  KSTOPBTF_WF7=1
	  WRITE(99,*)
      WRITE(99,'(" WAY 7!!!! BTF VALUE (Water-Fish) USED IS 0.0!!!! CALCULATED DOSE VALUE WILL BE 0.0 TOO!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  BTF_WF7_SUB=BTF_WF7
	  ENDIF
!
	  IF(CFISH(i,j,k).EQ.0.0)THEN
	  FISH=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  FISH=CFISH(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÕES
!----------------------------------------------------------------------------------------------------------------------------------------------------------------
      IF(KEYCONC_F7.EQV..FALSE.) THEN
	  CFISH(i,j,k)=CWATER*BTF_WF7
!
      IF((SD_CWATER.NE.0.0).OR.(SD_BTF_WF7.NE.0.0))THEN
	  SD_CFISH(i,j,k)=CFISH(i,j,k)*SQRT((SD_CWATER/WATER)**2+(SD_BTF_WF7/BTF_WF7_SUB)**2)
	  ELSE
	  SD_CFISH(i,j,k)=0.0
	  ENDIF
!
	  IF(CFISH(i,j,k).EQ.0.0)THEN
	  FISH=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  FISH=CFISH(i,j,k)
	  ENDIF
!
	  ENDIF
!----------------------------------------------------------------------------------------------------------------------------------------------------------------
!
!
      AD=(CFISH(i,j,k)*R_IRfish*FI*EF*ED)/(BW*AT)		 	
!
!
      AD2=AD
!
      SD_AD=AD2*SQRT((SD_CFISH(i,j,k)/FISH)**2+(SD_IRfish/R_IRfish)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
      RETURN
	  END
!
!
!*********************************************************************************
!*********************************************************************************
!
!
      SUBROUTINE DOSEWAY8(CHEMICAL,KSTOPBTF_WAVE8,KSTOPBTF_SAVE8,KSTOPfw8,NCHEM,NTIME,NLOCAL,BTF_SAVE8,SD_BTF_SAVE8,BTF_WAVE8,SD_BTF_WAVE8,i,j,k,CSOIL,SD_CSOIL,CWATER,SD_CWATER,&
	  CAVE,SD_CAVE,EF,SD_EF,ED,SD_ED,AT,SD_AT,BW,SD_BW,KEYCONC_AVE8,Fa,SD_Fa,Fp,SD_Fp,R_IRsolo,SD_IRsolo,R_IRagua,SD_IRagua,&
	  R_IRave,SD_IRave,FI,SD_FI,fw,SD_fw,AD,SD_AD)
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL KEYCONC_AVE8
	  CHARACTER(LEN=50)  :: CHEMICAL(500)	
      DIMENSION CWATER(NCHEM,NTIME,NLOCAL),CSOIL(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL),SD_CSOIL(NCHEM,NTIME,NLOCAL)
	  DIMENSION CAVE1(NCHEM,NTIME,NLOCAL),CAVE2(NCHEM,NTIME,NLOCAL),CAVE(NCHEM,NTIME,NLOCAL)
	  DIMENSION SD_CAVE1(NCHEM,NTIME,NLOCAL),SD_CAVE2(NCHEM,NTIME,NLOCAL),SD_CAVE(NCHEM,NTIME,NLOCAL)
!
!     PROTEÇÕES
!
	  IF(CWATER(i,j,k).EQ.0.0)THEN
	  WATER=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  WATER=CWATER(i,j,k)
	  ENDIF
	  IF(BTF_WAVE8.EQ.0.0)THEN
	  BTF_WAVE8_SUB=1.0
	  IF((KSTOPBTF_WAVE8.EQ.0).AND.(KEYCONC_AVE8.EQV..FALSE.))THEN
	  KSTOPBTF_WAVE8=1
	  WRITE(99,*)
      WRITE(99,'(" WAY 8!!!! BTF (Water-Bird) VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  BTF_WAVE8_SUB=BTF_WAVE8
	  ENDIF
!
	  IF(CSOIL(i,j,k).EQ.0.0)THEN
	  SOIL=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  SOIL=CSOIL(i,j,k)
	  ENDIF
	  IF(BTF_SAVE8.EQ.0.0)THEN
	  BTF_SAVE8_SUB=1.0
	  IF((KSTOPBTF_SAVE8.EQ.0).AND.(KEYCONC_AVE8.EQV..FALSE.))THEN
	  KSTOPBTF_SAVE8=1
	  WRITE(99,*)
      WRITE(99,'(" WAY 8!!!! BTF (Soil-Bird) VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  BTF_SAVE8_SUB=BTF_SAVE8
	  ENDIF
!
    IF(fw.EQ.0.0)THEN
	fw_SUB=1.0
	IF(KSTOPfw8.EQ.0)THEN
	KSTOPfw8=1
	WRITE(99,*)
    WRITE(99,'(" WAY 8!!!! fw VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	fw_SUB=fw
	ENDIF
!
	  IF(CAVE(i,j,k).EQ.0.0)THEN
	  AVE=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  AVE=CAVE(i,j,k)
	  ENDIF
!
!*************************************************************************************************************************************************************
      IF(KEYCONC_AVE8.EQV..FALSE.) THEN
	  CAVE1(i,j,k)=CWATER(i,j,k)*R_IRagua*fw*BTF_WAVE8
	  SD_CAVE1(i,j,k)=CAVE1(i,j,k)*SQRT((SD_CWATER(i,j,k)/WATER)**2+(SD_IRagua/R_IRagua)**2+(SD_fw/fw_SUB)**2+(SD_BTF_WAVE8/BTF_WAVE8_SUB)**2)
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
	  CAVE2(i,j,k)=CSOIL(i,j,k)*BTF_SAVE8*R_IRsolo*Fa*Fp
	  SD_CAVE2(i,j,k)=CAVE2(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_IRsolo/R_IRsolo)**2+(SD_Fa/Fa)**2+(SD_Fp/Fp)**2+(SD_BTF_SAVE8/BTF_SAVE8_SUB)**2)
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
	  CAVE(i,j,k)=CAVE1(i,j,k)+CAVE2(i,j,k)
	  IF((SD_CAVE1(i,j,k).NE.0.0).OR.(SD_CAVE2(i,j,k).NE.0.0))THEN
	  SD_CAVE(i,j,k)=SQRT(SD_CAVE1(i,j,k)**2+SD_CAVE2(i,j,k)**2)
	  ELSE
	  SD_CAVE(i,j,k)=0.0
	  ENDIF
!
	  IF(CAVE(i,j,k).EQ.0.0)THEN
	  AVE=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  AVE=CAVE(i,j,k)
	  ENDIF
!
	  ENDIF
!*************************************************************************************************************************************************************
!
!
      AD=(CAVE(i,j,k)*R_IRave*FI*EF*ED)/(BW*AT)
!
      AD2=AD
!
      SD_AD=AD2*SQRT((SD_CAVE(i,j,k)/AVE)**2+(SD_IRave/R_IRave)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
	  RETURN
	  END
!
!
!
!*********************************************************************************
!*********************************************************************************
!
      SUBROUTINE DOSEWAY9(CHEMICAL,KSTOPBTF_WEGG9,KSTOPBTF_SEGG9,KSTOPfw9,NCHEM,NTIME,NLOCAL,BTF_SEGG9,SD_BTF_SEGG9,BTF_WEGG9,SD_BTF_WEGG9,i,j,k,CSOIL,SD_CSOIL,CWATER,SD_CWATER,&
	  CEGG,SD_CEGG,EF,SD_EF,ED,SD_ED,AT,SD_AT,BW,SD_BW,KEYCONC_EGG9,Fa,SD_Fa,Fp,SD_Fp,R_IRsolo,SD_IRsolo,R_IRagua,SD_IRagua,&
	  R_IRegg,SD_IRegg,FI,SD_FI,fw,SD_fw,AD,SD_AD)
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL KEYCONC_EGG9
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
      DIMENSION CWATER(NCHEM,NTIME,NLOCAL),CSOIL(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL),SD_CSOIL(NCHEM,NTIME,NLOCAL)
	  DIMENSION CEGG1(NCHEM,NTIME,NLOCAL),CEGG2(NCHEM,NTIME,NLOCAL),CEGG(NCHEM,NTIME,NLOCAL)
	  DIMENSION SD_CEGG1(NCHEM,NTIME,NLOCAL),SD_CEGG2(NCHEM,NTIME,NLOCAL),SD_CEGG(NCHEM,NTIME,NLOCAL)
!
!     PROTEÇÃO
!
!
	  IF(CWATER(i,j,k).EQ.0.0)THEN
	  WATER=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  WATER=CWATER(i,j,k)
	  ENDIF
	  IF(BTF_WEGG9.EQ.0.0)THEN
	  BTF_WEGG9_SUB=1.0
	  IF((KSTOPBTF_WEGG9.EQ.0).AND.(KEYCONC_EGG9.EQV..FALSE.))THEN
	  KSTOPBTF_WEGG9=1
	  WRITE(99,*)
      WRITE(99,'(" WAY 9!!!! BTF (Water-Egg) VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  BTF_WEGG9_SUB=BTF_WEGG9
	  ENDIF
!
	  IF(CSOIL(i,j,k).EQ.0.0)THEN
	  SOIL=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  SOIL=CSOIL(i,j,k)
	  ENDIF
	  IF(BTF_SEGG9.EQ.0.0)THEN
	  BTF_SEGG9_SUB=1.0
	  IF((KSTOPBTF_SEGG9.EQ.0).AND.(KEYCONC_EGG9.EQV..FALSE.))THEN
	  KSTOPBTF_SEGG9=1
	  WRITE(99,*)
      WRITE(99,'(" WAY 9!!!! BTF (Soil-Egg) VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  BTF_SEGG9_SUB=BTF_SEGG9
	  ENDIF
!
    IF(fw.EQ.0.0)THEN
	fw_SUB=1.0
	IF(KSTOPfw9.EQ.0)THEN
	KSTOPfw9=1
	WRITE(99,*)
    WRITE(99,'(" WAY 9!!!! fw VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	ENDIF
	ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	fw_SUB=fw
	ENDIF
!
	  IF(CEGG(i,j,k).EQ.0.0)THEN
	  EGG=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  EGG=CEGG(i,j,k)
	  ENDIF
!
!
!     FIM PROTEÇÕES
!****************************************************************************************************************************************************************************
      IF(KEYCONC_EGG9.EQV..FALSE.) THEN
!
	  CEGG1(i,j,k)=CWATER(i,j,k)*R_IRagua*fw*BTF_WEGG9
	  SD_CEGG1(i,j,k)=CEGG1(i,j,k)*SQRT((SD_CWATER(i,j,k)/WATER)**2+(SD_IRagua/R_IRagua)**2+(SD_fw/fw_SUB)**2+(SD_BTF_WEGG9/BTF_WEGG9_SUB)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  CEGG2(i,j,k)=CSOIL(i,j,k)*BTF_SEGG9*R_IRsolo*Fa*Fp
	  SD_CEGG2(i,j,k)=CEGG2(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_IRsolo/R_IRsolo)**2+(SD_Fa/Fa)**2+(SD_Fp/Fp)**2+(SD_BTF_SEGG9/BTF_SEGG9_SUB)**2)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  CEGG(i,j,k)=CEGG1(i,j,k)+CEGG2(i,j,k)
	  IF((SD_CEGG1(i,j,k).NE.0.0).OR.(SD_CEGG2(i,j,k).NE.0.0))THEN
	  SD_CEGG(i,j,k)=SQRT(SD_CEGG1(i,j,k)**2+SD_CEGG2(i,j,k)**2)
	  ELSE
	  SD_CEGG(i,j,k)=0.0
	  ENDIF
!
	  IF(CEGG(i,j,k).EQ.0.0)THEN
	  EGG=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  EGG=CEGG(i,j,k)
	  ENDIF
!
	  ENDIF
!****************************************************************************************************************************************************************************
!
!      WRITE(*,*) CEGG1(i,j,k),CEGG2(i,j,k)
!
      AD=(CEGG(i,j,k)*R_IRegg*FI*EF*ED)/(BW*AT)
!
      AD2=AD
!
      SD_AD=AD2*SQRT((SD_CEGG(i,j,k)/EGG)**2+(SD_IRegg/R_IRegg)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
	  RETURN
	  END
!
!
!*********************************************************************************
!*********************************************************************************
!
      SUBROUTINE DOSEWAY10(CHEMICAL,KSTOPBTF_SG10,NCHEM,NTIME,NLOCAL,BTF_SG10,SD_BTF_SG10,i,j,k,CSOIL,SD_CSOIL,&
	  CGRAIN,SD_CGRAIN,EF,SD_EF,ED,SD_ED,AT,SD_AT,BW,SD_BW,KEYCONC_G10,R_IRgrain,SD_IRgrain,FI,SD_FI,AD,SD_AD)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	  LOGICAL KEYCONC_G10
	  CHARACTER(LEN=50)  :: CHEMICAL(500)
      DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),SD_CSOIL(NCHEM,NTIME,NLOCAL),CGRAIN(NCHEM,NTIME,NLOCAL),SD_CGRAIN(NCHEM,NTIME,NLOCAL)
!
!     PROTEÇÕES
!
!
	  IF(CSOIL(i,j,k).EQ.0.0)THEN
	  SOIL=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  SOIL=CSOIL(i,j,k)
	  ENDIF
	  IF(BTF_SG10.EQ.0.0)THEN
	  BTF_SG10_SUB=1.0
	  IF((KSTOPBTF_SG10.EQ.0).AND.(KEYCONC_G10.EQV..FALSE.))THEN
	  KSTOPBTF_SG10=1
	  WRITE(99,*)
	  WRITE(99,'(" WAY 10!!!! BTF (Soil-Grain) VALUE USED IS 0.0!!!! WARNING!!! CHEMICAL=",A30)')CHEMICAL(i)
	  ENDIF
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  BTF_SG10_SUB=BTF_SG10
	  ENDIF
!
	  IF(CGRAIN(i,j,k).EQ.0.0)THEN
	  GRAIN=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  GRAIN=CGRAIN(i,j,k)
	  ENDIF
!
!     FIM PROTEÇÃO
!************************************************************************************************************************************************************
      IF(KEYCONC_G10.EQV..FALSE.) THEN
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
	  CGRAIN(i,j,k)=CSOIL(i,j,k)*BTF_SG10*(1.00-0.85)
      IF((SD_CSOIL(i,j,k).NE.0.0).OR.(SD_BTF_SG10.NE.0.0))THEN
	  SD_CGRAIN(i,j,k)=CGRAIN(i,j,k)*SQRT((SD_CSOIL(i,j,k)/SOIL)**2+(SD_BTF_SG10/BTF_SG10_SUB)**2)
	  ELSE
	  SD_CGRAIN(i,j,k)=0.0
	  ENDIF
!
	  IF(CGRAIN(i,j,k).EQ.0.0)THEN
	  GRAIN=1.0
	  ELSE				    ! PROTEÇÃO CONTRA DIVISÃO DE UM NUMERO POR 0.0
	  GRAIN=CGRAIN(i,j,k)
	  ENDIF
!
	  ENDIF
!************************************************************************************************************************************************************
!
      AD=(CGRAIN(i,j,k)*R_IRgrain*FI*EF*ED)/(BW*AT)
!
      AD2=AD
!
      SD_AD=AD2*SQRT((SD_CGRAIN(i,j,k)/GRAIN)**2+(SD_IRgrain/R_IRgrain)**2+(SD_FI/FI)**2+(SD_EF/EF)**2+(SD_ED/ED)**2+(SD_BW/BW)**2+(SD_AT/AT)**2)
!
	  RETURN
	  END
!
!
!*********************************************************************************
!*********************************************************************************
!
!                 SUBROTINA PARA AVALIAÇÃO DE RISCO
!
!      sub-rotina que calcula os valores de HIag (HI agregado)
!   
!
       SUBROUTINE RISCOag(NIDADE,NVIAS,CHEMICAL,POLLUTANT,TYPE_POLLUTANT,TYPE_CHEMICAL,NPOL,NCHEM,NTIME,NLOCAL,HQ,SD_HQ,CR,SD_CR,HIag,SD_HIag,CRag,SD_CRag,HIag_ac,SD_HIag_ac,&
	   CRag_ac,SD_CRag_ac,HIag_tot,SD_HIag_tot,CRag_tot,SD_CRag_tot,HQ_tot,CR_tot,ED,SCENAR,INICIO)
!
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
       INTEGER SCENAR
!
	   CHARACTER(LEN=50)  :: CHEMICAL(500),POLLUTANT(500)
	   CHARACTER(LEN=14)  :: TYPE_CHEMICAL(2)
	   CHARACTER(LEN=20)  :: TYPE_POLLUTANT(500)
!
       DIMENSION ED(3,NIDADE)
!
       DIMENSION HIag(NIDADE,NCHEM,0:NTIME,NLOCAL),CRag(NIDADE,NCHEM,NTIME,NLOCAL)
	   DIMENSION SD_HIag(NIDADE,NCHEM,NTIME,NLOCAL),SD_HIag_TEMPORANEO(NIDADE,NCHEM,NTIME,NLOCAL),SD_CRag(NIDADE,NCHEM,NTIME,NLOCAL)
	   DIMENSION SD_CRag_TEMPORANEO(NIDADE,NCHEM,NTIME,NLOCAL)
	   DIMENSION HQ(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL),CR(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL)
	   DIMENSION SD_HQ(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL),SD_CR(NIDADE,NVIAS,NCHEM,NTIME,NLOCAL)
	   DIMENSION HIag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL),CRag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL)
	   DIMENSION SD_HIag_ac_TEMPORARIO(NIDADE,NCHEM,0:NTIME,NLOCAL),SD_HIag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL)
	   DIMENSION SD_CRag_ac_TEMPORARIO(NIDADE,NCHEM,0:NTIME,NLOCAL),SD_CRag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL)
	   DIMENSION HIag_tot(NIDADE,NCHEM,NLOCAL),CRag_tot(NIDADE,NCHEM,NLOCAL),SD_HIag_tot(NIDADE,NCHEM,NLOCAL),SD_CRag_tot(NIDADE,NCHEM,NLOCAL)
!
       DIMENSION HQ_ac(NIDADE,NVIAS,NCHEM,0:NTIME,NLOCAL),CRvias_ac(NIDADE,NVIAS,NCHEM,0:NTIME,NLOCAL)
	   DIMENSION HQ_cum_temp(NIDADE,NVIAS,0:NTIME,NLOCAL),CRvias_cum(NIDADE,NVIAS,0:NTIME,NLOCAL)
	   DIMENSION HQ_tot(NIDADE,NVIAS,NTIME,NLOCAL),CR_tot(NIDADE,NVIAS,NTIME,NLOCAL)
!
!
!
!	  INICIALIZATION BLOCK
!
      DO l=1,NIDADE
      DO i=1,NCHEM
      DO k=1,NLOCAL 
	  CRag_tot(l,i,k)=0.0
	  SD_CRag_tot(l,i,k)=0.0
	  HIag_tot(l,i,k)=0.0
	  SD_HIag_tot(l,i,k)=0.0
	  HIag(l,i,0,k)=0.0
      DO j=1,NTIME
      HIag(l,i,j,k)=0.0
	  SD_HIag_TEMPORANEO(l,i,j,k)=0.0
	  SD_HIag(l,i,j,k)=0.0
	  HIag_ac(l,i,j,k)=0.0
	  SD_HIag_ac(l,i,j,k)=0.0
      CRag(l,i,j,k)=0.0
	  SD_CRag_TEMPORANEO(l,i,j,k)=0.0
	  SD_CRag(l,i,j,k)=0.0
	  CRag_ac(l,i,j,k)=0.0
	  SD_CRag_ac(l,i,j,k)=0.0
!
      DO n=1,NVIAS
	  HQ_ac(l,n,i,j,k)=0.0
	  CRvias_ac(l,n,i,j,k)=0.0
!
      HQ_cum_temp(l,n,j,k)=0.0
	  CRvias_cum(l,n,j,k)=0.0
	  HQ_tot(l,n,j,k)=0.0
	  CR_tot(l,n,j,k)=0.0
!
      ENDDO	 ! ENDDO "n"
!
	  ENDDO	 ! ENDDO "j"
!
	  ENDDO	 ! ENDDO "k"
	  ENDDO	 ! ENDDO "i"
	  ENDDO	 ! ENDDO "l"
!
!
!    CICLO PARA METAIS NÃO CANCERIGENOS
!
!     DO N=1,19
!	 DO I=1,NCHEM
!	 DO J=1,NTIME
!    WRITE(*,*) (SD_HQ(n,i,j,k), K=1,NLOCAL)
!	 ENDDO
!	 ENDDO
!	 WRITE(*,'(" FIM N ",I3)') N
!	 ENDDO
!
!
		 LLCHEM=NCHEM
		 NEW_POL=NPOL
!
      do l=INICIO,NIDADE
!
!------------------
	  IF(SCENAR.EQ.1)THEN
	  NTIMEexp=ED(1,l)
	  ELSEIF(SCENAR.EQ.2)THEN
	  NTIMEexp=ED(2,l)
	  ELSEIF(SCENAR.EQ.3)THEN
	  NTIMEexp=ED(3,l)
	  ENDIF
!------------------

!
	  DO i=1,LLCHEM    
!
!
      DO j=1,NTIMEexp
      DO k=1,NLOCAL 
!
      DO n=1,NVIAS
!
!
      IF(HQ(l,n,i,j,k).gt.0.0) then
      HIag(l,i,j,k)= HIag(l,i,j,k) + HQ(l,n,i,j,k)
!	  WRITE(*,*)  SD_HQ(l,n,i,j,k)
	  SD_HIag_TEMPORANEO(l,i,j,k)=SD_HIag_TEMPORANEO(l,i,j,k)+SD_HQ(l,n,i,j,k)**2
	  SD_HIag(l,i,j,k)=SQRT(SD_HIag_TEMPORANEO(l,i,j,k))
      ELSE
	  ENDIF
!
!	  WRITE(*,*) SD_HIag(l,i,j,k)
!
	  IF(CR(l,n,i,j,k).gt.0.0) then
      CRag(l,i,j,k)= CRag(l,i,j,k) + CR(l,n,i,j,k)
	  SD_CRag_TEMPORANEO(l,i,j,k)=SD_CRag_TEMPORANEO(l,i,j,k)+SD_CR(l,n,i,j,k)**2
	  SD_CRag(l,i,j,k)=SQRT(SD_CRag_TEMPORANEO(l,i,j,k))
	  ELSE
	  ENDIF
!
!
      ENDDO	  !FIM DO CICLO "n"
!
      ENDDO	  !FIM DO CICLO "k"
      ENDDO	  !FIM DO CICLO "j"
      ENDDO	  !FIM DO CICLO "i"
!
!     
!							CICLO DE SOMA ACUMULADO TEMPORAL
!
!
 	  DO I=1,LLCHEM
!
      DO ii=1,NEW_POL
!								                  
	  IF (CHEMICAL(i).EQ.POLLUTANT(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
!
!      WRITE(*,'(" TYPE_POLLUTANT= ", A14)') TYPE_POLLUTANT(iii)
!	  WRITE(*,'(" TYPE_CHEMICAL2= ", A14)') TYPE_CHEMICAL(2)
!	  WRITE(*,'(" TYPE_CHEMICAL1= ", A14)') TYPE_CHEMICAL(1)
!
	  DO K=1,NLOCAL
      HIag_ac(l,i,0,k)=0.0
	  SD_HIag_ac_TEMPORARIO(l,i,0,k)=0.0
	  CRag_ac(l,i,0,k)=0.0
	  SD_CRag_ac_TEMPORARIO(l,i,0,k)=0.0
!
      DO n=1,NVIAS
      HQ_ac(l,n,i,0,k)=0.0
	  CRvias_ac(l,n,i,0,k)=0.0
	  ENDDO
!
	  DO J=1,NTIMEexp
!
	  IF(TYPE_POLLUTANT(iii).EQ.TYPE_CHEMICAL(2))THEN
	  HIag_ac(l,i,j,k)=HIag_ac(l,i,j-1,k)+HIag(l,i,j,k)
	  SD_HIag_ac_TEMPORARIO(l,i,j,k)=SD_HIag_ac_TEMPORARIO(l,i,j-1,k)+SD_HIag(l,i,j,k)**2
      SD_HIag_ac(l,i,j,k)=SQRT(SD_HIag_ac_TEMPORARIO(l,i,j,k))
!
!
	  WRITE(*,*)
	  WRITE(*,'(''*******************************************************************************************************************'')')
	  WRITE(*,'(''    WARNING!!!!!!! WARNING!!!!!!! WARNING!!!!!!! WARNING!!!!!!! WARNING!!!!!!! WARNING!!!!!!! WARNING!!!!!!!  '')')
	  WRITE(*,'(''*******************************************************************************************************************'')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' The inclusion of carcinogenic and non-carcinogenic risk assessment for chemical species that are cumulative in '')')
	  WRITE(*,'('' the human body is not available yet! '')')
	  WRITE(*,*)
 	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' CODE WILL STOP!!!! '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
!
      STOP
!
!     NOVO!
      DO n=1,NVIAS
      HQ_ac(l,n,i,j,k)=HQ_ac(l,n,i,j-1,k)+HQ(l,n,i,j,k)
	  ENDDO
!
       ELSEIF(TYPE_POLLUTANT(iii).EQ.TYPE_CHEMICAL(1))THEN
!
      HIag_ac(l,i,j,k)=HIag(l,i,j,k)
	  SD_HIag_ac(l,i,j,k)=SD_HIag(l,i,j,k)
!
!     NOVO
      DO n=1,NVIAS
      HQ_ac(l,n,i,j,k)=HQ(l,n,i,j,k)
	  ENDDO
!
	  ENDIF
!
	  CRag_ac(l,i,j,k)=CRag_ac(l,i,j-1,k)+CRag(l,i,j,k)
	  SD_CRag_ac_TEMPORARIO(l,i,j,k)=SD_CRag_ac_TEMPORARIO(l,i,j-1,k)+SD_CRag(l,i,j,k)**2
	  SD_CRag_ac(l,i,j,k)=SQRT(SD_CRag_ac_TEMPORARIO(l,i,j,k))
!
!     NOVO
      DO n=1,NVIAS
      CRvias_ac(l,n,i,j,k)=CRvias_ac(l,n,i,j-1,k)+CR(l,n,i,j,k)
	  ENDDO
!
!
	  ENDDO	  !FIM DO CICLO "J"
!
!
 	  HIag_tot(l,i,k)=HIag_ac(l,i,NTIMEexp,k)
	  SD_HIag_tot(l,i,k)=SD_HIag_ac(l,i,NTIMEexp,k)
!	  WRITE(*,*)HIag_tot(i,k)
!	  WRITE(*,*)SD_HIag_tot(i,k)
!
!
	  CRag_tot(l,i,k)=CRag_ac(l,i,NTIMEexp,k)
	  SD_CRag_tot(l,i,k)=SD_CRag_ac(l,i,NTIMEexp,k)
!
	  ENDDO	  !FIM DO CICLO "K"
	  ENDDO   !FIM DO CICLO "I"
!  
!	 DO I=1,NCHEM
!	 DO J=1,NTIME
!	  WRITE(*,*) (HIag_tot(i,k), K=1,NLOCAL)
!     WRITE(*,*) (HIag_ac(i,j,k), K=1,NLOCAL)
!	 ENDDO
!	 ENDDO 
!
      DO n=1,NVIAS
      DO k=1,NLOCAL
	  DO j=1,NTIMEexp
!
		 KKCHEM=NCHEM   
!
	  DO i=1,KKCHEM
!
      HQ_cum_temp(l,n,j,k)=HQ_cum_temp(l,n,j,k)+HQ_ac(l,n,i,j,k)
!
      CRvias_cum(l,n,j,k)=CRvias_cum(l,n,j,k)+CRvias_ac(l,n,i,j,k)
!
      ENDDO ! FIM DO i
!
      HQ_tot(l,n,j,k)=HQ_cum_temp(l,n,j,k)
	  CR_tot(l,n,j,k)=CRvias_cum(l,n,j,k)
!
      ENDDO ! FIM DO j
!
!
!
      ENDDO	 ! FIM DO k
	  ENDDO	 ! FIM DO n
!
	  ENDDO	 ! FIM DO l
!
      RETURN
	  END
!
!
!*********************************************************************************
!*********************************************************************************
!
!                 SUBROTINA PARA AVALIAÇÃO DE RISCO	CUMULATIVO
!
!      sub-rotina que calcula os valores de CRcum e HIcum (CR cumulativo)

       SUBROUTINE RISCOcum(NIDADE,NCHEM,NTIME,NLOCAL,HIag_ac,SD_HIag_ac,CRag,SD_CRag,HIcum_SUM,SD_HIcum_SUM,CRcum,SD_CRcum,&
	   HIcum_ac,SD_HIcum_ac,CRcum_ac,SD_CRcum_ac,HIcum_tot,SD_HIcum_tot,CRcum_tot,SD_CRcum_tot,ED,SCENAR,INICIO)
!
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
       INTEGER SCENAR
!
       DIMENSION ED(3,NIDADE)
!
       DIMENSION HIag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL),CRag(NIDADE,NCHEM,NTIME,NLOCAL)
       DIMENSION SD_HIag_ac(NIDADE,NCHEM,0:NTIME,NLOCAL),SD_CRag(NIDADE,NCHEM,NTIME,NLOCAL)
       DIMENSION HIcum(NIDADE,0:NTIME,NLOCAL),CRcum(NIDADE,NTIME,NLOCAL),SD_HIcum(NIDADE,NTIME,NLOCAL),SD_CRcum(NIDADE,NTIME,NLOCAL)
       DIMENSION HIcum_SUM(NIDADE,0:NTIME,NLOCAL),SD_HIcum_SUM(NIDADE,NTIME,NLOCAL),SD_HIcum_TEMP(NIDADE,NTIME,NLOCAL),SD_CRcum_TEMP(NIDADE,NTIME,NLOCAL)
	   DIMENSION HIcum_ac(NIDADE,0:NTIME, NLOCAL),CRcum_ac(NIDADE,0:NTIME, NLOCAL),HIcum_tot(NIDADE,NLOCAL),CRcum_tot(NIDADE,NLOCAL) 
	   DIMENSION SD_HIcum_ac(NIDADE,0:NTIME, NLOCAL),SD_CRcum_ac(NIDADE,0:NTIME, NLOCAL),SD_HIcum_tot(NIDADE,NLOCAL),SD_CRcum_tot(NIDADE,NLOCAL)
	   DIMENSION SD_HIcum_ac_TEMP(NIDADE,0:NTIME, NLOCAL),SD_CRcum_ac_TEMP(NIDADE,0:NTIME, NLOCAL) 
!
!       INICIALIZATION
!
      DO l=1,NIDADE
      DO j=1,NTIME
      DO k=1,NLOCAL
      HIcum_SUM(l,j,k)=0.0
	  HIcum_SUM(l,0,k)=0.0
	  SD_HIcum_SUM(l,j,k)=0.0
      CRcum(l,j,k)=0.0
	  SD_CRcum(l,j,k)=0.0
	  SD_CRcum_TEMP(l,j,k)=0.0
	  CRcum_ac(l,j,k)=0.0
	  SD_CRcum_ac(l,j,k)=0.0
	  SD_CRcum_ac_TEMP(l,j,k)=0.0
	  HIcum_ac(l,j,k)=0.0
	  SD_HIcum_ac(l,j,k)=0.0
	  SD_HIcum_ac_TEMP(l,j,k)=0.0
	  HIcum(l,j,k)=0.0
	  HIcum(l,0,k)=0.0
	  SD_HIcum(l,j,k)=0.0
	  SD_HIcum_TEMP(l,j,k)=0.0
	  ENDDO						
	  ENDDO
	  ENDDO
!
      DO l=1,NIDADE
      DO k=1,NLOCAL
	  CRcum_tot(l,k)=0.0
	  SD_CRcum_tot(l,k)=0.0
  	  HIcum_tot(l,k)=0.0
	  SD_HIcum_tot(l,k)=0.0
	  ENDDO
	  ENDDO
!
!
!    CICLO PARA METAIS NÃO CANCERIGENOS
!
!
!
	  JKCHEM=NCHEM
!
      DO l=INICIO,NIDADE
!
!------------------
	  IF(SCENAR.EQ.1)THEN
	  NTIMEexp=ED(1,l)
	  ELSEIF(SCENAR.EQ.2)THEN
	  NTIMEexp=ED(2,l)
	  ELSEIF(SCENAR.EQ.3)THEN
	  NTIMEexp=ED(3,l)
	  ENDIF
!------------------
!
      DO j=1,NTIMEexp
      DO k=1,NLOCAL 
      DO i=1,JKCHEM
!
!
      HIcum(l,j,k)= HIcum(l,j,k) + HIag_ac(l,i,j,k)
	  SD_HIcum_TEMP(l,j,k)= SD_HIcum_TEMP(l,j,k) + SD_HIag_ac(l,i,j,k)**2
	  SD_HIcum(l,j,k)=SQRT(SD_HIcum_TEMP(l,j,k))
!
      CRcum(l,j,k)=CRcum(l,j,k)+CRag(l,i,j,k)
	  SD_CRcum_TEMP(l,j,k)=SD_CRcum_TEMP(l,j,k)+SD_CRag(l,i,j,k)**2
	  SD_CRcum(l,j,k)=SQRT(SD_CRcum_TEMP(l,j,k))
!
      ENDDO
      ENDDO
      ENDDO
!
!      do j=1,ntime
!      write(*,*) (SD_CRcum(j,k), k=1,nlocal)
!	  enddo
!
      DO j=1,NTIMEexp
      DO k=1,NLOCAL 
      HIcum_SUM(l,j,k)= HIcum(l,j,k)
	  SD_HIcum_SUM(l,j,k)= SD_HIcum(l,j,k)
	  ENDDO
	  ENDDO
!
!
!
	  DO K=1,NLOCAL
	  HIcum_ac(l,0,k)=0.0
	  SD_HIcum_ac_TEMP(l,0,k)=0.0
	  CRcum_ac(l,0,k)=0.0
	  SD_CRcum_ac_TEMP(l,0,k)=0.0
!
!
	  DO J=1,NTIMEexp
!	  
      HIcum_ac(l,j,k)=HIcum_SUM(l,j,k)
	  SD_HIcum_ac(l,j,k)=SD_HIcum_SUM(l,j,k)
!
!	  HIcum_ac(j,k)= HIcum_ac(j-1,k)+HIcum_SUM(j,k)
!	  SD_HIcum_ac_TEMP(j,k)= SD_HIcum_ac_TEMP(j-1,k)+SD_HIcum_SUM(j,k)**2
!	  SD_HIcum_ac(j,k)=SQRT(SD_HIcum_ac_TEMP(j,k))
!
!
	  CRcum_ac(l,j,k)=CRcum_ac(l,j-1,k)+CRcum(l,j,k)
	  SD_CRcum_ac_TEMP(l,j,k)=SD_CRcum_ac_TEMP(l,j-1,k)+SD_CRcum(l,j,k)**2
	  SD_CRcum_ac(l,j,k)=SQRT(SD_CRcum_ac_TEMP(l,j,k))
!
	  ENDDO
  	  HIcum_tot(l,k)=HIcum_ac(l,NTIMEexp,k)
	  SD_HIcum_tot(l,k)=SD_HIcum_ac(l,NTIMEexp,k)
	  CRcum_tot(l,k)=CRcum_ac(l,NTIMEexp,k)
	  SD_CRcum_tot(l,k)=SD_CRcum_ac(l,NTIMEexp,k)
      ENDDO
!
      ENDDO  ! fim DO l=1,NIDADE
!
!	  do k=1,nlocal
!      do j=1,ed(scenar,2)
!      WRITE(*,*) HIcum(2,j,k)
!	  enddo
!	  write(*,*)
!	  enddo
!
      RETURN
	  END
!
!
!*********************************************************************************
!*********************************************************************************
!
      SUBROUTINE REDISTRI(iii,l,j,NVIAS,NIDADE,EF_INI,SD_EF_INI,BW_INI,SD_BW_INI,IR_INI,SD_IR_INI,ET_INI,SD_ET_INI,SA_INI,SD_SA_INI,&
	  EV_INI,SD_EV_INI,AF_INI,SD_AF_INI,EF,SD_EF,BW,SD_BW,IR,SD_IR,ET,SD_ET,SA,SD_SA,EV,SD_EV,AF,SD_AF,MUTAGENIC,ADAF,&
	  AT_INI,SD_AT_INI,AT,SD_AT) 
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      LOGICAL MUTAGENIC
!
      REAL*8 IR,IR_INI
!
	  DIMENSION IR_INI(NIDADE,20),SD_IR_INI(NIDADE,20),ET_INI(NIDADE,3),SD_ET_INI(NIDADE,3),SA_INI(NIDADE,2),SD_SA_INI(NIDADE,2)
	  DIMENSION EV_INI(NIDADE,2),SD_EV_INI(NIDADE,2),AF_INI(NIDADE),SD_AF_INI(NIDADE)
!!
	  DIMENSION EF_INI(NVIAS,NIDADE), BW_INI(NIDADE),AT_INI(2,NVIAS,NIDADE)
	  DIMENSION SD_EF_INI(NVIAS,NIDADE),SD_BW_INI(NIDADE),SD_AT_INI(2,NVIAS,NIDADE)
!
	  DIMENSION IR(20),SD_IR(20),ET(3),SD_ET(3),SA(2),SD_SA(2),EV(2),SD_EV(2),EF(NVIAS),SD_EF(NVIAS),MUTAGENIC(500,3),ADAF(3)
	  DIMENSION AT(2,NVIAS),SD_AT(2,NVIAS)
!
!
!
      DO IO=1,3
	  ADAF(IO)=1.0
	  ENDDO
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 1
!
      IF((l.EQ.1).AND.(j.EQ.1))THEN
	  DO I=1,20
	  IR(I)=IR_INI(1,I)
	  SD_IR(I)=SD_IR_INI(1,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(1,K)
	  SD_ET(K)=SD_ET_INI(1,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,1)
	  SD_EF(K)=SD_EF_INI(K,1)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,1)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,1)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(1,II)
	  SD_SA(II)=SD_SA_INI(1,II)
	  EV(II)=EV_INI(1,II)
	  SD_EV(II)=SD_EV_INI(1,II)
	  ENDDO 
	  AF=AF_INI(1)
	  SD_AF=SD_AF_INI(1)
	  BW=BW_INI(1)
	  SD_BW=SD_BW_INI(1)
!
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=10.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.1).AND.(j.EQ.2))THEN
	  DO I=1,20
	  IR(I)=IR_INI(2,I)
	  SD_IR(I)=SD_IR_INI(2,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(2,K)
	  SD_ET(K)=SD_ET_INI(2,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,2)
	  SD_EF(K)=SD_EF_INI(K,2)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,2)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,2)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(2,II)
	  SD_SA(II)=SD_SA_INI(2,II)
	  EV(II)=EV_INI(2,II)
	  SD_EV(II)=SD_EV_INI(2,II)
	  ENDDO 
	  AF=AF_INI(2)
	  SD_AF=SD_AF_INI(2)
	  BW=BW_INI(2)
	  SD_BW=SD_BW_INI(2)
 !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
!
	  ELSEIF((l.EQ.1).AND.(j.GE.3).AND.(j.LT.6))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(3,I)
	  SD_IR(I)=SD_IR_INI(3,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(3,K)
	  SD_ET(K)=SD_ET_INI(3,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,3)
	  SD_EF(K)=SD_EF_INI(K,3)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,3)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,3)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(3,II)
	  SD_SA(II)=SD_SA_INI(3,II)
	  EV(II)=EV_INI(3,II)
	  SD_EV(II)=SD_EV_INI(3,II)
	  ENDDO 
	  AF=AF_INI(3)
	  SD_AF=SD_AF_INI(3)
	  BW=BW_INI(3)
	  SD_BW=SD_BW_INI(3)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
!
	  ELSEIF((l.EQ.1).AND.(j.GE.6).AND.(j.LT.11))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(4,I)
	  SD_IR(I)=SD_IR_INI(4,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(4,K)
	  SD_ET(K)=SD_ET_INI(4,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,4)
	  SD_EF(K)=SD_EF_INI(K,4)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,4)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,4)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(4,II)
	  SD_SA(II)=SD_SA_INI(4,II)
	  EV(II)=EV_INI(4,II)
	  SD_EV(II)=SD_EV_INI(4,II)
	  ENDDO 
	  AF=AF_INI(4)
	  SD_AF=SD_AF_INI(4)
	  BW=BW_INI(4)
	  SD_BW=SD_BW_INI(4)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
!
	  ELSEIF((l.EQ.1).AND.(j.GE.11).AND.(j.LT.16))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(5,I)
	  SD_IR(I)=SD_IR_INI(5,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(5,K)
	  SD_ET(K)=SD_ET_INI(5,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,5)
	  SD_EF(K)=SD_EF_INI(K,5)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,5)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,5)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(5,II)
	  SD_SA(II)=SD_SA_INI(5,II)
	  EV(II)=EV_INI(5,II)
	  SD_EV(II)=SD_EV_INI(5,II)
	  ENDDO 
	  AF=AF_INI(5)
	  SD_AF=SD_AF_INI(5)
	  BW=BW_INI(5)
	  SD_BW=SD_BW_INI(5)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
!
	  ELSEIF((l.EQ.1).AND.(j.GE.16).AND.(j.LT.18))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(6,I)
	  SD_IR(I)=SD_IR_INI(6,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(6,K)
	  SD_ET(K)=SD_ET_INI(6,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,6)
	  SD_EF(K)=SD_EF_INI(K,6)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,6)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,6)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(6,II)
	  SD_SA(II)=SD_SA_INI(6,II)
	  EV(II)=EV_INI(6,II)
	  SD_EV(II)=SD_EV_INI(6,II)
	  ENDDO 
	  AF=AF_INI(6)
	  SD_AF=SD_AF_INI(6)
	  BW=BW_INI(6)
	  SD_BW=SD_BW_INI(6)
    !
      DO IO=1,3
      IF((MUTAGENIC(iii,IO).EQV..TRUE.).AND.(j.EQ.16))THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
!
	  ELSEIF((l.EQ.1).AND.(j.GE.18).AND.(j.LT.21))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
!
!
	  ELSEIF((l.EQ.1).AND.(j.GE.21).AND.(j.LT.65))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.1).AND.(j.GE.65))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
!
      ENDIF
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 2
!
      IF((l.EQ.2).AND.(j.EQ.1))THEN
	  DO I=1,20
	  IR(I)=IR_INI(2,I)
	  SD_IR(I)=SD_IR_INI(2,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(2,K)
	  SD_ET(K)=SD_ET_INI(2,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,2)
	  SD_EF(K)=SD_EF_INI(K,2)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,2)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,2)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(2,II)
	  SD_SA(II)=SD_SA_INI(2,II)
	  EV(II)=EV_INI(2,II)
	  SD_EV(II)=SD_EV_INI(2,II)
	  ENDDO 
	  AF=AF_INI(2)
	  SD_AF=SD_AF_INI(2)
	  BW=BW_INI(2)
	  SD_BW=SD_BW_INI(2)
 !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.2).AND.(j.GE.2).AND.(j.LT.5))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(3,I)
	  SD_IR(I)=SD_IR_INI(3,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(3,K)
	  SD_ET(K)=SD_ET_INI(3,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,3)
	  SD_EF(K)=SD_EF_INI(K,3)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,3)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,3)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(3,II)
	  SD_SA(II)=SD_SA_INI(3,II)
	  EV(II)=EV_INI(3,II)
	  SD_EV(II)=SD_EV_INI(3,II)
	  ENDDO 
	  AF=AF_INI(3)
	  SD_AF=SD_AF_INI(3)
	  BW=BW_INI(3)
	  SD_BW=SD_BW_INI(3)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.2).AND.(j.GE.5).AND.(j.LT.10))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(4,I)
	  SD_IR(I)=SD_IR_INI(4,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(4,K)
	  SD_ET(K)=SD_ET_INI(4,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,4)
	  SD_EF(K)=SD_EF_INI(K,4)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,4)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,4)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(4,II)
	  SD_SA(II)=SD_SA_INI(4,II)
	  EV(II)=EV_INI(4,II)
	  SD_EV(II)=SD_EV_INI(4,II)
	  ENDDO 
	  AF=AF_INI(4)
	  SD_AF=SD_AF_INI(4)
	  BW=BW_INI(4)
	  SD_BW=SD_BW_INI(4)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.2).AND.(j.GE.10).AND.(j.LT.15))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(5,I)
	  SD_IR(I)=SD_IR_INI(5,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(5,K)
	  SD_ET(K)=SD_ET_INI(5,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,5)
	  SD_EF(K)=SD_EF_INI(K,5)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,5)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,5)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(5,II)
	  SD_SA(II)=SD_SA_INI(5,II)
	  EV(II)=EV_INI(5,II)
	  SD_EV(II)=SD_EV_INI(5,II)
	  ENDDO 
	  AF=AF_INI(5)
	  SD_AF=SD_AF_INI(5)
	  BW=BW_INI(5)
	  SD_BW=SD_BW_INI(5)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.2).AND.(j.GE.15).AND.(j.LT.17))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(6,I)
	  SD_IR(I)=SD_IR_INI(6,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(6,K)
	  SD_ET(K)=SD_ET_INI(6,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,6)
	  SD_EF(K)=SD_EF_INI(K,6)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,6)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,6)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(6,II)
	  SD_SA(II)=SD_SA_INI(6,II)
	  EV(II)=EV_INI(6,II)
	  SD_EV(II)=SD_EV_INI(6,II)
	  ENDDO 
	  AF=AF_INI(6)
	  SD_AF=SD_AF_INI(6)
	  BW=BW_INI(6)
	  SD_BW=SD_BW_INI(6)
    !
      DO IO=1,3
      IF((MUTAGENIC(iii,IO).EQV..TRUE.).AND.(j.EQ.15))THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.2).AND.(j.GE.17).AND.(j.LT.20))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
      !
!
	  ELSEIF((l.EQ.2).AND.(j.GE.20).AND.(j.LT.64))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.2).AND.(j.GE.64))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 3
 !
!
	  IF((l.EQ.3).AND.(j.GE.1).AND.(j.LT.4))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(3,I)
	  SD_IR(I)=SD_IR_INI(3,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(3,K)
	  SD_ET(K)=SD_ET_INI(3,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,3)
	  SD_EF(K)=SD_EF_INI(K,3)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,3)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,3)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(3,II)
	  SD_SA(II)=SD_SA_INI(3,II)
	  EV(II)=EV_INI(3,II)
	  SD_EV(II)=SD_EV_INI(3,II)
	  ENDDO 
	  AF=AF_INI(3)
	  SD_AF=SD_AF_INI(3)
	  BW=BW_INI(3)
	  SD_BW=SD_BW_INI(3)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.3).AND.(j.GE.4).AND.(j.LT.9))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(4,I)
	  SD_IR(I)=SD_IR_INI(4,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(4,K)
	  SD_ET(K)=SD_ET_INI(4,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,4)
	  SD_EF(K)=SD_EF_INI(K,4)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,4)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,4)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(4,II)
	  SD_SA(II)=SD_SA_INI(4,II)
	  EV(II)=EV_INI(4,II)
	  SD_EV(II)=SD_EV_INI(4,II)
	  ENDDO 
	  AF=AF_INI(4)
	  SD_AF=SD_AF_INI(4)
	  BW=BW_INI(4)
	  SD_BW=SD_BW_INI(4)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.3).AND.(j.GE.9).AND.(j.LT.14))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(5,I)
	  SD_IR(I)=SD_IR_INI(5,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(5,K)
	  SD_ET(K)=SD_ET_INI(5,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,5)
	  SD_EF(K)=SD_EF_INI(K,5)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,5)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,5)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(5,II)
	  SD_SA(II)=SD_SA_INI(5,II)
	  EV(II)=EV_INI(5,II)
	  SD_EV(II)=SD_EV_INI(5,II)
	  ENDDO 
	  AF=AF_INI(5)
	  SD_AF=SD_AF_INI(5)
	  BW=BW_INI(5)
	  SD_BW=SD_BW_INI(5)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.3).AND.(j.GE.14).AND.(j.LT.16))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(6,I)
	  SD_IR(I)=SD_IR_INI(6,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(6,K)
	  SD_ET(K)=SD_ET_INI(6,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,6)
	  SD_EF(K)=SD_EF_INI(K,6)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,6)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,6)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(6,II)
	  SD_SA(II)=SD_SA_INI(6,II)
	  EV(II)=EV_INI(6,II)
	  SD_EV(II)=SD_EV_INI(6,II)
	  ENDDO 
	  AF=AF_INI(6)
	  SD_AF=SD_AF_INI(6)
	  BW=BW_INI(6)
	  SD_BW=SD_BW_INI(6)
    !
      DO IO=1,3
      IF((MUTAGENIC(iii,IO).EQV..TRUE.).AND. (j.EQ.14))THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.3).AND.(j.GE.16).AND.(j.LT.19))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
      !
!
	  ELSEIF((l.EQ.3).AND.(j.GE.19).AND.(j.LT.63))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.3).AND.(j.GE.63))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 4
  !
!
	  IF((l.EQ.4).AND.(j.GE.1).AND.(j.LT.6))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(4,I)
	  SD_IR(I)=SD_IR_INI(4,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(4,K)
	  SD_ET(K)=SD_ET_INI(4,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,4)
	  SD_EF(K)=SD_EF_INI(K,4)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,4)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,4)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(4,II)
	  SD_SA(II)=SD_SA_INI(4,II)
	  EV(II)=EV_INI(4,II)
	  SD_EV(II)=SD_EV_INI(4,II)
	  ENDDO 
	  AF=AF_INI(4)
	  SD_AF=SD_AF_INI(4)
	  BW=BW_INI(4)
	  SD_BW=SD_BW_INI(4)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.4).AND.(j.GE.6).AND.(j.LT.11))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(5,I)
	  SD_IR(I)=SD_IR_INI(5,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(5,K)
	  SD_ET(K)=SD_ET_INI(5,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,5)
	  SD_EF(K)=SD_EF_INI(K,5)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,5)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,5)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(5,II)
	  SD_SA(II)=SD_SA_INI(5,II)
	  EV(II)=EV_INI(5,II)
	  SD_EV(II)=SD_EV_INI(5,II)
	  ENDDO 
	  AF=AF_INI(5)
	  SD_AF=SD_AF_INI(5)
	  BW=BW_INI(5)
	  SD_BW=SD_BW_INI(5)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.4).AND.(j.GE.11).AND.(j.LT.13))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(6,I)
	  SD_IR(I)=SD_IR_INI(6,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(6,K)
	  SD_ET(K)=SD_ET_INI(6,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,6)
	  SD_EF(K)=SD_EF_INI(K,6)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,6)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,6)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(6,II)
	  SD_SA(II)=SD_SA_INI(6,II)
	  EV(II)=EV_INI(6,II)
	  SD_EV(II)=SD_EV_INI(6,II)
	  ENDDO 
	  AF=AF_INI(6)
	  SD_AF=SD_AF_INI(6)
	  BW=BW_INI(6)
	  SD_BW=SD_BW_INI(6)
    !
      DO IO=1,3
      IF((MUTAGENIC(iii,IO).EQV..TRUE.).AND.(j.EQ.11))THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.4).AND.(j.GE.13).AND.(j.LT.16))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
      !
!
	  ELSEIF((l.EQ.4).AND.(j.GE.16).AND.(j.LT.60))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.4).AND.(j.GE.60))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 5
  !
!
	  IF((l.EQ.5).AND.(j.GE.1).AND.(j.LT.6))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(5,I)
	  SD_IR(I)=SD_IR_INI(5,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(5,K)
	  SD_ET(K)=SD_ET_INI(5,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,5)
	  SD_EF(K)=SD_EF_INI(K,5)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,5)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,5)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(5,II)
	  SD_SA(II)=SD_SA_INI(5,II)
	  EV(II)=EV_INI(5,II)
	  SD_EV(II)=SD_EV_INI(5,II)
	  ENDDO 
	  AF=AF_INI(5)
	  SD_AF=SD_AF_INI(5)
	  BW=BW_INI(5)
	  SD_BW=SD_BW_INI(5)
  !
      DO IO=1,3
      IF(MUTAGENIC(iii,IO).EQV..TRUE.)THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.5).AND.(j.GE.6).AND.(j.LT.8))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(6,I)
	  SD_IR(I)=SD_IR_INI(6,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(6,K)
	  SD_ET(K)=SD_ET_INI(6,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,6)
	  SD_EF(K)=SD_EF_INI(K,6)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,6)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,6)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(6,II)
	  SD_SA(II)=SD_SA_INI(6,II)
	  EV(II)=EV_INI(6,II)
	  SD_EV(II)=SD_EV_INI(6,II)
	  ENDDO 
	  AF=AF_INI(6)
	  SD_AF=SD_AF_INI(6)
	  BW=BW_INI(6)
	  SD_BW=SD_BW_INI(6)
    !
      DO IO=1,3
      IF((MUTAGENIC(iii,IO).EQV..TRUE.).AND.(j.EQ.6))THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.5).AND.(j.GE.8).AND.(j.LT.11))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
      !
!
	  ELSEIF((l.EQ.5).AND.(j.GE.11).AND.(j.LT.55))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.5).AND.(j.GE.55))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 6

  !
!
	  IF((l.EQ.6).AND.(j.GE.1).AND.(j.LT.3))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(6,I)
	  SD_IR(I)=SD_IR_INI(6,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(6,K)
	  SD_ET(K)=SD_ET_INI(6,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,6)
	  SD_EF(K)=SD_EF_INI(K,6)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,6)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,6)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(6,II)
	  SD_SA(II)=SD_SA_INI(6,II)
	  EV(II)=EV_INI(6,II)
	  SD_EV(II)=SD_EV_INI(6,II)
	  ENDDO 
	  AF=AF_INI(6)
	  SD_AF=SD_AF_INI(6)
	  BW=BW_INI(6)
	  SD_BW=SD_BW_INI(6)
    !
      DO IO=1,3
      IF((MUTAGENIC(iii,IO).EQV..TRUE.).AND.(j.EQ.1))THEN
	  ADAF(IO)=3.0
	  ENDIF
	  ENDDO
!
	  ELSEIF((l.EQ.6).AND.(j.GE.3).AND.(j.LT.6))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
      !
!
	  ELSEIF((l.EQ.6).AND.(j.GE.6).AND.(j.LT.50))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.6).AND.(j.GE.50))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!
!     ----------------------------------------------------------------------
!     QUANDO l = 7
!
    !
!
	  IF((l.EQ.7).AND.(j.GE.1).AND.(j.LT.4))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(7,I)
	  SD_IR(I)=SD_IR_INI(7,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(7,K)
	  SD_ET(K)=SD_ET_INI(7,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,7)
	  SD_EF(K)=SD_EF_INI(K,7)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,7)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,7)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(7,II)
	  SD_SA(II)=SD_SA_INI(7,II)
	  EV(II)=EV_INI(7,II)
	  SD_EV(II)=SD_EV_INI(7,II)
	  ENDDO 
	  AF=AF_INI(7)
	  SD_AF=SD_AF_INI(7)
	  BW=BW_INI(7)
	  SD_BW=SD_BW_INI(7)
      !
!
	  ELSEIF((l.EQ.7).AND.(j.GE.4).AND.(j.LT.48))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO 
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.7).AND.(j.GE.48))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!     ----------------------------------------------------------------------
!     QUANDO l = 8
!
      !
!
	  IF((l.EQ.8).AND.(j.GE.1).AND.(j.LT.45))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(8,I)
	  SD_IR(I)=SD_IR_INI(8,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(8,K)
	  SD_ET(K)=SD_ET_INI(8,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,8)
	  SD_EF(K)=SD_EF_INI(K,8)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,8)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,8)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(8,II)
	  SD_SA(II)=SD_SA_INI(8,II)
	  EV(II)=EV_INI(8,II)
	  SD_EV(II)=SD_EV_INI(8,II)
	  ENDDO
	  AF=AF_INI(8)
	  SD_AF=SD_AF_INI(8)
	  BW=BW_INI(8)
	  SD_BW=SD_BW_INI(8)
      !
!
	  ELSEIF((l.EQ.8).AND.(j.GE.45))THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
	  IF(l.EQ.9)THEN		 
	  DO I=1,20
	  IR(I)=IR_INI(9,I)
	  SD_IR(I)=SD_IR_INI(9,I)
	  ENDDO
	  DO K=1,3
	  ET(K)=ET_INI(9,K)
	  SD_ET(K)=SD_ET_INI(9,K)
	  ENDDO
	  DO K=1,NVIAS
	  EF(K)=EF_INI(K,9)
	  SD_EF(K)=SD_EF_INI(K,9)
	  DO IXI=1,2 
	  AT(IXI,K)=AT_INI(IXI,K,9)
	  SD_AT(IXI,K)=SD_AT_INI(IXI,K,9)
	  ENDDO   !fim do(IXI)
	  ENDDO
	  DO II=1,2
	  SA(II)=SA_INI(9,II)
	  SD_SA(II)=SD_SA_INI(9,II)
	  EV(II)=EV_INI(9,II)
	  SD_EV(II)=SD_EV_INI(9,II)
	  ENDDO 
	  AF=AF_INI(9)
	  SD_AF=SD_AF_INI(9)
	  BW=BW_INI(9)
	  SD_BW=SD_BW_INI(9)
!
      ENDIF
!
!
!
	  RETURN
	  END
!
!*********************************************************************************************************************************************************
!*********************************************************************************************************************************************************
!
      SUBROUTINE RADIO(VARIASAO,NDURATION,NCHEM,NTIME,NLOCAL,NTYPECONC,CHEMICAL,RAD_POL,AMOLAR,VIDAMEIA,KEY_SD)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
	  LOGICAL KEYCONC,KEY_SD
	   CHARACTER(LEN=50)  :: CHEMICAL(500),RAD_POL(4)
	CHARACTER aspas*1
	DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),CWATER(NCHEM,NTIME,NLOCAL),CWATERDER(NCHEM,NTIME,NLOCAL),CWATEROTHER(NCHEM,NTIME,NLOCAL)
	DIMENSION CPAR(NCHEM,NTIME,NLOCAL),CSTEAM(NCHEM,NTIME,NLOCAL),CFRUIT(NCHEM,NTIME,NLOCAL),CLEAVES(NCHEM,NTIME,NLOCAL)
	DIMENSION CBEEF(NCHEM,NTIME,NLOCAL),CMILK(NCHEM,NTIME,NLOCAL),CAVE(NCHEM,NTIME,NLOCAL)
	DIMENSION CEGG(NCHEM,NTIME,NLOCAL),CFISH(NCHEM,NTIME,NLOCAL),CGRAIN(NCHEM,NTIME,NLOCAL),CSEDIMENT(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL),SD_CWATERDER(NCHEM,NTIME,NLOCAL),SD_CWATEROTHER(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CPAR(NCHEM,NTIME,NLOCAL),SD_CSTEAM(NCHEM,NTIME,NLOCAL),SD_CFRUIT(NCHEM,NTIME,NLOCAL),SD_CLEAVES(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CBEEF(NCHEM,NTIME,NLOCAL),SD_CMILK(NCHEM,NTIME,NLOCAL),SD_CAVE(NCHEM,NTIME,NLOCAL),SD_CEGG(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CFISH(NCHEM,NTIME,NLOCAL),SD_CGRAIN(NCHEM,NTIME,NLOCAL),SD_CSEDIMENT(NCHEM,NTIME,NLOCAL)
	DIMENSION KEYCONC(NCHEM,NTYPECONC),TIMESP(0:NTIME)  ! NTYPECONC SÃO OS 14 TIPOS DE CONCENTRAÇOES ACIMA !!!!
!
    DIMENSION AMOLAR(4),VIDAMEIA(4),CONC_RADIO(3,NTIME,NLOCAL),SD_CONC_RADIO(3,NTIME,NLOCAL)
	DIMENSION ATIVI_ESP(3,0:NTIME,NLOCAL),SD_ATIVI_ESP(3,NTIME,NLOCAL),ATIVI_EQU(0:NTIME,NLOCAL),SD_ATIVI_EQU(NTIME,NLOCAL)
	DIMENSION DOSE_ABS(NTIME,NLOCAL),SD_DOSE_ABS(NTIME,NLOCAL),DOSE_IN(NTIME,NLOCAL),SD_DOSE_IN(NTIME,NLOCAL)
	DIMENSION DOSE_OUT(NTIME,NLOCAL),SD_DOSE_OUT(NTIME,NLOCAL),CARC_IN(NTIME,NLOCAL),SD_CARC_IN(NTIME,NLOCAL)
	DIMENSION CARC_OUT(NTIME,NLOCAL),SD_CARC_OUT(NTIME,NLOCAL),CARC_IN_ac(0:NTIME,NLOCAL),CARC_OUT_ac(0:NTIME,NLOCAL)
	DIMENSION CARC_TOT_ac(NTIME,NLOCAL),SD_CARC_TOT_ac(NTIME,NLOCAL)
	DIMENSION SD_CARC_IN_ac_TEMP(0:NTIME,NLOCAL),SD_CARC_OUT_ac_TEMP(0:NTIME,NLOCAL),SD_CARC_IN_ac(NTIME,NLOCAL),SD_CARC_OUT_ac(NTIME,NLOCAL)
	DIMENSION HIND_EXT(0:NTIME,NLOCAL),HIND_INT(NTIME,NLOCAL),SD_HIND_EXT(NTIME,NLOCAL),SD_HIND_INT(NTIME,NLOCAL)
	DIMENSION KALL(NLOCAL),JALL(NLOCAL),IALL(NLOCAL)
!
!
      CALL CONCENTRATION(VARIASAO,NDURATION,NCHEM,NTYPECONC,CHEMICAL,&
	  CSOIL,CWATER,CPAR,CSTEAM,CFRUIT,CLEAVES,CBEEF,CMILK,CAVE,CEGG,CFISH,CGRAIN,CWATERDER,CWATEROTHER,CSEDIMENT,&
	  KEYCONC,NTIME,NLOCAL,TIMESP,&
      SD_CSOIL,SD_CWATER,SD_CPAR,SD_CSTEAM,SD_CFRUIT,SD_CLEAVES,SD_CBEEF,SD_CMILK,SD_CAVE,SD_CEGG,SD_CFISH,SD_CGRAIN,SD_CWATERDER,SD_CWATEROTHER,SD_CSEDIMENT)
!
!	  
	  DO k=1,NLOCAL
	  KALL(k)=0
	  JALL(k)=0
      IALL(k)=0
      DO i=1,3
	  ATIVI_ESP(i,0,k)=0.0
      ATIVI_EQU(0,k)=0.0
	  CARC_IN_ac(0,k)=0.0
	  CARC_OUT_ac(0,k)=0.0
	  HIND_EXT(0,k)=0.0
	  SD_CARC_IN_ac_TEMP(0,k)=0.0
	  SD_CARC_OUT_ac_TEMP(0,k)=0.0
	  DO j=1,NTIME
      CONC_RADIO(i,j,k)=0.0
	  SD_CONC_RADIO(i,j,k)=0.0
	  ATIVI_ESP(i,j,k)=0.0
	  SD_ATIVI_ESP(i,j,k)=0.0
!
	  ATIVI_EQU(j,k)=0.0
	  SD_ATIVI_EQU(j,k)=0.0
!
	  DOSE_ABS(j,k)=0.0
	  SD_DOSE_ABS(j,k)=0.0
!
	  DOSE_IN(j,k)=0.0
	  SD_DOSE_IN(j,k)=0.0
	  DOSE_OUT(j,k)=0.0
	  SD_DOSE_OUT(j,k)=0.0
!
	  CARC_IN(j,k)=0.0
	  SD_CARC_IN(j,k)=0.0
	  CARC_OUT(j,k)=0.0
	  SD_CARC_OUT(j,k)=0.0
!
	  CARC_IN_ac(j,k)=0.0
	  CARC_OUT_ac(j,k)=0.0
	  SD_CARC_IN_ac_TEMP(j,k)=0.0
	  SD_CARC_OUT_ac_TEMP(j,k)=0.0
	  SD_CARC_IN_ac(j,k)=0.0
	  SD_CARC_OUT_ac(j,k)=0.0
	  CARC_TOT_ac(j,k)=0.0
	  SD_CARC_TOT_ac(j,k)=0.0
!
      HIND_EXT(j,k)=0.0
	  HIND_INT(j,k)=0.0
	  SD_HIND_EXT(j,k)=0.0
	  SD_HIND_INT(j,k)=0.0
	  ENDDO
	  ENDDO
	  ENDDO
!
!
      ii=0
      iii=0
!
      DO i=1,NCHEM	  !CICLO POR ESPÉCIE QUÍMICA
!  
!
!								                  
	  IF (CHEMICAL(i).EQ.RAD_POL(3)) THEN
	  ii=1
	  ELSEIF(CHEMICAL(i).EQ.RAD_POL(4)) THEN
	  iii=2
	  ENDIF 
!
      kol=ii+iii
!
      DO k=1,NLOCAL		!CICLO POR LOCAL
!
      DO j=1,NTIME		!CICLO POR TEMPO
!
      IF(CHEMICAL(i).EQ.RAD_POL(1))THEN				 ! CONC_RADIO(1,j,k) = concentração do K-40
      CONC_RADIO(1,j,k)=CSOIL(i,j,k)
	  SD_CONC_RADIO(1,j,k)=SD_CSOIL(i,j,k)
      ELSEIF(CHEMICAL(i).EQ.RAD_POL(2))THEN			 ! CONC_RADIO(2,j,k) = concentração do Th-232
	  CONC_RADIO(2,j,k)=CSOIL(i,j,k)
	  SD_CONC_RADIO(2,j,k)=SD_CSOIL(i,j,k)
      ELSEIF(CHEMICAL(i).EQ.RAD_POL(3))THEN			 ! CONC_RADIO(3,j,k) = concentração do U-238 ou Ra-226
      CONC_RADIO(3,j,k)=CSOIL(i,j,k)
	  SD_CONC_RADIO(3,j,k)=SD_CSOIL(i,j,k)
      ELSEIF(CHEMICAL(i).EQ.RAD_POL(4))THEN
!
	  IF(kol.LT.3)THEN
      CONC_RADIO(3,j,k)=CSOIL(i,j,k)
	  SD_CONC_RADIO(3,j,k)=SD_CSOIL(i,j,k)
	  ENDIF
!
	  ENDIF
!
      ENDDO	   !FIM DO CICLO POR TEMPO
      ENDDO    !FIM DO CICLO POR LOCAL
      ENDDO   ! FIM CICLO POR ESPÉCIE QUÍMICA
!
!-----------------------------------------------------------
!	  CALCULO DA ATIVIDADE ESPECÍFICA
!-----------------------------------------------------------
!     
      DO i=1,3		! 1 representa K-40, 2 representa Th-232, 3 representa U-238 ou Ra-226
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
      IF(i.LT.3)THEN	  ! IF 1
      ATIVI_ESP(i,j,k)=ALOG(2.0)*(CONC_RADIO(i,j,k)*6.022E+23/AMOLAR(i))/VIDAMEIA(i)
!
	  IF(CONC_RADIO(i,j,k).NE.0.0)THEN	 ! IF 2
      SD_ATIVI_ESP(i,j,k)=ATIVI_ESP(i,j,k)*(SD_CONC_RADIO(i,j,k)/CONC_RADIO(i,j,k))
	  ELSE
	  SD_ATIVI_ESP(i,j,k)=0.0
	  ENDIF					             ! FIM IF 2
!
	  ELSE
!
      IF(kol.EQ.2)THEN					 ! IF 3
      ATIVI_ESP(i,j,k)=ALOG(2.0)*(CONC_RADIO(i,j,k)*6.022E+23/AMOLAR(4))/VIDAMEIA(4)
      ELSE
      ATIVI_ESP(i,j,k)=ALOG(2.0)*(CONC_RADIO(i,j,k)*6.022E+23/AMOLAR(3))/VIDAMEIA(3)
      ENDIF                              ! FIM IF 3

	  IF(CONC_RADIO(i,j,k).NE.0.0)THEN	 ! IF 4
      SD_ATIVI_ESP(i,j,k)=ATIVI_ESP(i,j,k)*(SD_CONC_RADIO(i,j,k)/CONC_RADIO(i,j,k))
	  ELSE
	  SD_ATIVI_ESP(i,j,k)=0.0
	  ENDIF								 ! FIM IF 4
!
      ENDIF	                             ! FIM IF 1
!
!
      ENDDO	   ! FIM DO CICLO POR LOCAL
	  ENDDO	   ! FIM DO CICLO POR TEMPO
      ENDDO    ! FIM DO CICLO POR RADIONUCLIDEOS
!    
!
!-----------------------------------------------------------
!	  CALCULO DA ATIVIDADE EQUIVALENTE
!-----------------------------------------------------------
!	  
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
      ATIVI_EQU(j,k)=0.077*ATIVI_ESP(1,j,k)+1.43*ATIVI_ESP(2,j,k)+ATIVI_ESP(3,j,k)
!
      IF((SD_ATIVI_ESP(1,j,k).NE.0.0).OR.(SD_ATIVI_ESP(2,j,k).NE.0.0).OR.(SD_ATIVI_ESP(3,j,k).NE.0.0))THEN
	  SD_ATIVI_EQU(j,k)=SQRT((0.077*SD_ATIVI_ESP(1,j,k))**2+(1.43*SD_ATIVI_ESP(2,j,k))**2+SD_ATIVI_ESP(3,j,k)**2)
	  ELSE
	  SD_ATIVI_EQU(j,k)=0.0
	  ENDIF
!
      ENDDO
	  ENDDO
!     
!-----------------------------------------------------------
!	  CALCULO DA DOSE ABSORVIDA
!-----------------------------------------------------------
!	  
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
      DOSE_ABS(j,k)=0.0417*ATIVI_ESP(1,j,k)+0.621*ATIVI_ESP(2,j,k)+0.462*ATIVI_ESP(3,j,k)
!
      IF((SD_ATIVI_ESP(1,j,k).NE.0.0).OR.(SD_ATIVI_ESP(2,j,k).NE.0.0).OR.(SD_ATIVI_ESP(3,j,k).NE.0.0))THEN
	  SD_DOSE_ABS(j,k)=SQRT((0.0417*SD_ATIVI_ESP(1,j,k))**2+(0.621*SD_ATIVI_ESP(2,j,k))**2+(0.462*SD_ATIVI_ESP(3,j,k))**2)
	  ELSE
	  SD_DOSE_ABS(j,k)=0.0
	  ENDIF
!
      ENDDO
	  ENDDO
!    
!-----------------------------------------------------------
!	  CALCULO DAS DOSES EFETIVAS INDOOR E OUTDOOR
!-----------------------------------------------------------
!	  
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
      DOSE_IN(j,k)=DOSE_ABS(j,k)*8760.0*0.8*7.0E-7
	  DOSE_OUT(j,k)=DOSE_ABS(j,k)*8760.0*0.2*7.0E-7
!
	  IF(DOSE_ABS(j,k).NE.0.0)THEN	 
      SD_DOSE_IN(j,k)=DOSE_IN(j,k)*(SD_DOSE_ABS(j,k)/DOSE_ABS(j,k))
      SD_DOSE_OUT(j,k)=DOSE_OUT(j,k)*(SD_DOSE_ABS(j,k)/DOSE_ABS(j,k))
	  ELSE
      SD_DOSE_IN(j,k)=0.0
      SD_DOSE_OUT(j,k)=0.0
	  ENDIF	
!
      ENDDO
	  ENDDO
!    
!-----------------------------------------------------------
!	  CALCULO DO RISCO CARCINOGENICO
!-----------------------------------------------------------
!	  
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
      CARC_IN(j,k)=DOSE_IN(j,k)/1000*78.0*0.05
	  CARC_OUT(j,k)=DOSE_OUT(j,k)/1000*78.0*0.05
!
	  IF(DOSE_IN(j,k).NE.0.0)THEN	 
      SD_CARC_IN(j,k)=CARC_IN(j,k)*(SD_DOSE_IN(j,k)/DOSE_IN(j,k))
      SD_CARC_OUT(j,k)=CARC_OUT(j,k)*(SD_DOSE_OUT(j,k)/DOSE_OUT(j,k))
	  ELSE
      SD_CARC_IN(j,k)=0.0
      SD_CARC_OUT(j,k)=0.0
	  ENDIF	
!
      ENDDO
	  ENDDO
!
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
	  CARC_IN_ac(j,k)=CARC_IN_ac(j-1,k)+CARC_IN(j,k)
!
	  SD_CARC_IN_ac_TEMP(j,k)=SD_CARC_IN_ac_TEMP(j-1,k)+SD_CARC_IN(j,k)**2
	  IF(SD_CARC_IN_ac_TEMP(j,k).NE.0.0)THEN
	  SD_CARC_IN_ac(j,k)=SQRT(SD_CARC_IN_ac_TEMP(j,k))	
	  ELSE
	  SD_CARC_IN_ac(j,k)=0.0
	  ENDIF
!
	  CARC_OUT_ac(j,k)=CARC_OUT_ac(j-1,k)+CARC_OUT(j,k)
	  SD_CARC_OUT_ac_TEMP(j,k)=SD_CARC_OUT_ac_TEMP(j-1,k)+SD_CARC_OUT(j,k)**2
	  IF(SD_CARC_OUT_ac_TEMP(j,k).NE.0.0)THEN
	  SD_CARC_OUT_ac(j,k)=SQRT(SD_CARC_OUT_ac_TEMP(j,k))	
	  ELSE
	  SD_CARC_OUT_ac(j,k)=0.0
	  ENDIF
!
      CARC_TOT_ac(j,k)=CARC_IN_ac(j,k)+CARC_OUT_ac(j,k)
	  IF((SD_CARC_OUT_ac(j,k).NE.0.0).OR.(SD_CARC_IN_ac(j,k).NE.0.0))THEN
	  SD_CARC_TOT_ac(j,k)=SQRT(SD_CARC_OUT_ac(j,k)**2+SD_CARC_IN_ac(j,k)**2)
	  ELSE
	  SD_CARC_OUT_ac(j,k)=0.0
	  ENDIF
!
      ENDDO
	  ENDDO
!     
!-----------------------------------------------------------
!	  CALCULO DO INDICE DE PERIGO EXTERNO E INTERNO
!-----------------------------------------------------------
!	  
      DO j=1,NTIME
	  DO k=1,NLOCAL
!
      HIND_EXT(j,k)=0.00020*ATIVI_ESP(1,j,k)+0.00386*ATIVI_ESP(2,j,k)+0.00270*ATIVI_ESP(3,j,k)
      HIND_INT(j,k)=0.00020*ATIVI_ESP(1,j,k)+0.00386*ATIVI_ESP(2,j,k)+0.00540*ATIVI_ESP(3,j,k)
!
      IF((SD_ATIVI_ESP(1,j,k).NE.0.0).OR.(SD_ATIVI_ESP(2,j,k).NE.0.0).OR.(SD_ATIVI_ESP(3,j,k).NE.0.0))THEN
	  SD_HIND_EXT(j,k)=SQRT((0.00020*SD_ATIVI_ESP(1,j,k))**2+(0.00386*SD_ATIVI_ESP(2,j,k))**2+(0.00270*SD_ATIVI_ESP(3,j,k))**2)
	  SD_HIND_INT(j,k)=SQRT((0.00020*SD_ATIVI_ESP(1,j,k))**2+(0.00386*SD_ATIVI_ESP(2,j,k))**2+(0.00540*SD_ATIVI_ESP(3,j,k))**2)
	  ELSE
	  SD_HIND_EXT(j,k)=0.0
	  SD_HIND_INT(j,k)=0.0
	  ENDIF
!
      ENDDO
	  ENDDO
!  
!
!-----------------------------------------------------------
!	  IMPRESSÃO DOS DADOS DE SAÍDA
!-----------------------------------------------------------
!
!
      aspas = '"'
!
	  open (UNIT=63, file='Results\Radiological Risk.json')	
!
	  WRITE(63,'("{")')
!  
	  WRITE(63,'(A1,"Specific Activities",A1,": [")')aspas,aspas
!
      IF(NTIME.EQ.1)THEN
      NCOMESO=1
	  ELSE
	  NCOMESO=2
	  ENDIF						
!
      DO k=1,NLOCAL
	  DO j=NCOMESO,NTIME
!
      IF((ATIVI_ESP(1,j,k).NE.ATIVI_ESP(1,j-1,k)).OR.(ATIVI_ESP(2,j,k).NE.ATIVI_ESP(2,j-1,k)).OR.(ATIVI_ESP(3,j,k).NE.ATIVI_ESP(3,j-1,k)))THEN
      KALL(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
!	  IF(j.EQ.1)THEN
!
      IF(KALL(k).EQ.1)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"K-40 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(1,j,k)
      write(63,'(A1,"K-40 error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ATIVI_ESP(1,j,k)
      write(63,'(A1,"Th-232 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(2,j,k)
      write(63,'(A1,"Th-232 error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ATIVI_ESP(2,j,k)
	  IF(kol.EQ.2)THEN
      write(63,'(A1,"Ra-226 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"Ra-226 error",A1,":",1x,ES12.5) )') aspas,aspas,SD_ATIVI_ESP(3,j,k)
	  ELSE
      write(63,'(A1,"U-238 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"U-238 error",A1,":",1x,ES12.5) )') aspas,aspas,SD_ATIVI_ESP(3,j,k)
      ENDIF
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
!
      ELSEIF((KALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"K-40 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(1,j,k)
      write(63,'(A1,"K-40 error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ATIVI_ESP(1,j,k)
      write(63,'(A1,"Th-232 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(2,j,k)
      write(63,'(A1,"Th-232 error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ATIVI_ESP(2,j,k)
	  IF(kol.EQ.2)THEN
      write(63,'(A1,"Ra-226 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"Ra-226 error",A1,":",1x,ES12.5) )') aspas,aspas,SD_ATIVI_ESP(3,j,k)
	  ELSE
      write(63,'(A1,"U-238 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"U-238 error",A1,":",1x,ES12.5) )') aspas,aspas,SD_ATIVI_ESP(3,j,k)
      ENDIF
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
      ENDIF
!
	  ELSE
!
!
      IF(KALL(k).EQ.1)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"K-40 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(1,j,k)
      write(63,'(A1,"K-40 error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Th-232 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(2,j,k)
      write(63,'(A1,"Th-232 error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
	  IF(kol.EQ.2)THEN
      write(63,'(A1,"Ra-226 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"Ra-226 error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
	  ELSE
      write(63,'(A1,"U-238 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"U-238 error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      ENDIF
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
      ELSEIF((KALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"K-40 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(1,j,k)
      write(63,'(A1,"K-40 error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Th-232 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(2,j,k)
      write(63,'(A1,"Th-232 error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
	  IF(kol.EQ.2)THEN
      write(63,'(A1,"Ra-226 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"Ra-226 error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
	  ELSE
      write(63,'(A1,"U-238 value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_ESP(3,j,k)
      write(63,'(A1,"U-238 error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      ENDIF
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ENDIF
!
!
	  ENDIF 
!
      ENDDO		! fim DO j=1,NTIME
!
	  ENDDO		! fim DO k=1,NLOCAL
!
!
	  WRITE(63,'("],")')
!
	  WRITE(63,'(A1,"Radium equivalent activities, Absorbed doses rates and Annual effective doses rates",A1,": [")')aspas,aspas
!
      DO k=1,NLOCAL
	  DO j=NCOMESO,NTIME
!
      IF(ATIVI_EQU(j,k).NE.ATIVI_EQU(j-1,k))THEN
      JALL(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
!
      IF(JALL(k).EQ.1)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"Ra-Eq value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_EQU(j,k)
      write(63,'(A1,"Ra-Eq error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ATIVI_EQU(j,k)
      write(63,'(A1,"AD value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_ABS(j,k)
      write(63,'(A1,"AD error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_DOSE_ABS(j,k)
      write(63,'(A1,"Heff-Out value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_OUT(j,k)
      write(63,'(A1,"Heff-Out error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_DOSE_OUT(j,k)
      write(63,'(A1,"Heff-In value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_IN(j,k)
      write(63,'(A1,"Heff-In error",A1,":",1x,ES12.5) )') aspas,aspas,SD_DOSE_IN(j,k)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ELSEIF((JALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Ra-Eq value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_EQU(j,k)
      write(63,'(A1,"Ra-Eq error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_ATIVI_EQU(j,k)
      write(63,'(A1,"AD value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_ABS(j,k)
      write(63,'(A1,"AD error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_DOSE_ABS(j,k)
      write(63,'(A1,"Heff-Out value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_OUT(j,k)
      write(63,'(A1,"Heff-Out error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_DOSE_OUT(j,k)
      write(63,'(A1,"Heff-In value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_IN(j,k)
      write(63,'(A1,"Heff-In error",A1,":",1x,ES12.5) )') aspas,aspas,SD_DOSE_IN(j,k)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
      ENDIF
!
	  ELSE
!
      IF(JALL(k).EQ.1)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"Ra-Eq value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_EQU(j,k)
      write(63,'(A1,"Ra-Eq error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"AD value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_ABS(j,k)
      write(63,'(A1,"AD error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Heff-Out value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_OUT(j,k)
      write(63,'(A1,"Heff-Out error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Heff-In value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_IN(j,k)
      write(63,'(A1,"Heff-In error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ELSEIF((JALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Ra-Eq value",A1,":",1x,ES12.5,",") )') aspas,aspas,ATIVI_EQU(j,k)
      write(63,'(A1,"Ra-Eq error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"AD value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_ABS(j,k)
      write(63,'(A1,"AD error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Heff-Out value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_OUT(j,k)
      write(63,'(A1,"Heff-Out error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Heff-In value",A1,":",1x,ES12.5,",") )') aspas,aspas,DOSE_IN(j,k)
      write(63,'(A1,"Heff-In error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ENDIF
!
	  ENDIF 
!
      ENDDO		! fim DO j=1,NTIME
!
	  ENDDO		! fim DO k=1,NLOCAL
!
!
	  WRITE(63,'("],")')
!
	  WRITE(63,'(A1,"Excess Lifetime Cancer Risk",A1,": [")')aspas,aspas
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"Indoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_IN(j,k)
      write(63,'(A1,"Indoor error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_CARC_IN(j,k)
      write(63,'(A1,"Accumulated indoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_IN_ac(j,k)
      write(63,'(A1,"Accumulated indoor error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_CARC_IN_ac(j,k)
      write(63,'(A1,"Outdoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_OUT(j,k)
      write(63,'(A1,"Outdoor error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_CARC_OUT(j,k)
      write(63,'(A1,"Accumulated outdoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_OUT_ac(j,k)
      write(63,'(A1,"Accumulated outdoor error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_CARC_OUT_ac(j,k)
      write(63,'(A1,"Accumulated total ELCR",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_TOT_ac(j,k)
      write(63,'(A1,"Accumulated total ELCR",A1,":",1x,ES12.5) )') aspas,aspas,SD_CARC_TOT_ac(j,k)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ELSE
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"Indoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_IN(j,k)
      write(63,'(A1,"Indoor error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Accumulated indoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_IN_ac(j,k)
      write(63,'(A1,"Accumulated indoor error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Outdoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_OUT(j,k)
      write(63,'(A1,"Outdoor error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Accumulated outdoor value",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_OUT_ac(j,k)
      write(63,'(A1,"Accumulated outdoor error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Accumulated total ELCR",A1,":",1x,ES12.5,",") )') aspas,aspas,CARC_TOT_ac(j,k)
      write(63,'(A1,"Accumulated total ELCR",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
	  ENDIF 
!
      ENDDO		! fim DO j=1,NTIME
!
!
	  ENDDO		! fim DO k=1,NLOCAL
!
!
	  WRITE(63,'("],")')
!
	  WRITE(63,'(A1,"External and internal Hazard Indexes",A1,": [")')aspas,aspas
!
!
      DO k=1,NLOCAL
	  DO j=NCOMESO,NTIME
!
      IF(HIND_EXT(j,k).NE.HIND_EXT(j-1,k))THEN
      IALL(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
	  IF(KEY_SD.EQV..TRUE.)THEN
!
      IF(IALL(k).EQ.1)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"External HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_EXT(j,k)
      write(63,'(A1,"External HI error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_HIND_EXT(j,k)
      write(63,'(A1,"Internal HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_INT(j,k)
      write(63,'(A1,"Internal HI error",A1,":",1x,ES12.5) )') aspas,aspas,SD_HIND_INT(j,k)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ELSEIF((IALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"External HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_EXT(j,k)
      write(63,'(A1,"External HI error",A1,":",1x,ES12.5,",") )') aspas,aspas,SD_HIND_EXT(j,k)
      write(63,'(A1,"Internal HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_INT(j,k)
      write(63,'(A1,"Internal HI error",A1,":",1x,ES12.5) )') aspas,aspas,SD_HIND_INT(j,k)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
      ENDIF
!
	  ELSE
!
      IF(IALL(k).EQ.1)THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(63,'(A1,"External HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_EXT(j,k)
      write(63,'(A1,"External HI error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Internal HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_INT(j,k)
      write(63,'(A1,"Internal HI error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ELSEIF((IALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(63,'("{")')
      write(63,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(63,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"External HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_EXT(j,k)
      write(63,'(A1,"External HI error",A1,":",1x,A1,"null",A1,",") )') aspas,aspas,aspas,aspas
      write(63,'(A1,"Internal HI value",A1,":",1x,ES12.5,",") )') aspas,aspas,HIND_INT(j,k)
      write(63,'(A1,"Internal HI error",A1,":",1x,A1,"null",A1) )') aspas,aspas,aspas,aspas
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(63,'("}")')
	  ELSE
	  WRITE(63,'("},")')
	  ENDIF
!
	  ENDIF
!
	  ENDIF 
!
      ENDDO		! fim DO j=1,NTIME
!
!
	  ENDDO		! fim DO k=1,NLOCAL
!
!
	  WRITE(63,'("]")')
!
	  WRITE(63,'("}")')
!
	  close(63)
!
	  RETURN
	  END
!
!
!
!*********************************************************************************************************************************************************
!*********************************************************************************************************************************************************
!
!
      SUBROUTINE INFORMATION(KKINFORMATION)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	   CHARACTER(LEN=10)  :: KKINFORMATION(6)
!
      WRITE(*,*)
	  WRITE(*,'("* ")')
	  WRITE(*,'("** ")')
	  WRITE(*,'("*** ")')
	  WRITE(*,'("****   See Information.out file created in the main folder ")')
	  WRITE(*,'("*** ")')
	  WRITE(*,'("** ")')
	  WRITE(*,'("* ")')
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*)
!
	  open (UNIT=39, file='Results\Information.txt')	
!
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
	  WRITE(39,'("                                                               INFORMATIONS ")')
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
!
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,'(" This file provides more information about the parameters present in the input sheets: ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Scenary ")')
	  WRITE(39,'(" - Concentration ")')
	  WRITE(39,'(" - Datachemical ")')
	  WRITE(39,'(" - Dataexp ")')
	  WRITE(39,'(" - Dataecological ")')
	  WRITE(39,*)
	  WRITE(39,*)
!
      IF((KKINFORMATION(1).EQ.'y').OR.(KKINFORMATION(1).EQ.'Y').OR.(KKINFORMATION(1).EQ.'yes').OR.(KKINFORMATION(1).EQ.'Yes').OR.(KKINFORMATION(1).EQ.'YES'))THEN
!
	  WRITE(39,'("                                         ******************************************************* ")')
	  WRITE(39,'("                                            INFORMATION ABOUT PARAMETERS OF THE SCENARY SHEET ")')
	  WRITE(39,'("                                         ******************************************************* ")')
	  WRITE(39,*)
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - SCENERY key = AGRICULTURAL   ---> The calculation will be performed in the AGRICULTURAL scenario ")')
	  WRITE(39,'("                                     The exposure routes choice must be made in section --- SHEME TO BE USED IF KEY SCENERY = 1 (AGRICULTURAL) --- ")')
	  WRITE(39,*)
	  WRITE(39,'(" - SCENERY key = INDUSTRIAL     ---> The calculation will be performed in the INDUSTRIAL scenario ")')
	  WRITE(39,'("                                     The exposure routes choice must be made in section --- SHEME TO BE USED IF KEY SCENERY = 2 (INDUSTRIAL) --- ")')
	  WRITE(39,*)
	  WRITE(39,'(" - SCENERY key = RESIDENTIAL    ---> The calculation will be performed in the RESIDENTIAL scenario ")')
	  WRITE(39,'("                                     The exposure routes choice must be made in section --- SHEME TO BE USED IF KEY SCENERY = 3 (RESIDENTIAL) --- ")')
	  WRITE(39,*)
	  WRITE(39,'(" - SCENERY key = IN NATURA      ---> The calculation will be performed in the IN NATURA scenario ")')
	  WRITE(39,'("                                     In this case, the human health risk will not be realized, as it is assumed that the risk assessment ")')
	  WRITE(39,'("                                     is being realized in an uninhabited area. ")')
	  WRITE(39,*)
	  WRITE(39,'("                                     Examples of possible scenarios IN NATURA: ")')
	  WRITE(39,'("                                     - Assessment of a mine in an uninhabited area  ")')
	  WRITE(39,'("                                     - Evaluation of a dump in an uninhabited area  ")')
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Num. of Chemicals ---> Inform the chemical species number that will be used in the risk assessments ")')
	  WRITE(39,*)
	  WRITE(39,'("                          If concentrations for 3 different chemical species are provided in the concentration sheet ")')
	  WRITE(39,'("                          then the value of --- Num. of Chemicals --- must be 3 ")')
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Num. of Times     ---> Inform the number of Times for which chemical species concentrations will be provided in the concentration sheet ")')
	  WRITE(39,*)
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Num. of Sites     ---> Inform the number of locations for which risk assessments will be realize   ")')
	  WRITE(39,*)
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Uncertainties key ")')
	  WRITE(39,*)
	  WRITE(39,'("        .TRUE.       ---> The uncertainties of all output parameters (HQ, CR, ...) will be calculated ")')
	  WRITE(39,'("        .FALSE.      ---> The uncertainties of all output parameters (HQ, CR, ...) will NOT be calculated ")')
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Key of NON-RADIOLOGICAL risk to human health ")')
	  WRITE(39,*)
	  WRITE(39,'("        .TRUE.       ---> The NON-RADIOLOGICAL human health risk assessment will be realized ")')
	  WRITE(39,'("        .FALSE.      ---> The NON-RADIOLOGICAL human health risk assessment will NOT be realized ")')
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Key of RADIOLOGICAL risk to human health ")')
	  WRITE(39,*)
	  WRITE(39,'("        .TRUE.       ---> The RADIOLOGICAL human health risk assessment will be realized ")')
	  WRITE(39,'("        .FALSE.      ---> The RADIOLOGICAL human health risk assessment will NOT be realized ")')
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Key of ECOLOGICAL risk to human health ")')
	  WRITE(39,*)
	  WRITE(39,'("        .TRUE.       ---> The ECOLOGICAL human health risk assessment will be realized ")')
	  WRITE(39,'("        .FALSE.      ---> The ECOLOGICAL human health risk assessment will NOT be realized ")')
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,'(" Example of enabling or disabling exposure pathway: ")')
	  WRITE(39,*)
	  WRITE(39,'(" WAY 1 = .TRUE.      = SOIL ----> HUMAN             The exposure route of accidental soil ingestion ")')
	  WRITE(39,'("                                                    will be considered for the risk calculation ")')
	  WRITE(39,*)
	  WRITE(39,'(" WAY 1 = .FALSE.     = SOIL ----> HUMAN             The exposure route of accidental soil ingestion ")')
	  WRITE(39,'("                                                    will NOT be considered for the risk calculation ")')
	  WRITE(39,*)
	  WRITE(39,'("-------------------- ")')
	  WRITE(39,'(" - Exposure duration assessment = Acute       ---> Human health risks will be calculated according to the ACUTE exposure pattern ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Exposure duration assessment = Subchronic  ---> Human health risks will be calculated according to the SUBCHRONIC exposure pattern ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Exposure duration assessment = Chronic     ---> Human health risks will be calculated according to the CHRONIC exposure pattern ")')
	  WRITE(39,*)
	  WRITE(39,'("                                                   How to assess the exposure pattern? ")')
	  WRITE(39,*)
	  WRITE(39,'("                                                   See flowchart on page 16 of the U.S. EPA 2009 reference: ")')
	  WRITE(39,'("                                                   U.S. EPA. Risk assessment guidance for superfund volume I: human health evaluation manual ")')
	  WRITE(39,'("                                                   Part F, Supplemental guidance for inhalation risk assessment, Washington DC, 2009. ")')
	  WRITE(39,*)
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
	  WRITE(39,*)
	  WRITE(39,*)
!
      ENDIF	    ! FIM DO IF KKINFORMATION(1)
!
!
      IF((KKINFORMATION(2).EQ.'y').OR.(KKINFORMATION(2).EQ.'Y').OR.(KKINFORMATION(2).EQ.'yes').OR.(KKINFORMATION(2).EQ.'Yes').OR.(KKINFORMATION(2).EQ.'YES'))THEN
!
	  WRITE(39,'("                                      ************************************************************* ")')
	  WRITE(39,'("                                         INFORMATION ABOUT PARAMETERS OF THE CONCENTRATION SHEET ")')
	  WRITE(39,'("                                      ************************************************************* ")')
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,'(" Parameters explanation through the example: ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Concentrations and respective uncertainties of the chemical species Benzo[a]pyrene and Lead ")')
	  WRITE(39,'(" - Concentrations and respective uncertainties in soil and drinking water matrices ")')
	  WRITE(39,'(" - Concentrations and respective uncertainties for 3 different locations ")')
	  WRITE(39,'(" - Concentrations and respective uncertainties in 2 different years (chronic risk) or 2 different events (acute risk) ")')
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,'("-------------------------------------------------------------------------------------")')
	  WRITE(39,'("                             Input Concentration Data")')
	  WRITE(39,'("-------------------------------------------------------------------------------------")')
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,'("                                   Chemical species 1")')
	  WRITE(39,'("                                   Benzo[a]pyrene           ---> Provide the name of the chemical species (the name must be provided as outlined in the HHRISK guide)")')
	  WRITE(39,'("                                   Matrix")')
	  WRITE(39,'("                                   SOIL                     ---> Provide the concentration matrix name. The possible names are: SOIL; DRINKING_WATER; BATH_WATER; OTHER_WATERS; SEDIMENTS; PARTICULATE; STEAM; FRUIT; LEAVES; MEAT; MILK; BIRD; EGG; FISH; GRAIN ")')
	  WRITE(39,'("Line 1: Values time 1        --->  10.47    20.62    4.68   ---> Each column represents a location. Column 1 = Location 1, Column 2 = Location 2, Column 3 = Location 3 .... ")')
	  WRITE(39,'("Line 2: Uncertainties time 1 --->   1.65     3.92    0.55    ")')
	  WRITE(39,'("Line 3: Values time 2        --->  11.62    34.09    8.68   ---> Odd lines are for providing concentration values at different times and locations ")')
	  WRITE(39,'("Line 4: Uncertainties time 2 --->   2.00     4.12    2.92   ---> Pairs lines are to provide the uncertainties of the concentration values at different times and locations ")')
	  WRITE(39,'("  Matrix")')
	  WRITE(39,'("  DRINKING_WATER ")')
	  WRITE(39,'("  1.00      4.78      12.57  ")')
	  WRITE(39,'("  0.08      0.42       3.11  ")')
	  WRITE(39,'("  2.13      5.86      21.95  ")')
	  WRITE(39,'("  0.44      0.75       5.21  ")')
	  WRITE(39,'("                                                           <----- There must be only 1 space between one chemical species and another")')
	  WRITE(39,'("  Chemical species 2")')
	  WRITE(39,'("  Pb      ")')
	  WRITE(39,'("  Matrix")')
	  WRITE(39,'("  SOIL ")')
	  WRITE(39,'("  2.00      9.66      23.65  ")')
	  WRITE(39,'("  0.16      1.07       4.22  ")')
	  WRITE(39,'("  4.26     10.43      42.63  ")')
	  WRITE(39,'("  0.88      1.50       5.33  ")')
	  WRITE(39,'("  Matrix")')
	  WRITE(39,'("  DRINKING_WATER ")')
	  WRITE(39,'("  4.00      6.77      12.57  ")')
	  WRITE(39,'("  0.54      1.00       3.11  ")')
	  WRITE(39,'("  7.33     15.86      16.19  ")')
	  WRITE(39,'("  0.99      2.57       2.37  ")')
	  WRITE(39,*)
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
!
	  WRITE(39,*)
	  WRITE(39,*)
!
      ENDIF	    ! FIM DO IF KKINFORMATION(2)
!
      IF((KKINFORMATION(3).EQ.'y').OR.(KKINFORMATION(3).EQ.'Y').OR.(KKINFORMATION(3).EQ.'yes').OR.(KKINFORMATION(3).EQ.'Yes').OR.(KKINFORMATION(3).EQ.'YES'))THEN
!
	  WRITE(39,'("                                       ************************************************************ ")')
	  WRITE(39,'("                                          INFORMATION ABOUT PARAMETERS OF THE DATACHEMICAL SHEET ")')
	  WRITE(39,'("                                       ************************************************************ ")')
	  WRITE(39,*)
	  WRITE(39,'("--------------  ")')
	  WRITE(39,'(" - NPOL               ---> Provide the number of chemical species present in the database ")')
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - KEY_SF = Disable   ---> The Slope Factors uncertainties will be those provided manually by the user ")')
	  WRITE(39,*)
	  WRITE(39,'(" - KEY_SF = Active    ---> The Slope Factors uncertainties will be calculated using the lognormal ")')
	  WRITE(39,'("                           statistical distribution as described by SASSI et al. (2007) ")')
	  WRITE(39,*)
	  WRITE(39,'("                           SASSI, G. et al. Quantitative estimation of uncertainty in human risk analysis.  ")')
	  WRITE(39,'("                           Journal of hazardous materials, v. 145, n. 1-2, p. 296-304, 2007.  ")')
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - First parameter is the name of the chemical species ---> For example: Cd, Fe, Pb, Benzo[a]pyrene, 9-Nitrophenanthrene, ... ")')
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - Second parameter defines if the chemical species present mutagenic mode of action when ingested or inhaled ")')
	  WRITE(39,*)
	  WRITE(39,'("    .TRUE.      ---> The chemical species present a mutagenic mode of action ")')
	  WRITE(39,'("    .FALSE.     ---> The chemical species NOT present a mutagenic mode of action ")')
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - Values and respective uncertainties of Relative bioavailability adjustment - BAF - ")')
	  WRITE(39,*)
	  WRITE(39,'("   The fraction of the chemical species absorbed by the gastrointestinal system ")')
	  WRITE(39,'("    1st  value = Soil ingestion ")')
	  WRITE(39,'("    2nd  value = Water ingestion ")')
	  WRITE(39,*)
	  WRITE(39,'("   The fraction of the chemical species absorbed by lungs ")')
	  WRITE(39,'("    3rd  value = Inhalation of particulate matter ")')
	  WRITE(39,'("    4th  value = Steam inhalation ")')
	  WRITE(39,*)
	  WRITE(39,'("   The fraction of the chemical species absorbed by skin ")')
	  WRITE(39,'("    5th  value = Dermal contact with soil  ---> This value must be 1.000 because the ABS value present in the dose calculation  ")')
	  WRITE(39,'("                                                by dermal contact with soil already does the BAF function ")')
	  WRITE(39,'("    6th  value = Dermal contact with water ---> This value must be 1.000 because the PC value present in the dose calculation ")')
	  WRITE(39,'("                                                by dermal contact with water already does the BAF function ")')
	  WRITE(39,*)
	  WRITE(39,'("   The fraction of the chemical species absorbed by the gastrointestinal system ")')
	  WRITE(39,'("    7th  value = Vegetables ingestion ")')
	  WRITE(39,'("    8th  value = Fruits ingestion ")')
	  WRITE(39,'("    9th  value = Bovine meat ingestion ")')
	  WRITE(39,'("   10th  value = Bovine milk ingestion ")')
	  WRITE(39,'("   11th  value = Bird meat ingestion ")')
	  WRITE(39,'("   12th  value = Eggs ingestion ")')
	  WRITE(39,'("   13th  value = Fish ingestion ")')
	  WRITE(39,'("   14th  value = Grains ingestion ")')
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - Values and respective uncertainties of Reference Doses - RfD ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st column value = RfD for ingestion of soil or food ")')
	  WRITE(39,'("    2nd column value = RfD for inhalation of particulate matter ")')
	  WRITE(39,'("    3rd column value = RfD for dermal contact with soil ")')
	  WRITE(39,'("    4th column value = RfD for ingestion of water ")')
	  WRITE(39,'("    5th column value = RfD for vapor inhalation ")')
	  WRITE(39,'("    6th column value = RfD for dermal contact with water  ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st line value = Chronic RfD values ")')
	  WRITE(39,'("    2nd line value = Subchronic RfD values ")')
	  WRITE(39,'("    3rd line value = Acute RfD values ")')
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - Values and respective uncertainties of Cancer Slope Factor - SF ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st column value = SF for ingestion of soil or food ")')
	  WRITE(39,'("    2nd column value = SF for inhalation of particulate matter ")')
	  WRITE(39,'("    3rd column value = SF for dermal contact with soil ")')
	  WRITE(39,'("    4th column value = SF for ingestion of water ")')
	  WRITE(39,'("    5th column value = SF for vapor inhalation ")')
	  WRITE(39,'("    6th column value = SF for dermal contact with water  ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st line value = Chronic SF values ")')
	  WRITE(39,'("    2nd line value = Subchronic SF values ")')
	  WRITE(39,'("    3rd line value = Acute SF values ")')
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - Varied parameters specific to each chemical species ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st  value = PC = Dermal permeability of the chemical species - Dermal contact with water")')
	  WRITE(39,'("    2nd  value = ABS = Absorption factor of the chemical species - Dermal contact with soil ")')
	  WRITE(39,'("    3rd  value = Kd = Soil-water partition coefficient - parameter used together with BTF ")')
	  WRITE(39,'("    4th  value = fw = Daily fraction of consumed water that is contaminated - parameter used together with BTF ")')	 
	  WRITE(39,*)
      WRITE(39,'("-------------- ")')
	  WRITE(39,'(" - Values and respective uncertainties of Biotransfer Factors - BTF ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st  value = Soil ---> Fruit ")')
	  WRITE(39,'("    2nd  value = Soil ---> Feed plants ")')
	  WRITE(39,'("    3rd  value = Feed plants ---> Bovine meat ")')
	  WRITE(39,'("    4th  value = Soil ---> Bovine meat ")')
	  WRITE(39,'("    5th  value = Feed plants ---> Bovine milk ")')
	  WRITE(39,'("    6th  value = Soil ---> Bovine milk ")')
	  WRITE(39,'("    7th  value = Soil ---> Vegetables ")')
	  WRITE(39,'("    8th  value = Water ---> Bovine meat ")')
	  WRITE(39,'("    9th  value = Water ---> Bovine milk  ")')
	  WRITE(39,'("   10th  value = Water ---> Bovine fish ")')
	  WRITE(39,'("   11th  value = Soil ---> Bird meat ")')
	  WRITE(39,'("   12th  value = Water ---> Bird meat  ")')
	  WRITE(39,'("   13th  value = Soil ---> Egg ")')
	  WRITE(39,'("   14th  value = Water ---> Egg ")')
	  WRITE(39,'("   15th  value = Soil ---> Grains   ")')
	  WRITE(39,*)
	  WRITE(39,*)
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
!
	  WRITE(39,*)
	  WRITE(39,*)
!
      ENDIF	    ! FIM DO IF KKINFORMATION(3)
!
      IF((KKINFORMATION(4).EQ.'y').OR.(KKINFORMATION(4).EQ.'Y').OR.(KKINFORMATION(4).EQ.'yes').OR.(KKINFORMATION(4).EQ.'Yes').OR.(KKINFORMATION(4).EQ.'YES'))THEN
!
	  WRITE(39,'("                                         ******************************************************* ")')
	  WRITE(39,'("                                            INFORMATION ABOUT PARAMETERS OF THE DATAEXP SHEET ")')
	  WRITE(39,'("                                         ******************************************************* ")')
	  WRITE(39,*)
	  WRITE(39,'("---------------------------  ")')
	  WRITE(39,'(" - Exposure duration  ED     ---> Defines how many years a person has been exposed to the contaminant ")')
	  WRITE(39,'("                                  These parameters will only be used for chronic risks ")')
	  WRITE(39,*)
	  WRITE(39,'("                                  Different ED values can be provided for each age group and scenario ")')
      WRITE(39,'("--------------------------- ")')
	  WRITE(39,'(" - Averaging time       AT   ---> Period over which exposure is averaged ")')
	  WRITE(39,'("   Carionogenic risks             AT for Carcinogenic effects should be provided in total days - Life expectancy x 365 days ")')
	  WRITE(39,'("                                  U.S. EPA (2011) considers the average life expectancy of 78 years. Therefore, carcinogenic AT must be 28470 days ")')
	  WRITE(39,*)
	  WRITE(39,'("                                  U.S. EPA, Exposure Factors Handbook. EPA/600/R-090/052F. 2011. ")')
	  WRITE(39,*)
	  WRITE(39,'("                                  Different AT - Carionogenic risks - values can be provided for each age group and scenario ")')
      WRITE(39,'("--------------------------- ")')
	  WRITE(39,'(" - Body weight          BW   ---> Body weight in kilograms for different age groups ")')
	  WRITE(39,*)
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
!
	  WRITE(39,*)
	  WRITE(39,*)
!
      ENDIF	    ! FIM DO IF KKINFORMATION(4)
!
!
      IF((KKINFORMATION(5).EQ.'y').OR.(KKINFORMATION(5).EQ.'Y').OR.(KKINFORMATION(5).EQ.'yes').OR.(KKINFORMATION(5).EQ.'Yes').OR.(KKINFORMATION(5).EQ.'YES'))THEN
!
	  WRITE(39,'("                                      ************************************************************** ")')
	  WRITE(39,'("                                         INFORMATION ABOUT PARAMETERS OF THE DATAECOLOGICAL SHEET ")')
	  WRITE(39,'("                                      ************************************************************** ")')
	  WRITE(39,*)
	  WRITE(39,'("----------------------------  ")')
	  WRITE(39,'(" - Total chemical Num.         ---> Provide the number of chemical species present in the database ")')
	  WRITE(39,*)
	  WRITE(39,'("---------------------------- ")')
	  WRITE(39,'(" - Water reference value = National   ---> The National water reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Water reference value = U.S EPA    ---> The U.S.EPA water reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Water reference value = WHO        ---> The WHO water reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'("---------------------------- ")')
	  WRITE(39,'(" - Soil reference value = Regional    ---> The Regional soil reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'("                          =*= Res.    ---> It will be used if the SCENARIO used is Residential  - key SCENARIO = 3 - ")')
	  WRITE(39,'("                          =*= Agr.    ---> It will be used if the SCENARIO used is Agricultural - key SCENARIO = 1 - ")')
	  WRITE(39,'("                          =*= Ind.    ---> It will be used if the SCENARIO used is Industrial   - key SCENARIO = 2 - ")')
	  WRITE(39,'("                          =*= In Nat  ---> It will be used if the SCENARIO used is In Nature    - key SCENARIO = 4 - ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Soil reference value = National    ---> The National soil reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Soil reference value = U.S. EPA    ---> The U.S.EPA soil reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - Sediment reference value = Regional  ---> The Regional sediment reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Sediment reference value = National  ---> The National sediment reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'(" - Sediment reference value = U.S. EPA  ---> The U.S.EPA sediment reference values will be used for the ecological indexes calculation. ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - First parameter is the name of the chemical species ---> For example: Cd, Fe, Pb, Benzo[a]pyrene, 9-Nitrophenanthrene, ... ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - Water Reference values for the chemical specie ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st  value = Water reference value established by the country where the ecological risk assessment is realized ")')
	  WRITE(39,'("    2nd  value = Water reference value established by the U.S. Environmental Protection Agency - U.S.EPA - ")')
	  WRITE(39,'("    3rd  value = Water reference value established by the World Health Organization - WHO - ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - Soil Reference values for the chemical specie ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st  value = Residential Soil Reference Value established by the Region or State where the ecological risk assessment is realized ")')
	  WRITE(39,'("    2nd  value = Agricultural Soil Reference Value established by the Region or State where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    3rd  value = Industrial Soil Reference Value established by the Region or State where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    4th  value = Prevention Soil Reference Value established by the Region or State where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    5th  value = Residential Soil Reference Value established by the country where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    6th  value = Agricultural Soil Reference Value established by the country where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    7th  value = Industrial Soil Reference Value established by the country where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    8th  value = Prevention Soil Reference Value established by the country where the ecological risk assessment is realized  ")')
	  WRITE(39,'("    9th  value = Residential Soil Reference Value established by the U.S. Environmental Protection Agency - U.S.EPA - ")')
	  WRITE(39,'("   10th  value = Agricultural Soil Reference Value established by the U.S. Environmental Protection Agency - U.S.EPA - ")')
	  WRITE(39,'("   11th  value = Industrial Soil Reference Value established by the U.S. Environmental Protection Agency - U.S.EPA - ")')
	  WRITE(39,'("   12th  value = Prevention Soil Reference Value established by the U.S. Environmental Protection Agency - U.S.EPA - ")')
	  WRITE(39,*)
	  WRITE(39,'("   13th  value = Reference values of the chemical species in the earth s crust ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - Sediment Reference values for the chemical specie ")')
	  WRITE(39,*)
	  WRITE(39,'("    1st  value = Sediment reference value established by the Region or State where the ecological risk assessment is realized ")')
	  WRITE(39,'("    2nd  value = Sediment reference value established by the country where the ecological risk assessment is realized ")')
	  WRITE(39,'("    3rd  value = Sediment reference value established by the U.S. Environmental Protection Agency - U.S.EPA - ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - OTHER PARAMETERS: ")')
	  WRITE(39,*)
	  WRITE(39,'("    =*= TEL        ---> Threshold Effect Level.   ")')
	  WRITE(39,'("    =*= PEL        ---> Probable Effect Level.   ")')
	  WRITE(39,'("    =*= ERM        ---> Effects Range-Maximum. ")')
	  WRITE(39,'("    =*= Tir        ---> Toxic response factor.   ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,'(" - PARAMETERS FOR CALCULATING TOXIC PRESSURE: ")')
	  WRITE(39,*)
	  WRITE(39,'("    The Alpha and Beta parameters are related to the chemical species sensitivity distribution curve (SSD).    ")')
	  WRITE(39,*)
	  WRITE(39,'("    =*= Alpha   --->  Corresponds to the average value of this distribution and is related to HC50.   ")')
	  WRITE(39,'("    =*= Beta    --->  Related to the width of the distribution itself. ")')
	  WRITE(39,*)
	  WRITE(39,'("------------------------------- ")')
	  WRITE(39,*)
	  WRITE(39,'("******************************************************************************************************************************************************************************")')
!
	  WRITE(39,*)
	  WRITE(39,*)
!
      ENDIF	    ! FIM DO IF KKINFORMATION(5)
!
!	  
      close(39)
!
      stop
!
 	  RETURN
	  END
!
!
!
!*********************************************************************************************************************************************************
!*********************************************************************************************************************************************************
!
!
      SUBROUTINE NOMINATION(CHEMICAL,NCHEM,NOMES)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
	   CHARACTER(LEN=50)  :: CHEMICAL(500),NOMES(500) 
!
!
!
      DO I=1,NCHEM
!
      IF(CHEMICAL(I).EQ.'p-Benzoquinone')THEN
	  NOMES(I)='p-BQ'
	  ELSEIF(CHEMICAL(I).EQ.'Acetophenone')THEN
	  NOMES(I)='APhe'
	  ELSEIF(CHEMICAL(I).EQ.'Naphthalene')THEN
	  NOMES(I)='Nap'
	  ELSEIF(CHEMICAL(I).EQ.'Naphthoquinone')THEN
	  NOMES(I)='Nap-Q'
	  ELSEIF(CHEMICAL(I).EQ.'Acenaphthylene')THEN
	  NOMES(I)='Acy'
	  ELSEIF(CHEMICAL(I).EQ.'Acenaphthene')THEN
	  NOMES(I)='Ace'
	  ELSEIF(CHEMICAL(I).EQ.'Fluorene')THEN
	  NOMES(I)='Flu'
	  ELSEIF(CHEMICAL(I).EQ.'2-Nitrobiphenyl')THEN
	  NOMES(I)='2-NBi'
	  ELSEIF(CHEMICAL(I).EQ.'Phenanthrene')THEN
	  NOMES(I)='Phe'
	  ELSEIF(CHEMICAL(I).EQ.'Anthracene')THEN
	  NOMES(I)='Ant'
	  ELSEIF(CHEMICAL(I).EQ.'5-Nitroacenaphthene')THEN
	  NOMES(I)='5-NAce'
	  ELSEIF(CHEMICAL(I).EQ.'Fluoranthene')THEN
	  NOMES(I)='Flt'
	  ELSEIF(CHEMICAL(I).EQ.'2-Nitrofluorene')THEN
	  NOMES(I)='2-NFlu'
	  ELSEIF(CHEMICAL(I).EQ.'Pyrene')THEN
	  NOMES(I)='Pyr'
	  ELSEIF(CHEMICAL(I).EQ.'Phenanthrenequinone')THEN
	  NOMES(I)='PheQ'
	  ELSEIF(CHEMICAL(I).EQ.'Retene')THEN
	  NOMES(I)='Ret'
	  ELSEIF(CHEMICAL(I).EQ.'9-Nitrophenanthrene')THEN
	  NOMES(I)='9-NPhe'
	  ELSEIF(CHEMICAL(I).EQ.'9-Nitroantracene')THEN
	  NOMES(I)='9-NAnt'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[a]fluorenone')THEN
	  NOMES(I)='B[a]F-one'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[a]anthracene')THEN
	  NOMES(I)='B[a]A'
	  ELSEIF(CHEMICAL(I).EQ.'Chrysene')THEN
	  NOMES(I)='Chr'
	  ELSEIF(CHEMICAL(I).EQ.'1-Nitropyrene')THEN
	  NOMES(I)='1-NPyr'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[b]fluoranthene')THEN
	  NOMES(I)='B[b]F'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[k]fluoranthene')THEN
	  NOMES(I)='B[k]F'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[e]pyrene')THEN
	  NOMES(I)='B[e]P'
	  ELSEIF(CHEMICAL(I).EQ.'Anthraquinone')THEN
	  NOMES(I)='AntQ'
	  ELSEIF(CHEMICAL(I).EQ.'3-Nitrobenzanthrone')THEN
	  NOMES(I)='3-NBAnt'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[a]pyrene')THEN
	  NOMES(I)='B[a]P'
	  ELSEIF(CHEMICAL(I).EQ.'6H-Benzo[cd]pyren-6-one')THEN
	  NOMES(I)='6H-BcdP'
	  ELSEIF(CHEMICAL(I).EQ.'Indeno[123-cd]pyrene')THEN
	  NOMES(I)='IcdP'
	  ELSEIF(CHEMICAL(I).EQ.'6-Nitrobenzo[a]pyrene')THEN
	  NOMES(I)='6-NBaP'
	  ELSEIF(CHEMICAL(I).EQ.'Dibenz[ah]anthracene')THEN
	  NOMES(I)='DahA'
	  ELSEIF(CHEMICAL(I).EQ.'Benzo[ghi]perylene')THEN
	  NOMES(I)='BghiP'
	  ELSEIF(CHEMICAL(I).EQ.'6-Nitrochrysene')THEN
	  NOMES(I)='6-NChr'
	  ELSE
	  NOMES(I)=CHEMICAL(I)
!
      ENDIF
!
      ENDDO  ! FIM DO(I)
!
 	  RETURN
	  END
!
!
!
!*********************************************************************************************************************************************************
!*********************************************************************************************************************************************************
!
      SUBROUTINE ECOLOGICAL(SCENAR,VARIASAO,NDURATION,NCHEM,NTIME,NLOCAL,NTYPECONC,CHEMICAL,SCENARIES)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
    INTEGER SCENAR
	REAL MUL_RWind,MUL_CFwater,MPIwat,IPIT_water,MUL_CFsoil,mCd_soil,IPIT_soil,IGEOsoil,MAX_PIsoil,MAX_PIsed
	REAL MUL_CFsed,mCd_sed,IPIT_sed,IGEOsed,mPELq,mERMq,Kd_MPI,MUL_RISKwat,MUL_RISKsoil,MUL_RISKsed,IRjFIN
!
	  LOGICAL KEYCONC
	   CHARACTER(LEN=50)  :: CHEMICAL(500),SPECIE(500),AGUA_REF(3),SOLO_REF(4),SED_REF(3),NOMES(500)
	   CHARACTER(LEN=14)  :: SCENARIES(4)
!
	CHARACTER EFREFsoil*2,EFREFsediment*2,aspas*1
	DIMENSION CSOIL(NCHEM,NTIME,NLOCAL),CWATER(NCHEM,NTIME,NLOCAL),CWATERDER(NCHEM,NTIME,NLOCAL),CWATEROTHER(NCHEM,NTIME,NLOCAL)
	DIMENSION CPAR(NCHEM,NTIME,NLOCAL),CSTEAM(NCHEM,NTIME,NLOCAL),CFRUIT(NCHEM,NTIME,NLOCAL),CLEAVES(NCHEM,NTIME,NLOCAL)
	DIMENSION CBEEF(NCHEM,NTIME,NLOCAL),CMILK(NCHEM,NTIME,NLOCAL),CAVE(NCHEM,NTIME,NLOCAL)
	DIMENSION CEGG(NCHEM,NTIME,NLOCAL),CFISH(NCHEM,NTIME,NLOCAL),CGRAIN(NCHEM,NTIME,NLOCAL),CSEDIMENT(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CSOIL(NCHEM,NTIME,NLOCAL),SD_CWATER(NCHEM,NTIME,NLOCAL),SD_CWATERDER(NCHEM,NTIME,NLOCAL),SD_CWATEROTHER(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CPAR(NCHEM,NTIME,NLOCAL),SD_CSTEAM(NCHEM,NTIME,NLOCAL),SD_CFRUIT(NCHEM,NTIME,NLOCAL),SD_CLEAVES(NCHEM,NTIME,NLOCAL)
		DIMENSION SD_CBEEF(NCHEM,NTIME,NLOCAL),SD_CMILK(NCHEM,NTIME,NLOCAL),SD_CAVE(NCHEM,NTIME,NLOCAL),SD_CEGG(NCHEM,NTIME,NLOCAL)
    	DIMENSION SD_CFISH(NCHEM,NTIME,NLOCAL),SD_CGRAIN(NCHEM,NTIME,NLOCAL),SD_CSEDIMENT(NCHEM,NTIME,NLOCAL)
	DIMENSION KEYCONC(NCHEM,NTYPECONC),TIMESP(0:NTIME)  ! NTYPECONC SÃO OS 14 TIPOS DE CONCENTRAÇOES ACIMA !!!!
!
      DIMENSION WATREF(0:500),CROSTA(500),SOILREF(0:500),SEDREF(0:500),TEL(500),PEL(500),ERM(500),TIR(500) 
	  DIMENSION ALFA(500),BETA(500),SEDREFNAC(500),SOILREFNAC(500),WATREFNAC(500)
!
      DIMENSION CFwater(NCHEM,0:NTIME,NLOCAL),RWind(NCHEM,NTIME,NLOCAL),MUL_CFwater(NTIME,NLOCAL),MPIwat(0:NTIME,NLOCAL)
	  DIMENSION	MUL_RWind(NTIME,NLOCAL),RWcomb(0:NTIME,NLOCAL),PASS_RWind(NCHEM,NTIME,NLOCAL),CFwaterNAC(NCHEM,NTIME,NLOCAL)
	  DIMENSION S_CFwaterNAC_water(NTIME,NLOCAL),IPIT_water(0:NTIME,NLOCAL)
	  DIMENSION	CFsoil(NCHEM,0:NTIME,NLOCAL),MUL_CFsoil(NTIME,NLOCAL),PLIsoil(0:NTIME,NLOCAL),SOMA_mCd_soil(NTIME,NLOCAL)
	  DIMENSION mCd_soil(NTIME,NLOCAL),TEMP_PERIsoil(NCHEM,NTIME,NLOCAL),PERIsoil(NTIME,NLOCAL),CFsoilNAC(NCHEM,NTIME,NLOCAL)
	  DIMENSION S_CFsoilNAC_soil(NTIME,NLOCAL),IPIT_soil(NTIME,NLOCAL),IGEOsoil(NCHEM,NTIME,NLOCAL),PIsoil(NCHEM,NTIME,NLOCAL)
	  DIMENSION EFCSOIL(0:NTIME,NLOCAL),EFCSEDIMENT(0:NTIME,NLOCAL),EFsoil(NCHEM,NTIME,NLOCAL),PINEWsoil(NTIME,NLOCAL)
	  DIMENSION P1_PINEWsoil(NTIME,NLOCAL),P2_PINEWsoil(NTIME,NLOCAL),COPIA_PIsoil(NCHEM,NTIME,NLOCAL),MAX_PIsoil(NTIME,NLOCAL)
!
	  DIMENSION	CFsed(NCHEM,0:NTIME,NLOCAL),MUL_CFsed(NTIME,NLOCAL),PLIsed(0:NTIME,NLOCAL),SOMA_mCd_sed(NTIME,NLOCAL)
	  DIMENSION mCd_sed(NTIME,NLOCAL),TEMP_PERIsed(NCHEM,NTIME,NLOCAL),PERIsed(NTIME,NLOCAL),CFsedNAC(NCHEM,NTIME,NLOCAL)
	  DIMENSION S_CFsedNAC_sed(NTIME,NLOCAL),IPIT_sed(NTIME,NLOCAL),IGEOsed(NCHEM,NTIME,NLOCAL),PIsed(NCHEM,NTIME,NLOCAL)
	  DIMENSION EFsed(NCHEM,NTIME,NLOCAL),PINEWsed(NTIME,NLOCAL),P1_PINEWsed(NTIME,NLOCAL),P2_PINEWsed(NTIME,NLOCAL)
	  DIMENSION COPIA_PIsed(NCHEM,NTIME,NLOCAL),MAX_PIsed(NTIME,NLOCAL),mPELq(NTIME,NLOCAL),TEMP_mPELq(NCHEM,NTIME,NLOCAL)
	  DIMENSION SOMA_mPELq(NTIME,NLOCAL),mERMq(NTIME,NLOCAL),TEMP_mERMq(NCHEM,NTIME,NLOCAL),SOMA_mERMq(NTIME,NLOCAL)
	  DIMENSION TEMP_TRI(NCHEM,NTIME,NLOCAL),TRI(NTIME,NLOCAL),Kd_MPI(0:NTIME,NLOCAL),ALFABETAwat(NCHEM,NTIME,NLOCAL)
	  DIMENSION ALFABETAsoil(NCHEM,NTIME,NLOCAL),ALFABETAsed(NCHEM,NTIME,NLOCAL),TPIwat(NCHEM,NTIME,NLOCAL)
	  DIMENSION	TPIsoil(NCHEM,NTIME,NLOCAL),TPIsed(NCHEM,NTIME,NLOCAL),ALFABETAwatREF(NCHEM,NTIME,NLOCAL),TPBGIwat(NCHEM,NTIME,NLOCAL)
	  DIMENSION ALFABETAsoilREF(NCHEM,NTIME,NLOCAL),TPBGIsoil(NCHEM,NTIME,NLOCAL),ALFABETAsedREF(NCHEM,NTIME,NLOCAL),TPBGIsed(NCHEM,NTIME,NLOCAL)
	  DIMENSION TPICORwat(NCHEM,NTIME,NLOCAL),TPICORsoil(NCHEM,NTIME,NLOCAL),TPICORsed(NCHEM,NTIME,NLOCAL)
	  DIMENSION MUL_RISKwat(NTIME,NLOCAL),MUL_RISKsoil(NTIME,NLOCAL),MUL_RISKsed(NTIME,NLOCAL),SUB_RISKwat(NCHEM,NTIME,NLOCAL)
	  DIMENSION	SUB_RISKsoil(NCHEM,NTIME,NLOCAL),SUB_RISKsed(NCHEM,NTIME,NLOCAL),RISKwat(0:NTIME,NLOCAL),RISKsoil(0:NTIME,NLOCAL),RISKsed(0:NTIME,NLOCAL)
      DIMENSION IRjFIN(0:NTIME,NLOCAL)
!
      DIMENSION IALL(NLOCAL),JALL1(NLOCAL),JALL2(NLOCAL),JALL3(NLOCAL),JALL4(NLOCAL)
!
      CALL CONCENTRATION(VARIASAO,NDURATION,NCHEM,NTYPECONC,CHEMICAL,&
	  CSOIL,CWATER,CPAR,CSTEAM,CFRUIT,CLEAVES,CBEEF,CMILK,CAVE,CEGG,CFISH,CGRAIN,CWATERDER,CWATEROTHER,CSEDIMENT,&
	  KEYCONC,NTIME,NLOCAL,TIMESP,&
      SD_CSOIL,SD_CWATER,SD_CPAR,SD_CSTEAM,SD_CFRUIT,SD_CLEAVES,SD_CBEEF,SD_CMILK,SD_CAVE,SD_CEGG,SD_CFISH,SD_CGRAIN,SD_CWATERDER,SD_CWATEROTHER,SD_CSEDIMENT)
!
!
      CALL ECODATA(SCENAR,NSPECIE,SPECIE,WATREF,CROSTA,SOILREF,SEDREF,TEL,PEL,ERM,TIR,ALFA,BETA,SEDREFNAC,SOILREFNAC,WATREFNAC,KEY_WATER,KEY_SOIL,KEY_SEDIMENT)
!
!
!
	  DO J=0,NTIME
	  DO K=1,NLOCAL
      EFCSOIL(J,K)=0.0
	  EFCSEDIMENT(J,K)=0.0
	  ENDDO
	  ENDDO
! 	   
!
      nnn=0
	  mmm=0
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA								                  
	  IF ((CHEMICAL(i).EQ.SPECIE(16)).AND.(KEYCONC(i,1).EQV..TRUE.)) THEN	   ! VALOR DE REFERENCIA SENDO O Al
	  EFCROSTAsolo=CROSTA(15)
	  EFREFsoil=CHEMICAL(i)
	  nnn=1
!
      DO J=1,NTIME
	  DO K=1,NLOCAL
	  EFCSOIL(J,K)=CSOIL(i,J,K)
	  ENDDO
	  ENDDO
	  GOTO 3
	  ENDIF 
 	  ENDDO
!
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA								                  
	  IF ((CHEMICAL(i).EQ.SPECIE(7)).AND.(KEYCONC(i,1).EQV..TRUE.)) THEN	   ! VALOR DE REFERENCIA SENDO O Fe
	  EFCROSTAsolo=CROSTA(7)
	  EFREFsoil=CHEMICAL(i)
	  nnn=1
!
      DO J=1,NTIME
	  DO K=1,NLOCAL
	  EFCSOIL(J,K)=CSOIL(i,J,K)
	  ENDDO
	  ENDDO
	  GOTO 3
	  ENDIF 
 	  ENDDO
!
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA								                  
	  IF ((CHEMICAL(i).EQ.SPECIE(9)).AND.(KEYCONC(i,1).EQV..TRUE.)) THEN	   ! VALOR DE REFERENCIA SENDO O Mn
	  EFCROSTAsolo=CROSTA(9)
	  EFREFsoil=CHEMICAL(i)
	  nnn=1
!
      DO J=1,NTIME
	  DO K=1,NLOCAL
	  EFCSOIL(J,K)=CSOIL(i,J,K)
	  ENDDO
	  ENDDO
	  GOTO 3
	  ENDIF 
 	  ENDDO
!
3     CONTINUE !----------------------------------------------------------------------------------------------------------------
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA								                  
	  IF ((CHEMICAL(i).EQ.SPECIE(16)).AND.(KEYCONC(i,15).EQV..TRUE.)) THEN	   ! VALOR DE REFERENCIA SENDO O Al
	  EFCROSTAsed=CROSTA(15)
	  EFREFsediment=CHEMICAL(i)
	  mmm=1
!
      DO J=1,NTIME
	  DO K=1,NLOCAL
	  EFCSEDIMENT(J,K)=CSEDIMENT(i,J,K)
	  ENDDO
	  ENDDO
	  GOTO 22
	  ENDIF 
 	  ENDDO
!
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA								                  
	  IF ((CHEMICAL(i).EQ.SPECIE(7)).AND.(KEYCONC(i,15).EQV..TRUE.)) THEN	   ! VALOR DE REFERENCIA SENDO O Fe
	  EFCROSTAsed=CROSTA(7)
	  EFREFsediment=CHEMICAL(i)
	  mmm=1
!
      DO J=1,NTIME
	  DO K=1,NLOCAL
	  EFCSEDIMENT(J,K)=CSEDIMENT(i,J,K)
	  ENDDO
	  ENDDO
	  GOTO 22
	  ENDIF 
 	  ENDDO
!
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA								                  
	  IF ((CHEMICAL(i).EQ.SPECIE(9)).AND.(KEYCONC(i,15).EQV..TRUE.)) THEN	   ! VALOR DE REFERENCIA SENDO O Mn
	  EFCROSTAsed=CROSTA(9)
	  EFREFsediment=CHEMICAL(i)
	  mmm=1
!
      DO J=1,NTIME
	  DO K=1,NLOCAL
	  EFCSEDIMENT(J,K)=CSEDIMENT(i,J,K)
	  ENDDO
	  ENDDO
	  GOTO 22
	  ENDIF 
 	  ENDDO
!
22    CONTINUE
!

	  IF (nnn.eq.0) THEN
!
      EFCROSTAsolo=0.0
!
      EFREFsoil='  '   
!
      WRITE(*,*)
	  WRITE(*,'("  WARNING!!! Soil concentrations of reference elements (Fe, Al or Mn) ")')
	  WRITE(*,'("  were not provided in the concentration sheet of the input file.")')
	  WRITE(*,'("  The soil enrichment factor (EF) calculation will not be performed!")')
	  WRITE(*,*)   
	  WRITE(*,'("  The soil EF value will be zero  ")')
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! Soil concentrations of reference elements (Fe, Al or Mn) were not provided in the concentration sheet of the input file.")')
	  WRITE(99,'(" The soil enrichment factor (EF) calculation will not be performed!  ----> The soil EF value will be zero")')
      WRITE(99,*)
!
	  ENDIF
!
	  IF (mmm.eq.0) THEN
!
      EFCROSTAsed=0.0
!
	  EFREFsediment='  '   
!
      WRITE(*,*)
	  WRITE(*,'("  WARNING!!! Sediment concentrations of reference elements (Fe, Al or Mn) ")')
	  WRITE(*,'("  were not provided in the concentration sheet of the input file.")')
	  WRITE(*,'("  The sediment enrichment factor (EF) calculation will not be performed!")')
	  WRITE(*,*)   
	  WRITE(*,'("  The sediment EF value will be zero  ")')
      WRITE(*,*)
!      
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! Sediment concentrations of reference elements (Fe, Al or Mn) were not provided in the concentration sheet of the input file.")')
	  WRITE(99,'(" The soil enrichment factor (EF) calculation will not be performed!  ----> The soil EF value will be zero")')
      WRITE(99,*)

!
	  ENDIF

!	
!
	  DO J=1,NTIME
	  DO K=1,NLOCAL
	  IALL(k)=0
	  JALL1(k)=0
	  JALL2(k)=0
	  JALL3(k)=0
	  JALL4(k)=0
	  MUL_CFwater(J,K)=1.0
	  MPIwat(J,K)=0.0
	  MPIwat(0,K)=0.0
	  S_CFwaterNAC_water(J,K)=0.0
	  MUL_RWind(J,K)=1.0
	  RWcomb(J,K)=0.0
	  RWcomb(0,K)=0.0
	  IPIT_water(J,K)=0.0
	  IPIT_water(0,K)=0.0
!
	  MUL_CFsoil(J,K)=1.0
	  PLIsoil(J,K)=0.0
	  PLIsoil(0,K)=0.0
	  SOMA_mCd_soil(J,K)=0.0
	  mCd_soil(J,K)=0.0
	  PERIsoil(J,K)=0.0
	  S_CFsoilNAC_soil(J,K)=0.0
	  IPIT_soil(J,K)=0.0
	  PINEWsoil(J,K)=0.0
	  P1_PINEWsoil(J,K)=0.0
	  P2_PINEWsoil(J,K)=0.0
	  MAX_PIsoil(J,K)=0.0
!
	  MUL_CFsed(J,K)=1.0
	  PLIsed(J,K)=0.0
	  PLIsed(0,K)=0.0
	  SOMA_mCd_sed(J,K)=0.0
	  mCd_sed(J,K)=0.0
	  PERIsed(J,K)=0.0
	  S_CFsedNAC_sed(J,K)=0.0
	  IPIT_sed(J,K)=0.0
	  PINEWsed(J,K)=0.0
	  P1_PINEWsed(J,K)=0.0
	  P2_PINEWsed(J,K)=0.0
	  MAX_PIsed(J,K)=0.0
	  SOMA_mPELq=0.0
	  mPELq(J,K)=0.0
	  SOMA_mERMq=0.0
	  mERMq(J,K)=0.0
	  TRI(J,K)=0.0
	  Kd_MPI(J,K)=0.0
	  Kd_MPI(0,K)=0.0
	  MUL_RISKwat(J,K)=1.0	  
	  MUL_RISKsoil(J,K)=1.0
	  MUL_RISKsed(J,K)=1.0
	  RISKwat(J,K)=0.0
	  RISKwat(0,K)=0.0
	  RISKsoil(J,K)=0.0
	  RISKsoil(0,K)=0.0
	  RISKsed(J,K)=0.0
	  RISKsed(0,K)=0.0
	  IRjFIN(J,K)=0.0
	  IRjFIN(0,K)=0.0
!
!
	  DO I=1,NCHEM
	  CFwater(I,J,K)=0.0
	  CFwater(I,0,K)=0.0
	  CFwaterNAC(I,J,K)=0.0
	  RWind(I,J,K)=0.0
	  PASS_RWind(I,J,K)=0.0
!
	  CFsoil(I,J,K)=0.0
	  CFsoil(I,0,K)=0.0
	  TEMP_PERIsoil(I,J,K)=0.0
	  CFsoilNAC(I,J,K)=0.0
	  IGEOsoil(I,J,K)=0.0
	  PIsoil(I,J,K)=0.0
	  EFsoil(I,J,K)=0.0
	  COPIA_PIsoil(I,J,K)=0.0
!
	  CFsed(I,J,K)=0.0
	  CFsed(I,0,K)=0.0
	  TEMP_PERIsed(I,J,K)=0.0
	  CFsedNAC(I,J,K)=0.0
	  IGEOsed(I,J,K)=0.0
	  PIsed(I,J,K)=0.0
	  EFsed(I,J,K)=0.0
	  COPIA_PIsed(I,J,K)=0.0
	  TEMP_mPELq(I,J,K)=0.0
	  TEMP_mERMq(I,J,K)=0.0
	  TEMP_TRI(I,J,K)=0.0
	  ALFABETAwat(I,J,K)=0.0
	  ALFABETAsoil(I,J,K)=0.0
	  ALFABETAsed(I,J,K)=0.0
	  TPIwat(I,J,K)=0.0
	  TPIsoil(I,J,K)=0.0
	  TPIsed(I,J,K)=0.0
	  ALFABETAwatREF(I,J,K)=0.0 
	  TPBGIwat(I,J,K)=0.0
	  ALFABETAsoilREF(I,J,K)=0.0 
	  TPBGIsoil(I,J,K)=0.0
	  ALFABETAsedREF(I,J,K)=0.0 
	  TPBGIsed(I,J,K)=0.0
	  TPICORwat(I,J,K)=0.0
	  TPICORsoil(I,J,K)=0.0
	  TPICORsed(I,J,K)=0.0
	  SUB_RISKwat(I,J,K)=0.0 
	  SUB_RISKsoil(I,J,K)=0.0
	  SUB_RISKsed(I,J,K)=0.0		 
!
	  ENDDO
	  ENDDO
	  ENDDO
!
!
	  NDIVISAO=0
	  NCOMPARTIMENTO1=0
	  NCOMPARTIMENTO2=0
	  NCOMPARTIMENTO3=0
!
!
      DO I=1,NCHEM
	  IF(KEYCONC(I,14).EQV..TRUE.)THEN
      NCOMPARTIMENTO1=1
	  GOTO 88
	  ENDIF
	  ENDDO
!
88    CONTINUE
!
      DO I=1,NCHEM
	  IF(KEYCONC(I,1).EQV..TRUE.)THEN
      NCOMPARTIMENTO2=1
	  GOTO 89
	  ENDIF
	  ENDDO 
!
89    CONTINUE
!
      DO I=1,NCHEM
	  IF(KEYCONC(I,15).EQV..TRUE.)THEN
      NCOMPARTIMENTO3=1
	  GOTO 90
	  ENDIF
	  ENDDO
!
90    CONTINUE
!
!
      NDIVISAO=NCOMPARTIMENTO1+NCOMPARTIMENTO2+NCOMPARTIMENTO3
!
!
      DO J=1,NTIME			! CICLO POR TEMPO
!
      DO K=1,NLOCAL			! CICLO POR LOCAL
!
      NRAIZW=0
	  NRAIZW2=0
	  NRAIZSO1=0
	  NRAIZSO2=0
	  NRAIZSO3=0
	  NRAIZSE1=0
	  NRAIZSE2=0
	  NPELTOT=0
	  NERMTOT=0
!
      DO I=1,NCHEM		  ! CICLO POR ESPÉCIE QUÍMICA
!
!
      iii=0
!
      DO ii=1,NSPECIE
!								                  
	  IF (CHEMICAL(i).EQ.SPECIE(ii)) THEN
	  iii=ii
	  ELSE
	  ENDIF 
	  ENDDO
	  IF (iii.eq.0) THEN
!
      IF((J.EQ.1).AND.(K.EQ.1))THEN
!
      WRITE(*,*)
	  WRITE(*,'("  WARNING!!! CHEMICAL NOT EXIST IN Dataecological DATABASE  --->  ",A30)') CHEMICAL(i)
	  WRITE(*,'("  The index values will all be zero for this Chemical species  ")')
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! CHEMICAL NOT EXIST IN Dataecological DATABASE  --->  ",A30)') CHEMICAL(i)
	  WRITE(99,'(" The index values will all be zero for this Chemical species  ")')
      WRITE(99,*)
!
      ENDIF
!
	  GOTO 14
	  ENDIF
!
!
!****************************************************************************************************************************************
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO PARA ÁGUA	  
!****************************************************************************************************************************************
!
      IF(KEYCONC(i,14).EQV..TRUE.)THEN
!
!
      IF(WATREF(iii).GT.0.0)THEN
      CFwater(I,J,K)=CWATEROTHER(I,J,K)/WATREF(iii)
	  ELSE
	  CFwater(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      RWind(I,J,K)=	CFwater(I,J,K)/(1+CFwater(I,J,K))
!
!-------------------------------------------------------------------
      IF(CFwater(I,J,K).NE.0.0)THEN
      MUL_CFwater(J,K)=MUL_CFwater(J,K)*CFwater(I,J,K)
	  NSTARTW=1
	  ELSE
	  MUL_CFwater(J,K)=MUL_CFwater(J,K)
	  NSTARTW=0
	  ENDIF
!
      NRAIZW=NRAIZW+NSTARTW
!
      PASS_RWind(I,J,K)=1.0-RWind(I,J,K)
!
	  MUL_RWind(J,K)=MUL_RWind(J,K)*PASS_RWind(I,J,K)
!     
!-------------------------------------------------------------------
      IF(WATREFNAC(iii).GT.0.0)THEN
      CFwaterNAC(I,J,K)=CWATEROTHER(I,J,K)/WATREFNAC(iii)
	  NSTARTW2=1
	  ELSE
	  CFwaterNAC(I,J,K)=0.0
	  NSTARTW2=0
	  ENDIF   
!
      NRAIZW2=NRAIZW2+NSTARTW2
!
      S_CFwaterNAC_water(J,K)=S_CFwaterNAC_water(J,K)+CFwaterNAC(I,J,K)
!-------------------------------------------------------------------
!
      IF((BETA(iii).GT.0.0).AND.(ALFA(iii).GT.0.0).AND.(CWATEROTHER(I,J,K).GT.0.0))THEN
      ALFABETAwat(I,J,K)=(LOG10(CWATEROTHER(I,J,K))-ALFA(iii))/BETA(iii)
!
      TPIwat(I,J,K)=1/(1+EXP(-ALFABETAwat(I,J,K)))
!
      ELSE
	  TPIwat(I,J,K)=0.0
	  ENDIF
!
!-------------------------------------------------------------------
!
      IF((BETA(iii).GT.0.0).AND.(ALFA(iii).GT.0.0).AND.(WATREF(iii).GT.0.0))THEN
      ALFABETAwatREF(I,J,K)=(LOG10(WATREF(iii))-ALFA(iii))/BETA(iii)
!
      TPBGIwat(I,J,K)=1/(1+EXP(-ALFABETAwatREF(I,J,K)))
!
      ELSE
	  TPBGIwat(I,J,K)=0.0
	  ENDIF
!
!-------------------------------------------------------------------
!
      TPICORwat(I,J,K)=(TPIwat(I,J,K)-TPBGIwat(I,J,K))/(1-TPBGIwat(I,J,K))
!
!-------------------------------------------------------------------
!
      SUB_RISKwat(I,J,K)=1.0-TPICORwat(I,J,K)
!
	  MUL_RISKwat(J,K)=MUL_RISKwat(J,K)*SUB_RISKwat(I,J,K)
!
!-------------------------------------------------------------------
!
      ELSEIF((KEYCONC(i,14).EQV..FALSE.).AND.(J.EQ.1).AND.(K.EQ.1))THEN
!
!
      WRITE(*,*)
	  WRITE(*,'("  WARNING!!! There are no values of contaminant concentrations for OTHER_WATER matrix")')
	  WRITE(*,'("  THE VALUES OF THE WATER INDICES WILL NOT BE CALCULATED FOR ----->  ",A30)')  CHEMICAL(i)
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! There are no values of contaminant concentrations for OTHER_WATER matrix")')    
	  WRITE(99,'(" THE VALUES OF THE WATER INDICES WILL NOT BE CALCULATED FOR ----->  ",A30)')  CHEMICAL(i)
      WRITE(99,*)
!
      ENDIF
!
!
!****************************************************************************************************************************************
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO PARA SOLO	  
!****************************************************************************************************************************************
!
      IF(KEYCONC(i,1).EQV..TRUE.)THEN
!
!
      IF(SOILREF(iii).GT.0.0)THEN
      CFsoil(I,J,K)=CSOIL(I,J,K)/SOILREF(iii)
	  ELSE
	  CFsoil(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      IF(CFsoil(I,J,K).NE.0.0)THEN
      MUL_CFsoil(J,K)=MUL_CFsoil(J,K)*CFsoil(I,J,K)
	  NSTARTSO1=1
	  ELSE
	  MUL_CFsoil(J,K)=MUL_CFsoil(J,K)
	  NSTARTSO1=0
	  ENDIF
!
      NRAIZSO1=NRAIZSO1+NSTARTSO1
!-------------------------------------------------------------------
!
      SOMA_mCd_soil(J,K)=SOMA_mCd_soil(J,K)+CFsoil(I,J,K)
!
!-------------------------------------------------------------------
      TEMP_PERIsoil(I,J,K)=CFsoil(I,J,K)*TIR(iii)
!
      PERIsoil(J,K)=PERIsoil(J,K)+TEMP_PERIsoil(I,J,K)
!-------------------------------------------------------------------
      IF(SOILREFNAC(iii).GT.0.0)THEN
      CFsoilNAC(I,J,K)=CSOIL(I,J,K)/SOILREFNAC(iii)
	  NSTARTSO2=1
	  ELSE
	  CFsoilNAC(I,J,K)=0.0
	  NSTARTSO2=0
	  ENDIF   
!
      NRAIZSO2=NRAIZSO2+NSTARTSO2
!
      S_CFsoilNAC_soil(J,K)=S_CFsoilNAC_soil(J,K)+CFsoilNAC(I,J,K)
!-------------------------------------------------------------------
      IF((CROSTA(iii).GT.0.0).AND.(CSOIL(I,J,K).GT.0.0))THEN
      IGEOsoil(I,J,K)=LOG10(CSOIL(I,J,K)/(1.5*CROSTA(iii)))/LOG10(2.0)
	  ELSE
	  IGEOsoil(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      IF(CROSTA(iii).GT.0.0)THEN
      PIsoil(I,J,K)=CSOIL(I,J,K)/CROSTA(iii)
	  NSTARTSO3=1
	  ELSE
	  PIsoil(I,J,K)=0.0
	  NSTARTSO3=0
	  ENDIF 
!
      NRAIZSO3=NRAIZSO3+NSTARTSO3
!
!-------------------------------------------------------------------
!
      IF((EFCSOIL(J,K).GT.0.0).AND.(EFCROSTAsolo.GT.0.0).AND.(CROSTA(iii).GT.0.0))THEN
      EFsoil(I,J,K)=(CSOIL(I,J,K)/EFCSOIL(J,K))/(CROSTA(iii)/EFCROSTAsolo)
	  ELSE
	  EFsoil(I,J,K)=0.0
	  ENDIF
!             [
!-------------------------------------------------------------------
!
      P1_PINEWsoil(J,K)=P1_PINEWsoil(J,K)+PIsoil(I,J,K)
!
      COPIA_PIsoil(I,J,K)=PIsoil(I,J,K)
!-------------------------------------------------------------------
! 
      IF((BETA(iii).GT.0.0).AND.(ALFA(iii).GT.0.0).AND.(CSOIL(I,J,K).GT.0.0))THEN
!
      ALFABETAsoil(I,J,K)=(LOG10(CSOIL(I,J,K))-ALFA(iii))/BETA(iii)
!
      TPIsoil(I,J,K)=1/(1+EXP(-ALFABETAsoil(I,J,K)))
!
      ELSE 
	  TPIsoil(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
!
      IF((BETA(iii).GT.0.0).AND.(ALFA(iii).GT.0.0).AND.(SOILREF(iii).GT.0.0))THEN
      ALFABETAsoilREF(I,J,K)=(LOG10(SOILREF(iii))-ALFA(iii))/BETA(iii)
!
      TPBGIsoil(I,J,K)=1/(1+EXP(-ALFABETAsoilREF(I,J,K)))
!
      ELSE
	  TPBGIsoil(I,J,K)=0.0
	  ENDIF
!
!-------------------------------------------------------------------
!
      TPICORsoil(I,J,K)=(TPIsoil(I,J,K)-TPBGIsoil(I,J,K))/(1-TPBGIsoil(I,J,K))
!
!-------------------------------------------------------------------
!
      SUB_RISKsoil(I,J,K)=1.0-TPICORsoil(I,J,K)
!
	  MUL_RISKsoil(J,K)=MUL_RISKsoil(J,K)*SUB_RISKsoil(I,J,K)
!
!-------------------------------------------------------------------

!
!
      ELSEIF((KEYCONC(i,1).EQV..FALSE.).AND.(J.EQ.1).AND.(K.EQ.1))THEN
!
!
      WRITE(*,*)
	  WRITE(*,'("  WARNING!!! There are no values of contaminant concentrations for SOIL matrix  ")')   
	  WRITE(*,'("  THE VALUES OF THE SOIL INDICES WILL NOT BE CALCULATED FOR ----->  ",A30)')  CHEMICAL(i)
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! There are no values of contaminant concentrations for SOIL matrix  ")')  
	  WRITE(99,'(" THE VALUES OF THE SOIL INDICES WILL NOT BE CALCULATED FOR ----->  ",A30)')  CHEMICAL(i)
      WRITE(99,*)
!
      ENDIF
!
!
!****************************************************************************************************************************************
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO PARA SEDIMENTOS	  
!****************************************************************************************************************************************
!
      IF(KEYCONC(i,15).EQV..TRUE.)THEN
!
!
      IF(SEDREF(iii).GT.0.0)THEN
      CFsed(I,J,K)=CSEDIMENT(I,J,K)/SEDREF(iii)
	  ELSE
	  CFsed(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      IF(CFsed(I,J,K).NE.0.0)THEN
      MUL_CFsed(J,K)=MUL_CFsed(J,K)*CFsed(I,J,K)
	  NSTARTSE1=1
	  ELSE
	  MUL_CFsed(J,K)=MUL_CFsed(J,K)
	  NSTARTSE1=0
	  ENDIF
!
      NRAIZSE1=NRAIZSE1+NSTARTSE1
!-------------------------------------------------------------------
!
      SOMA_mCd_sed(J,K)=SOMA_mCd_sed(J,K)+CFsed(I,J,K)
!
!-------------------------------------------------------------------
      TEMP_PERIsed(I,J,K)=CFsed(I,J,K)*TIR(iii)
!
      PERIsed(J,K)=PERIsed(J,K)+TEMP_PERIsed(I,J,K)
!-------------------------------------------------------------------
      IF(SEDREFNAC(iii).GT.0.0)THEN
      CFsedNAC(I,J,K)=CSEDIMENT(I,J,K)/SEDREFNAC(iii)
	  NSTARTSE2=1
	  ELSE
	  CFsedNAC(I,J,K)=0.0
	  NSTARTSE2=0
	  ENDIF   
!
      NRAIZSE2=NRAIZSE2+NSTARTSE2
!
      S_CFsedNAC_sed(J,K)=S_CFsedNAC_sed(J,K)+CFsedNAC(I,J,K)
!-------------------------------------------------------------------
      IF((SEDREF(iii).GT.0.0).AND.(CSEDIMENT(I,J,K).GT.0.0))THEN
      IGEOsed(I,J,K)=LOG10(CSEDIMENT(I,J,K)/(1.5*SEDREF(iii)))/LOG10(2.0)
	  ELSE
	  IGEOsed(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      PIsed(I,J,K)=CFsed(I,J,K)
!
!-------------------------------------------------------------------
!
      IF((EFCSEDIMENT(J,K).GT.0.0).AND.(EFCROSTAsed.GT.0.0).AND.(CROSTA(iii).GT.0.0))THEN
      EFsed(I,J,K)=(CSEDIMENT(I,J,K)/EFCSEDIMENT(J,K))/(CROSTA(iii)/EFCROSTAsed)
	  ELSE
	  EFsed(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
!
      P1_PINEWsed(J,K)=P1_PINEWsed(J,K)+PIsed(I,J,K)
!
      COPIA_PIsed(I,J,K)=PIsed(I,J,K)
!-------------------------------------------------------------------
      IF(PEL(iii).GT.0.0)THEN
      TEMP_mPELq(I,J,K)=CSEDIMENT(I,J,K)/PEL(iii)
	  NPEL=1
	  ELSE
	  TEMP_mPELq(I,J,K)=0.0
	  NPEL=0
	  ENDIF
!
      NPELTOT=NPELTOT+NPEL
!
      SOMA_mPELq(J,K)=SOMA_mPELq(J,K)+TEMP_mPELq(I,J,K)
!-------------------------------------------------------------------
      IF(ERM(iii).GT.0.0)THEN
      TEMP_mERMq(I,J,K)=CSEDIMENT(I,J,K)/ERM(iii)
	  NERM=1
	  ELSE
	  TEMP_mERMq(I,J,K)=0.0
	  NERM=0
	  ENDIF
!
      NERMTOT=NERMTOT+NERM
!
      SOMA_mERMq(J,K)=SOMA_mERMq(J,K)+TEMP_mERMq(I,J,K)
!-------------------------------------------------------------------
!
      IF((TEL(iii).GT.0.0).AND.(PEL(iii).GT.0.0).AND.(CSEDIMENT(I,J,K).GT.0.0))THEN
      TEMP_TRI(I,J,K)=SQRT(((CSEDIMENT(I,J,K)/TEL(iii))**2+(CSEDIMENT(I,J,K)/PEL(iii))**2)/2)
	  ELSE
	  TEMP_TRI(I,J,K)=0.0
	  ENDIF
!
      TRI(J,K)=TRI(J,K)+TEMP_TRI(I,J,K)
!
!-------------------------------------------------------------------
!
      IF((BETA(iii).GT.0.0).AND.(ALFA(iii).GT.0.0).AND.(CSEDIMENT(I,J,K).GT.0.0))THEN
!
      ALFABETAsed(I,J,K)=(LOG10(CSEDIMENT(I,J,K))-ALFA(iii))/BETA(iii)
!
	  TPIsed(I,J,K)=1/(1+EXP(-ALFABETAsed(I,J,K)))
!
      ELSE 
	  TPIsed(I,J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
!
      IF((BETA(iii).GT.0.0).AND.(ALFA(iii).GT.0.0).AND.(SEDREF(iii).GT.0.0))THEN
      ALFABETAsedREF(I,J,K)=(LOG10(SEDREF(iii))-ALFA(iii))/BETA(iii)
!
      TPBGIsed(I,J,K)=1/(1+EXP(-ALFABETAsedREF(I,J,K)))
!
      ELSE
	  TPBGIsed(I,J,K)=0.0
	  ENDIF
!
!-------------------------------------------------------------------
!
      TPICORsed(I,J,K)=(TPIsed(I,J,K)-TPBGIsed(I,J,K))/(1-TPBGIsed(I,J,K))
!
!-------------------------------------------------------------------
!
      SUB_RISKsed(I,J,K)=1.0-TPICORsed(I,J,K)
!
	  MUL_RISKsed(J,K)=MUL_RISKsed(J,K)*SUB_RISKsed(I,J,K)
!
!-------------------------------------------------------------------
!	  
!
      ELSEIF((KEYCONC(i,15).EQV..FALSE.).AND.(J.EQ.1).AND.(K.EQ.1))THEN
!
!
      WRITE(*,*)
	  WRITE(*,'("  WARNING!!! There are no values of contaminant concentrations for SEDIMENT matrix  ")')    
	  WRITE(*,'("  THE VALUES OF THE SEDIMENT INDICES WILL NOT BE CALCULATED FOR ----->  ",A30)')  CHEMICAL(i)
      WRITE(*,*)
!
      WRITE(99,*)
	  WRITE(99,'(" WARNING!!! There are no values of contaminant concentrations for SEDIMENT matrix  ")')   
	  WRITE(99,'(" THE VALUES OF THE SEDIMENT INDICES WILL NOT BE CALCULATED FOR ----->  ",A30)')  CHEMICAL(i)
      WRITE(99,*)
!
      ENDIF
!
!
14    CONTINUE
!
      ENDDO    ! ACABA O CICLO POR ESPÉCIE QUÍMICA (i)
!
!****************************************************************************************************************************************
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO PARA ÁGUA	  
!****************************************************************************************************************************************
!
!
      IF(NCOMPARTIMENTO1.EQ.1)THEN		  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      IF((MUL_CFwater(J,K).GT.0.0).AND.(NRAIZW.GT.0))THEN
      MPIwat(J,K)=MUL_CFwater(J,K)**(1.0/NRAIZW)
	  ELSE
	  MPIwat(J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      IF(NRAIZW2.GT.0)THEN
      IPIT_water(J,K)=S_CFwaterNAC_water(J,K)/NRAIZW2
      ENDIF
!-------------------------------------------------------------------
!
!
      RWcomb(J,K)=1.0-MUL_RWind(J,K)
!
!-------------------------------------------------------------------
!
      RISKwat(J,K)=1.0-MUL_RISKwat(J,K)
!
!-------------------------------------------------------------------
	  ENDIF						 	          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!****************************************************************************************************************************************
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO PARA SOLO	  
!****************************************************************************************************************************************
!		   			 
!
      IF(NCOMPARTIMENTO2.EQ.1)THEN		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      IF((MUL_CFsoil(J,K).GT.0.0).AND.(NRAIZSO1.GT.0))THEN
      PLIsoil(J,K)=MUL_CFsoil(J,K)**(1.0/NRAIZSO1)
	  ELSE
	  PLIsoil(J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      IF(NRAIZSO1.GT.0)THEN
      mCd_soil(J,K)=SOMA_mCd_soil(J,K)/NRAIZSO1
	  ENDIF
!-------------------------------------------------------------------
      IF(NRAIZSO2.GT.0)THEN
      IPIT_soil(J,K)=S_CFsoilNAC_soil(J,K)/NRAIZSO2
      ENDIF
!-------------------------------------------------------------------
      IF(NRAIZSO3.GT.0)THEN
      P2_PINEWsoil(J,K)=P1_PINEWsoil(J,K)/NRAIZSO3
	  ENDIF
!
!
      MAX_PIsoil(J,K)=MAXVAL(COPIA_PIsoil)
!
      DO IKK=1,NCHEM
	  COPIA_PIsoil(IKK,J,K)=0.0
	  ENDDO
!
      IF((P2_PINEWsoil(J,K).GT.0.0).OR.(MAX_PIsoil(J,K).GT.0.0))THEN
      PINEWsoil(J,K)=SQRT((P2_PINEWsoil(J,K)**2+MAX_PIsoil(J,K)**2)/2.0)
	  ENDIF
!
!-------------------------------------------------------------------
!
      RISKsoil(J,K)=1.0-MUL_RISKsoil(J,K)
!
!-------------------------------------------------------------------  
!
      ENDIF					    		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!****************************************************************************************************************************************
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO PARA SEDIMENTO	  
!****************************************************************************************************************************************
!
      IF(NCOMPARTIMENTO3.EQ.1)THEN		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		   			 
!
!
      IF((MUL_CFsed(J,K).GT.0.0).AND.(NRAIZSE1.GT.0))THEN
      PLIsed(J,K)=MUL_CFsed(J,K)**(1.0/NRAIZSE1)
	  ELSE
	  PLIsed(J,K)=0.0
	  ENDIF
!-------------------------------------------------------------------
      IF(NRAIZSE1.GT.0)THEN
      mCd_sed(J,K)=SOMA_mCd_sed(J,K)/NRAIZSE1
	  ENDIF
!-------------------------------------------------------------------
      IF(NRAIZSE2.GT.0)THEN
      IPIT_sed(J,K)=S_CFsedNAC_sed(J,K)/NRAIZSE2
	  ENDIF
!-------------------------------------------------------------------
      IF(NRAIZSE2.GT.0)THEN
      P2_PINEWsed(J,K)=P1_PINEWsed(J,K)/NRAIZSE2
	  ENDIF
!
      MAX_PIsed(J,K)=MAXVAL(COPIA_PIsed)
!
      DO IKK=1,NCHEM
	  COPIA_PIsed(IKK,J,K)=0.0
	  ENDDO
!
      IF((P2_PINEWsed(J,K).GT.0.0).OR.(MAX_PIsed(J,K).GT.0.0))THEN
      PINEWsed(J,K)=SQRT((P2_PINEWsed(J,K)**2+MAX_PIsed(J,K)**2)/2.0)
	  ENDIF
!
!-------------------------------------------------------------------
      IF(NPELTOT.GT.0)THEN
      mPELq(J,K)=SOMA_mPELq(J,K)/NPELTOT
	  ENDIF
!-------------------------------------------------------------------
      IF(NERMTOT.GT.0)THEN
      mERMq(J,K)=SOMA_mERMq(J,K)/NERMTOT
	  ENDIF
!-------------------------------------------------------------------
!
      RISKsed(J,K)=1.0-MUL_RISKsed(J,K)
!
!-------------------------------------------------------------------
!
      ENDIF					    		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!
!****************************************************************************************************************************************
!****************************************************************************************************************************************
!
!   CALCULOS DA AVALIAÇÃO DE RISCO ECOLÓGICO JUNÇÃO DE COMPARTIMENTOS	  
!
!****************************************************************************************************************************************
!****************************************************************************************************************************************
!
!
      IF((MPIwat(J,K).GT.0.0).AND.(PLIsed(J,K).GT.0.0))THEN
      Kd_MPI(J,K)=LOG10(PLIsed(J,K)/MPIwat(J,K))
	  ENDIF
!
!-------------------------------------------------------------------
!
!	   
!
!
      IF((NDIVISAO.GT.0).AND.(RISKwat(J,K).NE.1.0).AND.(RISKsoil(J,K).NE.1.0).AND.(RISKsed(J,K).NE.1.0))THEN
      IRjFIN(J,K)=1.0-(10.0**((LOG10(1.0-RISKwat(J,K))+LOG10(1.0-RISKsoil(J,K))+LOG10(1.0-RISKsed(J,K)))/NDIVISAO))
	  ENDIF
!
!
!
!
!****************************************************************************************************************************************
!
      ENDDO		! ACABA O CICLO POR LOCAL (K)
!
      ENDDO		! ACABA O CICLO POR TEMPO (J)
!
!
      AGUA_REF(1)='NATIONAL'
	  AGUA_REF(2)='U.S.EPA'
	  AGUA_REF(3)='WHO'
!											  
	  SOLO_REF(1)='REGIONAL'
	  SOLO_REF(2)='NATIONAL'
	  SOLO_REF(3)='U.S.EPA'
	  SOLO_REF(4)='IN NATURA'
!
	  SED_REF(1)='REGIONAL'
	  SED_REF(2)='NATIONAL'
	  SED_REF(3)='U.S.EPA'
!  
      aspas = '"'
!
      Mun_comp=0
!
      IF(NCOMPARTIMENTO1.EQ.1)THEN	
	  Mun_comp=1
	  ENDIF	
      IF(NCOMPARTIMENTO2.EQ.1)THEN	
	  Mun_comp=2
	  ENDIF
      IF(NCOMPARTIMENTO3.EQ.1)THEN	
	  Mun_comp=3
	  ENDIF  			 
!
      CALL NOMINATION(CHEMICAL,NCHEM,NOMES)
!
!-----------------------------------------------------------
!	  IMPRESSÃO DOS DADOS DE SAÍDA
!-----------------------------------------------------------
!
	  open (UNIT=68, file='Results\Individual.json')	
!
	  WRITE(68,'("{")')
!
	  WRITE(68,'(A1,"Selected scenario and reference values",A1,": [")')aspas,aspas
!
      WRITE(68,'("{")')
!
      write(68,'(A1,"Scenary",A1,":",1x,A1,A12,A1,",") )') aspas,aspas,aspas,SCENARIES(SCENAR),aspas
      write(68,'(A1,"Water reference value",A1,":",1x,A1,A8,A1,",") )') aspas,aspas,aspas,AGUA_REF(KEY_WATER),aspas
      write(68,'(A1,"Soil reference value",A1,":",1x,A1,A9,A1,",") )') aspas,aspas,aspas,SOLO_REF(KEY_SOIL),aspas
      write(68,'(A1,"Sediment reference value",A1,":",1x,A1,A9,A1) )') aspas,aspas,aspas,SED_REF(KEY_SEDIMENT),aspas
!
	  WRITE(68,'("}")')
!
	  WRITE(68,'("],")')
!
!
	  DO i=1,NCHEM
!
	  WRITE(68,'(A1,A10,A1,": [")')aspas,NOMES(i),aspas
!
      ijj=0
!
      DO ij=1,NSPECIE
!								                  
	  IF (CHEMICAL(i).EQ.SPECIE(ij)) THEN
	  ijj=ij
	  ELSE
	  ENDIF 
      ENDDO
!
!
      IF(NTIME.EQ.1)THEN
      NINICIO=1
	  ELSE
	  NINICIO=2
	  ENDIF
!
      DO k=1,NLOCAL
	  DO j=NINICIO,NTIME
!
      IF((CFsed(I,J,K).NE.CFsed(I,J-1,K)).OR.(CFsoil(I,J,K).NE.CFsoil(I,J-1,K)).OR.(CFwater(I,J,K).NE.CFwater(I,J-1,K)).OR.(EFCsoil(J,K).NE.EFCsoil(J-1,K)).OR.(EFCsediment(J,K).NE.EFCsediment(J-1,K)))THEN
      IALL(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
!
!	  IF(j.EQ.1)THEN
!
      IF(IALL(k).EQ.1)THEN
!
      WRITE(68,'("{")')
      write(68,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(68,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(68,'(A1,"CFwater",A1,":",1x,ES12.5,",") )') aspas,aspas,CFwater(I,J,K)
      write(68,'(A1,"CFsoil",A1,":",1x,ES12.5,",") )') aspas,aspas,CFsoil(I,J,K)
	  write(68,'(A1,"Igeo soil",A1,":",1x,ES12.5,",") )') aspas,aspas,IGEOsoil(I,J,K)
      write(68,'(A1,"EF soil - ",A2,A1,":",1x,ES12.5,",") )') aspas,EFREFsoil,aspas,EFsoil(I,J,K)
      write(68,'(A1,"PI soil",A1,":",1x,ES12.5,",") )') aspas,aspas,PIsoil(I,J,K)
      write(68,'(A1,"Igeo sediment",A1,":",1x,ES12.5,",") )') aspas,aspas,IGEOsed(I,J,K)
      write(68,'(A1,"EF sediment -",A2,A1,":",1x,ES12.5,",") )') aspas,EFREFsediment,aspas,EFsed(I,J,K)
      write(68,'(A1,"PI sediment",A1,":",1x,ES12.5) )') aspas,aspas,PIsed(I,J,K)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(68,'("}")')
	  ELSE
	  WRITE(68,'("},")')
	  ENDIF
!
	  ELSEIF((IALL(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(68,'("{")')
      write(68,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(68,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(68,'(A1,"CFwater",A1,":",1x,ES12.5,",") )') aspas,aspas,CFwater(I,J,K)
      write(68,'(A1,"CFsoil",A1,":",1x,ES12.5,",") )') aspas,aspas,CFsoil(I,J,K)
	  write(68,'(A1,"Igeo soil",A1,":",1x,ES12.5,",") )') aspas,aspas,IGEOsoil(I,J,K)
      write(68,'(A1,"EF soil -",A2,A1,":",1x,ES12.5,",") )') aspas,EFREFsoil,aspas,EFsoil(I,J,K)
      write(68,'(A1,"PI soil",A1,":",1x,ES12.5,",") )') aspas,aspas,PIsoil(I,J,K)
      write(68,'(A1,"Igeo sediment",A1,":",1x,ES12.5,",") )') aspas,aspas,IGEOsed(I,J,K)
      write(68,'(A1,"EF sediment -",A2,A1,":",1x,ES12.5,",") )') aspas,EFREFsediment,aspas,EFsed(I,J,K)
      write(68,'(A1,"PI sediment",A1,":",1x,ES12.5) )') aspas,aspas,PIsed(I,J,K)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(68,'("}")')
	  ELSE
	  WRITE(68,'("},")')
	  ENDIF
!
	  ENDIF
!
!
      ENDDO		! fim DO j=1,NTIME
!
	  ENDDO		! fim DO k=1,NLOCAL
!
      IF(i.EQ.NCHEM)THEN
	  WRITE(68,'("]")')
	  ELSE
	  WRITE(68,'("],")')
	  ENDIF
!
	  ENDDO	  ! FIM DO CICLO "i"
!	       
	  WRITE(68,'("}")')  
!
      CLOSE(68)
!
!
!-----------------------------------------------------------
!	  IMPRESSÃO DOS DADOS DE SAÍDA
!-----------------------------------------------------------
!
!
	  open (UNIT=70, file='Results\Combined.json')	
!
      WRITE(70,'("{")')	 
!
!-----------------------------------------------------------
!	  COMPARTIMENTO AGUA
!-----------------------------------------------------------
!
      IF(NCOMPARTIMENTO1.EQ.1)THEN		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		   			 
!
	  WRITE(70,'(A1,"Water compartment",A1,": [")')aspas,aspas
!
!
      IF(NTIME.EQ.1)THEN
      NINICIO1=1
	  ELSE
	  NINICIO1=2
	  ENDIF
!
      DO k=1,NLOCAL
	  DO j=NINICIO1,NTIME
!
      IF((MPIwat(J,K).NE.MPIwat(J-1,K)).OR.(RWcomb(J,K).NE.RWcomb(J-1,K)).OR.(IPIT_water(J,K).NE.IPIT_water(J-1,K)))THEN
      JALL1(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
!
!	  IF(j.EQ.1)THEN
!
      IF(JALL1(k).EQ.1)THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(70,'(A1,"MPI",A1,":",1x,ES12.5,",") )') aspas,aspas,MPIwat(J,K)
      write(70,'(A1,"Rw-comb",A1,":",1x,ES12.5,",") )') aspas,aspas,RWcomb(J,K)
	  write(70,'(A1,"IPIth",A1,":",1x,ES12.5) )') aspas,aspas,IPIT_water(J,K)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
      ELSEIF((JALL1(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(70,'("{")')
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(70,'(A1,"MPI",A1,":",1x,ES12.5,",") )') aspas,aspas,MPIwat(J,K)
      write(70,'(A1,"Rw-comb",A1,":",1x,ES12.5,",") )') aspas,aspas,RWcomb(J,K)
	  write(70,'(A1,"IPIth",A1,":",1x,ES12.5) )') aspas,aspas,IPIT_water(J,K)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
	  ENDIF
! 
      ENDDO		! fim DO j=1,NTIME
!
	  ENDDO		! fim DO k=1,NLOCAL
!
      IF((Mun_comp.EQ.1).and.(NDIVISAO.EQ.0))THEN
	  WRITE(70,'("]")')
	  ELSE
	  WRITE(70,'("],")')
	  ENDIF
!
      ENDIF						  		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------
!	  COMPARTIMENTO SOLO
!-----------------------------------------------------------
!
      IF(NCOMPARTIMENTO2.EQ.1)THEN		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		   			 
!
	  WRITE(70,'(A1,"Soil compartment",A1,": [")')aspas,aspas
!
!
      IF(NTIME.EQ.1)THEN
      NINICIO2=1
	  ELSE
	  NINICIO2=2
	  ENDIF
!
      DO k=1,NLOCAL
	  DO j=NINICIO2,NTIME
!
      IF(PLIsoil(J,K).NE.PLIsoil(J-1,K))THEN
      JALL2(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
!
!	  IF(j.EQ.1)THEN
!
      IF(JALL2(k).EQ.1)THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(70,'(A1,"PLI",A1,":",1x,ES12.5,",") )') aspas,aspas,PLIsoil(J,K)
      write(70,'(A1,"mCd",A1,":",1x,ES12.5,",") )') aspas,aspas,mCd_soil(J,K)
      write(70,'(A1,"PERI",A1,":",1x,ES12.5,",") )') aspas,aspas,PERIsoil(J,K)
      write(70,'(A1,"IPIT",A1,":",1x,ES12.5,",") )') aspas,aspas,IPIT_soil(J,K)
	  write(70,'(A1,"PInem",A1,":",1x,ES12.5) )') aspas,aspas,PINEWsoil(J,K)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
	  ELSEIF((JALL2(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(70,'(A1,"PLI",A1,":",1x,ES12.5,",") )') aspas,aspas,PLIsoil(J,K)
      write(70,'(A1,"mCd",A1,":",1x,ES12.5,",") )') aspas,aspas,mCd_soil(J,K)
      write(70,'(A1,"PERI",A1,":",1x,ES12.5,",") )') aspas,aspas,PERIsoil(J,K)
      write(70,'(A1,"IPIT",A1,":",1x,ES12.5,",") )') aspas,aspas,IPIT_soil(J,K)
	  write(70,'(A1,"PInem",A1,":",1x,ES12.5) )') aspas,aspas,PINEWsoil(J,K)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
	  ENDIF
!
      ENDDO		! fim DO j=1,NTIME
!
	  ENDDO		! fim DO k=1,NLOCAL
!
      IF((Mun_comp.EQ.2).and.(NDIVISAO.EQ.0))THEN
	  WRITE(70,'("]")')
	  ELSE
	  WRITE(70,'("],")')
	  ENDIF
!
      ENDIF						  		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------
!	  COMPARTIMENTO SEDIMENTO
!-----------------------------------------------------------
!
      IF(NCOMPARTIMENTO3.EQ.1)THEN		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		   			 
!
	  WRITE(70,'(A1,"Sediment compartment",A1,": [")')aspas,aspas
!
      IF(NTIME.EQ.1)THEN
      NINICIO3=1
	  ELSE
	  NINICIO3=2
	  ENDIF
!
      DO k=1,NLOCAL
	  DO j=NINICIO3,NTIME
!
      IF(PLIsed(J,K).NE.PLIsed(J-1,K))THEN
      JALL3(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
!
!	  IF(j.EQ.1)THEN
!
      IF(JALL3(k).EQ.1)THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(70,'(A1,"PLI",A1,":",1x,ES12.5,",") )') aspas,aspas,PLIsed(J,K)
      write(70,'(A1,"mCd",A1,":",1x,ES12.5,",") )') aspas,aspas,mCd_sed(J,K)
      write(70,'(A1,"PERI",A1,":",1x,ES12.5,",") )') aspas,aspas,PERIsed(J,K)
      write(70,'(A1,"IPIT",A1,":",1x,ES12.5,",") )') aspas,aspas,IPIT_sed(J,K)
	  write(70,'(A1,"PInem",A1,":",1x,ES12.5,",") )') aspas,aspas,PINEWsed(J,K)
	  write(70,'(A1,"m-PEL-q",A1,":",1x,ES12.5,",") )') aspas,aspas,mPELq(J,K)
      write(70,'(A1,"m-ERM-q",A1,":",1x,ES12.5,",") )') aspas,aspas,mERMq(J,K)
      write(70,'(A1,"TRI",A1,":",1x,ES12.5) )') aspas,aspas,TRI(J,K)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
      ELSEIF((JALL3(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(70,'(A1,"PLI",A1,":",1x,ES12.5,",") )') aspas,aspas,PLIsed(J,K)
      write(70,'(A1,"mCd",A1,":",1x,ES12.5,",") )') aspas,aspas,mCd_sed(J,K)
      write(70,'(A1,"PERI",A1,":",1x,ES12.5,",") )') aspas,aspas,PERIsed(J,K)
      write(70,'(A1,"IPIT",A1,":",1x,ES12.5,",") )') aspas,aspas,IPIT_sed(J,K)
	  write(70,'(A1,"PInem",A1,":",1x,ES12.5,",") )') aspas,aspas,PINEWsed(J,K)
	  write(70,'(A1,"m-PEL-q",A1,":",1x,ES12.5,",") )') aspas,aspas,mPELq(J,K)
      write(70,'(A1,"m-ERM-q",A1,":",1x,ES12.5,",") )') aspas,aspas,mERMq(J,K)
      write(70,'(A1,"TRI",A1,":",1x,ES12.5) )') aspas,aspas,TRI(J,K)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
	  ENDIF
! 
!
      ENDDO		! fim DO j=1,NTIME
!
	  ENDDO		! fim DO k=1,NLOCAL
!
      IF(NDIVISAO.GT.0)THEN
	  WRITE(70,'("],")')
	  ELSE
	  WRITE(70,'("]")')
	  ENDIF
!
      ENDIF						  		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!-----------------------------------------------------------
!	  COMPARTIMENTO MISTO
!-----------------------------------------------------------
!
      IF(NDIVISAO.GT.0)THEN		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		   			 
!
	  WRITE(70,'(A1,"Chemical Line of Evidence",A1,": [")')aspas,aspas
!
!
      IF(NTIME.EQ.1)THEN
      NINICIO4=1
	  ELSE
	  NINICIO4=2
	  ENDIF
!
      DO k=1,NLOCAL
	  DO j=NINICIO4,NTIME
!
      IF((Kd_MPI(J,K).NE.Kd_MPI(J-1,K)).OR.(RISKwat(J,K).NE.RISKwat(J-1,K)).OR.(RISKsoil(J,K).NE.RISKsoil(J-1,K)).OR.(RISKsed(J,K).NE.RISKsed(J-1,K)).OR.(IRjFIN(J,K).NE.IRjFIN(J-1,K)))THEN
      JALL4(k)=1
	  ENDIF
!
      ENDDO
	  ENDDO
!
!
      DO k=1,NLOCAL
	  DO j=1,NTIME
!
      IF(JALL4(k).EQ.1)THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,I3,",") )') aspas,aspas,j
      write(70,'(A1,"Kd-MPI",A1,":",1x,ES12.5,",") )') aspas,aspas,Kd_MPI(J,K)
      write(70,'(A1,"Risk ChemLoE water",A1,":",1x,ES12.5,",") )') aspas,aspas,RISKwat(J,K)
      write(70,'(A1,"Risk ChemLoE soil",A1,":",1x,ES12.5,",") )') aspas,aspas,RISKsoil(J,K)
      write(70,'(A1,"Risk ChemLoE sediment",A1,":",1x,ES12.5,",") )') aspas,aspas,RISKsed(J,K)
	  write(70,'(A1,"IR",A1,":",1x,ES12.5) )') aspas,aspas,IRjFIN(J,K)
      IF((J.EQ.NTIME).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
	  ELSEIF((JALL4(k).NE.1).AND.(j.EQ.1))THEN
!
      WRITE(70,'("{")')	 
      write(70,'(A1,"Local",A1,":",1x,I3,",") )') aspas,aspas,K
      write(70,'(A1,"Time",A1,":",1x,A1,"All",A1,",") )') aspas,aspas,aspas,aspas
      write(70,'(A1,"Kd-MPI",A1,":",1x,ES12.5,",") )') aspas,aspas,Kd_MPI(J,K)
      write(70,'(A1,"Risk ChemLoE water",A1,":",1x,ES12.5,",") )') aspas,aspas,RISKwat(J,K)
      write(70,'(A1,"Risk ChemLoE soil",A1,":",1x,ES12.5,",") )') aspas,aspas,RISKsoil(J,K)
      write(70,'(A1,"Risk ChemLoE sediment",A1,":",1x,ES12.5,",") )') aspas,aspas,RISKsed(J,K)
	  write(70,'(A1,"IR",A1,":",1x,ES12.5) )') aspas,aspas,IRjFIN(J,K)
      IF((J.EQ.1).and.(K.EQ.NLOCAL))THEN
	  WRITE(70,'("}")')
	  ELSE
	  WRITE(70,'("},")')
	  ENDIF
!
	  ENDIF
!
      ENDDO		! fim DO j=1,NTIME
!
!
	  ENDDO		! fim DO k=1,NLOCAL
!
	  WRITE(70,'("]")')
!
      ENDIF						  		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	  WRITE(70,'("}")')
!
!
      CLOSE(70)
!
!
 	  RETURN
	  END
!
!
!*********************************************************************************************************************************************************
!*********************************************************************************************************************************************************
!

!
      SUBROUTINE ECODATA (SCENAR,NSPECIE,SPECIE,WATREF,CROSTA,SOILREF,SEDREF,TEL,PEL,ERM,TIR,ALFA,BETA,SEDREFNAC,SOILREFNAC,WATREFNAC,KEY_WATER,KEY_SOIL,KEY_SEDIMENT)  
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      INTEGER SCENAR
!
	   CHARACTER(LEN=50)  :: SPECIE(500)
!
!
      DIMENSION WATREF_TEMP(500,3),SOILREF_TEMP(500,3,4),SEDREF_TEMP(500,3) ! O 3 REPRESENTA O TIPO DE REFERENCIA (1 = NACIONAL, 2 = USEPA, 3 = WHO)
      DIMENSION WATREF(0:500),CROSTA(500),SOILREF(0:500),SEDREF(0:500),TEL(500),PEL(500),ERM(500),TIR(500) 
	  DIMENSION ALFA(500),BETA(500),SEDREFNAC(500),SOILREFNAC(500),WATREFNAC(500)
!
!
!
      OPEN(UNIT=9,STATUS='OLD',FILE='dataecological.prn')     
!
!	
      READ(9,*)
      READ(9,*)
      READ(9,*)  ! READ USADOS PARA COLOCAR COMENTARIOS E FAZER ABERTURA DO FILE
      READ(9,*)
      READ(9,*)
      READ(9,*)
      READ(9,*)
      READ(9,*)
      READ(9,*)
      READ(9,*)
      READ(9,*) NSPECIE,KEY_WATER,KEY_SOIL,KEY_SEDIMENT
      READ(9,*)
!  
      IF((KEY_WATER.GT.3).OR.(KEY_WATER.LT.1))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE KEYS VALUES: Water reference value               MUST BE National, U.S. EPA OR WHO '')')
	  WRITE(*,*)
	  WRITE(*,'('' For more information, enable key information --- Help --- in the Scenary Sheet '')')
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF
!
      IF(((KEY_SOIL.GT.3).OR.(KEY_SOIL.LT.1)).OR.((KEY_SEDIMENT.GT.3).OR.(KEY_SEDIMENT.LT.1)))THEN
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,'('' THE KEYS VALUES: Soil reference value               MUST BE Regional, National or U.S. EPA '')')
	  WRITE(*,'(''                  Sediment reference value '')')
	  WRITE(*,*)
	  WRITE(*,'('' For more information, enable key information --- Help --- in the Scenary Sheet '')')
	  WRITE(*,*)
	  WRITE(*,'('' THE CODE WILL STOP '')')
	  WRITE(*,*)
	  WRITE(*,*)
	  WRITE(*,*)
	  stop
      ENDIF

!
      WATREF(0)=0.0
	  SOILREF(0)=0.0
	  SEDREF(0)=0.0
!
	  do i=1,500
	  WATREF(i)=0.0
	  SOILREF(i)=0.0
	  CROSTA(i)=0.0
	  SEDREF(i)=0.0
	  TEL(i)=0.0
	  PEL(i)=0.0
	  ERM(i)=0.0
	  TIR(i)=0.0
	  ALFA(i)=0.0
	  BETA(i)=0.0
      DO LOK=1,3
	  WATREF_TEMP(i,LOK)=0.0
	  SEDREF_TEMP(i,LOK)=0.0
	  DO IO=1,4
	  SOILREF_TEMP(i,LOK,IO)=0.0
	  ENDDO
	  ENDDO	
	  enddo
!
!
!
      DO i=1,NSPECIE
!
!
      READ(9,*)
	  READ(9,*)
	  READ(9,*) SPECIE(i)
!
      READ(9,*)
      READ(9,*)
	  READ(9,*)

	  READ(9,*)	WATREF_TEMP(i,1),WATREF_TEMP(i,2),WATREF_TEMP(i,3)	   ! 1 = NACIONAL, 2 = USEPA, 3 = WHO
!
      READ(9,*)
      READ(9,*)
	  READ(9,*)
	  READ(9,*)
!
	  READ(9,*)	SOILREF_TEMP(i,1,1),SOILREF_TEMP(i,1,2),SOILREF_TEMP(i,1,3),SOILREF_TEMP(i,1,4),SOILREF_TEMP(i,2,1),SOILREF_TEMP(i,2,2),SOILREF_TEMP(i,2,3),SOILREF_TEMP(i,2,4),SOILREF_TEMP(i,3,1),SOILREF_TEMP(i,3,2),SOILREF_TEMP(i,3,3),SOILREF_TEMP(i,3,4),CROSTA(i)	   ! 1 = REGIONAL, 2 = NACIONAL, 3 = USEPA   ! 1 = RESIDENCIAL, 2 = AGRÍCOLA, 3 = VALOR DE PREVENÇÃO
!
      READ(9,*)
      READ(9,*)
	  READ(9,*)
!
	  READ(9,*)	SEDREF_TEMP(i,1),SEDREF_TEMP(i,2),SEDREF_TEMP(i,3)   ! 1 = REGIONAL, 2 = NACIONAL, 3 = USEPA   ! 1 = RESIDENCIAL, 2 = AGRÍCOLA, 3 = VALOR DE PREVENÇÃO
!
      READ(9,*)
	  READ(9,*)
!
	  READ(9,*)	TEL(i),PEL(i),ERM(i),TIR(i),ALFA(i),BETA(i)   
!
!
!
      IF(KEY_WATER.EQ.1)THEN
	  WATREF(i)=WATREF_TEMP(i,1)
	  ELSEIF(KEY_WATER.EQ.2)THEN
	  WATREF(i)=WATREF_TEMP(i,2)
	  ELSEIF(KEY_WATER.EQ.3)THEN
	  WATREF(i)=WATREF_TEMP(i,3)
	  ENDIF
!
      IF(KEY_SOIL.EQ.1)THEN
!
      IF(SCENAR.EQ.3)THEN
	  SOILREF(i)=SOILREF_TEMP(i,1,1)
	  ELSEIF(SCENAR.EQ.1)THEN
	  SOILREF(i)=SOILREF_TEMP(i,1,2)
	  ELSEIF(SCENAR.EQ.4)THEN
	  SOILREF(i)=SOILREF_TEMP(i,1,4)
	  ELSEIF(SCENAR.EQ.2)THEN
	  SOILREF(i)=SOILREF_TEMP(i,1,3)
	  ENDIF
!
      ELSEIF(KEY_SOIL.EQ.2)THEN
!
      IF(SCENAR.EQ.3)THEN
	  SOILREF(i)=SOILREF_TEMP(i,2,1)
	  ELSEIF(SCENAR.EQ.1)THEN
	  SOILREF(i)=SOILREF_TEMP(i,2,2)
	  ELSEIF(SCENAR.EQ.4)THEN
	  SOILREF(i)=SOILREF_TEMP(i,2,4)
	  ELSEIF(SCENAR.EQ.2)THEN
	  SOILREF(i)=SOILREF_TEMP(i,2,3)
	  ENDIF
!
      ELSEIF(KEY_SOIL.EQ.3)THEN
!
      IF(SCENAR.EQ.3)THEN
	  SOILREF(i)=SOILREF_TEMP(i,3,1)
	  ELSEIF(SCENAR.EQ.1)THEN
	  SOILREF(i)=SOILREF_TEMP(i,3,2)
	  ELSEIF(SCENAR.EQ.4)THEN
	  SOILREF(i)=SOILREF_TEMP(i,3,4)
	  ELSEIF(SCENAR.EQ.2)THEN
	  SOILREF(i)=SOILREF_TEMP(i,3,3)
	  ENDIF
!
      ENDIF    ! ACABA O IF DE KEY_SOIL
!
!
      IF(KEY_SEDIMENT.EQ.1)THEN
	  SEDREF(i)=SEDREF_TEMP(i,1)
	  ELSEIF(KEY_SEDIMENT.EQ.2)THEN
	  SEDREF(i)=SEDREF_TEMP(i,2)
	  ELSEIF(KEY_SEDIMENT.EQ.3)THEN
	  SEDREF(i)=SEDREF_TEMP(i,3)
	  ENDIF
!
!
      IF(SCENAR.EQ.3)THEN
	  SOILREFNAC(i)=SOILREF_TEMP(i,2,1)
	  ELSEIF(SCENAR.EQ.1)THEN
	  SOILREFNAC(i)=SOILREF_TEMP(i,2,2)
	  ELSEIF(SCENAR.EQ.4)THEN
	  SOILREFNAC(i)=SOILREF_TEMP(i,2,4)
	  ELSEIF(SCENAR.EQ.2)THEN
	  SOILREFNAC(i)=SOILREF_TEMP(i,2,3)
	  ENDIF
!
	  SEDREFNAC(i)=SEDREF_TEMP(i,2)
!
      WATREFNAC(i)=WATREF_TEMP(i,1)
!
      ENDDO   ! FIM DO CICLO i
!
!
 	  RETURN
	  END

