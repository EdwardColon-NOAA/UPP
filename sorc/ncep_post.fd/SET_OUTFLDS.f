!> @file
!> @ brief set_outflds() reads post xml control file.
!>
!> @author J. Wang NCEP/EMC @date 2012-01-27

!> This routine reads the control file in xml format specifying
!> field(s) to post, and save all the field information in
!> a datatype array PSET.
!>
!> @param[in] KTH
!> @param[in] TH
!> @param[in] KPV
!> @param[in] PV
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2012-01-27 | Jun Wang | Initial
!> 2015-03-10 | Lin Gan  | Replace XML file with flat file implementation
!> 2019-10-30 | Bo Cui   | Removw "GOTO" Statement
!>
!> @author J. Wang NCEP/EMC @date 2012-01-27
      SUBROUTINE SET_OUTFLDS(kth,th,kpv,pv)
!

!
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
       use xml_perl_data,  only: paramset,post_avblflds
       use grib2_module,   only: num_pset,pset,nrecout,first_grbtbl,grib_info_init
       use lookup_mod,     only: ITB,JTB,ITBQ,JTBQ
       use ctlblk_mod,     only: npset, me, fld_info
       use rqstfld_mod,    only: mxfld, iget, ritehd, lvlsxml, datset, ident, &
                                 iavblfld, nfld, lvls
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     DECLARE VARIABLES.
!     
      integer, intent(in) :: KTH,KPV
      real, intent(in)    :: th(kth),pv(kpv)
!
      integer L,IFLD,MFLD,IAVBL,IREC,I,J
      CHARACTER*50 AVBLGRB_NAME
      logical :: FOUND_FLD
!
!******************************************************************************
!     START READCNTRL_XML HERE.
!     
!      IF(ME==0)THEN
!        WRITE(6,*)'READCNTRL_XML:  POSTING FCST HR ',IFHR,' FROM ',         &
!             IHRST,'UTC ',SDAT(1),'-',SDAT(2),'-',SDAT(3),' RUN'
!      ENDIF
!     
!     INITIALIZE VARIABLES.
!        ARRAY IGET IS THE "GET FIELD" FLAG ARRAY.
!
      DO IFLD=1,MXFLD
        IGET(IFLD) = -1
      enddo
!
!     SET FLAG TO OPEN NEW OUTPUT FILE
!
      LVLS   = 0
      RITEHD = .TRUE.

! allocate(lvlsxml(MXLVL,num_post_afld))

!$omp parallel do private(i,j)
      DO J=1,size(LVLSXML,2)
        DO I=1,size(LVLSXML,1)
          LVLSXML(I,J) = 0
        ENDDO
      ENDDO
!
      pset   = paramset(npset)
      datset = pset%datset
      if (me==0)print *,'in SET_OUTFLDS, num_pset=',num_pset,'datset=',trim(pset%datset),'npset=',npset
! 
!     NOW READ WHICH FIELDS ON 
!     WHICH LEVELS TO INTERPOLATE TO THE OUTPUT GRID.  THE
!     CHARACTER STRING "DONE" MARKS THE END OF THE OUTPUT
!     FIELD SPECIFICATIONS.
!
      call grib_info_init()
      MFLD = size(pset%param)

! LinGan set post_avblflds to current working paramset
! This is required for flat file solution to work for nmm

      post_avblflds%param =>paramset(npset)%param
      if (me==0) then
        write(*,*)'Size of pset is: ',MFLD
        write(*,*)'datset is: ',datset
        write(*,*)'MXFLD is: ',MXFLD
        write(*,*)'size of lvlsxml: ',size(lvlsxml)
        write(*,*)'size of post_avblflds param',size(post_avblflds%param)
      endif
      if(size(post_avblflds%param) <= 0) then
         write(0,*)'WRONG: post available fields not ready!!!'
         return
      endif
!
      IFLD = 0
      irec = 0
      DO I=1, MFLD
         
!     SEE IF REQUESTED FIELD IS AVAILABLE.  IF NOT, 
!     WRITE MESSAGE TO 6 AND DECREMENT FIELD 
!     COUNTER BY ONE.  THEN READ NEXT REQUESTED FIELD.
!     
!     GET POST AVAILBLE FIELD INDEX NUMBER FOR EACH OUTPUT FIELDS IN PSET
!

         FOUND_FLD = .false.

!        write(*,*)'cntfile,i=',i,'fld shortname=',trim(pset%param(i)%shortname)
!        write(*,*)'size(post_avblflds%param)=',size(post_avblflds%param)

         IFLD = IFLD + 1

!     segmentation fault occurred on nmm i=112

         IAVBL          = post_avblflds%param(i)%post_avblfldidx
         IGET(IAVBL)    = IFLD
         IDENT(IFLD)    = IAVBL
         IAVBLFLD(IFLD) = I
         FOUND_FLD      = .true.
         call set_lvlsxml(pset%param(i),ifld,irec,kpv,pv,kth,th)

      ENDDO

!     
!     ALL DONE FOUNDING REQUESTED FIELDS FOR current OUTPUT GRID.
!     SET NFLD TO TOTAL NUMBER OF REQUESTED OUTPUT FIELDS THAT 
!     ARE AVAILABLE., SET NRECOUT to total number of OUTPUT records 
!     NOTE: here NFLD i s total number of fields found in post_avblfld_table,
!           while nrecoutis the total number of grib messages that go 
!           into the output file. One fieldmay contain many different levels,
!           which each different level will be counted as one record
!
      NFLD    = IFLD
      NRECOUT = IREC
! Meng 04/19/18, add three fields for continous bucket
!      NRECOUT = IREC + 3
      allocate(fld_info(NRECOUT+100))
      do i=1,nrecout
        fld_info(i)%ifld     = 0
        fld_info(i)%lvl      = 0
        fld_info(i)%lvl1     = 0
        fld_info(i)%lvl2     = 0
        fld_info(i)%ntrange  = 0
        fld_info(i)%tinvstat = 0
      enddo
      if(me==0)write(*,*)'in readxml. nfld=',nfld,'nrecout=',nrecout
!
! skip creating ipv files if kth=0 and no isobaric fields are requested in ctl file      
!     if(kth == 0 .and. iget(013) <= 0) go to 999
!     
!     ECHO OUTPUT FIELDS/LEVELS TO 6.
!
!      IF(ME==0)THEN
!        WRITE(6,*)'BELOW ARE FIELD/LEVEL/SMOOTHING ',       &
!             'SPECIFICATIONS.,NFLD=',NFLD,'MXLVL=',MXLVL,'nrecout=',nrecout
!      ENDIF
!      DO 50 IFLD = 1,NFLD
!        IF(ME==0)THEN
!         i=IAVBLFLD(IFLD)
!         write(*,*)'readxml,ifld=',ifld,'iget(',IDENT(ifld),')=',iget(ident(ifld)),'iavbl=',IAVBLFLD(iget(ident(ifld))),'postvar=',trim(pset%param(i)%pname),  &
!             trim(pset%param(i)%fixed_sfc1_type),'lvls=',LVLS(:,ifld)
!         if(size(pset%param(i)%level)>0) then
!           WRITE(*,*) pset%param(i)%level
!         endif
!        ENDIF
! 50   CONTINUE
!     
!     END OF ROUTINE.
!     
 999  CONTINUE

       if(me==0)print *,'end of read_postcntrl_xml'
      RETURN
      END
