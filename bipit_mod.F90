!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module bipit_mod

!BOP
! !MODULE: bipit_mod
!
! !DESCRIPTION: started as copy of iage, being adjusted to a 
!               mixed layer residence time
!
! !REVISION HISTORY:
!  SVN:$Id: bipit_mod.F90 85763 2017-06-20 01:19:28Z altuntas@ucar.edu $
! 8/20/18 updating from 1 tracer, tau_ml, to 2, nutrient
! 8/30/18 adding particle (phyt-detritus), parti

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use blocks, only: nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km
   use domain, only: nblocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use kinds_mod
   use constants, only: c0, c1, c2, c10, char_blank, delim_fmt
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
   use passive_tracer_tools, only: ind_name_pair, rest_read_tracer_block
   use passive_tracer_tools, only: file_read_tracer_block, tracer_read
   use grid, only: dzr, zw, dz, KMT
   use vmix_kpp, only: HMXL
   use forcing_shf, only: SHF_QSW
   use sw_absorption, only: sw_absorb_frac
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: bipit_tracer_cnt,        &
             bipit_init,              &
             bipit_set_interior,      &
             bipit_reset

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      bipit_tracer_cnt = 3

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      tauml_ind = 1, &    ! bipit index
      nutri_ind = 2, &
      parti_ind = 3

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(bipit_tracer_cnt) :: &
      ind_name_table = (/ & 
       ind_name_pair(tauml_ind, 'TAUML'), &
       ind_name_pair(nutri_ind, 'NUTRI'), &
       ind_name_pair(parti_ind, 'PARTI') /)

!-----------------------------------------------------------------------
!  tavg ids for non-standard tavg variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_BIPIT_RESET_TEND       ! tavg id for surface reset tendency of BIPIT

!EOC
!*****************************************************************************
!--------------------------------------------------------------------------
! Additional module variables
!--------------------------------------------------------------------------
   integer(int_kind) :: &
          kn = 40                      ! layer to reset at and below; was 31
   integer(int_kind) :: &
          kL = 21                      ! layer to grow particles above
          
   real(r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      PLASTFLUX ! flux of P, sinking from above !added March 2019
contains

!*****************************************************************************
!BOP
! !IROUTINE: bipit_init
! !INTERFACE:

 subroutine bipit_init(bipit_ind_begin, init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize bipit tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  7/19/2018, accumulating for k<kn
!  7/27/2018, accumulates for zw(k)<HMXL

! !USES:

   use broadcast, only: broadcast_scalar
   use prognostic, only: curtime, oldtime, tracer_field
   use grid, only: KMT, n_topo_smooth, fill_points

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      bipit_ind_begin         ! starting index of bipit tracers in global tracer array
                             ! passed through to rest_read_tracer_block

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(bipit_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,bipit_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'bipit_mod:bipit_init'

   character(char_len) :: &
      init_bipit_option,           & ! option for initialization of bipit
      init_bipit_init_file,        & ! filename for option 'file'
      init_bipit_init_file_fmt       ! file format for option 'file'

   logical(log_kind) :: &
      lnml_found             ! Was bipit_nml found ?

   integer(int_kind) :: &
      n,                   & ! index for looping over tracers
      k,                   & ! index for looping over depth levels
      iblock,              & ! index for looping over blocks
      nml_error              ! namelist i/o error flag

!     l,                   & ! index for looping over time levels

   type(tracer_read), dimension(bipit_tracer_cnt) :: &
      tracer_init_ext        ! namelist variable for initializing tracers

   namelist /bipit_nml/ &
      init_bipit_option, init_bipit_init_file, tracer_init_ext, &
      init_bipit_init_file_fmt

   character (char_len) ::  &
      bipit_restart_filename  ! modified file name for restart file
!----------------------------------------------------------------------
!    set kn if needed
!---------------------------------------------------------------------   
!   kn = 4      ! reset levels 4:k; not needed here yet
!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(tauml_ind)%short_name = 'TAUML'
   tracer_d_module(tauml_ind)%long_name  = 'Mixed Layer Residence Time'
   tracer_d_module(tauml_ind)%units      = 'years'
   tracer_d_module(tauml_ind)%tend_units = 'years/s'
   tracer_d_module(tauml_ind)%flux_units = 'cm years/s'

   tracer_d_module(nutri_ind)%short_name = 'NUTRI'
   tracer_d_module(nutri_ind)%long_name  = 'Generalized nutrient'
   tracer_d_module(nutri_ind)%units      = 'mmol'
   tracer_d_module(nutri_ind)%tend_units = 'mmol/s'
   tracer_d_module(nutri_ind)%flux_units = 'cm mmol/s'
   
   tracer_d_module(parti_ind)%short_name = 'PARTI'
   tracer_d_module(parti_ind)%long_name  = 'Generalized phytoplankton/detritus'
   tracer_d_module(parti_ind)%units      = 'mmol'
   tracer_d_module(parti_ind)%tend_units = 'mmol/s'
   tracer_d_module(parti_ind)%flux_units = 'cm mmol/s'
!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_bipit_option = 'unknown'
   init_bipit_init_file = 'unknown'
   init_bipit_init_file_fmt = 'bin'

   do n = 1,bipit_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then  
         nml_error = -1
      else
         nml_error =  1      
      endif
      do while (nml_error > 0)
         read(nml_in, nml=bipit_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'init_bipit_option : ' // init_bipit_option)
      call document(subname, 'init_bipit_init_file : ' // init_bipit_init_file)
      call document(subname, 'bipit_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_bipit_option , master_task)
   call broadcast_scalar(init_bipit_init_file, master_task)
   call broadcast_scalar(init_bipit_init_file_fmt, master_task)

   do n = 1,bipit_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

!-----------------------------------------------------------------------
!  initialize tracers
!-----------------------------------------------------------------------

   select case (init_bipit_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup', 'ccsm_hybrid')
      TRACER_MODULE(:,:,:,tauml_ind,:,:) = c0
      TRACER_MODULE(:,:,1:5,nutri_ind,:,:) = c1
      TRACER_MODULE(:,:,6:10,nutri_ind,:,:) = 3 
      TRACER_MODULE(:,:,11:15,nutri_ind,:,:) = 5
      TRACER_MODULE(:,:,16:kL,nutri_ind,:,:) = c10
      TRACER_MODULE(:,:,kL+1:km,nutri_ind,:,:) = 20
      TRACER_MODULE(:,:,:,parti_ind,:,:) = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d BIPIT set to all zeros for TAUML and PARTI and a vertical structure for NUTRI' 
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif
       
   case ('restart', 'ccsm_continue', 'ccsm_branch' )

      bipit_restart_filename = char_blank

      if (init_bipit_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read bipit from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         bipit_restart_filename = read_restart_filename
         init_bipit_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         bipit_restart_filename = trim(init_bipit_init_file)

      endif

      call rest_read_tracer_block(bipit_ind_begin,          &
                                  init_bipit_init_file_fmt, &
                                  bipit_restart_filename,   &
                                  tracer_d_module,         &
                                  TRACER_MODULE)

   case ('file')
      call document(subname, 'bipit being read from separate file')

      call file_read_tracer_block(init_bipit_init_file_fmt, &
                                  init_bipit_init_file,     &
                                  tracer_d_module,         &
                                  ind_name_table,          &
                                  tracer_init_ext,         &
                                  TRACER_MODULE)
 
      if (n_topo_smooth > 0) then
   !      do n=1,bipit_tracer_cnt
         do k=1,km
            call fill_points(k,TRACER_MODULE(:,:,k,1,curtime,:), errorCode)

            if (errorCode /= POP_Success) then
               call POP_ErrorSet(errorCode, &
                  'bipit_init: error in fill_points')
               return
            endif
         enddo
   !      end do
      endif

   case default
      call document(subname, 'unknown init_bipit_option = ', init_bipit_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,bipit_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------

   call define_tavg_field(tavg_BIPIT_RESET_TEND, 'BIPIT_RESET_TEND',3,  & !2->3 on 7/20/18
                          long_name='reset tendency of BIPIT', &
                          units='years/s', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!EOC

 end subroutine bipit_init

!***********************************************************************
!BOP
! !IROUTINE: bipit_set_interior
! !INTERFACE:

 subroutine bipit_set_interior(k, DTRACER_MODULE, bid, TRACER)

! !DESCRIPTION:
!  set interior source/sink term for bipit tracer
!
! !REVISION HISTORY:
!  7/19/2018 updated source to be 0 for k>=kn, rather than for k=1
!  7/20      c1/c2 rather than c1/seconds_per_year, because was getting -inf at surface
!  7/23 above did not fix inf; trying source for only top 2 layers to see if layer 3 is then finite
!  7/25 above did not fix; trying just layer 2 as source in case layer 1 has issues due to variable volume
!  7/27 above problems due to robert filter for 1-degree pop; using avgfit instead fixes it
!       updating to accumulate for z<HMXL
!  8/20 when using HMXL, HBLT, or HMXL_DR, the smoothness of those depths (related to physics) 
!       prevents the problems with advection
!       Now adding NUTRI
!  8/29 Nutri working without light, adding light effects
! !USES:

   use time_management, only: seconds_in_year

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: &
      k                   ! vertical level index

   integer(int_kind), intent(in) :: bid

   real(r8), dimension(nx_block,ny_block,km,bipit_tracer_cnt), intent(in) :: &
      TRACER
! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,bipit_tracer_cnt), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink term

!  local variables
   !real(r8), dimension(nx_block,ny_block) :: & !now declared at module level
   !   PLASTFLUX ! flux of P, sinking from above !added March 2019
   real(r8), dimension(nx_block,ny_block) :: &
      DTAUML_LOC != TRACER_MODULE(:,:,:,tauml_ind)      ! tauml tracer
   real(r8), dimension(nx_block,ny_block) :: &
      DPARTI_LOC != TRACER_MODULE(:,:,:,tauml_ind)      ! parti tracer
   real(r8), dimension(nx_block,ny_block) :: &
      DNUTRI_LOC != TRACER_MODULE(:,:,:,nutri_ind)      ! nutri tracer
   real(r8), dimension(nx_block,ny_block,km) :: &
      NUTRI_LOC ! = TRACER(:,:,:,nutri_ind)      ! nutri tracer
   real(r8), dimension(nx_block,ny_block,km) :: &
      PARTI_LOC ! = TRACER(:,:,:,nutri_ind)      ! parti tracer
   real(r8), dimension(km) :: &
      sw_absorb
   real(r8), dimension(nx_block,ny_block,km) :: &
      sw_3
   real(r8), dimension(nx_block,ny_block) :: &
      sw_mlmean
   real(r8), dimension(nx_block,ny_block) :: &
      MUN       !N/(N+k)
   real(r8), dimension(nx_block,ny_block) :: &
      MUL       !(1-exp(-0.05*sw_absorb*SHF))
   real(r8), parameter :: &
      WSINK = 2.5*100.0/86400.0 !m/day * 100cm/m / 86400s/day
   real(r8), parameter :: &
      ALPHAI = 0.035 
   real(r8), parameter :: &
      MUMAX = 0.37/86400.0 !doubling/day  / 86400s/day
   real(r8), parameter :: &
      KN = 3.2 
   real(r8), parameter :: &
      TAUDECAY = 50*86400.0 !decay from P to nothing (lifetime of P+lifetime of D)
   integer(int_kind) :: &
      ki
!EOP
!BOC
!  
!----------------------------------------------------------------------
!    set kn
!---------------------------------------------------------------------   
!   kn = 10      ! restore levels 10:km
!-----------------------------------------------------------------------
   !WSINK=1.0*100/86400 !m/day * 100cm/m / 86400s/day
   sw_absorb = c0 
    do ki=1,kL
         call sw_absorb_frac(zw(ki),sw_absorb(ki))
      where (zw(ki) < HMXL(:,:,bid)) 
        sw_3(:,:,ki)=sw_absorb(ki)*dz(ki)
   !    sw_3(:,:,ki) = 1-exp(-ALPHAI*10000*sw_absorb(k)*SHF_QSW(:,:,bid))
      elsewhere
         sw_3(:,:,ki) = c0
      end where
    end do
    sw_mlmean=sum(sw_3,3)/HMXL(:,:,bid)
 
    where (zw(k) < HMXL(:,:,bid)) 
       DTAUML_LOC = c1 / seconds_in_year
    end where

    where (zw(k-1) < HMXL(:,:,bid) .AND. zw(k) > HMXL(:,:,bid)) 
          DTAUML_LOC  = (HMXL(:,:,bid)-zw(k-1))*dzr(k)*c1/seconds_in_year
    end where

    where (zw(k-1) > HMXL(:,:,bid)) 
          DTAUML_LOC = c0
    end where
    
    NUTRI_LOC=TRACER(:,:,:,nutri_ind)
    PARTI_LOC=TRACER(:,:,:,parti_ind)
    DPARTI_LOC=c0
    DNUTRI_LOC=c0
    if (k<kL) then
       call sw_absorb_frac(zw(k),sw_absorb(k))
       where (zw(k) < HMXL(:,:,bid))
       MUL = 1-exp(-ALPHAI*10000*sw_mlmean*SHF_QSW(:,:,bid))
    !   MUL = sw_mlmean
       elsewhere
       MUL = 1-exp(-ALPHAI*10000*sw_absorb(k)*SHF_QSW(:,:,bid))
       end where
     where (NUTRI_LOC(:,:,k)>0)
       MUN = NUTRI_LOC(:,:,k)/(NUTRI_LOC(:,:,k)+KN)
       DNUTRI_LOC = -MUMAX*MUL*MUN
       DPARTI_LOC = -DNUTRI_LOC !- PARTI_LOC(:,:,k)/(30.0*86400)
     end where
    endif
   !sinking
    if (k.eq.1) then
       where (PARTI_LOC(:,:,k)>0)
          PLASTFLUX(:,:,bid)= WSINK*PARTI_LOC(:,:,k)  !W=1000cm/day=10m/day
          DPARTI_LOC=DPARTI_LOC - PLASTFLUX(:,:,bid)*dzr(k) !W=1000cm/day=10m/day
       elsewhere
          PLASTFLUX(:,:,bid)=0.0
       end where
    else
       where (PARTI_LOC(:,:,k)>0)
          DPARTI_LOC=DPARTI_LOC + PLASTFLUX(:,:,bid)*dzr(k)  
          PLASTFLUX(:,:,bid)=WSINK*PARTI_LOC(:,:,k)  
          DPARTI_LOC=DPARTI_LOC - PLASTFLUX(:,:,bid)*dzr(k)                
       elsewhere
          PLASTFLUX(:,:,bid)=0.0
       end where
    endif
   !decay
    where (PARTI_LOC(:,:,k)>0)
       DPARTI_LOC = DPARTI_LOC - PARTI_LOC(:,:,k)/TAUDECAY
    end where
    DTRACER_MODULE(:,:,tauml_ind)=DTAUML_LOC
    DTRACER_MODULE(:,:,nutri_ind)=DNUTRI_LOC 
    DTRACER_MODULE(:,:,parti_ind)=DPARTI_LOC 
!-----------------------------------------------------------------------
!EOC

 end subroutine bipit_set_interior

!***********************************************************************
!BOP
! !IROUTINE: bipit_reset
! !INTERFACE:

 subroutine bipit_reset(TRACER_MODULE, bid)

! !DESCRIPTION:
!  reset surface value for bipit tracer
!
! !REVISION HISTORY:
!  7/19/2018 reset Tracer for kn:km, rather than just level 1
!             did not update tavg_BIPIT_RESET_TEND accordingly yet
!  7/10 changed WORK loop; the deep values are 0 as desired, but 
!               much of surface is -inf
! 7/27/2018 WORK loop is not being run
!           changing reset to set to 0 below HMXL and reduce for cells partly below
! 8/20/2018 resetting NUTRI to 10 for depth level 10 and below

! !USES:

   use time_management, only: mix_pass, c2dtt
   use prognostic, only: PSURF, newtime
   use constants, only: grav
   use grid, only: dz

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: bid

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,km,bipit_tracer_cnt), intent(inout) :: &
      TRACER_MODULE      ! bipit tracer

!EOP
!BOC
!----------------------------------------------------------------------
!    set kn
!---------------------------------------------------------------------   
!   kn = 10      ! reset levels 10:km
!----------------------------------------------------------------------
!  local variables
!----------------------------------------------------------------------
   real (r8), dimension(nx_block,ny_block) :: WORK

   real(r8), dimension(nx_block,ny_block,km) :: &
      TAUML_LOC != TRACER_MODULE(:,:,:,tauml_ind)      ! tauml tracer
   real(r8), dimension(nx_block,ny_block,km) :: &
      NUTRI_LOC != TRACER_MODULE(:,:,:,nutri_ind)      ! nutri tracer
   real(r8), dimension(nx_block,ny_block,km) :: &
      PARTI_LOC != TRACER_MODULE(:,:,:,nutri_ind)      ! nutri tracer
   integer(int_kind) :: &
      k                   ! vertical level index
!-----------------------------------------------------------------------

!   if (mix_pass /= 1) then
!      if (accumulate_tavg_now(tavg_BIPIT_RESET_TEND)) then
!         do k=kn,km
!            WORK = -TRACER_MODULE(:,:,k,bid) !/ c2dtt(k)
!         !   WORK = WORK * (c1 + PSURF(:,:,newtime,bid)/grav/dz(k))
!            call accumulate_tavg_field(WORK,tavg_BIPIT_RESET_TEND,bid,k)
!         enddo
!      endif
!   endif

!   TRACER_MODULE(:,:,kn:km,:) = c0
  TAUML_LOC = TRACER_MODULE(:,:,:,tauml_ind)
  NUTRI_LOC = TRACER_MODULE(:,:,:,nutri_ind)
  PARTI_LOC = TRACER_MODULE(:,:,:,parti_ind)
  do k=1,km
    where (zw(k-1) < HMXL(:,:,bid) .AND. zw(k) > HMXL(:,:,bid)) 
          TAUML_LOC(:,:,k)  = (HMXL(:,:,bid)-zw(k-1))*dzr(k)*TAUML_LOC(:,:,k)
    end where

    where (zw(k-1) > HMXL(:,:,bid)) 
          TAUML_LOC(:,:,k) = c0
    end where
  enddo

  NUTRI_LOC(:,:,kn:km) = 18
  where (NUTRI_LOC<0)
    NUTRI_LOC=c0
  end where
  where (PARTI_LOC<0)
    PARTI_LOC=c0
  end where
  TRACER_MODULE(:,:,:,tauml_ind)=TAUML_LOC
  TRACER_MODULE(:,:,:,nutri_ind)=NUTRI_LOC
  TRACER_MODULE(:,:,:,parti_ind)=PARTI_LOC
!-----------------------------------------------------------------------
!EOC

 end subroutine bipit_reset

!***********************************************************************

end module bipit_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
