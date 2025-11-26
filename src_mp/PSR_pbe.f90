!**********************************************************************************************

subroutine psr_pbe()

  !**********************************************************************************************
  !
  ! Perfectly stirred reactor for the homogeneous PBE
  ! Stelios Rigopoulos (06/05/2017)
  ! Modified 25/06/2020
  !
  !**********************************************************************************************
  use pbe_mod, only: growth_function, jet_cl_model, amb_temp, amb_rho, current_temp, current_rho, p_water, tau_g, current_XH2O, dv, m, v, v_m
  use pbe_mod, only: Loss_Sw, Production_Sw, Smw, p_sat_liq, p_sat_ice, amb_p, G_mixing_line, Smw_time_series, step_update, activation_logical, consumption_logical
  use pbe_mod, only: plume_cooling_rate, T_time_series, dt_input, ratio_max
  use pbe_mod, only: kappa_bins, part_rho_bins, v0_bins, ni_new, inps_distribution_logical, v, nuclei_logical, activation_logical_bins, S_vc_bins, ni_type, Loss_Sw_bins, m_source_bins, m_source_pbe
  use ice_microphys_mod
  
  implicit none
  
  double precision, allocatable :: ni(:)
  !real, allocatable :: Smw_time_series(:)
  
  double precision moment(0:1)
  double precision int_time,tin,current_time,meansize,dt,particle_mass
  
  integer i,i_step,n_steps,iflag,flowflag,nin,i_write,n_write,i_writesp
  integer agg_kernel_update,n_pbe_grid
  
  !**********************************************************************************************
  
  ! Initialisation
  
  
  ! Initialise PBE 
  call pbe_read(n_pbe_grid)
  allocate(ni(n_pbe_grid))
  allocate(ni_type(n_pbe_grid))  
  call pbe_grid()
  call pbe_init(ni)   
  
  !Initialize inps_distribution_logical to false if then inps_distribution_logical is set to .true. in ice.in, then pbe_grid and pbe_init will be overwritten 
  inps_distribution_logical = .false.


  ! Read PSR input data
  open(30,file='psr/psr.in')
  do i=1,2
    read(30,*)
  end do
  read(30,*) int_time
  read(30,*) dt
  dt_input = dt
  read(30,*) agg_kernel_update
  read(30,*) i_writesp
  read(30,*) n_write
  close(30)
  
  if (growth_function>=4) then
     do i=1,m
      ni(i) = ni(i) / dv(i)
      ni_type(i) = 0.
     end do
     call pbe_ice_read()
     current_temp = amb_temp
     current_rho = amb_rho
     !Open output file
     open(999,file='pbe/ice_jet_temperature.out')
     allocate(Loss_Sw_bins(n_pbe_grid))
     allocate(m_source_bins(n_pbe_grid))
     allocate(m_source_pbe(n_pbe_grid))
     m_source_pbe(:) = 0.
  endif

  if (inps_distribution_logical) then 
    write(*,*) 'Re-initialize PBE grid and operating conditions based on user input file named ice_nucleating_particles.in'    
    call pbe_file_init()
    allocate(activation_logical_bins(n_pbe_grid))
    allocate(S_vc_bins(n_pbe_grid))
    do i=1,m
      ni(i) = ni_new(i) / dv(i)
      !write(*,*)'v_m',v_m(i)
      !v(i) = v0_bins(i)
      S_vc_bins(i) = 0.
      if (nuclei_logical(i)) then
        activation_logical_bins(i) = .false.
      else
        activation_logical_bins(i) = .true.
      endif
    end do
  else
    allocate(nuclei_logical(n_pbe_grid))
    do i=1,m
      nuclei_logical(i) = .false. 
    end do
  endif


  ! Initialise PSR integration
  n_steps = int_time/dt
  current_time= 0.D0
  i_write = 0

  !Allocate array with supersaturation
  allocate(Smw_time_series(1:n_steps*1000))
  allocate(T_time_series(1:n_steps*1000))
  step_update = 0
  
  !Initialize logical for droplet activation
  activation_logical = .false.

  ! Update the thermodynamic state if ice growth and if jet centerline model is activated
  if ((growth_function>=4).and.(jet_cl_model>0)) then
    if (jet_cl_model==1) then
      call pbe_ice_update(current_time, dt, current_temp, current_rho)
    elseif (jet_cl_model==2) then
      write(*,*) 'Using LES data...'
      call pbe_ice_update_LES(current_time, dt, current_temp, current_rho, current_XH2O)
    endif
  endif
  
  !----------------------------------------------------------------------------------------------


  ! Integration
  
  ! Write initial moments
  call PBE_moments(ni,moment,meansize,particle_mass)
  i_step = 1
  !do i_step = 1,n_steps
  do while (current_time < int_time)  
    step_update = i_step
    !Write temperature and growth timescale (tau_g) to output
    if (growth_function>=4) then
      
      !Compute and update saturation ratio based on the input data
      call p_sat_liq_murphy_koop(p_sat_liq)
      call p_sat_ice_murphy_koop(p_sat_ice)
      
      !Append to array
      !call append_scalar(Smw_time_series, Smw)
      Smw_time_series(i_step) = Smw
      T_time_series(i_step) = current_temp
      
      !Supersaturation consumption
      if ((i_step > 1) .and. (consumption_logical) ) then
        ! if (Production_Sw>Loss_Sw .and. Smw_time_series(i_step)>1) then
        !   ! write(*,*) 'Loss_Sw',Loss_Sw
        !   ! write(*,*) 'Production_Sw',Production_Sw
        !   write(*,*) 'temperature', current_temp
        ! endif
        Smw_time_series(i_step) = Smw_time_series(i_step-1) + (Production_Sw - Loss_Sw)*dt 
        plume_cooling_rate = (T_time_series(i_step) - T_time_series(i_step-1))/dt
      endif   
      p_water = Smw_time_series(i_step) * p_sat_liq

      !Force mixing line to not go below ice saturation - we are forcing persisten contrails conditions for jet model = 1
      ! if ((p_water<p_sat_ice) .and. (activation_logical)) then
      !   if (jet_cl_model .eq. 1) then
      !     p_water = p_sat_ice
      !     Smw_time_series(i_step) = p_water/p_sat_liq        
      !   endif
      ! endif


      !Write to output file
      if (inps_distribution_logical) then 
        write(999,1001) current_time,current_temp,current_rho,p_water,tau_g,current_XH2O,Loss_Sw,Smw,Smw_time_series(i_step),maxval(S_vc_bins),moment(0),moment(1),meansize,amb_p,particle_mass,sum(m_source_pbe(:))
      else
        if (activation_logical) then
          write(999,1001) current_time,current_temp,current_rho,p_water,tau_g,current_XH2O,Loss_Sw,Smw,Smw_time_series(i_step),1.0,moment(0),moment(1),meansize,amb_p,particle_mass,sum(m_source_pbe(:))
        else
          write(999,1001) current_time,current_temp,current_rho,p_water,tau_g,current_XH2O,Loss_Sw,Smw,Smw_time_series(i_step),0.0,moment(0),moment(1),meansize,amb_p,particle_mass,sum(m_source_pbe(:))
        endif
      endif  
    endif
  
    ! The following should be done if the kernel should be updated at each time step due to e.g. 
    ! temperature dependency
    if (agg_kernel_update==1) then
      ! Insert here the expression for updating the kernel
      ! agg_kernel_const = 
      call PBE_agg_beta(2)
    end if
  
    ! Integrate
    call pbe_integ(ni,dt)

    !write(*,*) 'dt',dt
  
    ! Calculate moments
    call pbe_moments(ni,moment,meansize,particle_mass)

 
    !if (growth_function>=4) then
    !  if (moment(0)==0) then
    !    call pbe_init(ni)
    !    do i=1,m
    !      ni(i) = ni(i) / dv(i)
    !    end do
    !    call pbe_moments(ni,moment,meansize)
    !  !else 
    !    !write(*,*) 'moment(0) is no longer zero so ice growth starts'
    !  endif
    !endif
  
    ! Write moments
    current_time = current_time + dt
    i_step = i_step + 1
  
    ! Write PSD
    if ((i_write==n_write).or.(i_step==n_steps)) then
      i_write = 0
      call pbe_output(ni,i_writesp,i_step)
    end if
    i_write = i_write + 1
  
    ! Update the thermodynamic state if ice growth and if jet centerline model is activated
    if ((growth_function>=4).and.(jet_cl_model>0)) then
      if (jet_cl_model==1) then
        call pbe_ice_update(current_time, dt, current_temp, current_rho)
      elseif (jet_cl_model==2) then
        call pbe_ice_update_LES(current_time, dt, current_temp, current_rho, current_XH2O)
      endif
  
      !call saturation_ratio(Smw_time_series, k, dt, Smw, Loss_Sw)
    endif
  
  end do

  write(*,*) 'user specified integration time',int_time
  write(*,*) 'actual integration time',current_time
  !write(*,*) 'ratio_max',ratio_max
  
  !----------------------------------------------------------------------------------------------
  
  ! Deallocate arrays
  deallocate(ni)
  deallocate(nuclei_logical)
  if (inps_distribution_logical) then 
    deallocate(kappa_bins)
    deallocate(part_rho_bins)
    deallocate(v0_bins)
  endif
  call PBE_deallocate()
  
  !Close ice ouptut file
  if (growth_function>=4) then
    close(999)
  endif   
  
  1001 format(16E25.12E3)
  
  end subroutine psr_pbe
  
  !**********************************************************************************************