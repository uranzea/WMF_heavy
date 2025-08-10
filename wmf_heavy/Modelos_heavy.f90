!=======================================================================
!  MÓDULO: models
!  Descripción:
!    Implementa un conjunto de modelos hidrológicos distribuidos para 
!    simulación en cuencas: lluvia, escorrentía, sedimentos, 
!    deslizamientos e inundaciones.
!
!  Autores originales:
!    Nicolás Velásquez Girón  (nicolas.velasquezgiron@gmail.com)
!    Adaptado de código de Esneider Zapata Atehortua
!
!  Licencia:
!    Software Libre bajo GNU GPL v3 o posterior.
!=======================================================================

module models
  implicit none

  !--------------------------------------------------------------------
  ! VARIABLES GLOBALES DE LA MODELACIÓN
  !--------------------------------------------------------------------

  ! Coordenadas esquina inferior izquierda de la cuenca
  real :: xll, yll

  ! Valor que representa dato nulo
  real :: noData

  ! Tamaños de celda
  real :: dx        ! Largo del mapa leído
  real :: dxP       ! Largo proyectado del mapa

  ! Dimensiones de malla
  integer :: ncols, nrows          ! Número de columnas y filas
  integer :: nceldas               ! Número total de celdas en la cuenca
  integer :: verbose                ! Indicador nivel de salida por pantalla
  
  ! Rutas de archivo
  character(len=500) :: rute_speed
  character(len=500) :: rute_storage

  ! Identificadores y posiciones de eventos de lluvia
  integer, allocatable :: idEvento(:)
  integer, allocatable :: posEvento(:)

  ! Para separación de flujos por tipo de lluvia
  integer, allocatable :: posConv(:)
  integer, allocatable :: posStra(:)

  ! Propiedades geomorfológicas
  integer, allocatable :: drena(:,:)         ! Topología de la cuenca
  integer, allocatable :: unit_type(:,:)     ! Tipo de celda
  real, allocatable    :: hill_long(:,:), hill_slope(:,:)
  real, allocatable    :: stream_long(:,:), stream_slope(:,:)
  real, allocatable    :: stream_width(:,:)
  real, allocatable    :: elem_area(:,:)

  ! Propiedades físicas
  real, allocatable :: v_coef(:,:)           ! Velocidades verticales
  real, allocatable :: v_exp(:,:)            ! Exponentes verticales
  real, allocatable :: h_coef(:,:)           ! Velocidades horizontales
  real, allocatable :: h_exp(:,:)            ! Exponentes horizontales
  real, allocatable :: Max_capilar(:,:)
  real, allocatable :: Max_gravita(:,:)
  real, allocatable :: Max_aquifer(:,:)
  real, allocatable :: Retorned(:,:)
  real, allocatable :: EvpSerie(:)

  real :: retorno_gr   ! 0 = no retorno, 1 = sí retorno tanque 3→2
  real :: retorno_aq   ! Retorno máximo del tanque 4→3

  ! Variables de control y configuración del modelo
  real    :: dt                    ! Paso de tiempo de la modelación
  integer :: rain_first_point
  integer :: sim_sediments
  integer :: sim_slides
  integer :: sim_floods
  integer :: save_storage
  integer :: save_speed
  integer :: save_retorno
  integer :: save_vfluxes
  integer :: save_rc
  integer :: show_storage
  integer :: show_speed
  integer :: show_mean_speed
  integer :: show_mean_retorno
  integer :: show_area
  integer :: separate_fluxes
  integer :: separate_rain
  integer :: speed_type(3)
  integer, allocatable :: control(:,:)      
  integer, allocatable :: control_h(:,:)    
  integer, allocatable :: guarda_cond(:)
  integer, allocatable :: guarda_vfluxes(:)
  integer :: calc_niter

  ! Resultados globales
  real, allocatable :: Storage(:,:)
  real              :: storage_constant
  real, allocatable :: Speed_map(:,:)
  real, allocatable :: Mean_Rain(:,:)
  real, allocatable :: Acum_rain(:,:)
  real, allocatable :: Fluxes(:,:)
  real, allocatable :: Storage_conv(:,:)
  real, allocatable :: Storage_stra(:,:)
  real, allocatable :: mean_storage(:,:)
  real, allocatable :: mean_speed(:,:)
  real, allocatable :: mean_retorno(:)
  real, allocatable :: mean_vfluxes(:,:)
  real, allocatable :: vfluxes(:,:)
  real, allocatable :: rc_coef(:,:)

  ! Parámetros para sedimentos
  real :: sed_factor
  real :: wi(3), Qskr, G, diametro(3)
  real :: qlin_sed
  real :: ERO(3), EROt(3), DEP(3), DEPt(3)
  real, allocatable :: VolERO(:), VolDEPo(:)
  real, allocatable :: Vs(:,:), Vd(:,:), VSc(:,:), Vdc(:,:)
  real, allocatable :: Krus(:,:), Crus(:,:), Prus(:,:)
  real, allocatable :: PArLiAc(:,:)

  ! Variables para modelo de deslizamientos
  integer :: sl_GullieNoGullie
  real    :: sl_FS
  real, allocatable :: sl_RiskVector(:,:)
  real, allocatable :: sl_SlideOcurrence(:,:)
  integer, allocatable :: sl_SlideAcumulate(:,:)
  integer, allocatable :: sl_SlideNcellTime(:)
  real, allocatable    :: sl_Zcrit(:,:), sl_Zmin(:,:), sl_Zmax(:,:), sl_Bo(:,:)
  real, allocatable    :: sl_Zs(:,:)
  real                 :: sl_GammaW
  real, allocatable    :: sl_GammaS(:,:), sl_Cohesion(:,:)
  real, allocatable    :: sl_FrictionAngle(:,:), sl_RadSlope(:,:)

  ! Variables del modelo de inundaciones
  real, allocatable :: flood_Q(:,:), flood_Qsed(:,:)
  real, allocatable :: flood_h(:,:)
  integer, allocatable :: flood_flood(:,:)
  real, allocatable :: flood_speed(:,:), flood_ufr(:,:)
  real              :: flood_rdf
  real, allocatable :: flood_Cr(:,:)
  integer, allocatable :: flood_eval(:,:)
  real :: flood_area, flood_diff, flood_sec_tam
  real :: flood_AV
  real, allocatable :: flood_w(:,:), flood_d50(:,:)
  integer, allocatable :: flood_aquien(:,:)
  real, allocatable :: flood_hand(:,:), flood_loc_hand(:,:)
  real, allocatable :: flood_sections(:,:), flood_sec_cells(:,:)
  real :: flood_Cmax
  real, allocatable :: flood_slope(:,:)
  real :: flood_dw, flood_dsed, flood_umbral
  integer :: flood_max_iter
  real :: flood_step, flood_hmax
  real, allocatable :: flood_profundidad(:,:)

  !=====================================================================
  ! SUBRUTINA PRINCIPAL
  ! shia_v1
  !---------------------------------------------------------------------
  ! Ejecuta el modelo hidrológico de la cuenca:
  !   - Lectura de lluvia
  !   - Cálculo de escorrentía
  !   - Flujos verticales y horizontales
  !   - Procesos opcionales: sedimentos, deslizamientos, inundaciones
  !=====================================================================
  contains

  subroutine shia_v1( ruta_bin, ruta_hdr, calib, StoIn, HspeedIn, &
                      N_cel, N_cont, N_contH, N_reg, Q, Qsed, Qseparated, &
                      Hum, St1, St3, balance, speed, AreaControl, StoOut, &
                      ruta_storage, ruta_speed, ruta_vfluxes, ruta_binConv, &
                      ruta_binStra, ruta_hdrConv, ruta_hdrStra, Qsep_byrain, &
                      ruta_retorno, ruta_rc )

    integer, intent(in) :: N_cel, N_reg, N_cont, N_contH
    real,    intent(in) :: calib(11)
    character(len=500), intent(in) :: ruta_bin, ruta_hdr
    character(len=500), intent(in), optional :: ruta_storage
    character(len=500), intent(in), optional :: ruta_binConv, ruta_hdrConv
    character(len=500), intent(in), optional :: ruta_binStra, ruta_hdrStra
    character(len=500), intent(in), optional :: ruta_speed, ruta_retorno
    character(len=500), intent(in), optional :: ruta_vfluxes, ruta_rc
    real,    intent(in), optional :: StoIn(5, N_cel)
    real,    intent(in), optional :: HspeedIn(4, N_cel)
    real, intent(out) :: Hum(N_contH, N_reg)
    real, intent(out) :: St1(N_contH, N_reg), St3(N_contH, N_reg)
    real, intent(out) :: Q(N_cont, N_reg)
    real, intent(out) :: Qsed(N_cont, 3, N_reg)
    real, intent(out) :: Qseparated(N_cont, 3, N_reg)
    real, intent(out) :: Qsep_byrain(N_cont, 2, N_reg)
    real, intent(out) :: StoOut(5, N_cel), balance(N_reg)
    real, intent(out) :: speed(N_cont, N_reg)
    real, intent(out) :: AreaControl(N_cont, N_reg)

    real :: Rain(N_cel)
    integer :: RainInt(N_cel)
    integer :: Conv(N_cel), Stra(N_cel)
    integer :: Co, St, Res
    real    :: rain_sum
    integer :: celda, tiempo, drenaid
    integer :: control_cont, controlh_cont, i
    real    :: entradas, salidas, StoAtras
    real    :: m3_mmHill(N_cel), m3_mmRivers(N_cel)
    real    :: vflux(4), hflux(4), hflux_c(4), hflux_s(4)
    real    :: Ret, Ret_aq, Evp_loss
    real    :: QfluxesOut(3)
    real    :: vspeed(4, N_cel), hspeed(4, N_cel)
    real    :: section_area
    real    :: H(3, N_cel)
    real    :: Area_coef(nceldas)
    real    :: Vsal_sed(3)

    call rain_read_ascii_table( ruta_hdr, N_reg )

    if (allocated(Mean_Rain)) deallocate(Mean_Rain)
    allocate(Mean_Rain(1, N_reg))
    Mean_Rain = 0.0

    if (allocated(Acum_rain)) deallocate(Acum_rain)
    allocate(Acum_rain(1, N_cel))
    Acum_rain = 0.0

    m3_mmHill   = elem_area(1, :) / 1000.0
    m3_mmRivers = elem_area(1, :) / 1000.0

    Q = 0.0

    if (StoIn(1,1) .gt. 0) then
      StoOut = StoIn
    else
      StoOut = Storage + storage_constant
    end if

    entradas = 0.0
    salidas  = 0.0
    balance  = 0.0

    if (HspeedIn(1,1) .gt. 0) then
      hspeed = HspeedIn
    else
      do i = 1, 4
        vspeed(i, :) = v_coef(i, :) * calib(i) * dt
        if (speed_type(i) == 1) then
          hspeed(i, :) = h_coef(i, :) * calib(i + 4)
        else
          hspeed(i, :) = 0.0
        end if
      end do
    end if

    H(1, :) = Max_capilar(1, :) * calib(9)
    H(2, :) = Max_gravita(1, :) * calib(10)
    H(3, :) = Max_aquifer(1, :) * calib(11)
    if (retorno_gr == 1) then
      if (allocated(Retorned)) deallocate(Retorned)
      allocate(Retorned(1, N_cel))
      Retorned = 0.0

      if (show_mean_retorno == 1) then
        if (allocated(mean_retorno)) deallocate(mean_retorno)
        allocate(mean_retorno(N_reg))
        mean_retorno = 0.0
      end if
    end if

    if (save_speed == 1) then
      if (allocated(Speed_map)) deallocate(Speed_map)
      allocate(Speed_map(4, N_cel))
    end if

    if (sim_sediments == 1) then
      call sed_allocate(N_cel)
      Qsed = 0.0
    end if

    if (sim_slides == 1) then
      call slide_allocate(N_cel, N_reg)
    end if

    if (sim_floods == 1) then
      call flood_allocate(N_cel)
    end if

    if (separate_fluxes == 1) then
      if (allocated(Fluxes)) deallocate(Fluxes)
      allocate(Fluxes(3, N_cel))
      Fluxes = 0.0
    end if

    if (separate_rain == 1) then
      if (allocated(Storage_conv)) deallocate(Storage_conv)
      if (allocated(Storage_stra)) deallocate(Storage_stra)
      allocate(Storage_conv(5, N_cel), Storage_stra(5, N_cel))
      Storage_conv = 0.0
      Storage_stra = 0.0
      hflux_c = 0.0
      hflux_s = 0.0
      Conv    = 0
      Stra    = 0
      Qsep_byrain = 0.0
      call rain_read_ascii_table_separate(ruta_hdrConv, ruta_hdrStra, N_reg)
    end if

    if (show_storage == 1) then
      if (allocated(mean_storage)) deallocate(mean_storage)
      allocate(mean_storage(5, N_reg))
      mean_storage = 0.0
    end if

    if (show_mean_speed == 1) then
      if (allocated(mean_speed)) deallocate(mean_speed)
      allocate(mean_speed(4, N_reg))
      mean_speed = 0.0
    end if

    if (save_storage == 1) then
      if (.not. allocated(guarda_cond)) then
        allocate(guarda_cond(N_reg))
        guarda_cond = 0
      else
        if (size(guarda_cond) < N_reg) then
          deallocate(guarda_cond)
          allocate(guarda_cond(N_reg))
          guarda_cond = 0
        end if
      end if
    end if

    if (save_vfluxes == 1) then
      if (allocated(vfluxes)) deallocate(vfluxes)
      allocate(vfluxes(4, N_cel))

      if (.not. allocated(guarda_vfluxes)) then
        allocate(guarda_vfluxes(N_reg))
        guarda_vfluxes = 0
      else
        if (size(guarda_vfluxes) < N_reg) then
          deallocate(guarda_vfluxes)
          allocate(guarda_vfluxes(N_reg))
          guarda_vfluxes = 0
        end if
      end if

      if (allocated(mean_vfluxes)) deallocate(mean_vfluxes)
      allocate(mean_vfluxes(4, N_reg))
      mean_vfluxes = 0.0
    end if

    if (save_rc == 1) then
      if (allocated(rc_coef)) deallocate(rc_coef)
      allocate(rc_coef(2, N_cel))
      rc_coef = 0.0
    end if

    !=============================================================
    !             EJECUCIÓN DEL MODELO - BUCLE TEMPORAL
    !=============================================================
    do tiempo = 1, N_reg

      control_cont  = 2
      controlh_cont = 1

      StoAtras = sum(StoOut)

      if (posEvento(tiempo) == 1) then
        Rain = 0.0
      else
        call read_int_basin(ruta_bin, posEvento(tiempo), &
                            N_cel, RainInt, Res)
        Rain = RainInt / 1000.0

        if (separate_rain == 1) then
          call read_int_basin(ruta_binConv, posConv(tiempo), N_cel, Conv, Res)
          call read_int_basin(ruta_binStra, posStra(tiempo), N_cel, Stra, Res)
        end if
      end if

      rain_sum = 0.0
      Acum_rain(1, :) = Acum_rain(1, :) + Rain

      do celda = 1, N_cel

        drenaid = N_cel - drena(1, celda) + 1
        entradas = entradas + Rain(celda)
        rain_sum = rain_sum + Rain(celda)

        if (separate_rain == 1) then
          Co = Conv(celda) / 1000
          St = Stra(celda) / 1000
        end if

        vflux(1) = max(0.0, Rain(celda) - H(1, celda) + StoOut(1, celda))
        StoOut(1, celda) = StoOut(1, celda) + Rain(celda) - vflux(1)

        Evp_loss = min( EvpSerie(tiempo) * vspeed(1, celda) * &
                        (StoOut(1, celda) / H(1, celda))**0.6, &
                        StoOut(1, celda) )

        if (separate_rain == 1) then
          if (StoOut(1, celda) > 0.0) then
            Storage_conv(1, celda) = max(0.0, Storage_conv(1, celda) - Evp_loss * Storage_conv(1, celda) / StoOut(1, celda))
            Storage_stra(1, celda) = max(0.0, Storage_stra(1, celda) - Evp_loss * Storage_stra(1, celda) / StoOut(1, celda))
          end if
        end if

        StoOut(1, celda) = StoOut(1, celda) - Evp_loss

        do i = 1, 3
          vflux(i+1) = min(vflux(i), vspeed(i+1, celda))
          StoOut(i+1, celda) = StoOut(i+1, celda) + vflux(i) - vflux(i+1)
        end do

        if (retorno_aq > 0) then
          Ret_aq = max(0.0, StoOut(4, celda) - H(3, celda))
          StoOut(3, celda) = StoOut(3, celda) + Ret_aq
          StoOut(4, celda) = StoOut(4, celda) - Ret_aq
        end if

        if (retorno_gr > 0) then
          Ret = max(0.0, StoOut(3, celda) - H(2, celda))
          StoOut(2, celda) = StoOut(2, celda) + Ret
          StoOut(3, celda) = StoOut(3, celda) - Ret
          Retorned(1, celda) = Retorned(1, celda) + Ret
          vflux(2) = vflux(2) + Ret
          vflux(3) = vflux(3) - Ret
        end if

        if (save_vfluxes == 1) then
          do i = 1, 4
            vfluxes(i, celda) = vflux(i)
          end do
        end if

        if (save_rc == 1) then
          rc_coef(1, celda) = rc_coef(1, celda) + vflux(2) - vflux(3)
          rc_coef(2, celda) = rc_coef(2, celda) + Rain(celda)
        end if

        if (separate_rain == 1) then
          Storage_conv(1, celda) = Storage_conv(1, celda) + (Rain(celda) - vflux(1)) * Co
          Storage_stra(1, celda) = Storage_stra(1, celda) + (Rain(celda) - vflux(1)) * St

          do i = 1, 3
            Storage_conv(i+1, celda) = Storage_conv(i+1, celda) + (vflux(i) - vflux(i+1)) * Co
            Storage_stra(i+1, celda) = Storage_stra(i+1, celda) + (vflux(i) - vflux(i+1)) * St
          end do

          if (retorno_gr > 0) then
            Storage_conv(2, celda) = Storage_conv(2, celda) + Ret * Co
            Storage_conv(3, celda) = Storage_conv(3, celda) - Ret * Co
            Storage_stra(2, celda) = Storage_stra(2, celda) + Ret * St
            Storage_stra(3, celda) = Storage_stra(3, celda) - Ret * St
          end if
        end if

        salidas = salidas + vflux(4) + Evp_loss

        do i = 1, 3
          select case (speed_type(i))
          case (1)
            hflux(i) = (1 - hill_long(1, celda) / (hspeed(i, celda) * dt + hill_long(1, celda))) * StoOut(i+1, celda)
          case (2)
            call calc_speed( StoOut(i+1, celda) * m3_mmHill(celda), &
                             h_coef(i, celda) * calib(i+4), &
                             h_exp(i, celda), hill_long(1, celda), &
                             hspeed(i, celda), section_area )
            hflux(i) = min(section_area * hspeed(i, celda) * dt / m3_mmHill(celda), StoOut(i+1, celda))

            if (sim_sediments == 1 .and. i == 1) then
              qlin_sed = hflux(i) * m3_mmHill(celda) / (dt * dxp)
              call sed_hillslope( sed_factor, hspeed(1, celda), StoOut(2, celda), &
                                   hill_slope(1, celda), celda, drenaid, unit_type(1, celda) )
            end if
          end select

          if (separate_rain == 1) then
            if (StoOut(i+1, celda) > 0) then
              hflux_c(i) = hflux(i) * Storage_conv(i+1, celda) / StoOut(i+1, celda)
              hflux_s(i) = hflux(i) * Storage_stra(i+1, celda) / StoOut(i+1, celda)
            else
              hflux_c(i) = 0.0
              hflux_s(i) = 0.0
            end if

            Storage_conv(i+1, celda) = Storage_conv(i+1, celda) - hflux_c(i)
            Storage_stra(i+1, celda) = Storage_stra(i+1, celda) - hflux_s(i)
          end if

          StoOut(i+1, celda) = StoOut(i+1, celda) - hflux(i)
        end do
      end do ! Fin loop celdas

      Mean_Rain(1, tiempo) = rain_sum / N_cel

      if (save_storage == 1) then
        if (guarda_cond(tiempo) > 0) then
          call write_float_basin(ruta_storage, StoOut, guarda_cond(tiempo), N_cel, 5)
        end if
      end if

      if (save_vfluxes == 1) then
        mean_vfluxes(:, tiempo) = sum(vfluxes, dim=2) / N_cel
        if (guarda_vfluxes(tiempo) > 0) then
          call write_float_basin(ruta_vfluxes, vfluxes, guarda_vfluxes(tiempo), N_cel, 4)
        end if
      end if

      if (save_speed == 1) then
        call write_float_basin(ruta_speed, hspeed, tiempo, N_cel, 4)
      end if

      if (save_retorno == 1) then
        call write_float_basin(ruta_retorno, Retorned, tiempo, N_cel, 1)
      end if

      if (show_storage == 1) then
        mean_storage(:, tiempo) = sum(StoOut, dim=2) / N_cel
      end if

      if (show_mean_speed == 1) then
        mean_speed(:, tiempo) = sum(hspeed, dim=2) / N_cel
      end if

      if (show_mean_retorno == 1) then
        mean_retorno(tiempo) = sum(Retorned) / N_cel
      end if

      if (show_mean_retorno == 1 .or. save_retorno == 1) then
        Retorned = 0.0
      end if

      balance(tiempo) = sum(StoOut) - StoAtras - entradas + salidas
      entradas = 0.0
      salidas  = 0.0

      if (verbose == 1) then
        print *, real(tiempo) / N_reg
      end if

    end do ! Fin loop tiempo

    if (save_rc == 1) then
      where(rc_coef(1, :) > rc_coef(2, :)) rc_coef(1, :) = rc_coef(2, :)
      call write_float_basin(ruta_rc, rc_coef, 1, N_cel, 2)
    end if

  end subroutine shia_v1

  !-------------------------------------------------------------
  ! Leer vector de reales de archivo binario
  subroutine read_float_basin(ruta, record, N_cel, vect, Res)
    integer, intent(in)           :: record, N_cel
    character(len=*), intent(in)  :: ruta
    real, intent(out)             :: vect(N_cel)
    integer, intent(out)          :: Res

    open(10, file=ruta, form='unformatted', status='old', access='direct', RECL=4*N_cel)
    read(10, rec=record, iostat=Res) vect
    if (Res /= 0) print *, 'Error: Se ha tratado de leer un valor fuera del rango'
    close(10)
  end subroutine read_float_basin

  !-------------------------------------------------------------
  ! Leer matriz de reales (N_col x N_cel) de binario
  subroutine read_float_basin_Ncol(ruta, record, N_cel, N_col, vect, Res)
    integer, intent(in)           :: record, N_cel, N_col
    character(len=*), intent(in)  :: ruta
    real, intent(out)             :: vect(N_col, N_cel)
    integer, intent(out)          :: Res

    open(10, file=ruta, form='unformatted', status='old', access='direct', RECL=4*N_cel*N_col)
    read(10, rec=record, iostat=Res) vect
    if (Res /= 0) print *, 'Error: Se ha tratado de leer un valor fuera del rango'
    close(10)
  end subroutine read_float_basin_Ncol

  !-------------------------------------------------------------
  ! Leer vector de enteros de binario
  subroutine read_int_basin(ruta, record, N_cel, vect, Res)
    integer, intent(in)           :: record, N_cel
    character(len=*), intent(in)  :: ruta
    integer, intent(out)          :: vect(N_cel)
    integer, intent(out)          :: Res

    open(10, file=ruta, form='unformatted', status='old', access='direct', RECL=4*N_cel)
    read(10, rec=record, iostat=Res) vect
    if (Res /= 0) print *, 'Error: Se ha tratado de leer un valor fuera del rango'
    close(10)
  end subroutine read_int_basin

  !-------------------------------------------------------------
  ! Escribir matriz de reales en binario (N_col x N_cel)
  subroutine write_float_basin(ruta, vect, record, N_cel, N_col)
    integer, intent(in)          :: record, N_cel, N_col
    character(len=*), intent(in) :: ruta
    real, intent(in)             :: vect(N_col, N_cel)
    character(len=10)            :: estado

    estado = 'old'
    if (record == 1) estado = 'replace'

    if (record > 0) then
      open(10, file=ruta, form='unformatted', status=estado, access='direct', RECL=4*N_col*N_cel)
      write(10, rec=record) vect
      close(10)
    end if
  end subroutine write_float_basin

  !-------------------------------------------------------------
  ! Escribir matriz de enteros en binario (N_col x N_cel)
  subroutine write_int_basin(ruta, vect, record, N_cel, N_col)
    integer, intent(in)          :: record, N_cel, N_col
    character(len=*), intent(in) :: ruta
    integer, intent(in)          :: vect(N_col, N_cel)
    character(len=10)            :: estado

    estado = 'old'
    if (record == 1) estado = 'replace'

    open(10, file=ruta, form='unformatted', status=estado, access='direct', RECL=4*N_col*N_cel)
    write(10, rec=record) vect
    close(10)
  end subroutine write_int_basin
  !=============================================================
  ! SUBRUTINAS DE LECTURA E INTERPOLACIÓN DE LLUVIA
  !=============================================================

  !-------------------------------------------------------------
  ! Lectura de tabla de eventos de lluvia en ASCII (general)
  subroutine rain_read_ascii_table(ruta, Nintervals)
    character(len=*), intent(in) :: ruta
    integer, intent(in)          :: Nintervals
    character(len=20) :: oe
    integer          :: i, cont, Ntotal

    if (allocated(idEvento))   deallocate(idEvento)
    if (allocated(posEvento))  deallocate(posEvento)
    allocate(idEvento(Nintervals), posEvento(Nintervals))

    open(unit=10, file=ruta, status='old', action='read')
    do i = 1, 3
      read(10,*) oe, oe, oe, Ntotal
    end do
    read(10,*)
    read(10,*)
    read(10,*)

    if (Nintervals + rain_first_point <= Ntotal + 1) then
      if (rain_first_point > 1) then
        do i = 1, rain_first_point
          read(10,*)
        end do
      end if

      cont = 1
      do i = 1, Nintervals
        read(10,*) idEvento(cont), posEvento(cont), oe
        cont = cont + 1
      end do
    end if
    close(10)
  end subroutine rain_read_ascii_table

  !-------------------------------------------------------------
  ! Lectura tablas separadas de lluvia convectiva / estratiforme
  subroutine rain_read_ascii_table_separate(rutaConv, rutaStra, Nintervals)
    character(len=*), intent(in) :: rutaConv, rutaStra
    integer, intent(in)          :: Nintervals
    character(len=20) :: oe
    integer          :: i, cont, Ntotal

    if (allocated(posConv))  deallocate(posConv)
    if (allocated(posStra))  deallocate(posStra)
    allocate(posConv(Nintervals), posStra(Nintervals))

    ! Convectiva
    open(unit=10, file=rutaConv, status='old', action='read')
    do i = 1, 3
      read(10,*) oe, oe, oe, Ntotal
    end do
    read(10,*)
    read(10,*)
    read(10,*)
    if (Nintervals + rain_first_point <= Ntotal) then
      if (rain_first_point > 1) then
        do i = 1, rain_first_point
          read(10,*)
        end do
      end if
      cont = 1
      do i = 1, Nintervals
        read(10,*) oe, posConv(cont), oe
        cont = cont + 1
      end do
    end if
    close(10)

    ! Estratiforme
    open(unit=10, file=rutaStra, status='old', action='read')
    do i = 1, 3
      read(10,*) oe, oe, oe, Ntotal
    end do
    read(10,*)
    read(10,*)
    read(10,*)
    if (Nintervals + rain_first_point <= Ntotal) then
      if (rain_first_point > 1) then
        do i = 1, rain_first_point
          read(10,*)
        end do
      end if
      cont = 1
      do i = 1, Nintervals
        read(10,*) oe, posStra(cont), oe
        cont = cont + 1
      end do
    end if
    close(10)
  end subroutine rain_read_ascii_table_separate

  !-------------------------------------------------------------
  ! Prepara mapeo celda-triángulo para método TIN
  subroutine rain_pre_mit(tin_perte, xy_basin, TIN, coord, nceldas, ntin, ncoord)
    integer, intent(in) :: ntin, nceldas, TIN(3, ntin), ncoord
    real,    intent(in) :: coord(2, ncoord), xy_basin(2, nceldas)
    integer, intent(out):: tin_perte(1, nceldas)

    integer :: celda, triangulo
    real    :: p0x, p1x, p2x, p0y, p1y, p2y, px, py, Area, s, t

    do celda = 1, nceldas
      px = xy_basin(1, celda)
      py = xy_basin(2, celda)
      do triangulo = 1, ntin
        p0x = coord(1, TIN(1, triangulo))
        p1x = coord(1, TIN(2, triangulo))
        p2x = coord(1, TIN(3, triangulo))

        p0y = coord(2, TIN(1, triangulo))
        p1y = coord(2, TIN(2, triangulo))
        p2y = coord(2, TIN(3, triangulo))

        Area = abs(0.5 * (-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y))
        s = 1.0 / (2.0*Area) * (p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py)
        t = 1.0 / (2.0*Area) * (p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py)

        if (s > 0.0 .and. t > 0.0) then
          tin_perte(1, celda) = triangulo
        end if
      end do
    end do
  end subroutine rain_pre_mit

  !-------------------------------------------------------------
  ! Interpolación de lluvia método TIN
  subroutine rain_mit( xy_basin, coord, rain, tin, tin_perte, nceldas, nhills, ncoord, &
                       ntin, nreg, ruta, umbral, meanRain, posIds, maskVector )
    integer, intent(in)  :: nceldas, nreg, ncoord, nhills
    integer, intent(in)  :: ntin, tin(3, ntin), tin_perte(1, nceldas), maskVector(nceldas)
    real,    intent(in)  :: xy_basin(2, nceldas), coord(2, ncoord), rain(ncoord, nreg), umbral
    character(len=*), intent(in) :: ruta
    real,    intent(out) :: meanRain(nreg)
    integer, intent(out) :: posIds(nreg)

    real    :: ax(nceldas), ay(nceldas), bx(nceldas), by(nceldas)
    real    :: cx(nceldas), cy(nceldas), det1, det2, det3, det4, coef1(nceldas)
    real    :: Cel_x, Cel_y, az, bz, cz, coef2, campo(nceldas), campoHill(nhills)
    integer :: campoInt(nceldas), mascara(nceldas), campoIntHill(nhills)
    integer :: Cel_pert, celdas, tiempo, celdas_hills, cont

    if (sum(maskVector) == nceldas) then
      celdas_hills = 1
    else if (sum(maskVector) == nhills) then
      celdas_hills = 2
    end if

    campoInt = 0
    mascara  = 1

    if (celdas_hills == 2) then
      call basin_subbasin_map2subbasin(maskVector, campo, campoHill, nhills, nceldas, mascara, celdas_hills)
      campoIntHill = campoHill
      call write_int_basin(ruta, campoIntHill, 1, nhills, 1)
    else
      call write_int_basin(ruta, campoInt, 1, nceldas, 1)
    end if

    ax = coord(1, tin(1, tin_perte(1, :)))
    ay = coord(2, tin(1, tin_perte(1, :)))
    bx = coord(1, tin(2, tin_perte(1, :)))
    by = coord(2, tin(2, tin_perte(1, :)))
    cx = coord(1, tin(3, tin_perte(1, :)))
    cy = coord(2, tin(3, tin_perte(1, :)))
    coef1 = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)

    cont = 2
    do tiempo = 1, nreg
      do celdas = 1, nceldas
        Cel_pert = tin_perte(1, celdas)
        Cel_x    = xy_basin(1, celdas)
        Cel_y    = xy_basin(2, celdas)

        az = max(rain(tin(1, Cel_pert), tiempo), 0.0)
        bz = max(rain(tin(2, Cel_pert), tiempo), 0.0)
        cz = max(rain(tin(3, Cel_pert), tiempo), 0.0)

        det1 = (cx(celdas) - ax(celdas)) * (Cel_y - ay(celdas))
        det2 = (by(celdas) - ay(celdas)) * (Cel_x - ax(celdas))
        det3 = (cy(celdas) - ay(celdas)) * (Cel_x - ax(celdas))
        det4 = (Cel_y - ay(celdas)) * (bx(celdas) - ax(celdas))

        coef2 = det1*(bz - az) + det2*(cz - az) - det3*(bz - az) - det4*(cz - az)
        campo(celdas) = max(az - coef2 / coef1(celdas), 0.0)
      end do

      if (sum(campo) > umbral .and. count(campo > umbral) > 0) then
        meanRain(tiempo) = sum(campo) / count(campo > 0)
        if (celdas_hills == 1) then
          campoInt = campo * 1000
          call write_int_basin(ruta, campoInt, cont, nceldas, 1)
        else if (celdas_hills == 2) then
          call basin_subbasin_map2subbasin(maskVector, campo, campoHill, nhills, nceldas, mascara, celdas_hills)
          campoIntHill = campoHill * 1000
          call write_int_basin(ruta, campoIntHill, cont, nhills, 1)
        end if
        posIds(tiempo) = cont
        cont = cont + 1
      else
        meanRain(tiempo) = 0.0
        posIds(tiempo) = 1
      end if
    end do
  end subroutine rain_mit

  !-------------------------------------------------------------
  ! Interpolación de lluvia método IDW
  subroutine rain_idw( xy_basin, coord, rain, pp, nceldas, ncoord, nreg, nhills, ruta, umbral, &
                       meanRain, posIds, maskVector )
    integer, intent(in)  :: nceldas, ncoord, nreg, nhills
    integer, intent(in)  :: maskVector(nceldas)
    real,    intent(in)  :: xy_basin(2, nceldas), coord(2, ncoord)
    real,    intent(in)  :: rain(ncoord, nreg), pp, umbral
    character(len=*), intent(in) :: ruta
    real,    intent(out) :: meanRain(nreg)
    integer, intent(out) :: posIds(nreg)

    real    :: W(ncoord, nceldas), Wr, campo(nceldas), valor, campoHill(nhills)
    integer :: campoInt(nceldas), campoIntHill(nhills), mascara(nceldas)
    integer :: tiempo, celda, i, cont, celdas_hills

    if (sum(maskVector) == nceldas) then
      celdas_hills = 1
    else
      celdas_hills = 2
    end if

    campoInt = 0
    mascara  = 1

    if (celdas_hills == 2) then
      call basin_subbasin_map2subbasin(maskVector, campo, campoHill, nhills, nceldas, mascara, celdas_hills)
      campoIntHill = campoHill
      call write_int_basin(ruta, campoIntHill, 1, nhills, 1)
    else
      call write_int_basin(ruta, campoInt, 1, nceldas, 1)
    end if

    do i = 1, ncoord
      W(i, :) = 1.0 / (sqrt((coord(1, i) - xy_basin(1, :))**2 + (coord(2, i) - xy_basin(2, :))**2))**pp
    end do

    cont = 2
    do tiempo = 1, nreg
      do celda = 1, nceldas
        Wr = sum(W(:, celda) * rain(:, tiempo), mask = rain(:, tiempo) > 0.0)
        valor = max(Wr / sum(W(:, celda), mask = rain(:, tiempo) >= 0.0), 0.0)
        if (valor == valor - 1) then
          campo(celda) = 0.0
        else
          campo(celda) = valor
        end if
      end do

      if (sum(campo) > umbral .and. count(campo > umbral) > 0) then
        meanRain(tiempo) = sum(campo) / count(campo > 0)
        if (celdas_hills == 1) then
          campoInt = campo * 1000
          call write_int_basin(ruta, campoInt, cont, nceldas, 1)
        else
          call basin_subbasin_map2subbasin(maskVector, campo, campoHill, nhills, nceldas, mascara, celdas_hills)
          campoIntHill = campoHill * 1000
          call write_int_basin(ruta, campoIntHill, cont, nhills, 1)
        end if
        posIds(tiempo) = cont
        cont = cont + 1
      else
        meanRain(tiempo) = 0.0
        posIds(tiempo) = 1
      end if
    end do
  end subroutine rain_idw

  !=============================================================
  ! SOLUCIÓN DE ONDA CINEMÁTICA
  !=============================================================
  subroutine calc_speed(sm, coef, expo, elem_long, speed, area)
    real, intent(in)    :: sm, coef, expo, elem_long
    real, intent(inout) :: speed
    real, intent(out)   :: area
    real :: new_speed
    integer :: i

    do i = 1, calc_niter
      area       = sm / (elem_long + speed * dt)   ! [m2]
      new_speed  = coef * (area**expo)             ! [m/s]
      speed      = (2*new_speed + speed) / 3       ! Promedio suavizado
    end do
  end subroutine calc_speed

  !=============================================================
  ! MÓDULO DE SEDIMENTOS
  !=============================================================

  subroutine sed_allocate(N_cel)
    integer, intent(in) :: N_cel

    if (.not. allocated(Vd)) then
      allocate(Vd(3, N_cel))
      Vd = 0.0
    end if
    if (.not. allocated(Vs))   allocate(Vs(3, N_cel))
    if (.not. allocated(VolERO)) allocate(VolERO(N_cel))
    if (.not. allocated(VolDEPo)) allocate(VolDEPo(N_cel))
    if (.not. allocated(VSc))  allocate(VSc(3, N_cel))
    if (.not. allocated(Vdc))  allocate(Vdc(3, N_cel))

    Vs      = 0.0
    VSc     = 0.0
    Vdc     = 0.0
    VolERO  = 0.0
    VolDEPo = 0.0
    EROt    = 0.0
    DEPt    = 0.0
  end subroutine sed_allocate

  subroutine sed_hillslope(alfa, v2, S2, So, celda, drena_id, tipo)
    real, intent(in) :: v2, So, alfa, S2
    integer, intent(in) :: celda, drena_id, tipo

    integer :: i
    real :: qsSUStot, qsEROtot, qsSUS(3), qsERO(3)
    real :: SUStot, DEPtot, totXSScap, REScap
    real :: Qskr, cap, Adv, Te(3)

    DEP = 0.0
    qsSUStot = 0.0
    qsSUS    = 0.0
    qsERO    = 0.0
    SUStot   = 0.0
    DEPtot   = 0.0

    Qskr = alfa * (So ** 1.664) * (qlin_sed ** 2.035) * Krus(1, celda) * Crus(1, celda) * Prus(1, celda) * dxp * dt

    do i = 1, 3
      if (S2 / 1000 > wi(i) * dt) then
        Te(i) = wi(i) * dt * 1000 / S2
      else
        Te(i) = 1.0
      end if

      Vd(i, celda) = Vd(i, celda) + Te(i) * Vs(i, celda)
      Vs(i, celda) = Vs(i, celda) - Te(i) * Vs(i, celda)
      DEP(i)       = Te(i) * Vs(i, celda)

      SUStot = SUStot + Vs(i, celda)
      DEPtot = DEPtot + Vd(i, celda)
    end do

    do i = 1, 3
      if (Vs(i, celda) > 0.0) then
        if (Qskr < SUStot) then
          cap = Qskr * Vs(i, celda) / SUStot
          Adv = Vs(i, celda) * v2 * dt / (dxp + v2 * dt)
          qsSUS(i) = max(cap, Adv)
        else
          qsSUS(i) = Vs(i, celda)
        end if
      end if
      qsSUStot = qsSUStot + qsSUS(i)
      Vs(i, celda) = Vs(i, celda) - qsSUS(i)

      if (tipo == 1) then
        Vs(i, drena_id) = Vs(i, drena_id) + qsSUS(i)
      else
        VSc(i, celda)   = VSc(i, celda) + qsSUS(i)
      end if
    end do

    totXSScap = max(0.0, Qskr - qsSUStot)

    if (totXSScap > 0.0 .and. DEPtot > 0.0) then
      do i = 1, 3
        if (Vd(i, celda) > 0.0) then
          qsERO(i) = PArLiAc(i, celda) * totXSScap / 100.0
          if (tipo == 1) then
            Vs(i, drena_id) = Vs(i, drena_id) + qsERO(i)
          else
            VSc(i, celda)   = VSc(i, celda) + qsERO(i)
          end if
        end if
      end do
    end if

    EROt = EROt + qsERO
    VolERO(celda)  = VolERO(celda) + sum(qsERO)
    VolDEPo(celda) = VolDEPo(celda) + sum(DEP)
  end subroutine sed_hillslope

  subroutine sed_channel(S5, v5, Q5, So, area_sec, celda, drena_id, VolSal)
    real, intent(in) :: S5, v5, So, area_sec, Q5
    integer, intent(in) :: celda, drena_id
    real, intent(out) :: VolSal(3)

    integer :: i
    real, parameter :: grav=9.8
    real, parameter :: Gsed=2.65
    real :: Cw(3), Te(3), qsSUS(3), qsBM(3), XSScap, VolEH, AdvF, supply, Rh

    DEP = 0.0
    qsSUS = 0.0
    qsBM  = 0.0

    if (S5 > 0.0) then
      Rh = area_sec * 1000 / S5
    else
      Rh = 0.0
    end if

    do i = 1, 3
      Cw(i) = 0.05 * (Gsed / (Gsed - 1)) * (v5 * So) / sqrt((Gsed - 1) * grav * diametro(i) / 1000.0) * &
              sqrt((Rh * So) / ((Gsed - 1) * (diametro(i) / 1000.0)))

      if (S5 / 1000 > wi(i) * dt) then
        Te(i) = wi(i) * dt * 1000 / S5
      else
        Te(i) = 1.0
      end if

      Vdc(i, celda) = Vdc(i, celda) + Te(i) * VSc(i, celda)
      VSc(i, celda) = VSc(i, celda) - Te(i) * VSc(i, celda)
      DEP(i)        = DEP(i) + Te(i) * VSc(i, celda)
    end do

    do i = 1, 3
      supply = VSc(i, celda) + Vdc(i, celda)
      qsSUS(i) = 0.0
      if (supply > 0.0) then
        VolEH = Q5 * Cw(i) * dt / 2.65
        AdvF  = min(1.0, v5 * dt / (dxp + v5 * dt))
        qsSUS(i) = VSc(i, celda) * AdvF
        XSScap   = max(0.0, VolEH - qsSUS(i))
        qsBM(i)  = min(XSScap, Vdc(i, celda) * AdvF)
        DEP(i)   = max(DEP(i) - qsBM(i), 0.0)
        VSc(i, celda) = VSc(i, celda) - qsSUS(i)
        Vdc(i, celda) = Vdc(i, celda) - qsBM(i)
        if (drena_id /= 0) VSc(i, drena_id) = VSc(i, drena_id) + qsSUS(i) + qsBM(i)
        VolSal(i) = (qsSUS(i) + qsBM(i)) / dt
      end if
    end do

    VolDEPo(celda) = VolDEPo(celda) + sum(DEP)
  end subroutine sed_channel

  !=============================================================
  ! MÓDULO DE DESLIZAMIENTOS
  !=============================================================

  subroutine slide_allocate(N_cel, N_reg)
    integer, intent(in) :: N_cel, N_reg

    if (allocated(sl_Zmin))           deallocate(sl_Zmin)
    if (allocated(sl_Zmax))           deallocate(sl_Zmax)
    if (allocated(sl_Zcrit))          deallocate(sl_Zcrit)
    if (allocated(sl_Bo))             deallocate(sl_Bo)
    if (allocated(sl_SlideOcurrence)) deallocate(sl_SlideOcurrence)
    if (allocated(sl_RiskVector))     deallocate(sl_RiskVector)
    if (allocated(sl_SlideNcellTime)) deallocate(sl_SlideNcellTime)
    if (allocated(sl_SlideAcumulate)) deallocate(sl_SlideAcumulate)

    allocate(sl_Zmin(1, N_cel), sl_Zmax(1, N_cel), sl_Zcrit(1, N_cel), sl_Bo(1, N_cel))
    allocate(sl_SlideOcurrence(1, N_cel), sl_RiskVector(1, N_cel))
    allocate(sl_SlideNcellTime(N_reg), sl_SlideAcumulate(1, N_cel))

    sl_Zmin = sl_Cohesion / ((sl_GammaW * (COS(sl_RadSlope))**2 * TAN(sl_FrictionAngle)) + &
                             (sl_GammaS * (COS(sl_RadSlope))**2 * (TAN(sl_RadSlope) - TAN(sl_FrictionAngle))))

    sl_Zcrit = (sl_GammaS / sl_GammaW) * sl_Zs * (1.0 - (TAN(sl_RadSlope) / TAN(sl_FrictionAngle))) + &
               (sl_Cohesion / (sl_GammaW * (COS(sl_RadSlope))**2 * TAN(sl_FrictionAngle)))

    sl_Zmax = sl_Cohesion / ((sl_GammaS * (COS(sl_RadSlope))**2) * (TAN(sl_RadSlope) - TAN(sl_FrictionAngle)))

    sl_Bo   = ATAN(-TAN(sl_FrictionAngle * (sl_GammaW - sl_GammaS) / sl_GammaS))

    sl_RiskVector = 2
    where (unit_type == 3) sl_RiskVector = 0
    where (sl_RadSlope < sl_Bo) sl_RiskVector = 0
    where (sl_Zs < sl_Zmin .and. sl_RadSlope < sl_Bo) sl_RiskVector = 0
    where (sl_Zs > sl_Zmax .and. sl_RadSlope > sl_Bo) sl_RiskVector = 1

    sl_SlideOcurrence = 0
    sl_SlideNcellTime = 0
    sl_SlideAcumulate = 0
  end subroutine slide_allocate

  subroutine slide_ocurrence(N_cel, timeStep, cell, StorageT3, MaxStoT3)
    integer, intent(in) :: cell, N_cel, timeStep
    real,    intent(in) :: StorageT3, MaxStoT3
    real :: Zw, Num, Den

    if (sl_RiskVector(1, cell) == 2 .and. sl_SlideOcurrence(1, cell) == 0) then
      Zw = (StorageT3 * sl_Zs(1, cell)) / MaxStoT3
      if (Zw >= sl_Zcrit(1, cell)) then
        sl_SlideOcurrence(1, cell) = 1
        sl_SlideAcumulate(1, cell) = 1
        sl_SlideNcellTime(timeStep) = sl_SlideNcellTime(timeStep) + 1
        if (sl_GullieNoGullie == 1) call slide_hill2gullie(N_cel, cell, timeStep)
      else
        Num = sl_Cohesion(1, cell) + (sl_GammaS(1, cell) * sl_Zs(1, cell) - Zw * sl_GammaW) * &
              (cos(sl_RadSlope(1, cell)))**2 * TAN(sl_FrictionAngle(1, cell))
        Den = sl_GammaS(1, cell) * sl_Zs(1, cell) * sin(sl_RadSlope(1, cell)) * cos(sl_RadSlope(1, cell))
        if (Num / Den <= sl_FS) then
          sl_SlideOcurrence(1, cell) = 1
          sl_SlideAcumulate(1, cell) = 1
          sl_SlideNcellTime(timeStep) = sl_SlideNcellTime(timeStep) + 1
          if (sl_GullieNoGullie == 1) call slide_hill2gullie(N_cel, cell, timeStep)
        end if
      end if
    end if
  end subroutine slide_ocurrence

  subroutine slide_hill2gullie(N_cel, cell, timeSt)
    integer, intent(in) :: cell, N_cel, timeSt
    logical :: Flag
    integer :: localCell, direction

    if (unit_type(1, cell) <= 2) then
      Flag = .true.
      localCell = cell
      do while (Flag)
        direction = N_cel - drena(1, localCell) + 1
        if (unit_type(1, direction) > unit_type(1, localCell)) then
          Flag = .false.
        else
          sl_SlideOcurrence(1, direction) = 3
          sl_SlideAcumulate(1, direction) = sl_SlideAcumulate(1, direction) + 1
          sl_SlideNcellTime(timeSt) = sl_SlideNcellTime(timeSt) + 1
          localCell = direction
        end if
      end do
    end if
  end subroutine slide_hill2gullie

  !=============================================================
  ! MÓDULO DE INUNDACIONES
  !=============================================================

  subroutine flood_allocate(N_cel)
    integer, intent(in) :: N_cel

    if (.not. allocated(flood_Q))         allocate(flood_Q(1, N_cel))
    if (.not. allocated(flood_Qsed))      allocate(flood_Qsed(1, N_cel))
    if (.not. allocated(flood_speed))     allocate(flood_speed(1, N_cel))
    if (.not. allocated(flood_h))         allocate(flood_h(1, N_cel))
    if (.not. allocated(flood_ufr))       allocate(flood_ufr(1, N_cel))
    if (.not. allocated(flood_Cr))        allocate(flood_Cr(1, N_cel))
    if (.not. allocated(flood_eval))      allocate(flood_eval(1, N_cel))
    if (.not. allocated(flood_flood))     allocate(flood_flood(1, N_cel))
    if (.not. allocated(flood_loc_hand))  allocate(flood_loc_hand(1, N_cel))
    if (.not. allocated(flood_slope))     allocate(flood_slope(1, N_cel))

    flood_sec_tam = size(flood_sections(:, 1))

    if (.not. allocated(flood_profundidad)) &
      allocate(flood_profundidad(int(flood_sec_tam), N_cel))

    flood_eval  = 0
    flood_flood = 0
  end subroutine flood_allocate

  subroutine flood_params(celda)
    integer, intent(in) :: celda
    real :: d50, cr, dw, dsed, cmax

    d50  = flood_d50(1, celda)
    dw   = flood_dw
    dsed = flood_dsed
    cmax = flood_Cmax

    flood_h(1, celda) = flood_Q(1, celda) / (flood_speed(1, celda) * flood_w(1, celda))
    if (flood_h(1, celda) > flood_hmax) flood_h(1, celda) = flood_hmax

    flood_ufr(1, celda) = flood_speed(1, celda) / ((5.75) * log(flood_h(1, celda) / flood_d50(1, celda)) + 6.25)
    cr = flood_Cmax * (0.06 * flood_h(1, celda)) ** (0.2 / flood_ufr(1, celda))
    flood_Cr(1, celda) = cr

    flood_rdf = (1. / d50) * ((9.81 / 0.0128) * (cr + (1 - cr) * (dw / dsed)))**0.5 * ((cmax / cr)**(1. / 3) - 1.0)
    flood_Qsed(1, celda) = flood_Q(1, celda) * (1.0 + (flood_Cr(1, celda) / (1 - flood_Cr(1, celda))))
  end subroutine flood_params

  subroutine flood_find_hdiff(celda)
    integer, intent(in) :: celda

    flood_loc_hand(1, :) = 99999
    where(flood_aquien == celda) flood_loc_hand = flood_hand
    call QsortC(flood_loc_hand(1, :))
  end subroutine flood_find_hdiff

  subroutine flood_debris_flow(celda, areas, dif)
    integer, intent(in)  :: celda
    real, intent(out)    :: areas, dif

    integer :: i
    real :: Qpar, Qp, h, area

    h = 0.0
    area = 0.0
    areas = 0.0
    Qp = 0.0
    i = 1

    do while (flood_loc_hand(1, i) /= 0 .and. flood_loc_hand(1, i) /= 9999 .and. &
              Qp < flood_Qsed(1, celda) .and. i < flood_max_iter)

      dif   = flood_loc_hand(1, i)
      areas = areas + area
      area  = (dif - h) * i * dx
      Qpar  = ((2. / 5) * flood_rdf * (dif - h) ** (3. / 2) * flood_slope(1, celda) * (1. / 2)) * area
      h     = dif
      Qp    = Qp + Qpar
      i     = i + 1
    end do
  end subroutine flood_debris_flow

  subroutine flood_debris_flow2(celda, area, dif)
    integer, intent(in)  :: celda
    real, intent(out)    :: area, dif

    integer :: i
    real :: Qp, fondo, linea
    real :: diferencias(int(flood_sec_tam)), positions(int(flood_sec_tam))

    area = 0.0
    Qp   = 0.0
    i    = 1

    fondo = flood_sections((floor(flood_sec_tam) / 2) + 1, celda)

    do while (Qp < flood_Qsed(1, celda) .and. i * flood_step < 50)
      dif = i * flood_step
      linea = fondo + dif
      diferencias = linea - flood_sections(:, celda)
      area = sum(diferencias, mask = diferencias > 0.0) * dxp

      Qp = ((2. / 5) * flood_rdf * (i * flood_step) ** (3. / 2) * flood_slope(1, celda) * (1. / 2)) * area
      i = i + 1
    end do

    flood_profundidad(:, celda) = diferencias
    positions = 999999
    where(diferencias > 0.0) positions = flood_sec_cells(:, celda)
    call QsortC(positions)

    do i = 1, count(diferencias > 0.0)
      flood_flood(1, int(positions(i))) = 1
    end do
    flood_flood(1, celda) = 2
  end subroutine flood_debris_flow2

  !=============================================================
  ! UTILIDADES VARIAS
  !=============================================================

  subroutine basin_subbasin_map2subbasin(sub_pert, basin_var, subbasin_sum, &
                                         n_nodos, nceldas, cauce, sum_mean)
    integer, intent(in)          :: n_nodos, nceldas, sum_mean
    integer, intent(in)          :: sub_pert(nceldas)
    real,    intent(in)          :: basin_var(nceldas)
    integer, intent(in), optional:: cauce(nceldas)
    real,    intent(out)         :: subbasin_sum(n_nodos)

    integer :: i, posi, cont_valores
    real    :: suma_valores

    do i = 1, n_nodos
      posi = n_nodos - i + 1
      if (present(cauce)) then
        suma_valores = sum(basin_var, mask = sub_pert == posi .and. cauce == 1)
        cont_valores = count(sub_pert == posi .and. cauce == 1)
      else
        suma_valores = sum(basin_var, mask = sub_pert == posi)
        cont_valores = count(sub_pert == posi)
      end if

      if (cont_valores > 0) then
        if (sum_mean == 1) then
          subbasin_sum(i) = suma_valores / cont_valores
        else if (sum_mean == 2) then
          subbasin_sum(i) = suma_valores
        end if
      else
        subbasin_sum(i) = 0.0
      end if
    end do
  end subroutine basin_subbasin_map2subbasin

  recursive subroutine QsortC(A)
    real, intent(inout), dimension(:) :: A
    integer :: iq
    if (size(A) > 1) then
      call Partition(A, iq)
      call QsortC(A(:iq-1))
      call QsortC(A(iq:))
    end if
  end subroutine QsortC

  subroutine Partition(A, marker)
    real, intent(inout), dimension(:) :: A
    integer, intent(out)              :: marker
    integer :: i, j
    real    :: temp, x

    x = A(1)
    i = 0
    j = size(A) + 1

    do
      j = j - 1
      do
        if (A(j) <= x) exit
        j = j - 1
      end do

      i = i + 1
      do
        if (A(i) >= x) exit
        i = i + 1
      end do

      if (i < j) then
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
      else if (i == j) then
        marker = i + 1
        return
      else
        marker = i
        return
      end if
    end do
  end subroutine Partition

end module models

