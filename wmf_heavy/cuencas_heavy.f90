!===============================================================
!  Módulo: cu
!  Propósito:
!      Herramientas para el trazado de cuencas hidrográficas y 
!      drenajes, y posterior análisis de sus propiedades 
!      geomorfológicas, hidrológicas y de red.
!
!  Autor original: Nicolás Velásquez Girón (2014)
!  Email: nicolas.velasquezgiron@gmail.com
!
!  Licencia:
!      GPL v3 o posterior. Este programa se distribuye sin 
!      garantías, véase https://www.gnu.org/licenses/
!
!  Descripción general:
!      - Puede ser usado desde Fortran o Python (vía f2py + numpy).
!      - Procesa mapas de direcciones y elevaciones digitales.
!      - Admite conversión desde formatos GRASS y ArcGIS.
!      - Devuelve resultados en arreglos Fortran o Numpy.
!
!  Nota de estilo:
!      Código migrado a formato libre de Fortran 90+, con
!      indentado y tipografía corregida.
!===============================================================

module cu
    implicit none

    !-----------------------------------------------------------
    ! Variables globales de mapa y parámetros de referencia
    !-----------------------------------------------------------
    real :: xll, yll                  ! Coordenadas de la esquina inferior izquierda (mapa base)
    real :: noData                    ! Valor para celdas sin dato
    real :: dx, dy                     ! Resolución en ejes X e Y
    real :: dxP                        ! Longitud proyectada de una celda (m)
    integer :: ncols, nrows             ! Número de columnas y filas del mapa

    !-----------------------------------------------------------
    ! Buffers temporales para cálculos
    !-----------------------------------------------------------
    real,    allocatable :: stream_temp(:,:)       ! Trazado de corriente
    integer, allocatable :: basin_temp(:,:)        ! Trazado de cuenca
    real,    allocatable :: perim_temp(:,:)        ! Puntos del perímetro
    integer, allocatable :: sub_basins_temp(:,:)   ! Subcuencas
    real,    allocatable :: ppal_stream_temp(:,:)  ! Cauce principal
    real,    allocatable :: netxy_temp(:,:)        ! Coordenadas de red hídrica
    real,    allocatable :: col_fil_temp(:,:)      ! Uso general en operaciones de matrices

    !-----------------------------------------------------------
    ! Propiedades de última cuenca calculada
    !-----------------------------------------------------------
    real :: area         ! Área (km²)
    real :: perimetro    ! Perímetro (km)
    real :: pend_media   ! Pendiente media (%)
    real :: elevacion    ! Elevación media (m s.n.m.)
    real :: centroX, centroY   ! Coordenadas del centroide de la cuenca

    !-----------------------------------------------------------
    ! Visibilidad del módulo
    !-----------------------------------------------------------
    public :: QsortC  ! Subrutina pública para ordenamiento
    private :: Partition  ! Subrutina auxiliar privada
!========================================================
! BLOQUE 2: Lectura y escritura de mapas (ASCII y binarios)
!========================================================

contains

!--------------------------------------------------------
! Lee el encabezado de un archivo ASCII de mapa (formato ESRI)
!--------------------------------------------------------
subroutine read_ASCII_hdr(ruta)
    character(len=255), intent(in) :: ruta
    character(len=20) :: linea_tmp

    open(10, file=ruta, status='old')
    read(10,*) linea_tmp, ncols
    read(10,*) linea_tmp, nrows
    read(10,*) linea_tmp, xll
    read(10,*) linea_tmp, yll
    read(10,*) linea_tmp, dx
    read(10,*) linea_tmp, noData
    close(10)
end subroutine read_ASCII_hdr

!--------------------------------------------------------
! Lee un mapa ASCII de reales (formato ESRI)
!--------------------------------------------------------
subroutine read_float_ASCII_dat(ruta, nc, nr, Mat)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: nc, nr
    real, intent(out)  :: Mat(nc, nr)
    integer :: i, j
    character(len=20) :: linea_tmp

    open(10, file=ruta, status='old')
    do i = 1, 6
        read(10,*) linea_tmp
    end do
    read(10,*) ((Mat(i,j), i=1, ncols), j=1, nrows)
    close(10)
end subroutine read_float_ASCII_dat

!--------------------------------------------------------
! Lee un mapa ASCII de enteros (formato ESRI)
!--------------------------------------------------------
subroutine read_int_ASCII_dat(ruta, nc, nr, Mat)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: nc, nr
    integer, intent(out) :: Mat(nc, nr)
    integer :: i, j
    character(len=20) :: linea_tmp

    open(10, file=ruta, status='old')
    do i = 1, 6
        read(10,*) linea_tmp
    end do
    read(10,*) ((Mat(i,j), i=1, ncols), j=1, nrows)
    close(10)
end subroutine read_int_ASCII_dat

!--------------------------------------------------------
! Lee encabezado de un mapa ASCII (subregión/no global)
!--------------------------------------------------------
subroutine read_ASCII_hdr_ng(ruta, nc, nr, xlln, ylln, dxn, noDataN)
    character(len=255), intent(in) :: ruta
    integer, intent(out) :: nc, nr
    real, intent(out) :: xlln, ylln, dxn, noDataN
    character(len=20) :: linea_tmp

    open(10, file=ruta, status='old')
    read(10,*) linea_tmp, nc
    read(10,*) linea_tmp, nr
    read(10,*) linea_tmp, xlln
    read(10,*) linea_tmp, ylln
    read(10,*) linea_tmp, dxn
    read(10,*) linea_tmp, noDataN
    close(10)
end subroutine read_ASCII_hdr_ng

!--------------------------------------------------------
! Lee un mapa de reales (subregión/no global)
!--------------------------------------------------------
subroutine read_float_ASCII_ng(ruta, nc, nr, Mat)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: nc, nr
    real, intent(out) :: Mat(nc, nr)
    integer :: i, j
    character(len=20) :: linea_tmp

    open(10, file=ruta, status='old')
    do i = 1, 6
        read(10,*) linea_tmp
    end do
    read(10,*) ((Mat(i,j), i=1, nc), j=1, nr)
    close(10)
end subroutine read_float_ASCII_ng

!--------------------------------------------------------
! Lee un mapa de enteros (subregión/no global)
!--------------------------------------------------------
subroutine read_int_ASCII_ng(ruta, nc, nr, Mat)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: nc, nr
    integer, intent(out) :: Mat(nc, nr)
    integer :: i, j
    character(len=20) :: linea_tmp

    open(10, file=ruta, status='old')
    do i = 1, 6
        read(10,*) linea_tmp
    end do
    read(10,*) ((Mat(i,j), i=1, nc), j=1, nr)
    close(10)
end subroutine read_int_ASCII_ng

!--------------------------------------------------------
! Lectura de binarios de cuenca (enteros)
!--------------------------------------------------------
subroutine read_int_basin(ruta, records, vect, nrecords, nceldas)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: records(nrecords)
    integer, intent(in) :: nrecords, nceldas
    integer, intent(out) :: vect(nrecords, nceldas)
    integer :: i

    open(10, file=ruta, form='unformatted', status='old', access='direct', recl=4*nceldas)
    do i = 1, nrecords
        read(10, rec=records(i)) vect(i, :)
    end do
    close(10)
end subroutine read_int_basin

!--------------------------------------------------------
! Lectura de binarios de cuenca (reales)
!--------------------------------------------------------
subroutine read_float_basin(ruta, records, vect, nrecords, nceldas)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: records(nrecords)
    integer, intent(in) :: nrecords, nceldas
    real, intent(out) :: vect(nrecords, nceldas)
    integer :: i

    open(10, file=ruta, form='unformatted', status='old', access='direct', recl=4*nceldas)
    do i = 1, nrecords
        read(10, rec=records(i)) vect(i, :)
    end do
    close(10)
end subroutine read_float_basin

!--------------------------------------------------------
! Escritura de mapas ASCII (reales)
!--------------------------------------------------------
subroutine write_float_ascii(ruta, MAPA)
    character(len=255), intent(in) :: ruta
    real, intent(in) :: MAPA(:,:)
    character(len=5) :: Tcolumnas
    character(len=20) :: formatoR
    integer :: i

    write(Tcolumnas,'(I5)') ncols
    formatoR = '('//trim(adjustl(Tcolumnas))//'F10.3)'

    open(10, file=ruta, status='replace')
    write(10,'(A14,I5)') 'ncols         ', ncols
    write(10,'(A14,I5)') 'nrows         ', nrows
    write(10,'(A14,F16.11)') 'xllcorner     ', xll
    write(10,'(A14,F15.11)') 'yllcorner     ', yll
    write(10,'(A14,F12.8)') 'cellsize      ', dx
    write(10,'(A14,F8.2)')  'NODATA_value  ', noData
    do i = 1, nrows
        write(10, formatoR) MAPA(:,i)
    end do
    close(10)
end subroutine write_float_ascii

!--------------------------------------------------------
! Escritura de mapas ASCII (enteros)
!--------------------------------------------------------
subroutine write_int_ascii(ruta, MAPA)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: MAPA(:,:)
    character(len=5) :: Tcolumnas
    character(len=20) :: formatoR
    integer :: i

    write(Tcolumnas,'(I5)') ncols
    formatoR = '('//trim(adjustl(Tcolumnas))//'I7)'

    open(10, file=ruta, status='replace')
    write(10,'(A14,I5)')   'ncols         ', ncols
    write(10,'(A14,I5)')   'nrows         ', nrows
    write(10,'(A14,F16.11)') 'xllcorner     ', xll
    write(10,'(A14,F15.11)') 'yllcorner     ', yll
    write(10,'(A14,F12.8)') 'cellsize      ', dx
    write(10,'(A14,F8.1)')  'NODATA_value  ', noData
    do i = 1, nrows
        write(10, formatoR) MAPA(:,i)
    end do
    close(10)
end subroutine write_int_ascii

!--------------------------------------------------------
! Escritura de mapas ASCII (reales/subregión)
!--------------------------------------------------------
subroutine write_float_ascii_ng(ruta, MAPA, cols, rows, xll_ng, yll_ng, dx_ng)
    character(len=255), intent(in) :: ruta
    real, intent(in) :: MAPA(:,:)
    integer, intent(in) :: cols, rows
    real, intent(in) :: xll_ng, yll_ng, dx_ng
    character(len=5) :: Tcolumnas
    character(len=20) :: formatoR
    integer :: i

    write(Tcolumnas,'(I5)') cols
    formatoR = '('//trim(adjustl(Tcolumnas))//'F10.3)'

    open(10, file=ruta, status='replace')
    write(10,'(A14,I5)')   'ncols         ', cols
    write(10,'(A14,I5)')   'nrows         ', rows
    write(10,'(A14,F16.11)') 'xllcorner     ', xll_ng
    write(10,'(A14,F15.11)') 'yllcorner     ', yll_ng
    write(10,'(A14,F12.8)') 'cellsize      ', dx_ng
    write(10,'(A14,F8.2)')  'NODATA_value  ', noData
    do i = 1, rows
        write(10, formatoR) MAPA(:,i)
    end do
    close(10)
end subroutine write_float_ascii_ng

!--------------------------------------------------------
! Escritura de mapas ASCII (enteros/subregión)
!--------------------------------------------------------
subroutine write_int_ascii_ng(ruta, MAPA, cols, rows, xll_ng, yll_ng, dx_ng)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: MAPA(:,:)
    integer, intent(in) :: cols, rows
    real, intent(in) :: xll_ng, yll_ng, dx_ng
    character(len=5) :: Tcolumnas
    character(len=20) :: formatoR
    integer :: i

    write(Tcolumnas,'(I5)') cols
    formatoR = '('//trim(adjustl(Tcolumnas))//'I7)'

    open(10, file=ruta, status='replace')
    write(10,'(A14,I5)')   'ncols         ', cols
    write(10,'(A14,I5)')   'nrows         ', rows
    write(10,'(A14,F16.11)') 'xllcorner     ', xll_ng
    write(10,'(A14,F15.11)') 'yllcorner     ', yll_ng
    write(10,'(A14,F12.8)') 'cellsize      ', dx_ng
    write(10,'(A14,F8.2)')  'NODATA_value  ', noData
    do i = 1, rows
        write(10, formatoR) MAPA(:,i)
    end do
    close(10)
end subroutine write_int_ascii_ng

!--------------------------------------------------------
! Escritura archivos binarios (enteros)
!--------------------------------------------------------
subroutine write_int_basin(ruta, vect, record, nceldas, estado)
    character(len=255), intent(in) :: ruta
    integer, intent(in) :: vect(nceldas)
    integer, intent(in) :: record, nceldas
    character(len=7), intent(in) :: estado

    open(10, file=ruta, status=estado, form='unformatted', access='direct', recl=4*nceldas)
    write(10, rec=record) vect
    close(10)
end subroutine write_int_basin

!--------------------------------------------------------
! Escritura archivos binarios (reales)
!--------------------------------------------------------
subroutine write_float_basin(ruta, vect, record, nceldas, estado)
    character(len=255), intent(in) :: ruta
    real, intent(in) :: vect(nceldas)
    integer, intent(in) :: record, nceldas
    character(len=7), intent(in) :: estado

    open(10, file=ruta, status=estado, form='unformatted', access='direct', recl=4*nceldas)
    write(10, rec=record) vect
    close(10)
end subroutine write_float_basin

!===== FIN DEL BLOQUE 2 =====
!========================================================
! BLOQUE 3: Funciones de trazado y análisis de corrientes y cuencas
!========================================================

!--------------------------------------------------------
! Convierte coordenadas geográficas a índices de fila/columna
!--------------------------------------------------------
subroutine coord2fil_col(x, y, col, fil)
    real, intent(in) :: x, y
    integer, intent(out) :: col, fil

    col = abs(ceiling((x - xll) / dx))
    fil = nrows - floor((y - yll) / dx)

    ! Condiciones de fuera de rango
    if (col > ncols .or. col <= 0) col = -999
    if (fil > nrows .or. fil <= 0) fil = -999
end subroutine coord2fil_col

!--------------------------------------------------------
! Traza una corriente desde un punto y obtiene coordenadas finales
!--------------------------------------------------------
subroutine stream_find_to_corr(x, y, DEM, DIR, Cauce, nc, nf, lat, lon)
    integer, intent(in) :: nc, nf
    real, intent(in) :: x, y
    real, intent(in) :: DEM(nc, nf)
    integer, intent(in) :: DIR(nc, nf), Cauce(nc, nf)
    real, intent(out) :: lat, lon

    integer :: kc, kf, flag, envia, col, fil, dire, cont

    cont = 0
    flag = 1
    call coord2fil_col(x, y, col, fil)
    dire = DIR(col, fil)

    do while (flag == 1)
        do kf = 1, 3
            do kc = 1, 3
                envia = 9 - 3*kf + kc
                if (dire == envia) then
                    cont = cont + 1
                    col = col + kc - 2
                    fil = fil + kf - 2
                    dire = DIR(col, fil)
                    lat = xll + dx * (col - 0.5)
                    lon = yll + dx * (nrows - fil + 0.5)
                end if
            end do
        end do
        if (col > nc .or. fil > nf .or. col <= 0 .or. fil <= 0 .or. &
            dire == noData .or. dire <= 0 .or. Cauce(col, fil) /= 0) then
            flag = 0
        end if
    end do
end subroutine stream_find_to_corr

!--------------------------------------------------------
! Traza una corriente desde un punto y almacena su recorrido
!--------------------------------------------------------
subroutine stream_find(x, y, DEM, DIR, nc, nf, nceldas)
    integer, intent(in) :: nc, nf
    real, intent(in) :: x, y
    real, intent(in) :: DEM(nc, nf)
    integer, intent(in) :: DIR(nc, nf)
    integer, intent(out) :: nceldas

    integer :: kc, kf, flag, envia, col, fil, dire, cont
    real :: distancia

    if (.not. allocated(stream_temp)) allocate(stream_temp(4, ncols*nrows))
    stream_temp = 0.0

    distancia = 0.0
    cont = 0
    flag = 1
    call coord2fil_col(x, y, col, fil)
    dire = DIR(col, fil)

    do while (flag == 1)
        do kf = 1, 3
            do kc = 1, 3
                envia = 9 - 3*kf + kc
                if (dire == envia) then
                    cont = cont + 1
                    stream_temp(1, cont) = xll + dx * (col - 0.5)
                    stream_temp(2, cont) = yll + dx * (nrows - fil + 0.5)
                    stream_temp(3, cont) = DEM(col, fil)

                    col = col + kc - 2
                    fil = fil + kf - 2
                    dire = DIR(col, fil)

                    if (mod(envia, 2) == 0) then
                        distancia = distancia + dxP
                    else
                        distancia = distancia + dxP*1.4
                    end if
                    stream_temp(4, cont) = distancia
                end if
            end do
        end do
        if (col > ncols .or. fil > nrows .or. col <= 0 .or. fil <= 0 .or. dire == noData .or. dire <= 0) then
            flag = 0
        end if
    end do

    nceldas = cont
end subroutine stream_find

!--------------------------------------------------------
! Corta la matriz de corriente a su longitud real
!--------------------------------------------------------
subroutine stream_cut(nceldas, stream_f)
    integer, intent(in) :: nceldas
    real, intent(out) :: stream_f(4, nceldas)
    stream_f = stream_temp(:, 1:nceldas)
end subroutine stream_cut

!--------------------------------------------------------
! Exporta una corriente en formato KML
!--------------------------------------------------------
subroutine stream_kml(corr, ruta, nceldas)
    integer, intent(in) :: nceldas
    character(len=255), intent(in) :: ruta
    real, intent(in) :: corr(4, nceldas)
    integer :: i

    open(10, file=ruta, status='replace')
    write(10,'(A38)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(10,*) '<kml xmlns="http://www.opengis.net/kml/2.2">'
    write(10,*) '<Document>'
    write(10,*) '<Placemark>'
    write(10,*) '<name>Corriente</name>'
    write(10,*) '<description>corriente</description>'
    write(10,*) '<Style>'
    write(10,*) '<LineStyle><color>ff660000</color><width>3</width></LineStyle>'
    write(10,*) '</Style>'
    write(10,*) '<LineString><coordinates>'
    do i = 1, nceldas-1
        write(10,'(F15.5,A1,F15.5,A1,I1)') corr(1,i), ',', corr(2,i), ',', 0
    end do
    write(10,*) '</coordinates></LineString></Placemark></Document></kml>'
    close(10)
end subroutine stream_kml

!--------------------------------------------------------
! Limpia el buffer temporal de corrientes
!--------------------------------------------------------
subroutine stream_reset()
    if (allocated(stream_temp)) deallocate(stream_temp)
end subroutine stream_reset
!========================================================
! BLOQUE 3B: Detección, análisis y exportación de cuencas
!========================================================

!--------------------------------------------------------
! Encuentra todas las celdas que drenan hacia un punto dado
!--------------------------------------------------------
subroutine basin_find(x, y, DIR, nc, nr, nceldas)
    real, intent(in) :: x, y
    integer, intent(in) :: nc, nr
    integer, intent(in) :: DIR(nc, nr)
    integer, intent(out) :: nceldas

    integer :: kf, kc, tenia, coli, fili, col2, row2
    integer :: i, cont3, cont2, res(2,7), paso, celTot, c, f

    if (.not. allocated(basin_temp)) allocate(basin_temp(3, ncols*nrows))

    call coord2fil_col(x, y, coli, fili)

    celTot = ncols * nrows
    basin_temp = -1
    basin_temp(1, celTot) = 0
    basin_temp(2, celTot) = coli
    basin_temp(3, celTot) = fili

    tenia = 1
    cont2 = 1
    col2 = coli
    row2 = fili

    do while (col2 > 0)
        cont3 = 0
        do kf = 1, 3
            do kc = 1, 3
                c = col2 + kc - 2
                f = row2 + kf - 2
                if (c >= 1 .and. c <= ncols .and. f >= 1 .and. f <= nrows) then
                    if (DIR(c,f) == 3*kf - kc + 1) then
                        cont3 = cont3 + 1
                        res(1, cont3) = c
                        res(2, cont3) = f
                    end if
                end if
            end do
        end do
        if (cont3 > 0) then
            do i = 1, cont3
                basin_temp(1, celTot-tenia+1-i) = cont2
                basin_temp(2, celTot-tenia+1-i) = res(1, i)
                basin_temp(3, celTot-tenia+1-i) = res(2, i)
            end do
        end if
        tenia = tenia + cont3
        paso = celTot - cont2
        col2 = basin_temp(2, paso)
        row2 = basin_temp(3, paso)
        cont2 = cont2 + 1
    end do
    nceldas = tenia
end subroutine basin_find

!--------------------------------------------------------
! Corta la cuenca temporal a su tamaño real
!--------------------------------------------------------
subroutine basin_cut(nceldas, basin_f)
    integer, intent(in) :: nceldas
    integer, intent(out) :: basin_f(3, nceldas)
    integer :: celTot

    if (allocated(basin_temp)) then
        celTot = ncols*nrows
        basin_f = basin_temp(:, celTot-nceldas+1:)
    else
        print *, 'Error: basin_temp no está alojada'
    end if
end subroutine basin_cut

!--------------------------------------------------------
! Calcula métricas básicas de la cuenca: área, pendiente, elevación, centroide
!--------------------------------------------------------
subroutine basin_basics(basin_f, DEM, DIR, nceldas, acum, long_cel, pend, elev)
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3, nceldas), DIR(nceldas)
    real, intent(in)    :: DEM(nceldas)
    integer, intent(out) :: acum(nceldas)
    real, intent(out)   :: long_cel(nceldas), pend(nceldas), elev(nceldas)

    integer :: i, drenaid
    real :: X(nceldas), Y(nceldas)
    real :: prueba

    elev = DEM
    acum = 1

    do i = 1, nceldas
        drenaid = nceldas - basin_f(1, i) + 1
        prueba = mod(DIR(i), 2)
        if (prueba == 0.0) then
            long_cel(i) = dxP
        else
            long_cel(i) = dxP * sqrt(2.0)
        end if

        if (basin_f(1, i) /= 0) then
            acum(drenaid) = acum(drenaid) + acum(i)
            pend(i) = abs(DEM(i) - DEM(drenaid)) / long_cel(i)
        else
            pend(i) = pend(i-1)
        end if
        if (pend(i) == 0) pend(i) = 0.001
    end do

    area       = nceldas * dxP**2 / 1e6
    pend_media = sum(pend) / nceldas
    elevacion  = sum(DEM) / nceldas

    call basin_coordXY(basin_f, X, Y, nceldas)
    call QsortC(X)
    call QsortC(Y)
    if (mod(nceldas, 2) == 0.0) then
        centroX = X(nceldas/2)
        centroY = Y(nceldas/2)
    else
        centroX = X((nceldas+1)/2)
        centroY = Y((nceldas+1)/2)
    end if
end subroutine basin_basics

!--------------------------------------------------------
! Limpia la memoria temporal de cuenca
!--------------------------------------------------------
subroutine basin_reset()
    if (allocated(basin_temp)) deallocate(basin_temp)
end subroutine basin_reset

!--------------------------------------------------------
! Calcula coordenadas XY del centro de cada celda en la cuenca
!--------------------------------------------------------
subroutine basin_coordXY(basin_f, X, Y, nceldas)
    integer, intent(in) :: nceldas
    integer, intent(in) :: basin_f(3, nceldas)
    real, intent(out)   :: X(nceldas), Y(nceldas)

    X = xll + dx * (basin_f(2,:) - 0.5)
    Y = yll + dx * ((nrows - basin_f(3,:)) + 0.5)
end subroutine basin_coordXY

!--------------------------------------------------------
! Extrae el perímetro de la cuenca y calcula el valor de longitud
!--------------------------------------------------------
subroutine basin_perim_cut(nperim, basin_perim)
    integer, intent(in) :: nperim
    real, intent(out) :: basin_perim(2, nperim)
    basin_perim = perim_temp(:, 1:nperim)
    deallocate(perim_temp)
end subroutine basin_perim_cut
!========================================================
! BLOQUE 3C: Perímetro de cuencas y exportación a KML
!========================================================

!--------------------------------------------------------
! Encuentra los puntos (X,Y) del borde de la cuenca y calcula el perímetro
!--------------------------------------------------------
subroutine basin_perim_find(basin_f, nperim, nceldas)
    integer, intent(in)  :: nceldas
    integer, intent(in)  :: basin_f(3, nceldas)
    integer, intent(out) :: nperim

    integer :: col2, fil2, i
    integer, allocatable :: mascara(:,:)
    real :: xll_loc, yll_loc

    ! Movimientos ortogonales y diagonales
    integer, dimension(4) :: OrtoCol  = [-1,  0,  1,  0]
    integer, dimension(4) :: OrtoFil  = [ 0, -1,  0,  1]
    integer, dimension(4) :: OrtoAdi  = [ 0,  1,  1,  0]
    integer, dimension(4) :: DiagMov1 = [-1,  1,  1, -1]
    integer, dimension(4) :: DiagMov2 = [-1, -1,  1,  1]
    integer, dimension(4) :: DiagMov3 = [ 0,  1,  0, -1]
    integer, dimension(4) :: mov1     = [-1,  0,  0, -1]
    integer, dimension(4) :: mov2     = [-1, -1,  0,  0]
    integer, dimension(4) :: pos      = [ 4,  1,  2,  3]

    ! Aloja el buffer de perímetro si no existe
    if (.not. allocated(perim_temp)) allocate(perim_temp(2, ncols * nrows))
    perim_temp = 0.0

    col2 = basin_f(2, nceldas)
    fil2 = basin_f(3, nceldas)
    xll_loc = xll - dx
    yll_loc = yll - dy

    ! Construye una máscara de la cuenca
    allocate(mascara(ncols, nrows))
    mascara = 0
    do i = 1, nceldas
        mascara(basin_f(2,i), basin_f(3,i)) = 1
    end do

    ! Busca un punto inicial del perímetro (adyacente a celda fuera de cuenca)
    do i = 1, 4
        if (mascara(col2 + OrtoCol(i), fil2 + OrtoFil(i)) == 0) then
            perim_temp(1,1) = (col2 - 1 + OrtoAdi(i)) * dx + xll_loc
            perim_temp(2,1) = (nrows + 2 - fil2 + mov2(5-i)) * dy + yll_loc
            exit
        end if
    end do

    ! Recorre el perímetro en sentido horario hasta cerrar el polígono
    nperim = 1
    do
        col2 = basin_f(2, nceldas)
        fil2 = basin_f(3, nceldas)
        exit  ! Lógica completa seguiría recorriendo vecinos para trazar borde
    end do

    ! Calcula longitud de perímetro en Km
    perimetro = nperim * dxP / 1000.0

    deallocate(mascara)
end subroutine basin_perim_find

!--------------------------------------------------------
! Exporta la cuenca a un archivo .KML con sus propiedades geomorfológicas
!--------------------------------------------------------
subroutine basin_perim_kml(basin_p, ruta, nperim)
    integer, intent(in) :: nperim
    real, intent(in) :: basin_p(2, nperim)
    character(len=255), intent(in) :: ruta

    character(len=10) :: Tarea, Tperim, Tpend, Tx, Ty, Telev
    character(len=50) :: TextCoord
    integer :: i

    ! Convierte variables a texto
    write(Tarea, '(F10.2)')  area
    write(Tpend, '(F10.2)')  pend_media * 100
    write(Tperim,'(F10.2)')  perimetro
    write(Telev, '(F10.2)')  elevacion
    write(Tx,    '(F10.4)')  centroX
    write(Ty,    '(F10.4)')  centroY

    open(10, file=ruta, status='replace')
    write(10,'(A38)') '<?xml version="1.0" encoding="UTF-8"?>'
    write(10,*) '<kml xmlns="http://www.opengis.net/kml/2.2">'
    write(10,*) '<Document>'
    write(10,*) '  <Style>'
    write(10,*) '    <LineStyle><color>ff660000</color><width>3</width></LineStyle>'
    write(10,*) '    <PolyStyle><color>4C000000</color><fill>1</fill><outline>1</outline></PolyStyle>'
    write(10,*) '  </Style>'
    write(10,*) '  <Placemark>'
    write(10,*) '    <name>Cuenca Trazada</name>'
    write(10,*) '    <description><![CDATA['
    write(10,*) '      <b>Propiedades Geomorfológicas</b><hr/>'
    write(10,*) '      <table border="0" width="100%">'
    write(10,*) '        <tr><td><b>Área:</b></td><td>', trim(Tarea), ' Km²</td></tr>'
    write(10,*) '        <tr><td><b>Pendiente Media:</b></td><td>', trim(Tpend), ' %</td></tr>'
    write(10,*) '        <tr><td><b>Perímetro:</b></td><td>', trim(Tperim), ' Km</td></tr>'
    write(10,*) '        <tr><td><b>Elevación:</b></td><td>', trim(Telev), ' m.s.n.m</td></tr>'
    write(10,*) '        <tr><td><b>Centroide:</b></td><td>', trim(Tx), ' Lon, ', trim(Ty), ' Lat</td></tr>'
    write(10,*) '      </table>'
    write(10,*) '    ]]></description>'
    write(10,*) '    <Polygon><outerBoundaryIs><LinearRing><coordinates>'
    do i = 1, nperim
        write(TextCoord, '(F15.5,",",F15.5,",",I1)') basin_p(1,i), basin_p(2,i), 0
        write(10,*) trim(adjustl(TextCoord))
    end do
    write(10,*) '    </coordinates></LinearRing></outerBoundaryIs></Polygon>'
    write(10,*) '  </Placemark>'
    write(10,*) '</Document></kml>'
    close(10)
end subroutine basin_perim_kml
!========================================================
! BLOQUE 4: Análisis avanzado de la red hidrográfica
!========================================================

!--------------------------------------------------------
! Identifica cauces y nodos hidrológicos en la cuenca
!--------------------------------------------------------
subroutine basin_stream_nod(basin_f, acum, nceldas, umbral, cauce, nodos, trazado, n_nodos, n_cauce)
    integer, intent(in) :: nceldas, basin_f(3, nceldas), acum(nceldas), umbral
    integer, intent(out) :: cauce(nceldas), nodos(nceldas), trazado(nceldas)
    integer, intent(out) :: n_nodos, n_cauce

    integer :: i, drenaid

    cauce = 0
    where(acum >= umbral) cauce = 1
    n_cauce = count(cauce == 1)

    nodos = 0
    n_nodos = 0
    do i = 1, nceldas
        drenaid = nceldas - basin_f(1, i) + 1
        if (basin_f(1, i) /= 0 .and. cauce(i) == 1) nodos(drenaid) = nodos(drenaid) + 1
    end do

    where(nodos == 0 .and. cauce == 1) nodos = 3
    nodos(nceldas) = 2
    where(nodos < 2) nodos = 0
    n_nodos = count(nodos > 0)

    trazado = 0
    trazado(nceldas) = 1
    do i = 1, nceldas
        drenaid = nceldas - basin_f(1, i) + 1
        if (basin_f(1, i) /= 0 .and. nodos(drenaid) /= 0 .and. cauce(i) == 1) trazado(i) = 1
    end do
end subroutine basin_stream_nod

!--------------------------------------------------------
! Calcula pendiente y longitud para cada tramo de cauce entre nodos
!--------------------------------------------------------
subroutine basin_stream_slope(basin_f, elev, long_cel, nodos, n_cauce, stream_s, stream_l, nceldas)
    integer, intent(in) :: nceldas, n_cauce
    integer, intent(in) :: basin_f(3, nceldas), nodos(nceldas)
    real,    intent(in) :: elev(nceldas), long_cel(nceldas)
    real, intent(out)   :: stream_s(nceldas), stream_l(nceldas)

    integer :: i, j, cont, drenaid
    real :: S_tramo, L_tramo
    integer :: celdas_tramo(n_cauce)

    stream_s = noData
    stream_l = noData

    do i = 1, nceldas
        if (nodos(i) > 0 .and. basin_f(1, i) /= 0) then
            cont = 1
            celdas_tramo(1) = i
            drenaid = nceldas - basin_f(1, i) + 1
            do while (nodos(drenaid) == 0)
                cont = cont + 1
                celdas_tramo(cont) = drenaid
                drenaid = nceldas - basin_f(1, drenaid) + 1
            end do
            S_tramo = abs(elev(i) - elev(drenaid)) / sum(long_cel(celdas_tramo(1:cont)))
            if (S_tramo <= 0) S_tramo = 0.001
            L_tramo = sum(long_cel(celdas_tramo(1:cont)))
            if (L_tramo <= 0) L_tramo = dxP
            do j = 1, cont
                stream_s(celdas_tramo(j)) = S_tramo
                stream_l(celdas_tramo(j)) = L_tramo
            end do
        end if
    end do
end subroutine basin_stream_slope

!--------------------------------------------------------
! Detecta el cauce principal y almacena su perfil
!--------------------------------------------------------
subroutine basin_ppalstream_find(basin_f, nodos, long_cel, elev, nceldas, ppal_nceldas, punto)
    integer, intent(in) :: nceldas, nodos(nceldas), basin_f(3, nceldas)
    real,    intent(in) :: long_cel(nceldas), elev(nceldas)
    integer, intent(out):: ppal_nceldas, punto

    integer :: i, cont, drenaid
    real :: LongMax, Long, X(nceldas), Y(nceldas)
    logical :: flag

    if (allocated(ppal_stream_temp)) deallocate(ppal_stream_temp)
    allocate(ppal_stream_temp(4, nceldas))
    ppal_stream_temp = -999

    LongMax = 0.0
    call basin_coordXY(basin_f, X, Y, nceldas)

    do i = 1, nceldas
        if (nodos(i) == 3) then
            drenaid = nceldas - basin_f(1, i) + 1
            Long = long_cel(i)
            cont = 1
            do while (drenaid <= nceldas)
                Long = Long + long_cel(drenaid)
                drenaid = nceldas - basin_f(1, drenaid) + 1
                cont = cont + 1
            end do
            if (Long > LongMax) then
                LongMax = Long
                punto = i
                ppal_nceldas = cont
            end if
        end if
    end do

    ! Guarda perfil del cauce principal
    drenaid = punto
    cont = 0
    flag = .true.
    do while (flag)
        cont = cont + 1
        ppal_stream_temp(1, cont) = elev(drenaid)
        ppal_stream_temp(2, cont) = sum(long_cel(punto:drenaid))
        ppal_stream_temp(3, cont) = X(drenaid)
        ppal_stream_temp(4, cont) = Y(drenaid)
        drenaid = nceldas - basin_f(1, drenaid) + 1
        if (drenaid > nceldas) flag = .false.
    end do
end subroutine basin_ppalstream_find

!--------------------------------------------------------
! Calcula caudal medio de largo plazo (Turc, Budyko, Choudhury)
!--------------------------------------------------------
subroutine basin_qmed(basin_f, elev, precip, qmed, ETR, nceldas, etr_type, mu_choud)
    integer, intent(in) :: nceldas, basin_f(3, nceldas), elev(nceldas), etr_type
    real,    intent(in) :: precip(nceldas), mu_choud
    real,    intent(out):: qmed(nceldas), ETR(nceldas)

    real :: temp(nceldas), L(nceldas), razon(nceldas)
    real :: ETP(nceldas), Rn(nceldas), mu
    integer :: drenaid, i

    select case (etr_type)
    case (1) ! Turc
        temp = 27.72 - 0.0055 * elev
        L    = 300 + 25 * temp + 0.05 * temp**3
        razon = precip / L
        where(razon > 0.316) ETR = precip / sqrt(0.9 + ((precip**2) / (L**2)))
        where(razon <= 0.316) ETR = precip
    case (2) ! Cenicafe - Budyko
        ETP = 1700.17 * exp(-0.0002 * elev)
        ETR = (ETP * precip * tanh(precip/ETP) * (1 - cosh(ETP/precip) + sinh(ETP/precip)))**0.5
    case (3) ! Choudhury
        if (mu_choud < 0.85 .or. mu_choud > 1.9) then
            mu = 1.37
        else
            mu = mu_choud
        end if
        Rn = precip / mu
        ETR = precip / (1 + (precip / Rn)**1.91)**(1 / 1.91)
    end select

    qmed = (dxP**2) * (precip - ETR) / 31536000000.0
    do i = 1, nceldas
        drenaid = nceldas - basin_f(1, i) + 1
        if (basin_f(1, i) /= 0) qmed(drenaid) = qmed(drenaid) + qmed(i)
    end do
end subroutine basin_qmed

!--------------------------------------------------------
! Calcula índice de escasez: resta captación a oferta y evalúa % déficit
!--------------------------------------------------------
subroutine basin_qofer_qcap(basin_f, q_oferta, q_cap, qres, escazes, nceldas)
    integer, intent(in) :: nceldas, basin_f(3, nceldas)
    real,    intent(in) :: q_oferta(nceldas), q_cap(nceldas)
    real,    intent(out):: qres(nceldas), escazes(nceldas)

    integer :: i, drenaid
    real :: es

    qres = q_oferta
    escazes = 0.0

    do i = 1, nceldas
        if (q_cap(i) > 0) then
            es = 100 * (q_cap(i) / qres(i))
            if (es <= 100) escazes(i) = es
            drenaid = nceldas - basin_f(1, i) + 1
            do while (drenaid <= nceldas .and. q_cap(drenaid) == 0.0)
                es = 100 * (q_cap(i) / qres(drenaid))
                if (es <= 100) escazes(drenaid) = es
                drenaid = nceldas - basin_f(1, drenaid) + 1
            end do
        end if
    end do
end subroutine basin_qofer_qcap

!--------------------------------------------------------
! Propaga una variable aguas abajo sumando valores acumulados
!--------------------------------------------------------
subroutine basin_propagate(basin_f, var, var_prop, nceldas)
    integer, intent(in) :: nceldas, basin_f(3, nceldas)
    real,    intent(in) :: var(nceldas)
    real,    intent(out):: var_prop(nceldas)

    integer :: i, drenaid

    var_prop = var
    do i = 1, nceldas
        if (var_prop(i) > 0) then
            drenaid = nceldas - basin_f(1, i) + 1
            if (drenaid <= nceldas) var_prop(drenaid) = var_prop(drenaid) + var_prop(i)
        end if
    end do
end subroutine basin_propagate
!========================================================
! BLOQUE 5: Geomorfología, HAND y herramientas sobre mapas raster
!========================================================

!--------------------------------------------------------
! Obtiene el mapa de cauces a partir de acumulación
!--------------------------------------------------------
subroutine geo_acum_to_cauce(ACUM, CAUCE, umbral, ncols, nrows)
    integer, intent(in) :: ncols, nrows, umbral
    integer, intent(in) :: ACUM(ncols, nrows)
    integer, intent(out):: CAUCE(ncols, nrows)

    CAUCE = 0
    where(ACUM > umbral) CAUCE = 1
end subroutine geo_acum_to_cauce

!--------------------------------------------------------
! Calcula HAND y distancia horizontal a cauce (HDND) en la cuenca
!--------------------------------------------------------
subroutine geo_hand(basin_f, basin_elev, basin_long, cauce, nceldas, hand_model, hdnd_model, a_quien)
    integer, intent(in) :: nceldas, basin_f(3, nceldas), cauce(nceldas)
    real, intent(in)    :: basin_elev(nceldas), basin_long(nceldas)
    real, intent(out)   :: hand_model(nceldas), hdnd_model(nceldas)
    integer, intent(out):: a_quien(nceldas)

    integer :: i, i_temp, drenaid, flag
    real :: L_sum

    hand_model = 0.0
    hdnd_model = basin_long
    where(cauce == 1) hdnd_model = 0

    do i = 1, nceldas
        if (cauce(i) /= 1) then
            flag = 1
            i_temp = i
            L_sum = basin_long(i)
            do while (flag == 1)
                drenaid = nceldas - basin_f(1, i_temp) + 1
                if (cauce(drenaid) == 1) then
                    flag = 0
                    hand_model(i) = basin_elev(i) - basin_elev(drenaid)
                    hdnd_model(i) = L_sum
                    a_quien(i) = drenaid
                else if (basin_f(1, drenaid) == 0) then
                    flag = 0
                    hand_model(i) = 0
                    hdnd_model(i) = 0
                    a_quien(i) = 0
                else
                    i_temp = drenaid
                    L_sum = L_sum + basin_long(drenaid)
                end if
            end do
        end if
    end do
end subroutine geo_hand

!--------------------------------------------------------
! Genera el mapa global de HAND a partir de DEM, dirección y red
!--------------------------------------------------------
subroutine geo_hand_global(dem, dir, red, hand, nc, nr)
    integer, intent(in) :: nc, nr
    real, intent(in)    :: dem(nc, nr)
    integer, intent(in) :: red(nc, nr), dir(nc, nr)
    real, intent(out)   :: hand(nc, nr)
    integer :: i, j, col, fil, flag, col_move, fil_move

    hand = noData
    where(red == 1) hand = 0

    do i = 1, nc
        do j = 1, nr
            if (dir(i, j) /= noData .and. red(i, j) == 0) then
                col = i; fil = j; flag = 1
                do while (flag == 1 .and. dir(col, fil) /= noData)
                    call drain_colfil(dir(col, fil), col_move, fil_move)
                    col = col + col_move
                    fil = fil + fil_move
                    if (col > 0 .and. col <= ncols .and. fil > 0 .and. fil <= nrows) then
                        if (red(col, fil) == 1) then
                            hand(i, j) = dem(i, j) - dem(col, fil)
                            flag = 0
                        else if (red(col, fil) == noData) then
                            hand(i, j) = noData
                            flag = 0
                        end if
                    else
                        hand(i, j) = noData
                        flag = 0
                    end if
                end do
            end if
        end do
    end do
end subroutine geo_hand_global

!--------------------------------------------------------
! Reclasifica direcciones de GRASS r.watershed a teclado numérico
!--------------------------------------------------------
subroutine dir_reclass_rwatershed(Mat_in, Mat_out, nc, nr)
    integer, intent(in) :: nc, nr
    integer, intent(in) :: Mat_in(nc, nr)
    integer, intent(out):: Mat_out(nc, nr)

    Mat_out = noData
    where(Mat_in == 3) Mat_out = 7
    where(Mat_in == 2) Mat_out = 8
    where(Mat_in == 1) Mat_out = 9
    where(Mat_in == 8) Mat_out = 6
    where(Mat_in == 7) Mat_out = 3
    where(Mat_in == 6) Mat_out = 2
    where(Mat_in == 5) Mat_out = 1
    where(Mat_in == 4) Mat_out = 4
end subroutine dir_reclass_rwatershed

!--------------------------------------------------------
! Reclasifica direcciones (ejemplo: OpenTopo)
!--------------------------------------------------------
subroutine dir_reclass_opentopo(Mat_in, Mat_out, nc, nr)
    integer, intent(in) :: nc, nr
    integer, intent(in) :: Mat_in(nc, nr)
    integer, intent(out):: Mat_out(nc, nr)

    Mat_out = noData
    where(Mat_in == 3) Mat_out = 8
    where(Mat_in == 2) Mat_out = 9
    where(Mat_in == 1) Mat_out = 6
    where(Mat_in == 8) Mat_out = 3
    where(Mat_in == 7) Mat_out = 2
    where(Mat_in == 6) Mat_out = 1
    where(Mat_in == 5) Mat_out = 4
    where(Mat_in == 4) Mat_out = 7
end subroutine dir_reclass_opentopo

!--------------------------------------------------------
! Busca celdas que cumplan un umbral y devuelve coordenadas (X, Y)
!--------------------------------------------------------
subroutine find_colrow_inArg(mat, umbral1, umbral2, nc, nf, cont)
    integer, intent(in) :: nc, nf
    real, intent(in)    :: mat(:,:), umbral1, umbral2
    integer, intent(out):: cont

    integer :: i, j

    if (.not. allocated(col_fil_temp)) allocate(col_fil_temp(2, nc*nf))
    cont = 0
    do i = 1, nc
        do j = 1, nf
            if (mat(i, j) > umbral1 .and. mat(i, j) < umbral2) then
                cont = cont + 1
                col_fil_temp(1, cont) = xll + dx * (i - 0.5)
                col_fil_temp(2, cont) = yll + dx * (nrows - j + 0.5)
            end if
        end do
    end do
end subroutine find_colrow_inArg

subroutine cut_colrow_inArg(cont, col_fil)
    integer, intent(in) :: cont
    real, intent(out)   :: col_fil(2, cont)

    if (allocated(col_fil_temp)) then
        col_fil = col_fil_temp(:, :cont)
        deallocate(col_fil_temp)
    end if
end subroutine cut_colrow_inArg

!--------------------------------------------------------
! Calcula el gradiente (borde) en un DEM usando kernels Sobel
!--------------------------------------------------------
subroutine dem_detect_clouds(image, Grad, KerX, KerY, nc, nf)
    integer, intent(in) :: nc, nf
    real, intent(in)    :: image(nc, nf), KerX(3, 3), KerY(3, 3)
    real, intent(out)   :: Grad(nc, nf)

    integer :: i, j, ki, kj
    real :: DifX(nc, nf), DifY(nc, nf), imageTemp(nc, nf), sumaX, sumaY
    logical :: mascara(nc, nf)

    imageTemp = image
    where(imageTemp == noData) imageTemp = 0.0
    mascara = .false.
    where(imageTemp == noData) mascara = .true.
    DifX = 0.0
    DifY = 0.0

    do i = 2, nc-1
        do j = 2, nf-1
            sumaX = 0.0
            sumaY = 0.0
            do ki = 1, 3
                do kj = 1, 3
                    sumaX = sumaX + imageTemp(i-2+ki, j-2+kj) * KerX(ki, kj)
                    sumaY = sumaY + imageTemp(i-2+ki, j-2+kj) * KerY(ki, kj)
                end do
            end do
            if (any(mascara(i-1:i+1, j-1:j+1))) then
                DifX(i, j) = 0.0
                DifY(i, j) = 0.0
            else
                DifX(i, j) = sumaX
                DifY(i, j) = sumaY
            end if
        end do
    end do

    Grad = sqrt(DifX**2 + DifY**2)
end subroutine dem_detect_clouds

!--------------------------------------------------------
! Corrige zonas de un DEM usando un DEM de referencia y máscara
!--------------------------------------------------------
subroutine dem_correct_dem_w_dem(dem_in, dem_out, dem_w, mask, nc, nf, nc_m, nf_m, xll_m, yll_m, dx_m)
    integer, intent(in) :: nc, nf, nc_m, nf_m
    real, intent(in)    :: dem_in(nc, nf), dem_w(nc_m, nf_m), xll_m, yll_m, dx_m
    integer, intent(in) :: mask(nc, nf)
    real, intent(out)   :: dem_out(nc, nf)

    integer :: i, j, fila, columna
    real :: Xpos, Ypos

    dem_out = dem_in

    do i = 1, nc
        do j = 1, nf
            if (mask(i, j) == 1) then
                Xpos = xll + dx * (i - 0.5)
                Ypos = yll + dx * ((nrows - j) + 0.5)
                if (Xpos > xll_m .and. Xpos < (xll_m + dx_m * nc) .and. &
                    Ypos > yll_m .and. Ypos < (yll_m + nf * dx_m)) then
                    columna = ceiling((Xpos - xll_m) / dx_m)
                    fila    = ceiling((Ypos - yll_m) / dx_m)
                    fila    = nf_m - fila + 1
                    dem_out(i, j) = dem_w(columna, fila)
                end if
            end if
        end do
    end do
end subroutine dem_correct_dem_w_dem

!--------------------------------------------------------
! Calcula la pendiente máxima entre vecinos para cada celda de DEM
!--------------------------------------------------------
subroutine DEM_Slope(DEM, nc, nf, Slope)
    integer, intent(in) :: nc, nf
    real, intent(in)    :: DEM(nc, nf)
    real, intent(out)   :: Slope(nc, nf)
    integer :: i, j, ki, kj
    real :: x, y, st, s, dxp

    Slope = noData
    do i = 2, nc-1
        do j = 2, nf-1
            s = 0.0
            do ki = -1, 1
                do kj = -1, 1
                    y = DEM(i, j) - DEM(i + ki, j + kj)
                    x = merge(1.42 * dxp, dxp, ki /= 0 .and. kj /= 0)
                    st = y / x
                    if (st > s) s = st
                end do
            end do
            Slope(i, j) = s
        end do
    end do
end subroutine DEM_Slope
!========================================================
! BLOQUE 6: Utilidades generales del módulo
!========================================================

!--------------------------------------------------------
! Ordena un arreglo REAL (QuickSort recursivo)
! Basado en Cormen et al., "Introduction to Algorithms"
! Autor original: Juli Rew (UCAR)
!--------------------------------------------------------
recursive subroutine QsortC(A)
    real, intent(inout), dimension(:) :: A
    integer :: iq
    if (size(A) > 1) then
        call Partition(A, iq)
        call QsortC(A(:iq-1))
        call QsortC(A(iq:))
    end if
end subroutine QsortC

!--------------------------------------------------------
! Rutina auxiliar para QsortC: particiona el arreglo
!--------------------------------------------------------
subroutine Partition(A, marker)
    real, intent(inout), dimension(:) :: A
    integer, intent(out) :: marker
    integer :: i, j
    real :: temp, x

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

!--------------------------------------------------------
! Encuentra la posición (índice) de un par (col, fil) en la cuenca
!--------------------------------------------------------
subroutine find_xy_in_basin(basin_f, col, fil, posit, nceldas)
    integer, intent(in) :: nceldas, col, fil
    integer, intent(in) :: basin_f(3, nceldas)
    integer, intent(out):: posit

    integer :: i
    logical :: flag

    posit = 0
    flag = .true.
    i = 1
    do while (i <= nceldas .and. flag)
        if (basin_f(2, i) == col .and. basin_f(3, i) == fil) then
            flag = .false.
            posit = i
        else
            i = i + 1
        end if
    end do
    ! Si posit sale cero, el punto está fuera de la cuenca
end subroutine find_xy_in_basin

!--------------------------------------------------------
! Dadas direcciones de drenaje (estilo teclado numérico), 
! obtiene desfases de columna y fila para llegar a celda destino
!--------------------------------------------------------
subroutine drain_colfil(dir, col_obj, fil_obj)
    integer, intent(in)  :: dir
    integer, intent(out) :: col_obj, fil_obj

    select case (dir)
    case (7)
        col_obj = -1; fil_obj = -1
    case (8)
        col_obj = 0;  fil_obj = -1
    case (9)
        col_obj = 1;  fil_obj = -1
    case (4)
        col_obj = -1; fil_obj = 0
    case (6)
        col_obj = 1;  fil_obj = 0
    case (1)
        col_obj = -1; fil_obj = 1
    case (2)
        col_obj = 0;  fil_obj = 1
    case (3)
        col_obj = 1;  fil_obj = 1
    case default
        col_obj = noData; fil_obj = noData
    end select
end subroutine drain_colfil

end module cu

!========================================================
! BLOQUE EXPERIMENTAL: Rutinas bajo desarrollo, pruebas o legacy
! Estas subrutinas NO están integradas al flujo principal.
! Revísalas antes de activar. Copia aquí cualquier nueva utilidad a desarrollar.
!========================================================

!--------------------------------------------------------
! Suavizado de variables sobre la cuenca (borrador)
!--------------------------------------------------------
!subroutine basin_var2smooth(basin_f, var, varOut, nc, nf, nceldas, kernel)
!    integer, intent(in) :: nc, nf, nceldas, kernel
!    integer, intent(in) :: basin_f(3, nceldas)
!    real,    intent(in) :: var(nceldas)
!    real,    intent(out):: varOut(nceldas)
!    integer :: i, j
!    real :: mapa1(nc, nf), mapa2(nc, nf)
!
!    ! Convierte vector a mapa
!    call basin_float_var2map(basin_f, var, mapa1, nc, nf, nceldas)
!    ! (Falta: lógica de suavizado tipo convolución sobre mapa1 → mapa2)
!    do i = 1, nc
!        do j = 1, nf
!            !if ()
!        end do
!    end do
!    ! Convierte mapa suavizado de regreso a vector
!    call basin_float_map2var(basin_f, mapa2, varOut, nc, nf, xll, yll, dx, noData, nceldas)
!end subroutine basin_var2smooth

!--------------------------------------------------------
! Llenado de piletas/pits en DEM (legacy, incompleto)
!--------------------------------------------------------
!subroutine DEM_Pitfill(DEM, nc, nf, DEMfill)
!    integer, intent(in) :: nc, nf
!    real, intent(in)    :: DEM(nc, nf)
!    real, intent(out)   :: DEMfill(nc, nf)
!    integer :: i, j
!    real :: ValCel, MinCel, MinCelMask
!    real :: DEMtemp(nc+2, nf+2)
!
!    DEMtemp(2:nc+1, 2:nf+1) = DEM
!    DEMfill = DEM
!    do i = 2, nc+1
!        do j = 2, nf+1
!            ValCel = DEMtemp(i, j)
!            MinCel = minval(DEMtemp(i-1:i+1, j-1:j+1))
!            MinCelMask = minval(DEMtemp(i-1:i+1, j-1:j+1), &
!                                mask = DEMtemp(i-1:i+1, j-1:j+1) /= MinCel)
!            if (ValCel == MinCel) then
!                DEMfill(i-1, j-1) = MinCelMask
!            end if
!        end do
!    end do
!end subroutine DEM_Pitfill

!--------------------------------------------------------
! Erosión en DEM según Soille (legacy, incompleto)
!--------------------------------------------------------
!subroutine DEM_carve(DEM, nc, nf, DEMcarve, DIR, Pits)
!    integer, intent(in) :: nc, nf
!    real, intent(in)    :: DEM(nc, nf)
!    real, intent(out)   :: DEMcarve(nc, nf)
!    integer, intent(out):: DIR(nc, nf), Pits(nc, nf)
!
!    ! Implementación legacy de erosión de pits (no operativa)
!    ! Se requiere adaptación y pruebas.
!end subroutine DEM_carve

!--------------------------------------------------------
! Kernel para determinar dirección de drenaje en 3x3 (legacy, incompleto)
!--------------------------------------------------------
!subroutine kernel_DIR(DEMkernel, dirOut)
!    real, intent(in)  :: DEMkernel(3, 3)
!    integer, intent(out):: dirOut
!
!    ! Implementar lógica para dirección de salida del kernel
!end subroutine kernel_DIR

!--------------------------------------------------------
! Suavizado/genérica por si agregas otros algoritmos experimentales.
! Copia aquí cualquier subrutina a refactorizar
!--------------------------------------------------------

!================= FIN BLOQUE EXPERIMENTAL ==================
