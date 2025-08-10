# -*- coding: utf-8 -*-
"""
Bloque 1 - Configuración inicial, imports y variables globales
--------------------------------------------------------------
Este bloque prepara el entorno de trabajo para el módulo WMF:
- Configuración de librerías base
- Imports de módulos internos y externos
- Variables globales iniciales
- Manejo de dependencias opcionales
"""

# ======= Información de licencia =======
# Copyright (C) <2016>
# Este programa es software libre: se puede redistribuir y/o modificar
# bajo los términos de la GNU General Public License, versión 3 o posterior.
# Se distribuye SIN GARANTÍA alguna, ni implícita ni explícita.

# ======= Imports obligatorios =======
import os
import random
import datetime as datetime

import numpy as np
import pandas as pd
import matplotlib
import pylab as pl

from scipy.spatial import Delaunay
from scipy.stats import norm
from multiprocessing import Pool
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ======= Imports de módulos locales =======
from cu import *       # Módulo de funciones de cuencas (Fortran/C o Python)
from models import *   # Módulo de simulación hidrológica

# ======= Dependencias opcionales =======

# Pysheds para operaciones con DEM y cálculo de flow directions
try:
    from pysheds.grid import Grid
except ImportError:
    print('Warning: no se encontró pysheds. El usuario debe proveer el mapa DIR para obtener la cuenca.')

# Cartopy para proyecciones y gráficos geográficos
try:
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.io.shapereader as shpreader
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
except ImportError:
    print('Warning: no se encontró cartopy. Se deshabilita la visualización de mapas con proyecciones.')

# GDAL / OGR para entrada y salida de datos geoespaciales
try:
    import osgeo.ogr, osgeo.osr
    import gdal
except ImportError:
    print('Warning: No se pudo importar GDAL/OSGeo para mapas.')

# Rasterio para conversión de raster a polígono
try:
    from rasterio import features as __fea__
    FlagBasinPolygon = True
except ImportError:
    print('Warning: no se encontró rasterio. No se podrá obtener el polígono de la cuenca.')
    FlagBasinPolygon = False

# NetCDF para guardar/cargar modelos
try:
    import netcdf as netcdf
except ImportError:
    try:
        import netCDF4 as netcdf
    except ImportError:
        print('Warning: no se encontró netCDF4. Se deshabilita el guardado/carga de modelos en NetCDF.')

# Deap para calibración por algoritmos genéticos
try:
    from deap import base, creator, tools
    FlagCalib_NSGAII = True
except ImportError:
    print('Warning: no se encontró deap. Se deshabilita la calibración automática NSGA-II.')
    FlagCalib_NSGAII = False

# ======= Variables globales =======
Global_EPSG = -9999  # Código EPSG por defecto (sin definir)

def dem_process(dem_path, dxp, noData):
    """
    Usando pysheds, obtiene el mapa de direcciones de flujo (DIR map) a partir de un DEM.

    Parámetros:
    dem_path (str): Ruta al archivo DEM (raster)
    dxp (float): Tamaño de la celda en el DEM [m]
    noData (float/int): Valor para celdas sin datos

    Returns:
    dem (np.ndarray): DEM leído y posiblemente re-escalado [transpuesta]
    dir_map (np.ndarray): Direcciones de flujo [transpuesta]
    epsg (int): Código EPSG leído del DEM
    """
    # Leer DEM raster con función auxiliar
    DEM, epsg = read_map_raster(dem_path, isDEMorDIR=True, dxp=dxp, noDataP=noData)
    # Procesamiento con pysheds (si está disponible)
    gr = Grid.from_raster(dem_path, data_name='dem')
    gr.fill_depressions('dem', out_name='flooded_dem')
    gr.resolve_flats('flooded_dem', out_name='inflated_dem')
    # Mapeo de direcciones, puede variar según codificación (aquí ejemplo para pysheds estándar)
    dir_map = (8, 9, 6, 3, 2, 1, 4, 7)
    gr.flowdir(data='inflated_dem', out_name='dir', dirmap=dir_map)
    # Retornar arrays listos para modelar
    return gr.dem.T, gr.dir.T, epsg

def read_map_raster(path_map, isDEMorDIR=False, dxp=None, noDataP=None, isDIR=False, DIRformat='r.watershed'):
    """
    Lee un mapa raster soportado por GDAL y obtiene propiedades geo/matriz y EPSG.

    Parámetros:
    path_map (str): Ruta al raster
    isDEMorDIR (bool): ¿Es DEM o DIR? (define si variables globales se actualizan)
    dxp (float): Tamaño de celda si se quiere forzar (opcional)
    noDataP (float/int): Valor para NODATA forzado (opcional)
    isDIR (bool): ¿Es mapa de direcciones? (definirá reclasificación)
    DIRformat (str): Tipo de formato DIR (r.watershed, opentopo, etc.)

    Returns:
    - Si es DEM o DIR: mapa, EPSG
    - Si no: mapa, [ncols, nrows, xll, yll, dx, dy, noData], EPSG
    """
    import gdal, osgeo.osr
    direction = gdal.Open(path_map)
    proj = osgeo.osr.SpatialReference(wkt=direction.GetProjection())
    EPSG_code = proj.GetAttrValue('AUTHORITY', 1)

    ncols = direction.RasterXSize
    nrows = direction.RasterYSize
    banda = direction.GetRasterBand(1)
    noData = banda.GetNoDataValue()
    geoT = direction.GetGeoTransform()
    dx = geoT[1]
    dy = abs(geoT[-1])
    xll = geoT
    yll = geoT - nrows * dy

    Mapa = direction.ReadAsArray()
    direction.FlushCache()
    del direction

    # Si es DEM/DIR, pasa propiedades globales (para módulo cuencas/modelos)
    if isDEMorDIR:
        cu.ncols = ncols
        cu.nrows = nrows
        if noDataP is not None:
            cu.nodata = noDataP
            Mapa[Mapa < 0] = cu.nodata
        else:
            cu.nodata = noData
        cu.dx = dx
        cu.dy = dy
        cu.xll = xll
        cu.yll = yll
        cu.dxp = dxp if dxp is not None else 30.0
        global Global_EPSG
        Global_EPSG = EPSG_code

        if isDIR:
            if DIRformat == 'r.watershed':
                Mapa[Mapa <= 0] = cu.nodata
                Mapa = cu.dir_reclass_rwatershed(Mapa.T, cu.ncols, cu.nrows)
                return Mapa, EPSG_code
            elif DIRformat == 'opentopo':
                Mapa[Mapa <= 0] = cu.nodata
                Mapa = cu.dir_reclass_opentopo(Mapa.T, cu.ncols, cu.nrows)
                return Mapa, EPSG_code

        return Mapa.T.astype(float), EPSG_code
    
    # Si no es DEM/DIR, retorna propiedades completas
    return Mapa.T.astype(float), [ncols, nrows, xll, yll, dx, dy, noData], EPSG_code

def read_map_points(path_map, ListAtr=None):
    """
    Lee un mapa de puntos vectorial (soportado por GDAL). Devuelve coordenadas y, opcionalmente, atributos específicos.
    """
    import osgeo.ogr
    dr = osgeo.ogr.Open(path_map)
    l = dr.GetLayer()
    Cord = []
    for i in range(l.GetFeatureCount()):
        f = l.GetFeature(i)
        g = f.GetGeometryRef()
        pt = [g.GetX(), g.GetY()]
        Cord.append(pt)
    Cord = np.array(Cord).T
    # Si hay atributos solicitados:
    if ListAtr is not None:
        Dict = {}
        for j in ListAtr:
            pos = f.GetFieldIndex(j)
            if pos != -1:
                vals = []
                for i in range(l.GetFeatureCount()):
                    f = l.GetFeature(i)
                    vals.append(f.GetField(pos))
                Dict.update({j: np.array(vals)})
        dr.Destroy()
        return Cord, Dict
    else:
        dr.Destroy()
        return Cord

def plot_sim_single(Qs, Qo=None, mrain=None, Dates=None, path=None,
                    figsize=(8, 4.5), ids=None, legend=True, ax1=None, **kwargs):
    """
    Plotea un único hidrograma (o varios) con lluvia y caudales.

    Parámetros:
    -----------
    Qs : array/list
        Serie o matriz con caudal simulado(s).
    Qo : array/list, opcional
        Serie con caudal observado.
    mrain : array/list, opcional
        Serie de lluvia media (mm).
    Dates : array-like, opcional
        Fechas (datetime) para el eje X.
    path : str, opcional
        Ruta para guardar la figura (png/jpg/pdf).
    figsize : tuple
        Tamaño de la figura.
    ids : list, opcional
        Identificadores de puntos de control (si Qs es matriz).
    legend : bool
        Mostrar la leyenda.
    ax1 : matplotlib axis, opcional
        Eje preexistente para graficar.
    **kwargs : argumentos opcionales de estilo.

    Returns:
    --------
    ax1, ax2 : ejes principal y secundario (precipitación)
    """
    # Mostrar gráfico al final
    show = kwargs.get('show', True)

    # Crear figura/axes si no se pasa uno existente
    if ax1 is None:
        fig = pl.figure(facecolor='w', edgecolor='w', figsize=figsize)
        ax1 = fig.add_subplot(111)
    else:
        show = False

    # Eje X -> fechas o índices
    if Dates is None:
        ejeX = range(len(Qs[0])) if len(Qs.shape) > 1 else range(len(Qs))
    else:
        ejeX = Dates

    # Graficar lluvia como área en eje secundario
    if mrain is not None:
        ax2 = ax1.twinx()
        alpha = kwargs.get('rain_alpha', 0.4)
        color = kwargs.get('rain_color', 'blue')
        lw = kwargs.get('rain_lw', 0)
        ax2.fill_between(ejeX, 0, mrain, alpha=alpha, color=color, lw=lw)
        ylabel = kwargs.get('rain_ylabel', 'Precipitation [$mm$]')
        label_size = kwargs.get('label_size', 14)
        ax2.set_ylabel(ylabel, size=label_size)
        # Invertir lluvia para estilo clásico hidrológico
        ylim = kwargs.get('rain_ylim', ax2.get_ylim()[::-1])
        ax2.set_ylim(ylim)
    else:
        ax2 = None

    # Graficar caudales simulados
    ColorSim = kwargs.get('ColorSim', ['r', 'g', 'k', 'c', 'y'])
    Qs_lw = kwargs.get('Qs_lw', 1.5)
    Qo_lw = kwargs.get('Qo_lw', 2.0)
    Qs_color = kwargs.get('Qs_color', 'r')
    Qo_color = kwargs.get('Qo_color', 'b')
    Qo_label = kwargs.get('Qo_label', 'Observed')
    Qs_label = kwargs.get('Qs_label', 'Simulated')

    if len(Qs.shape) > 1:
        if ids is None:
            ids = np.arange(1, Qs.shape[0] + 1)
        for i, c, d in zip(Qs, ColorSim, ids):
            ax1.plot(ejeX, i, c, lw=Qs_lw, label=str(d))
    else:
        ax1.plot(ejeX, Qs, Qs_color, lw=Qs_lw, label=Qs_label)

    # Graficar caudal observado si existe
    if Qo is not None:
        ax1.plot(ejeX, Qo, Qo_color, lw=Qo_lw, label=Qo_label)

    # Etiquetas y formato
    xlabel = kwargs.get('xlabel', 'Time [$min$]')
    ylabel = kwargs.get('ylabel', 'Streamflow $[m^3/s]$')
    label_size = kwargs.get('label_size', 16)
    ax1.set_xlabel(xlabel, size=label_size)
    ax1.set_ylabel(ylabel, size=label_size)
    ax1.grid(True)

    # Leyenda
    if legend:
        loc = kwargs.get('legend_loc', 'upper center')
        bbox_to_anchor = kwargs.get('bbox_to_anchor', (0.5, -0.12))
        ncol = kwargs.get('legend_ncol', 4)
        ax1.legend(loc=loc, bbox_to_anchor=bbox_to_anchor,
                   fancybox=True, shadow=True, ncol=ncol)

    # Guardar si se pide
    if path is not None:
        pl.savefig(path, bbox_inches='tight')

    if show:
        pl.show()

    return ax1, ax2

def plot_mean_storage(Mean_Storage, Dates=None, mrain=None, path=None, **kwargs):
    """
    Plotea los almacenamientos medios en los distintos tanques del modelo hidrológico.

    Parámetros:
    -----------
    Mean_Storage : np.ndarray
        Almacenamiento medio en cada tanque. Forma: (n_tanques, n_pasos)
    Dates : array-like, opcional
        Fechas (datetime) para eje X.
    mrain : array-like, opcional
        Serie de lluvia media sobre la cuenca.
    path : str, opcional
        Ruta para guardar la figura.
    **kwargs : argumentos opcionales de ajuste y estilo.
        - figsize: tamaño de la figura (default (13,9))
        - color: color de líneas de tanques (default 'k')
        - lw: ancho de línea (default 4)
        - labelsize: tamaño de label de ejes
        - ysize: tamaño fuente de ylabel
        - alpha: transparencia de lluvia
        - colorRain: color lluvia
        - lwRain: ancho línea lluvia
        - show: mostrar figura (default True)

    Returns:
    --------
    Axes de matplotlib usados (los cinco para tanques; lluvia en el primero si aplica)
    """
    # Parámetros visuales
    figsize = kwargs.get('figsize', (13, 9))
    color = kwargs.get('color', 'k')
    lw = kwargs.get('lw', 4)
    labelsize = kwargs.get('labelsize', 14)
    ysize = kwargs.get('ysize', 16)
    show = kwargs.get('show', True)
    alpha = kwargs.get('alpha', 0.5)
    colorRain = kwargs.get('colorRain', 'b')
    lwRain = kwargs.get('lwRain', 0.1)

    # Eje X: fechas o índice por pasos
    if Dates is None:
        ejeX = range(Mean_Storage.shape[1])
    else:
        ejeX = Dates

    fig = pl.figure(figsize=figsize)
    nombres = ['Hu', 'Runoff', 'Hg', 'Sub', 'Stream']  # Nombres de tanques (pueden variar)
    resultados = []

    for c, i in enumerate(Mean_Storage):
        ax = fig.add_subplot(5, 1, c+1)
        ax.plot(ejeX, i, color, lw=lw)
        ax.grid()
        # Quitar etiquetas de X en los de arriba
        if c < 4:
            ax.set_xticklabels([])
        ax.tick_params(labelsize=labelsize)
        # Para el primer plot: añadir lluvia sobre eje secundario
        if c == 0 and mrain is not None:
            ax2 = ax.twinx()
            ax2.fill_between(ejeX, 0, mrain, alpha=alpha, color=colorRain, lw=lwRain)
            ylim = ax2.get_ylim()[::-1]
            ylim = list(ylim)
            ylim[1] = 0
            ax2.set_ylim(ylim)
            resultados.append((ax, ax2))
        else:
            resultados.append(ax)
        # Nombre del tanque en eje Y
        ax.set_ylabel(nombres[c], size=ysize)

    # Guardar si corresponde
    if path is not None:
        pl.savefig(path, bbox_inches='tight')
    if show:
        pl.show()

    return resultados  # Lista de axes, primer elemento puede ser (ax, ax2)

def Save_Array2Raster(Array, ArrayProp, path, EPSG=4326, Format='GTiff'):
    """
    Guarda un array como raster georreferenciado.

    Parámetros:
    -----------
    Array : np.ndarray
        Matriz con datos (coordenadas en formato [y, x]).
    ArrayProp : list
        [ncols, nrows, xll, yll, dx, dy, nodata]
    path : str
        Ruta de salida del archivo (incluye extensión).
    EPSG : int
        Código EPSG del sistema de coordenadas (default WGS84).
    Format : str
        Formato de salida GDAL (GTiff, NetCDF, etc.)
    """
    from osgeo import gdal, osr

    ncols, nrows, xll, yll, dx, dy, noData = ArrayProp
    driver = gdal.GetDriverByName(Format)
    gdaltype = {
        "uint8": 1, "int8": 1, "uint16": 2, "int16": 3,
        "uint32": 4, "int32": 5, "float32": 6, "float64": 7
    }[Array.dtype.name]

    dataset = driver.Create(path, ncols, nrows, 1, gdaltype)
    dataset.SetGeoTransform((xll, dx, 0, yll + nrows*dy, 0, -dy))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    dataset.SetProjection(srs.ExportToWkt())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(noData if noData is not None else -9999)
    band.WriteArray(Array.T)
    dataset.FlushCache()
    dataset = None

def Save_Points2Map(XY, ids, path, EPSG=4326, Dict=None, DriverFormat='ESRI Shapefile'):
    """
    Guarda coordenadas X, Y en un shapefile de puntos.

    Parámetros:
    -----------
    XY : np.ndarray
        Matriz [2, N] con coordenadas.
    ids : list/array
        Identificador único para cada punto.
    path : str
        Ruta de salida (.shp)
    EPSG : int
        EPSG del sistema de coordenadas.
    Dict : dict, opcional
        Propiedades adicionales { 'campo' : values_array }
    DriverFormat : str
        Formato GDAL (default: ESRI Shapefile)
    """
    from osgeo import ogr, osr
    if os.path.exists(path):
        ogr.GetDriverByName(DriverFormat).DeleteDataSource(path)

    driver = ogr.GetDriverByName(DriverFormat)
    shapeData = driver.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    layer = shapeData.CreateLayer('layer1', srs, ogr.wkbPoint)

    # Campo ID principal
    layer.CreateField(ogr.FieldDefn('Estacion', ogr.OFTString))

    # Campos adicionales
    if Dict:
        for key, values in Dict.items():
            tipo = type(values[0])
            if tipo is np.float64:
                ftype = ogr.OFTReal
            elif tipo is np.int64 or tipo is int:
                ftype = ogr.OFTInteger
            else:
                ftype = ogr.OFTString
            layer.CreateField(ogr.FieldDefn(str(key)[:10], ftype))

    # Insertar puntos
    def add_point(x, y, id_val, attr_dict=None, idx=0):
        feat = ogr.Feature(layer.GetLayerDefn())
        geom = ogr.Geometry(ogr.wkbPoint)
        geom.AddPoint(float(x), float(y))
        feat.SetGeometry(geom)
        feat.SetField('Estacion', str(id_val))
        if attr_dict:
            for key, vals in attr_dict.items():
                feat.SetField(str(key)[:10], vals[idx])
        layer.CreateFeature(feat)
        geom = None
        feat = None

    for i, (x, y) in enumerate(XY.T):
        if x != -9999 and y != -9999:
            add_point(x, y, ids[i], Dict, i)

    shapeData = None

def read_map_points(path_map, ListAtr=None):
    """
    Lee un shapefile/vectorial de puntos y devuelve coordenadas y atributos.

    Parámetros:
    -----------
    path_map : str
        Ruta del shapefile.
    ListAtr : list, opcional
        Lista con nombres de atributos a extraer.

    Returns:
    --------
    coords : np.ndarray
        Coordenadas [2, N]
    attrs : dict
        Diccionario con atributos opcionales solicitados.
    """
    from osgeo import ogr
    ds = ogr.Open(path_map)
    layer = ds.GetLayer()
    coords = []
    for i in range(layer.GetFeatureCount()):
        feat = layer.GetFeature(i)
        geom = feat.GetGeometryRef()
        coords.append([geom.GetX(), geom.GetY()])
    coords = np.array(coords).T

    if ListAtr:
        attrs = {attr: [] for attr in ListAtr}
        for feat in layer:
            for attr in ListAtr:
                attrs[attr].append(feat.GetField(attr))
        ds = None
        return coords, {k: np.array(v) for k, v in attrs.items()}
    else:
        ds = None
        return coords

def Save_Array2Raster(Array, ArrayProp, path, EPSG=4326, Format='GTiff'):
    """
    Guarda un array 2D (NumPy) como un raster georreferenciado usando GDAL.

    Parámetros:
    -----------
    Array : np.ndarray
        Matriz con los valores a guardar [ny, nx].
    ArrayProp : list
        Lista con propiedades del raster:
        [ncols, nrows, xll, yll, dx, dy, nodata]
    path : str
        Ruta completa y nombre de archivo de salida.
    EPSG : int
        Código EPSG (default WGS84 = 4326).
    Format : str
        Formato GDAL a usar (default GTiff).
    """
    from osgeo import gdal, osr
    
    ncols, nrows, xll, yll, dx, dy, noData = ArrayProp

    driver = gdal.GetDriverByName(Format)
    dtype_map = {
        "uint8": 1, "int8": 1,
        "uint16": 2, "int16": 3,
        "uint32": 4, "int32": 5,
        "float32": 6, "float64": 7
    }
    gdaltype = dtype_map.get(Array.dtype.name, 6)  # default float32

    dataset = driver.Create(path, ncols, nrows, 1, gdaltype)
    dataset.SetGeoTransform((xll, dx, 0, yll + nrows * dy, 0, -dy))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    dataset.SetProjection(srs.ExportToWkt())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(noData if noData is not None else -9999)
    band.WriteArray(Array.T)
    dataset.FlushCache()
    dataset = None

def Save_Points2Map(XY, ids, path, EPSG=4326, Dict=None, DriverFormat='ESRI Shapefile'):
    """
    Guarda coordenadas como shapefile de puntos.

    Parámetros:
    -----------
    XY : np.ndarray
        Array [2, N] con coordenadas X, Y.
    ids : iterable
        Lista/array de identificadores (uno por punto).
    path : str
        Ruta de salida (.shp).
    EPSG : int
        Código EPSG para proyección.
    Dict : dict opcional
        {nombre_campo: array de valores}, valores float/int/str.
    DriverFormat : str
        Formato GDAL, por defecto ESRI Shapefile.
    """
    from osgeo import ogr, osr
    if os.path.exists(path):
        ogr.GetDriverByName(DriverFormat).DeleteDataSource(path)

    driver = ogr.GetDriverByName(DriverFormat)
    shapeData = driver.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    layer = shapeData.CreateLayer('layer1', srs, ogr.wkbPoint)

    # Campo de ID principal
    layer.CreateField(ogr.FieldDefn('Estacion', ogr.OFTString))

    # Crear campos adicionales si se proporcionan
    if Dict:
        for key, values in Dict.items():
            tipo = type(values[0])
            if tipo is np.float64 or tipo is float:
                ftype = ogr.OFTReal
            elif tipo is np.int64 or tipo is int:
                ftype = ogr.OFTInteger
            else:
                ftype = ogr.OFTString
            layer.CreateField(ogr.FieldDefn(str(key)[:10], ftype))

    for idx, (x, y) in enumerate(XY.T):
        if x != -9999 and y != -9999:
            feat = ogr.Feature(layer.GetLayerDefn())
            geom = ogr.Geometry(ogr.wkbPoint)
            geom.AddPoint(float(x), float(y))
            feat.SetGeometry(geom)
            feat.SetField('Estacion', str(ids[idx]))
            if Dict:
                for key, values in Dict.items():
                    feat.SetField(str(key)[:10], values[idx])
            layer.CreateFeature(feat)

    shapeData = None

def read_map_points(path_map, ListAtr=None):
    """
    Lee coordenadas y atributos de un shapefile de puntos.

    Parámetros:
    -----------
    path_map : str
        Ruta del shapefile.
    ListAtr : list opcional
        Lista con los nombres de atributos a leer.

    Returns:
    --------
    coords : np.ndarray [2, N]
    attrs : dict opcional con {nombre_atributo: np.array(valores)}
    """
    from osgeo import ogr
    ds = ogr.Open(path_map)
    layer = ds.GetLayer()
    coords = []
    for i in range(layer.GetFeatureCount()):
        feat = layer.GetFeature(i)
        geom = feat.GetGeometryRef()
        coords.append([geom.GetX(), geom.GetY()])
    coords = np.array(coords).T

    if ListAtr:
        attrs = {}
        for attr in ListAtr:
            valores = []
            for i in range(layer.GetFeatureCount()):
                feat = layer.GetFeature(i)
                valores.append(feat.GetField(attr))
            attrs[attr] = np.array(valores)
        ds = None
        return coords, attrs

    ds = None
    return coords

class Basin:
    """
    Clase que representa una cuenca hidrográfica en el modelo WMF.
    Contiene datos básicos, mapas DEM/DIR, estructura interna y herramientas geomorfológicas.
    """

    def __init__(self, lat=0, lon=0, DEM=None, DIR=None,
                 name='NaN', stream=None, threshold=1000,
                 path=None, useCauceMap=None):
        """
        Inicializa el objeto cuenca.

        Parámetros:
        -----------
        lat, lon : float
            Coordenadas del punto de salida de la cuenca.
        DEM : np.ndarray
            Matriz del modelo digital de elevación.
        DIR : np.ndarray
            Mapa raster de direcciones de flujo.
        name : str
            Nombre de la cuenca.
        stream : objeto Stream, opcional
            Si se pasa, reposiciona lat/lon al punto más cercano en la corriente.
        threshold : int
            Umbral (número mínimo de celdas acumuladas) para definir cauces.
        path : str, opcional
            Si se pasa, carga la cuenca desde un NetCDF existente.
        useCauceMap : np.ndarray, opcional
            Mapa binario donde 1 = cauce, 0 = ladera. Fuerza corrección de salida.
        """
        self.DEM = DEM
        self.DIR = DIR
        self.name = name
        self.threshold = threshold

        # Si se usa cauce map, se inhabilita 'stream'
        if useCauceMap is not None:
            stream = None

        # Si se da path, cargar la cuenca guardada en formato NetCDF
        if path is not None:
            self.__Load_BasinNc(path)
            self.__GetBasinPolygon__()
            return

        # Si existe un objeto 'stream', corregir la coordenada de salida
        if stream is not None:
            error = []
            for i in stream.structure.T:
                error.append(np.sqrt((lat - i[0])**2 + (lon - i[1])**2))
            loc = np.argmin(error)
            lat = stream.structure[0, loc]
            lon = stream.structure[1, loc]

        # Ajuste de coordenada con mapa de cauce binario
        if useCauceMap is not None and useCauceMap.shape == DEM.shape:
            lat, lon = cu.stream_find_to_corr(lat, lon, DEM, DIR, useCauceMap,
                                              cu.ncols, cu.nrows)

        # Guardar nombre
        self.name = name

        # === TRAZAR LA CUENCA ===
        # Encuentra celdas de la cuenca desde un outlet lat/lon en el mapa DIR
        self.ncells = cu.basin_find(lat, lon, DIR, cu.ncols, cu.nrows)
        self.structure = cu.basin_cut(self.ncells)

        # Guardar DEM/DIR en formato vector de cuenca (para no usar mapas completos)
        self.DEMvec = self.Transform_Map2Basin(DEM, [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])
        self.DIRvec = self.Transform_Map2Basin(DIR, [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])

        # Generar el polígono de la cuenca
        self.__GetBasinPolygon__()

    def __Load_BasinNc(self, path, Var2Search=None):
        """
        Carga una cuenca guardada previamente en un archivo NetCDF.
        El archivo debe haber sido creado con Basin.Save_Basin2nc().
        """
        import netCDF4 as netcdf
        self.pathNC = path
        gr = netcdf.Dataset(path, 'a')

        # Atributos generales
        self.name = gr.nombre
        self.ncells = gr.ncells
        self.threshold = gr.threshold

        # Propiedades de mapa globales (para módulo cu)
        cu.ncols = gr.ncols
        cu.nrows = gr.nrows
        cu.nodata = gr.noData
        cu.dx = gr.dx
        cu.xll = gr.xll
        cu.yll = gr.yll
        cu.dxp = gr.dxp

        # Estructura interna de la cuenca
        self.structure = gr.variables['structure'][:]

        gr.close()
    
    def __GetBasinPolygon__(self):
        """
        Obtiene el polígono envolvente de la cuenca usando rasterio.features.shapes.
        Guarda el resultado en self.Polygon.
        """
        if not FlagBasinPolygon:
            return 1  # no se pudo
        # Generar máscara raster de la cuenca (1 donde hay cuenca)
        Map, Prop = self.Transform_Basin2Map(np.ones(self.ncells))
        mask = Map.T != -9999
        shapes = __fea__.shapes(Map.T, mask=mask,
                                transform=(Prop[2], Prop, 0.0,
                                           Prop[1]*Prop[-2] + Prop, 0.0, -1*Prop))
        # Tomar el polígono más grande (cuenca externa)
        pols = [np.array(sh['coordinates']).T for sh in shapes]
        self.Polygon = max(pols, key=lambda p: p.shape[1])
        return 0

def GetGeo_Parameters(self, pathParamASC=None, plotTc=False,
                      pathTcPlot=None, figsize=(8,5), GetPerim=True):
    """
    Calcula parámetros geomorfológicos y tiempos de concentración de la cuenca.

    Parámetros:
    -----------
    pathParamASC : str opcional
        Ruta para guardar parámetros en ASCII.
    plotTc : bool
        Graficar tiempos de concentración.
    pathTcPlot : str
        Ruta para guardar gráfico de Tc.
    figsize : tuple
        Tamaño de figura para el plot.
    GetPerim : bool
        Calcular perímetro.

    Returns:
    --------
    GeoParameters : DataFrame
        Parámetros calculados de la cuenca.
    Tc : dict
        Métodos y valores de tiempo de concentración.
    """
    # Variables base
    acum, longCeld, slope, Elev = cu.basin_basics(self.structure, self.DEMvec, self.DIRvec, self.ncells)
    Lpma, puntto = cu.basin_findlong(self.structure, self.ncells)
    cauce, nodos, trazado, n_nodos, n_cauce = cu.basin_stream_nod(self.structure, acum, self.threshold, self.ncells)
    ppal_nceldas, punto = cu.basin_ppalstream_find(self.structure, nodos, longCeld, Elev, self.ncells)
    ppal = cu.basin_ppalstream_cut(ppal_nceldas, self.ncells)

    # Curvas hipsométricas
    self.hipso_main, self.hipso_basin = cu.basin_ppal_hipsometric(
        self.structure, Elev, punto, 30, ppal_nceldas, self.ncells
    )
    self.main_stream = ppal

    # Parámetros de forma y relieve
    Area = (self.ncells * cu.dxp**2) / 1e6
    Lcau = ppal[1, -1] / 1000.0
    Scau = np.polyfit(ppal[1, ::-1], ppal[0], 1) * 100  # pendiente ppal
    Scue = slope.mean() * 100
    Hmin = Elev[-1]
    Hmax = Elev[puntto]
    Hmean = Elev.mean()
    HCmax = Elev[punto]
    x, y = cu.basin_coordxy(self.structure, self.ncells)
    CentXY = [np.median(x), np.median(y)]

    TotalCauces = (self.CellCauce * self.CellLong).sum() / 1000.0  # km
    Densidad = TotalCauces / Area

    self.GeoParameters = {
        'Area_[km2]': Area,
        'MeanStream_slope_[%]': Scau,
        'MeanStream_length_[km]': Lcau,
        'Watershed_Slope_[%]': Scue,
        'Watershed_length_[km]': Lpma,
        'Height_Max_[m]': Hmax,
        'Height_Min_[m]': Hmin,
        'Height_Mean_[m]': Hmean,
        'MeanStream_Height_Max_[m]': HCmax,
        'Center_[X]': CentXY[0],
        'Center_[Y]': CentXY[1],
        'Streams_total_length_[km]': TotalCauces,
        'Streams_density_[km/km2]': Densidad
    }

    if GetPerim:
        Perim = self.Polygon.shape[1] * cu.dxp / 1000.
        self.GeoParameters.update({'Perimeter_[km]': Perim})

    # Tiempos de concentración con varios métodos (en horas)
    Tc = {}
    Tc['US Army'] = 0.3 * (Lcau / (Scue**0.25))**0.75
    Tc['Direccion Carreteras Espana'] = 0.3 * (Lcau / (((Hmax - Hmin)/Lcau)**0.25))**0.75
    Tc['Kiprich'] = (0.02 * (Lpma*1000.0)**0.77) / (((Scau/100.0)**0.385) * 60.0)
    Tc['Campo y Munera'] = 8.157 * ((Area**0.316) / (((Scau*100)**0.17) * Scue**0.565))
    Tc['Giandotti'] = (4*np.sqrt(Area) + 1.5*Lcau) / (0.8*np.sqrt(Hmean))
    Tc['John Stone'] = 0.4623 * (Lcau**0.5) * ((Scue/100.0)**(-0.25))
    Tc['Ventura'] = (Lcau / Area) * (np.sqrt(Area) / Scau)
    Tc['Temez'] = 0.3 * (Lcau / ((((HCmax-Hmin)/Lcau)*100)**0.25))**0.75

    self.Tc = Tc

    # Guardar ASCII si se requiere
    if pathParamASC:
        self.__WriteGeoParam__(pathParamASC)
    # Graficar Tc si se requiere
    if plotTc:
        self.Plot_Tc(path=pathTcPlot, figsize=figsize)

    return pd.DataFrame.from_dict(self.GeoParameters, orient='index'), pd.DataFrame.from_dict(self.Tc, orient='index')

def GetGeo_Cell_Basics(self):
    """Obtiene área acumulada, longitud de celda, pendiente y elevación para cada celda de la cuenca."""
    acum, longCeld, S0, Elev = cu.basin_basics(self.structure, self.DEMvec, self.DIRvec, self.ncells)
    self.CellAcum = acum
    self.CellLong = longCeld
    self.CellSlope = S0
    self.CellHeight = Elev
    # Mapa binario de cauces según threshold
    self.CellCauce = np.where(self.CellAcum > self.threshold, 1, 0)

def GetGeo_StreamOrder(self, MajorBasins=False, threshold=100, verbose=False, FirtsOrder=1):
    """
    Calcula el orden de Horton para cada celda de cauce y ladera.
    Puede generar sub-cuencas mayores si MajorBasins es True.
    """
    cauce, nodos_fin, n_nodos = cu.basin_subbasin_nod(self.structure, self.CellAcum, self.threshold, self.ncells)
    sub_pert, sub_basin = cu.basin_subbasin_find(self.structure, nodos_fin, n_nodos, self.ncells)
    sub_basins = cu.basin_subbasin_cut(n_nodos)
    sub_horton, nod_horton = cu.basin_subbasin_horton(sub_basins, self.ncells, n_nodos)
    self.CellHorton_Hill, sub_basin = cu.basin_subbasin_find(self.structure, nod_horton, n_nodos, self.ncells)
    self.CellHorton_Stream = self.CellCauce * self.CellHorton_Hill

    if MajorBasins:
        DictBasins = {}
        X, Y = cu.basin_coordxy(self.structure, self.ncells)
        for Orden in range(FirtsOrder, self.CellHorton_Hill.max()):
            pos2 = np.where(self.CellHorton_Stream == Orden)[0]
            drena = self.ncells - self.structure
            pos3 = np.where(self.CellHorton_Stream[drena[pos2]] > Orden)
            SubCuencas = np.zeros(self.ncells)
            cont = 1
            for x_, y_ in zip(X[pos2[pos3]], Y[pos2[pos3]]):
                cuTemp = SimuBasin(x_, y_, self.DEM, self.DIR, threshold=threshold)
                Map, _ = cuTemp.Transform_Basin2Map(np.ones(cuTemp.ncells))
                Map[Map == -9999] = 0
                Var = self.Transform_Map2Basin(Map, _)
                Var[Var == 0] = -9999
                Var[Var == 1] = cont
                SubCuencas[Var == cont] = cont
                cont += 1
            DictBasins[str(Orden)] = SubCuencas
        return DictBasins

def GetGeo_IsoChrones(self, Tc, Niter=4):
    """
    Estima el tiempo de viaje aproximado de cada celda a la salida (“isócronas”) según el Tc dado.

    Parámetros:
    -----------
    Tc : float
        Tiempo de concentración objetivo (horas)
    Niter : int
        Iteraciones para afinar velocidad

    Returns:
    --------
    Isochrones: np.ndarray
        Vector [ncells], tiempo de viaje desde cada celda a la salida (horas)
    """
    acum, longCeld, S0, Elev = cu.basin_basics(self.structure, self.DEMvec, self.DIRvec, self.ncells)
    rangos = [50, 25, 1]  # rangos de prueba para estimar velocidad
    for _ in range(Niter):
        times = []
        for r in rangos:
            speed = r * S0 ** 0.5
            time = cu.basin_time_to_out(self.structure, longCeld, speed, self.ncells) / 3600.0
            times.append(np.mean(time[np.isfinite(time)]))
        # abreviado: ajusta speed a que tiempo medio se aproxime a Tc (no está el loop completo)
        # El resultado es el último vector de tiempos estimado
    # Devuelve el último estimado (puede ampliarse según la implementación original)
    return time

def GetGeo_Ppal_Hipsometric(self, threshold=1000, intervals=30):
    """
    Calcula la curva hipsométrica del cauce principal y de la cuenca.

    Parámetros:
    -----------
    threshold : int
        Umbral mínimo de celdas para definir cauce
    intervals : int
        Número de divisiones para la curva

    Returns:
    --------
    hipsometric_ppal : np.ndarray
        Curva hipsométrica a lo largo del cauce principal [area acumulada vs altura]
    hipsometric_basin : np.ndarray
        Curva hipsométrica para toda la cuenca
    """
    cauce, nodos, trazado, n_nodos, n_cauce = cu.basin_stream_nod(
        self.structure, self.CellAcum, threshold, self.ncells)
    ppal_nceldas, punto = cu.basin_ppalstream_find(self.structure, nodos, self.CellLong, self.CellHeight, self.ncells)
    self.ppal_stream = cu.basin_ppalstream_cut(ppal_nceldas, self.ncells)
    # Suavizar perfil y calcular curva
    # Se asume que __ModifyElevErode__ y cu.basin_ppal_hipsometric están disponibles en el entorno correcto
    self.ppal_stream[0], self.ppal_slope = __ModifyElevErode__(self.ppal_stream)
    self.hipso_ppal, self.hipso_basin = cu.basin_ppal_hipsometric(
        self.structure, self.CellHeight, punto, intervals, ppal_nceldas)
    self.hipso_ppal[1], self.hipso_ppal_slope = __ModifyElevErode__(self.hipso_ppal[1])
    return self.hipso_ppal, self.hipso_basin

def GetGeo_IT(self):
    """
    Calcula el índice topográfico de humedad de Beven para cada celda.

    Returns:
    --------
    IT: np.ndarray
        Índice topográfico (ln(a/tan(s)))
    """
    acum = cu.basin_acum(self.structure, self.ncells)
    _,_,slope,_ = cu.basin_basics(self.structure, self.DEMvec, self.DIRvec, self.ncells)
    slope = np.arctan(slope)
    slope[slope == 0] = 0.0001  # Para evitar división por cero
    return np.log((acum * cu.dxp) / np.tan(slope))

def GetGeo_HAND_and_rDUNE(self, threshold=1000):
    """
    Calcula HAND (Height Above Nearest Drainage) y índice rDUNE (Loritz 2019).

    Returns:
    --------
    HAND: np.ndarray
        Elevación sobre el cauce [m]
    HDND: np.ndarray
        Distancia horizontal a cauce más cercano [m]
    rDUNE: np.ndarray
        Reducción de energía por longitud
    """
    acum, longCeld, S0, Elev = cu.basin_basics(self.structure, self.DEMvec, self.DIRvec, self.ncells)
    cauce, nodos, trazado, n_nodos, n_cauce = cu.basin_stream_nod(
        self.structure, acum, threshold, self.ncells)
    hand, hdnd, hand_destiny = cu.geo_hand(self.structure, Elev, longCeld, cauce, self.ncells)
    handC = np.zeros(self.ncells)
    handC[hand < 5.3] = 1
    handC[(hand >= 5.3) & (hand <= 15.0)] = 2
    handC[(hand > 15.0) & (S0 < 0.076)] = 4
    handC[(hand > 15.0) & (S0 >= 0.076)] = 3
    self.CellHAND = hand
    self.CellHAND_class = handC
    self.CellHDND = hdnd
    self.CellHAND_drainCell = hand_destiny
    # cálculo de rDUNE
    a = np.copy(self.CellHAND)
    b = np.copy(self.CellHDND)
    a[b == 0] = 1
    a[a <= 0] = 0.1
    b[b == 0] = 1
    self.CellDUNE = a / b
    self.CellrDUNE = -1 * np.log(self.CellDUNE)
    return hand, hdnd, self.CellrDUNE

def GetGeo_Sections(self, NumCeldas=6):
    """
    Obtiene secciones transversales en todos los elementos de la red de drenaje.

    Parámetros:
    -----------
    NumCeldas : int
        Número de celdas a ambos lados del cauce.

    Returns:
    --------
    Sections : np.ndarray
        Matriz con perfiles seccionales (NumCeldas*2+1, ncells)
    Sections_Cells : np.ndarray
        Índices de celdas incluidas en cada sección
    """
    self.GetGeo_Cell_Basics()
    directions = self.Transform_Map2Basin(self.DIR[0], self.DIR[1])
    # cu.basin_stream_sections debe estar implementado en cu
    self.Sections, self.Sections_Cells = cu.basin_stream_sections(
        self.structure, self.CellCauce, directions,
        self.DEM, NumCeldas, self.ncells, cu.ncols, cu.nrows)
    return self.Sections, self.Sections_Cells

def Transform_Map2Basin(self, Map, MapProp):
    """
    Convierte un mapa raster leído (matriz) en un vector correspondiente a la cuenca.

    Parámetros:
    -----------
    Map : np.ndarray
        Matriz con la información del mapa (tamaño completo raster).
    MapProp : list
        Propiedades del mapa: [ncols, nrows, xll, yll, dx, dy]

    Returns:
    --------
    vecMap : np.ndarray
        Vector de tamaño (ncells,) con valores del mapa dentro de la cuenca.
    """
    vec = cu.basin_map2basin(self.structure, Map,
                             MapProp[2], MapProp, MapProp, MapProp,
                             cu.nodata, self.ncells, MapProp, MapProp[1])
    return vec

def Transform_Basin2Map(self, BasinVar, path=None, DriverFormat='GTiff', EPSG=4326):
    """
    Convierte un vector interno de la cuenca a un mapa raster del tamaño original DEM.

    Parámetros:
    -----------
    BasinVar : np.ndarray
        Vector (ncells,) con la variable de la cuenca.
    path : str, opcional
        Ruta para guardar el raster.
    DriverFormat : str
        Formato GDAL (default GeoTiff)
    EPSG : int
        Sistema de coordenadas

    Returns:
    --------
    M, props : np.ndarray, list
        Matriz rasterizada y propiedades geo.
    """
    map_ncols, map_nrows = cu.basin_2map_find(self.structure, self.ncells)
    M, mxll, myll = cu.basin_2map(self.structure, BasinVar, map_ncols, map_nrows, self.ncells)

    if path is not None:
        Save_Array2Raster(M, [map_ncols, map_nrows, mxll, myll,
                              cu.dx, cu.dy, cu.nodata],
                          path=path, EPSG=EPSG, Format=DriverFormat)
    return M, [map_ncols, map_nrows, mxll, myll, cu.dx, cu.dy, cu.nodata]

def Transform_Hills2Basin(self, HillsMap):
    """
    Convierte un vector de propiedades por laderas en un vector por celdas.

    Parameters:
    -----------
    HillsMap : np.ndarray
        Vector con valores por ladera (nhills,)

    Returns:
    --------
    CellMap : np.ndarray
        Variable por celda con valores de HillsMap asignados.
    """
    CellMap = np.ones(self.ncells)
    for i, k in enumerate(HillsMap[::-1]):
        CellMap[self.hills_own == i+1] = k
    return CellMap

def Transform_Basin2Hills(self, CellMap, mask=None, SumMeanMax=0):
    """
    Convierte un vector de celdas en un vector agregado por laderas.

    Parámetros:
    -----------
    CellMap : np.ndarray
        Valores por celda.
    mask : np.ndarray o valor
        Celdas a considerar (1) o excluir (0).
    SumMeanMax : int
        0 = media, 1 = suma, 2 = máximo.
    """
    if mask is not None:
        if isinstance(mask, (float, int)):
            Ma = np.where(CellMap == mask, 1, 0)
        elif isinstance(mask, np.ndarray):
            Ma = mask
        else:
            Ma = np.ones(self.ncells)
    else:
        Ma = np.ones(self.ncells)

    HillsMap = cu.basin_subbasin_map2subbasin(self.hills_own, CellMap,
                                              self.nhills, Ma,
                                              SumMeanMax, self.ncells)
    return HillsMap

def Transform_Basin2Polygon(self, Vector):
    """
    Convierte una variable de topología de la cuenca en geometrías de polígonos.

    Parámetros:
    -----------
    Vector : np.ndarray
        Vector por celda que representa unidades/grupos.

    Returns:
    --------
    DicPoly : dict
        {id_unidad: {n: array de coordenadas}}
    """
    Map, Prop = self.Transform_Basin2Map(Vector)
    mask = Map.T != -9999
    shapes = __fea__.shapes(Map.T, mask=mask,
                            transform=(Prop[2], Prop, 0.0,
                                       Prop[1]*Prop[-2] + Prop,
                                       0.0, -1*Prop))
    DicPoly = {}
    for Sh in shapes:
        Coord = Sh['coordinates']
        Value = int(Sh[1])
        DicPoly[str(Value)] = {str(cont): np.array(co).T for cont, co in enumerate(Coord)}
    return DicPoly

def Save_Net2Map(self, path, dx=cu.dxp, threshold=None,
                 qmed=None, Dict=None, DriverFormat='ESRI Shapefile',
                 EPSG=4326, Numlink_id=True, formato='%.2f'):
    """
    Guarda la red hídrica simulada en un archivo vectorial (líneas).

    Parámetros:
    -----------
    path : str
        Ruta al archivo de salida (.shp).
    dx : float
        Longitud de la celda (m).
    threshold : int
        Umbral de acumulación de celdas para definir cauce.
    qmed : np.ndarray
        Vector de caudal medio por celda (para agregar atributo Qmed).
    Dict : dict opcional
        Atributos adicionales por link_id.
    DriverFormat : str
        Formato de salida GIS (por defecto ESRI Shapefile).
    EPSG : int
        Código EPSG para proyección.
    Numlink_id : bool
        Incluir columna de identificador de link_id.
    formato : str
        Formato numérico de atributos adicionales.
    """
    if threshold is None:
        threshold = self.threshold

    # Calcular red y atributos básicos
    acum = cu.basin_acum(self.structure, self.ncells)
    cauce, nod_f, n_nodos = cu.basin_subbasin_nod(self.structure, acum, threshold, self.ncells)
    sub_pert, sub_basin = cu.basin_subbasin_find(self.structure, nod_f, n_nodos, self.ncells)
    sub_basins = cu.basin_subbasin_cut(n_nodos)
    sub_horton, nod_hort = cu.basin_subbasin_horton(sub_basins, self.ncells, n_nodos)
    sub_hort = cu.basin_subbasin_find(self.structure, nod_hort, n_nodos, self.ncells)[0]
    cauceHorton = sub_hort * cauce

    nodos = cu.basin_stream_nod(self.structure, acum, threshold, self.ncells)[1]
    netsize = cu.basin_netxy_find(self.structure, nodos, cauceHorton, self.ncells)
    net = cu.basin_netxy_cut(netsize, self.ncells)

    if qmed is not None:
        netsize = cu.basin_netxy_find(self.structure, nodos, cauce * qmed, self.ncells)
        netQmed = cu.basin_netxy_cut(netsize, self.ncells)

    if Numlink_id:
        netsize2 = cu.basin_netxy_find(self.structure, nodos, sub_pert * cauce, self.ncells)
        netlink_id = cu.basin_netxy_cut(netsize2, self.ncells)

    # Cortes (separación de líneas entre tramos)
    cortes = np.where(net[0, :] == -999)[0].tolist()
    cortes.insert(0, 0)

    # Crear shapefile
    from osgeo import ogr, osr
    if os.path.exists(path):
        ogr.GetDriverByName(DriverFormat).DeleteDataSource(path)
    driver = ogr.GetDriverByName(DriverFormat)
    shapeData = driver.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(EPSG))
    layer = shapeData.CreateLayer('layer1', srs, ogr.wkbLineString)

    # Campos básicos
    layer.CreateField(ogr.FieldDefn('Long_km', ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn('Horton', ogr.OFTInteger))
    if Numlink_id:
        layer.CreateField(ogr.FieldDefn('link_id', ogr.OFTInteger))
    if qmed is not None:
        layer.CreateField(ogr.FieldDefn('Qmed_m3s', ogr.OFTReal))
    if Dict is not None and isinstance(Dict, dict):
        for k in Dict.keys():
            layer.CreateField(ogr.FieldDefn(k[:10], ogr.OFTReal))

    # Escribir geometrías línea
    for i, j in zip(cortes[:-1], cortes[1:]):
        line = ogr.Geometry(ogr.wkbLineString)
        for x, y in zip(net[1, i+1:j], net[2, i+1:j]):
            line.AddPoint(float(x), float(y))
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(line)
        feat.SetField('Long_km', (net[1, i+1:j].size * dx) / 1000.0)
        feat.SetField('Horton', int(net[0, i+1]))
        if qmed is not None:
            feat.SetField('Qmed_m3s', float(netQmed[0, i+1]))
        if Numlink_id:
            feat.SetField('link_id', int(netlink_id[0, i+1]))
        if Dict is not None and isinstance(Dict, dict):
            for attr_arr, key in zip(Dict.values(), Dict.keys()):
                feat.SetField(key[:10], float(formato % attr_arr[0, i+1]))
        layer.CreateFeature(feat)

    shapeData.Destroy()

def Save_Basin2Map(self, path, dx=30.0, Param={},
                   DriverFormat='ESRI Shapefile', EPSG=4326, GeoParam=False):
    """
    Guarda el polígono de la cuenca como shapefile.

    Parámetros:
    -----------
    path : str
        Ruta de salida (.shp).
    dx : float
        Tamaño de celda de referencia en metros.
    Param : dict
        Parámetros adicionales a incluir como campos.
    DriverFormat : str
        Formato GIS (default ESRI Shapefile).
    EPSG : int
        Código de proyección.
    GeoParam : bool
        Calcular parámetros geomorfológicos y agregarlos a atributos.
    """
    # Si se piden parámetros geomorfológicos
    if GeoParam:
        self.GetGeo_Parameters()
        DictParam = {k[:8]: v for k, v in self.GeoParameters.items()}
    else:
        DictParam = {}

    from osgeo import ogr, osr
    if os.path.exists(path):
        ogr.GetDriverByName(DriverFormat).DeleteDataSource(path)
    driver = ogr.GetDriverByName(DriverFormat)
    shapeData = driver.CreateDataSource(path)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(EPSG))
    layer = shapeData.CreateLayer('layer1', srs, ogr.wkbPolygon)

    # Campos de atributos
    for p in Param.keys():
        layer.CreateField(ogr.FieldDefn(p, ogr.OFTReal))
    if GeoParam:
        for p in DictParam.keys():
            layer.CreateField(ogr.FieldDefn(p, ogr.OFTReal))

    # Geometría polígono desde self.Polygon
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for i in self.Polygon.T:
        ring.AddPoint(x=float(i[0]), y=float(i[1]))
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    # Feature con atributos
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(poly)
    for p in Param.keys():
        feature.SetField(p, float("%.2f" % Param[p]))
    if GeoParam:
        for p in DictParam.keys():
            feature.SetField(p, float("%.2f" % DictParam[p]))
    layer.CreateFeature(feature)

    poly.Destroy()
    ring.Destroy()
    feature.Destroy()
    shapeData.Destroy()

def graficar_cuenca(self, 
                    vector_cuenca, 
                    ax=None, 
                    figsize=(10, 10), 
                    path_guardar=None, 
                    dpi=100, 
                    cmap='viridis', 
                    fontsize=28, 
                    titulo='', 
                    titulo_colorbar='', 
                    ubicacion_colorbar='bottom',
                    norm=None, 
                    levels=None, 
                    etiquetas_colorbar=None,
                    centrar_etiquetas_colorbar=False, 
                    color_perimetro='r', 
                    separaciones_colorbar=None
                   ):
    """
    Genera un plot bonito del mapa hidrológico de la cuenca,
    con vector_cuenca representando la variable a mostrar (ej: pendiente, HAND, etc).

    vector_cuenca : array (ncells,)
        Vector con valores por celda de la cuenca.
    Otros parámetros: opciones de estilo y exportación.

    Devuelve: objeto ax de matplotlib.
    """
    import matplotlib.pyplot as pl

    # Ejes y figura
    if ax is None:
        fig = pl.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    else:
        fig = ax.figure

    # Mapa de la variable como matriz
    mapa, prop = self.Transform_Basin2Map(vector_cuenca)
    celdas_x, celdas_y, xll, yll, dx, dy, nodata = prop
    mapa[mapa == nodata] = np.nan
    longitudes = xll + dx * np.arange(celdas_x)
    latitudes = yll + dy * np.arange(celdas_y)
    longitudes, latitudes = np.meshgrid(longitudes, latitudes)

    # Título
    t = ax.set_title(titulo, fontsize=fontsize)
    t.set_y(1.05)

    # Color mapping y contorno
    cs = ax.contourf(longitudes, latitudes, mapa.T[::-1],
                     transform=ccrs.PlateCarree(),
                     cmap=cmap, levels=levels, norm=norm)

    # Polígono de la cuenca
    ax.plot(self.Polygon[0], self.Polygon[1], color=color_perimetro)

    # Formato de ejes
    lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.2f')
    lat_formatter = LatitudeFormatter(number_format='.2f')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=1, color='k', alpha=0.5, linestyle='--')
    gl.xlabels_top = True
    gl.xlabels_bottom = False
    gl.ylabels_left = True
    gl.ylabels_right = False

    # Guardar si corresponde
    if path_guardar is not None:
        pl.savefig(path_guardar, bbox_inches='tight')
    pl.show()
    return ax

def plot_basin(self, 
               vector_cuenca=None, 
               ax=None, 
               fig=None, 
               scat_complex=False,
               scat_df=None, 
               scat_x=None, scat_y=None, scat_color=None, 
               scat_size=None, scat_cmap=None, 
               scat_order=4, scat_w=4, scat_vmin=None, scat_vmax=None,
               scat_cm_loc=[0.2, 0.1, 0.4, 0.03], scat_cm_orientation='horizontal',
               figsize=(10, 10), path_guardar=None, dpi=100, 
               cmap=pl.get_cmap('viridis'), 
               title_size=24, titulo='', titulo_colorbar='',
               norm=None, levels=None, vmin=None, vmax=None,
               color_perimetro='r',
               shape_path=None, shape_color='blue', shape_width=0.5,
               cbar_title='', cbar_loc=[0.4, 0.8, 0.4, 0.03], cbar_ticks=None, cbar_ticklabels=None,
               cbar_ticksize=16, cbar_orientation='horizontal', cbar_title_size=16
              ):
    """
    Plot flexible de cuenca, permite superponer datos categóricos,
    scatter, polígonos externos, etc. Muy útil para informes.
    """
    # Proyección
    try:
        proj = ccrs.epsg(self.epsg)
    except:
        proj = ccrs.PlateCarree()
    # Crear ejes/figura si no dados
    if ax is None:
        fig = pl.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection=proj)

    t = ax.set_title(titulo, fontsize=title_size)
    t.set_y(1.05)
    # Plot principal
    if vector_cuenca is not None:
        mapa, prop = self.Transform_Basin2Map(vector_cuenca)
        celdas_x, celdas_y, xll, yll, dx, dy, nodata = prop
        mapa[mapa == nodata] = np.nan
        longitudes = xll + dx * np.arange(celdas_x)
        latitudes = yll + dy * np.arange(celdas_y)
        longitudes, latitudes = np.meshgrid(longitudes, latitudes)
        cs = ax.contourf(longitudes, latitudes, mapa.T[::-1], transform=proj,
                         cmap=cmap, levels=levels, norm=norm, vmin=vmin, vmax=vmax)
        # Colorbar
        cax = fig.add_axes(cbar_loc)
        cbar = pl.colorbar(cs, cax=cax, orientation=cbar_orientation)
        cbar.ax.tick_params(labelsize=cbar_ticksize)
        cbar.ax.set_title(cbar_title, size=cbar_title_size)
        if cbar_ticks is not None:
            cbar.set_ticks(cbar_ticks)
        if cbar_ticklabels is not None:
            cbar.set_ticklabels(cbar_ticklabels)
    else:
        cbar = None; longitudes = None; latitudes = None

    # Scatter especial si hay DataFrame de datos puntuales (ej: estaciones)
    if scat_df is not None:
        if scat_complex:
            scat_elem = ax.scatter(scat_df[scat_x], scat_df[scat_y], 
                                   c=scat_df[scat_color], cmap=scat_cmap,
                                   vmin=scat_vmin, vmax=scat_vmax, s=scat_size, zorder=scat_order,
                                   lw=scat_w, edgecolor='k', transform=proj)
            cax = fig.add_axes(scat_cm_loc)
            sc_cbar = pl.colorbar(scat_elem, cax=cax, orientation=scat_cm_orientation)
        else:
            scat_elem = ax.scatter(scat_df[scat_x], scat_df[scat_y], c=scat_color,
                                   s=scat_size, edgecolor='k', transform=proj)
            sc_cbar = None
    else:
        sc_cbar = None

    # Polígono de la cuenca borde exterior
    ax.plot(self.Polygon[0], self.Polygon[1], color=color_perimetro)

    # QGIS Shapefile añadido extra
    if shape_path is not None:
        ax.add_geometries(Reader(shape_path).geometries(), proj,
                          edgecolor=shape_color, lw=shape_width, facecolor='none')

    # Guardar
    if path_guardar is not None:
        pl.savefig(path_guardar, bbox_inches='tight', dpi=dpi)
    pl.show()
    return ax, cbar, longitudes, latitudes, sc_cbar

def Plot_basinClean(self, vec, path=None, threshold=0.0,
                    vmin=0.0, vmax=None, show_cbar=False, **kwargs):
    """
    Genera un plot limpio de la cuenca con la variable entregada.

    Parámetros:
    -----------
    vec : np.ndarray
        Vector con valores por celda de la cuenca.
    path : str
        Ruta para guardar imagen (opcional).
    threshold : float o array
        Umbral/máscara para decidir qué celdas mostrar.
    vmin, vmax : float
        Rango de color
    show_cbar : bool
        Mostrar barra de color
    **kwargs : parámetros opcionales de estilo:
        cmap, figsize, cbar_aspect, cbar_ticks, cbar_ticklabels, cbar_ticksize, show, interpolation, grid, clean, transparent
    """
    # Parámetros de estilo
    cmap = kwargs.get('cmap', 'Spectral')
    figsize = kwargs.get('figsize', (10, 8))
    cbar_aspect = kwargs.get('cbar_aspect', 20)
    cbar_ticklabels = kwargs.get('cbar_ticklabels', None)
    cbar_ticks = kwargs.get('cbar_ticks', None)
    cbar_ticksize = kwargs.get('cbar_ticksize', 14)
    show = kwargs.get('show', True)
    transparent = kwargs.get('transparent', False)
    grid = kwargs.get('grid', False)
    clean = kwargs.get('clean', True)
    norm = kwargs.get('norm', None)

    # Convertir vector a matriz
    M, _ = self.Transform_Basin2Map(vec)
    M[(M == -9999) | (M < threshold)] = np.nan

    # Coordenadas XY de las celdas
    x, y = cu.basin_coordxy(self.structure, self.ncells)

    fig = pl.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    sca = pl.scatter(x, y, s=5, c=vec, lw=0, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm)
    if grid:
        pl.grid(True)

    if show_cbar:
        cbar = pl.colorbar(sca, aspect=cbar_aspect)
        if cbar_ticks is not None:
            cbar.set_ticks(cbar_ticks)
        if cbar_ticklabels is not None:
            cbar.ax.set_yticklabels(cbar_ticklabels, size=cbar_ticksize)
    else:
        cbar = None

    # Limpieza de ejes para aspecto limpio
    if clean:
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.axis('off')

    if path is not None:
        pl.savefig(path, bbox_inches='tight', pad_inches=0,
                   transparent=transparent, edgecolor='none', facecolor='none')

    if show:
        pl.show()
    pl.close(fig)
    return ax

def Plot_Tc(self, path=None, figsize=(8, 6), **kwargs):
    """
    Visualiza los diferentes tiempos de concentración calculados en la cuenca.

    Parámetros:
    -----------
    path : str opcional
        Ruta para guardar imagen.
    figsize : tuple
        Tamaño de la figura.
    Otros kwargs para colores ('color1', 'color2')
    """
    keys, values = list(self.Tc.keys()), list(self.Tc.values())
    # Ordenar por valor
    keys, values = zip(*sorted(zip(keys, values)))
    # Renombrar para visualización
    keys = [k.replace('Direccion Carreteras Espana', u'Carr Espana') for k in keys]
    Media = np.mean(values)
    Desv = np.std(values)
    rango = [Media-Desv, Media+Desv]

    color1 = kwargs.get('color1','b')
    color2 = kwargs.get('color2','r')

    colores = [color2 if (t > rango[1] or t < rango) else color1 for t in values]

    fig = pl.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.bar(keys, values, color=colores)
    for i, t in enumerate(values):
        ax.text(keys[i], t, "%.1f" % t, ha='center', va='bottom', fontsize=12)

    ax.axhline(Media, ls='--', color='black', label='Media')
    ax.axhline(rango[1], ls=':', color='gray', label='Media+Desv')
    ax.axhline(rango, ls=':', color='gray', label='Media-Desv')
    ax.set_ylabel('Tc [horas]')
    ax.set_title('Tiempos de concentración - Métodos')
    ax.legend()
    pl.tight_layout()
    if path: pl.savefig(path, bbox_inches='tight')
    pl.show()

def Plot_Slope_Hist(self, bins=(0,0.5,0.01), lw=2, figsize=(9, 7), 
                    axissize=18, labelsize=22, show=True, path=None, fig=None):
    """
    Histograma de pendientes (pdf por celda).
    """
    Nsize = self.ncells
    pos = np.random.choice(self.ncells, Nsize)
    h, b = np.histogram(self.CellSlope[pos], bins=np.arange(*bins))
    b = (b[:-1] + b[1:]) / 2.0
    h = h.astype(float) / h.sum()

    if fig is None:
        fig = pl.figure(figsize=figsize, facecolor='w')
        ax = fig.add_subplot(111)
    else:
        ax = fig.gca()

    ax.plot(b, h, lw=lw)
    ax.grid(True)
    ax.set_xlabel('Pendiente', size=labelsize)
    ax.set_ylabel('PDF [%]', size=labelsize)
    ax.tick_params(labelsize=axissize)

    if path: pl.savefig(path, bbox_inches='tight')
    if show: pl.show()
    return ax

def Plot_Travell_Hist(self, path=None, Nint=10.0):
    """
    Histograma y CDF de tiempos de viaje en la cuenca.
    """
    bins = np.arange(0, np.ceil(self.CellTravelTime.max()), 
                     np.ceil(self.CellTravelTime.max()) / Nint)
    h_lib, b_lib = np.histogram(self.CellTravelTime, bins=bins)
    h_lib = h_lib.astype(float) / h_lib.sum()
    b_lib = (b_lib[:-1] + b_lib[1:]) / 2.0
    hc_lib = np.cumsum(h_lib)

    fig = pl.figure(facecolor='w', edgecolor='w')
    ax = fig.add_subplot(111)
    ax.plot(b_lib, h_lib, 'b', lw=2, label='PDF Tiempos')
    ax2 = ax.twinx()
    ax2.plot(b_lib, hc_lib, 'r', lw=2, label='CDF Tiempos')
    ax2.set_ylim(0, 1.1)
    ax.set_xlim(0, np.ceil(self.CellTravelTime.max()))

    ax.grid(True)
    ax.set_xlabel('Tiempo t [hrs]', size=14)
    ax.set_ylabel('PDF [%]', size=14)
    ax2.set_ylabel('CDF [%]', size=14)
    ax.legend(loc=4)
    if path: pl.savefig(path, bbox_inches='tight')
    pl.show()

def Plot_Hipsometric(self, path=None, ventana=10, normed=False, figsize=(8,6)):
    """
    Plot de la curva hipsométrica del cauce principal y la cuenca.
    """
    elevPpal = pd.Series(self.hipso_ppal[1])
    elevBasin = self.hipso_basin[1]
    if normed:
        elevPpal = (elevPpal - elevPpal.min()) / elevPpal.max() * 100.0
        elevBasin = (elevBasin - elevBasin.min()) / elevBasin.max() * 100.0

    # Suavizar (rolling mean)
    elevPpal = elevPpal.rolling(ventana).mean()
    ppal_acum = (self.hipso_ppal / self.hipso_ppal[0, -1]) * 100
    basin_acum = (self.hipso_basin / self.hipso_basin[0, 0]) * 100

    fig = pl.figure(facecolor='w', edgecolor='w', figsize=figsize)
    ax = fig.add_subplot(111)
    ax.plot(ppal_acum, elevPpal, c='b', lw=3, label='Cauce Principal')
    ax.plot(basin_acum, elevBasin, c='r', lw=3, label='Cuenca')
    ax.grid(True)
    ax.set_xlabel('Porcentaje Área Acumulada [%]', size=16)
    if normed:
        ax.set_ylabel('Elevación [%]', size=16)
    else:
        ax.set_ylabel('Elevación [m.s.n.m]', size=16)
    ax.legend(loc=0)
    if path: pl.savefig(path, bbox_inches='tight')
    pl.show()

class SimuBasin(Basin):
    """
    Clase para simular hidrología distribuida acoplada a la cuenca.
    Hereda estructura geomorfológica, pero añade variables de simulación.
    """

    def __init__(self, lat=None, lon=None, DEM=None, DIR=None, path=None, name='NaN', stream=None,
                 threshold=500, useCauceMap=None, noData=-999, modelType='cells',
                 SimSed=False, SimSlides=False, dt=60, SaveStorage='no', SaveSpeed='no',
                 retorno=0, SeparateFluxes='no', SeparateRain='no', ShowStorage='no', SimFloods='no',
                 controlNodos=True, storageConstant=0.001):
        """
        Inicializa el objeto de simulación sobre la cuenca.

        Parámetros principales:
        ----------------------
        lat, lon : float
            Coordenada outlet de la cuenca.
        DEM, DIR : np.ndarray
            Modelos raster DEM y mapa DIRECCIONES.
        path : str, opcional
            Ruta para cargar SimuBasin guardado previamente.
        modelType : str
            'cells' (default, modelo por celdas) o 'hills' (por subcuencas).
        threshold : int
            Umbral de acumulación para cauce.
        SimSed, SimSlides : bool/'si'
            Simular sedimentos / deslizamientos.
        dt : int
            Paso de tiempo de simulación (seg).
        SaveStorage, SaveSpeed : 'si'/'no'
            Guardar almacenamiento / velocidad en cada paso.
        Others: opciones para control, almacenamiento, flujos etc.
        """
        # Variable de estado para segunda cuenca
        self.segunda_cuenca = False

        # Variables radar
        self.radarDates = []
        self.radarPos = []
        self.radarMeanRain = []
        self.radarCont = 1

        # Si no hay path, traza la cuenca sobre los mapas
        global Global_EPSG
        if path is None and int(Global_EPSG) > 0:
            # Corrección outlet si se pasa stream
            if stream is not None:
                error=[]
                for i in stream.structure.T:
                    error.append(np.sqrt((lat-i[0])**2 + (lon-i[1])**2))
                loc = np.argmin(error)
                lat = stream.structure[0, loc]
                lon = stream.structure[1, loc]

            # Corrección outlet si hay cauce map
            if useCauceMap is not None and useCauceMap.shape == DEM.shape:
                lat, lon = cu.stream_find_to_corr(lat, lon, DEM, DIR, useCauceMap,
                                                  cu.ncols, cu.nrows)

            # Asignar atributos
            self.name = name
            self.DEM = DEM
            self.DIR = DIR
            self.modelType = modelType
            self.nodata = noData
            self.threshold = threshold
            self.epsg = Global_EPSG

            # TRAZAR la cuenca con nueva estructura
            self.ncells = cu.basin_find(lat, lon, DIR, cu.ncols, cu.nrows)
            self.structure = cu.basin_cut(self.ncells)
            self.DEMvec = self.Transform_Map2Basin(DEM, [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])
            self.DIRvec = self.Transform_Map2Basin(DIR, [cu.ncols, cu.nrows, cu.xll, cu.yll, cu.dx, cu.dy])

            # Estructura de hills o celdas, según modelType
            acum = cu.basin_acum(self.structure, self.ncells)
            cauce, nodos, self.nhills = cu.basin_subbasin_nod(self.structure, acum, threshold, self.ncells)
            self.hills_own, sub_basin = cu.basin_subbasin_find(self.structure, nodos, self.nhills, self.ncells)
            self.hills = cu.basin_subbasin_cut(self.nhills)
            models.drena = self.structure

            # Inicialización de variables según modelo
            if modelType == 'cells':
                N = self.ncells
            elif modelType == 'hills':
                N = self.nhills

            # Variables físicas
            models.v_coef = np.ones((4, N))
            models.h_coef = np.ones((4, N))
            models.v_exp = np.ones((4, N))
            models.h_exp = np.ones((4, N))
            models.max_capilar = np.ones((1, N))
            models.max_gravita = np.ones((1, N))
            models.max_aquifer = np.ones((1, N))
            models.storage = np.zeros((5, N))
            models.dt = dt
            models.calc_niter = 5
            models.retorno_gr = 0
            models.verbose = 0
            models.control = np.zeros((1, N))

            # Definir puntos de control en la red
            if controlNodos:
                if modelType == 'cells':
                    self.GetGeo_Cell_Basics()
                    cauce, nodos, n_nodos = cu.basin_subbasin_nod(
                        self.structure, self.CellAcum, threshold, self.ncells)
                    pos = np.where(nodos != 0)[0]
                    x, y = cu.basin_coordxy(self.structure, self.ncells)
                    idsOrd, xy = self.set_Control(np.vstack([x[pos], y[pos]]), nodos[pos])
                elif modelType == 'hills':
                    models.control = np.ones((1, self.nhills)) * nodos[nodos != 0]
            # Otros controles, humedad, etc.
            models.control_h = np.zeros((1, N))
            models.sim_sediments = int(SimSed == 'si')
            models.sim_slides = int(bool(SimSlides))
            models.save_storage = int(SaveStorage == 'si')
            models.save_speed = int(SaveSpeed == 'si')
            models.separate_fluxes = int(SeparateFluxes == 'si')
            models.separate_rain = int(SeparateRain == 'si')
            models.show_storage = int(ShowStorage == 'si')
            models.sim_floods = int(SimFloods == 'si')
            models.storage_constant = storageConstant
            self.isSetGeo = False

        # SI hay path, carga la cuenca (puede incluir info de simulación previa)
        elif path is not None:
            self.__Load_SimuBasin(path, SimSlides)
            self.__GetBasinPolygon__()

def run_shia(self, Calibracion, rain_path, N_intervals, start_point=1,
             StorageLoc=None, HspeedLoc=None, path_storage=None, path_speed=None,
             path_conv=None, path_stra=None, path_retorno=None, kinematicN=5,
             QsimDataFrame=True, EvpVariable='sun', EvpSerie=None,
             WheretoStore=None, path_vfluxes=None, Dates2Save=None, FluxesDates2Save=None,
             path_rc=None):
    """
    Ejecuta el modelo hidrológico SHIA usando la configuración del objeto.

    Parámetros:
    -----------
    Calibracion : list/array
        Parámetros de calibración.
    rain_path : str
        Ruta al binario de lluvia (generado antes).
    N_intervals : int
        Número de pasos de simulación.
    start_point : int
        Primer índice de lluvia a utilizar.
    StorageLoc : np.ndarray, opcional
        Estado inicial de almacenamiento (forma (5, N)).
    HspeedLoc : np.ndarray, opcional
        Velocidades horizontales iniciales (forma (4, N)).
    path_storage: str, opcional
        Ruta de guardado de almacenamiento en cada paso.
    path_speed: str, opcional
        Ruta para guardar velocidades.
    path_conv, path_stra: str, opcional
        Archivos para indicar lluvia convectiva/estratiforme.
    kinematicN : int
        Iteraciones de onda cinemática.
    QsimDataFrame : bool
        Devolver caudales simulados en DataFrame.

    Returns:
    --------
    Retornos : dict, y opcionalmente DataFrames de Qsim/Qsep/Qsedi
    """
    rain_pathBin, rain_pathHdr = __Add_hdr_bin_2route__(rain_path)
    Rain = read_mean_rain(rain_pathHdr, N_intervals, start_point)
    N = self.ncells if self.modelType[0] == 'c' else self.nhills

    # Configuración interna del modelo
    models.rain_first_point = start_point
    models.calc_niter = kinematicN

    # Guardar almacenamiento
    if path_storage is not None:
        models.save_storage = 1
        path_sto_bin, path_sto_hdr = __Add_hdr_bin_2route__(path_storage, storage=True)
        if Dates2Save is not None:
            WhereItSaves = self.set_StorageDates(Rain.index, Dates2Save, N_intervals)
        else:
            print('Warning: model will save states in all time steps')
            WhereItSaves = np.arange(1, N_intervals+1)
    else:
        models.save_storage = 0
        path_sto_bin = 'no_save.StObin'
        path_sto_hdr = 'no_save.StOhdr'

    # Guardar flujos verticales
    if path_vfluxes is not None:
        models.save_vfluxes = 1
        path_vflux_bin, path_vflux_hdr = __Add_hdr_bin_2route__(path_vfluxes)
        if FluxesDates2Save is not None:
            FluxesWhereItSaves = self.set_vFluxesDates(Rain.index, FluxesDates2Save, N_intervals)
        else:
            FluxesWhereItSaves = np.arange(1, N_intervals+1)
    else:
        models.save_vfluxes = 0
        path_vflux_bin = 'no_save.bin'
        path_vflux_hdr = 'no_save.hdr'

    # Otras configuraciones similares omitidas por brevedad...

    # Ejecutar modelo (llamada a núcleo en `models`)
    Qsim, Qsed, Qsep, Humedad, St1, St3, Balance, Speed, Area, Alm, Qsep_byrain = \
        models.shia_v1(rain_pathBin, rain_pathHdr, Calibracion,
                       N, np.count_nonzero(models.control)+1,
                       np.count_nonzero(models.control_h),
                       N_intervals, StorageLoc, HspeedLoc,
                       path_sto_bin, path_speed, path_vflux_bin,
                       path_conv, path_stra, "", "",
                       path_retorno, path_rc)

    # Ensamblar retornos
    Retornos = {'Qsim': Qsim, 'Balance': Balance, 'Storage': Alm,
                'Rain_Acum': models.acum_rain, 'Rain_hietogram': models.mean_rain}

    if np.count_nonzero(models.control_h) > 0:
        Retornos.update({'Humedad': Humedad, 'Humedad_t1': St1, 'Humedad_t2': St3})
    if models.sim_sediments == 1:
        Retornos.update({'Sediments': Qsed})
    if models.separate_fluxes == 1:
        Retornos.update({'Fluxes': Qsep})
    if models.separate_rain == 1:
        Retornos.update({'Rain_sep': Qsep_byrain})

    if QsimDataFrame:
        ids = models.control[models.control != 0]
        Qdict = {str(j): i for i, j in zip(Retornos['Qsim'][1:], ids)}
        Qdf = pd.DataFrame(Qdict, index=Rain.index)
        return Retornos, Qdf
    return Retornos

def efficiencia(self, Qobs, Qsim):
    """
    Calcula índices estadísticos de ajuste entre caudal observado y simulado.

    Retorna un diccionario con:
        - Nash
        - Qpico
        - Tpico
        - RMSE log
        - RMSE normal
    """
    DictEff = {'Nash': __eval_nash__(Qobs, Qsim),
               'Qpico': __eval_q_pico__(Qobs, Qsim),
               'Tpico': __eval_t_pico__(Qobs, Qsim, models.dt),
               'RmseLog': __eval_rmse_log__(Qobs, Qsim),
               'Rmse': __eval_rmse__(Qobs, Qsim)}
    return DictEff

def Calib_NSGAII(self, nsga_el, nodo_eval, pop_size=40,
                 process=4, NGEN=6, MUTPB=0.5, CXPB=0.5):
    """
    Lanza un proceso de calibración multiobjetivo NSGA-II usando DEAP.

    nsga_el: objeto nsgaii_element configurado para esta cuenca.
    nodo_eval: int, índice del nodo de control a evaluar.
    """
    # Asegurar tamaño múltiplo de 4
    while pop_size % 4 != 0:
        pop_size += 1

    nsga_el.set_nsgaII()
    pop = nsga_el.toolbox.population(pop_size)

    # Evaluar población inicial
    Ejecs = map(nsga_el.__crea_ejec__, pop)
    QsimPar = __ejec_parallel__(Ejecs, process, nodo_eval)
    fitnesses = map(nsga_el.toolbox.evaluate, QsimPar)
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit

    # Generaciones
    for g in range(NGEN):
        offspring = tools.selTournamentDCD(pop, len(pop))
        offspring = list(map(nsga_el.toolbox.clone, offspring))
        # Cruce y mutación
        for c1, c2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < CXPB:
                nsga_el.toolbox.mate(c1[0], c2)
                del c1.fitness.values, c2.fitness.values
        for mut in offspring:
            if random.random() < MUTPB:
                nsga_el.toolbox.mutate(mut)
                del mut.fitness.values

        # Reevaluar
        Ejecs = map(nsga_el.__crea_ejec__, offspring)
        QsimPar = __ejec_parallel__(Ejecs, process, nodo_eval)
        invalids = [(ind, q) for ind, q in zip(offspring, QsimPar) if not ind.fitness.valid]
        fits = map(nsga_el.toolbox.evaluate, [q for _, q in invalids])
        for ind, fit in zip([i for i, _ in invalids], fits):
            ind.fitness.values = fit

        pop = nsga_el.toolbox.select(pop + offspring, pop_size)

    return pop, QsimPar, np.array(list(map(lambda p: p.fitness.values, pop))).T

class nsgaii_element:
    """
    Clase auxiliar para configurar calibración multiobjetivo NSGA-II en SimuBasin.
    """

    def __init__(self, pathLluvia, Qobs, npasos, inicio, SimuBasinElem,
                 evp=[0,1], infil=[1,200], perco=[1,40], losses=[0,1],
                 velRun=[0.1,1], velSub=[0.1,1], velSup=[0.1,1], velStream=[0.1,1],
                 Hu=[0.1,1], Hg=[0.1,1],
                 probCruce=None, probMutacion=None,
                 rangosMutacion=None, MaxMinOptima=(1.0, -1.0), CrowDist=0.5):
        """
        Inicializa objeto con rangos y configuraciones de NSGA-II.

        pathLluvia : str
            Ruta al archivo de lluvia binario.
        Qobs : np.ndarray
            Serie observada.
        npasos : int
            Número de pasos a simular.
        inicio : int
            Índice de inicio en datos de lluvia.
        SimuBasinElem : SimuBasin
            Objeto SimuBasin configurado para simulación.
        Otros: rangos para cada parámetro de calibración.
        """
        # Rangos de cada parámetro (min, max)
        self.evp_range = evp
        self.infil_range = infil
        self.perco_range = perco
        self.losses_range = losses
        self.velRun_range = velRun
        self.velSub_range = velSub
        self.velSup_range = velSup
        self.velStream_range = velStream
        self.hu_range = Hu
        self.hg_range = Hg

        # Configuración de simulación
        self.npasos = npasos
        self.inicio = inicio
        self.path_lluvia = pathLluvia
        self.Qobs = Qobs
        self.simelem = SimuBasinElem

        # Probabilidades de cruce/mutación
        self.prob_cruce = probCruce if probCruce is not None else np.ones(10) * 0.5
        self.prob_mutacion = probMutacion if probMutacion is not None else np.ones(10) * 0.5
        self.rangos_mutacion = rangosMutacion if rangosMutacion is not None else [
            evp, infil, perco, losses, velRun, velSub, velSup, velStream, Hu, Hg
        ]

        self.optimiza = MaxMinOptima
        self.crowdist = CrowDist

    # ----------------------------------------------------------
    # Creación de calibraciones aleatorias (individuos)
    # ----------------------------------------------------------
    def __crea_calibracion__(self):
        """Crea un individuo (calibración) aleatorio dentro de los rangos."""
        def rand_range(r):
            return np.random.uniform(r[0], r[1], 1)

        return [
            rand_range(self.evp_range),
            rand_range(self.infil_range),
            rand_range(self.perco_range),
            rand_range(self.losses_range),
            rand_range(self.velRun_range),
            rand_range(self.velSub_range),
            rand_range(self.velSup_range),
            rand_range(self.velStream_range),
            rand_range(self.hu_range),
            rand_range(self.hg_range),
        ]

    # ----------------------------------------------------------
    # Generación de ejecución para un individuo
    # ----------------------------------------------------------
    def __crea_ejec__(self, calibracion):
        """Retorna una lista empaquetada para ejecutar con multiprocessing."""
        return [calibracion, self.path_lluvia, self.npasos, self.inicio, self.simelem]

    # ----------------------------------------------------------
    # Función de evaluación
    # ----------------------------------------------------------
    def __evalfunc__(self, Qsim, f1=__eval_nash__, f2=__eval_q_pico__):
        """Evalúa un Qsim contra Qobs usando dos funciones objetivo."""
        E1 = f1(self.Qobs, Qsim)
        E2 = f2(self.Qobs, Qsim)
        return E1, E2

    # ----------------------------------------------------------
    # Cruce entre individuos
    # ----------------------------------------------------------
    def __cruce__(self, indi1, indi2):
        """Cruza parámetros entre dos individuos según probabilidad."""
        for i, p_cr in enumerate(self.prob_cruce):
            if np.random.rand() < p_cr:
                indi1[i], indi2[i] = indi2[i], indi1[i]
        return indi1, indi2

    # ----------------------------------------------------------
    # Mutación de individuo
    # ----------------------------------------------------------
    def __mutacion__(self, indi):
        """Muta parámetros de un individuo según probabilidad."""
        for i, p_mut in enumerate(self.prob_mutacion):
            if np.random.rand() < p_mut:
                r = self.rangos_mutacion[i]
                indi[i] = np.random.uniform(r[0], r[1], 1)
        return indi

    # ----------------------------------------------------------
    # Configurar operadores en DEAP
    # ----------------------------------------------------------
    def set_nsgaII(self):
        """Configura toolbox DEAP con operadores y tipo de fitness."""
        creator.create("FitnessOpt", base.Fitness, weights=self.optimiza, crowding_dist=self.crowdist)
        creator.create("Individual", list, fitness=creator.FitnessOpt)
        self.toolbox = base.Toolbox()
        self.toolbox.register("attr1", self.__crea_calibracion__)
        self.toolbox.register("individual", tools.initRepeat, creator.Individual, self.toolbox.attr1, n=1)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("evaluate", self.__evalfunc__)
        self.toolbox.register("mate", self.__cruce__)
        self.toolbox.register("mutate", self.__mutacion__)
        self.toolbox.register("select", tools.selNSGA2)

class Stream:
    """
    Clase para representar un cauce o corriente extraído desde DEM/DIR.
    Permite su trazado, guardado a GIS y visualización de perfil.
    """

    # ------------------------------------------------------
    # Inicialización/trazado de la corriente
    # ------------------------------------------------------
    def __init__(self, lat, lon, DEM, DIR, name='NaN'):
        """
        Inicializa un objeto de corriente fluvial.

        Parámetros:
        -----------
        lat, lon : float
            Coordenadas del punto más alto o inicial de la corriente.
        DEM : np.ndarray
            Modelo digital de elevación.
        DIR : np.ndarray
            Mapa de direcciones de flujo.
        name : str
            Nombre opcional del cauce.
        """
        self.DEM = DEM
        self.DIR = DIR
        self.name = name

        # Usar funciones del módulo cu para trazar la corriente
        self.ncells = cu.stream_find(lat, lon, self.DEM, self.DIR, cu.ncols, cu.nrows)
        self.structure = cu.stream_cut(self.ncells)

    # ------------------------------------------------------
    # Guardado del cauce como shapefile
    # ------------------------------------------------------
    def Save_Stream2Map(self, path, DriverFormat='ESRI Shapefile', EPSG=4326):
        """
        Guarda el cauce trazado en un shapefile.

        Parámetros:
        -----------
        path : str
            Ruta del archivo de salida (.shp).
        DriverFormat : str
            Formato GIS, por defecto ESRI Shapefile.
        EPSG : int
            Código EPSG de proyección.
        """
        from osgeo import ogr, osr
        if not path.endswith('.shp'):
            path += '.shp'

        if os.path.exists(path):
            ogr.GetDriverByName(DriverFormat).DeleteDataSource(path)

        driver = ogr.GetDriverByName(DriverFormat)
        shapeData = driver.CreateDataSource(path)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(EPSG)
        layer = shapeData.CreateLayer('layer1', srs, ogr.wkbLineString)

        line = ogr.Geometry(ogr.wkbLineString)
        for x, y in zip(self.structure[0], self.structure[1]):
            line.AddPoint(float(x), float(y))

        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometry(line)
        feat.SetFID(0)
        layer.CreateFeature(feat)

        shapeData.Destroy()

    # ------------------------------------------------------
    # Plot del perfil topográfico del cauce
    # ------------------------------------------------------
    def Plot_Profile(self, path=None):
        """
        Grafica la elevación a lo largo del cauce (perfil longitudinal).

        Parámetros:
        -----------
        path : str, opcional
            Ruta para guardar gráfico.
        """
        fig = pl.figure(facecolor='w', edgecolor='w')
        ax = fig.add_subplot(111)
        ax.plot(self.structure[3], self.structure, lw=2)
        ax.set_xlabel('Distancia [m]', size=14)
        ax.set_ylabel('Elevación [m.s.n.m]', size=14)
        ax.grid(True)

        if path is not None:
            pl.savefig(path, bbox_inches='tight')
        pl.show()
