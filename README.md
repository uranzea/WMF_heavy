# Water Model Framework (WMF)

Water Model Framework (WMF) es un conjunto avanzado de herramientas para modelación hidrológica distribuida y análisis geomorfológico de cuencas hidrográficas. Este proyecto integra procesamiento de datos geoespaciales,modelos hidrológicos dinámicos, visualización avanzada y calibración multiobjetivo en un paquete Python flexible y extensible.

---

## Características principales

- Procesamiento de Modelos Digitales de Elevación (DEM) y mapas de direcciones de flujo (DIR).
- Delimitación automática y análisis geomorfológico de cuencas y subcuencas.
- Modelo hidrológico distribuido basado en SHIA con simulación dinámica.
- Soporte para simulación de sedimentos, deslizamientos e inundaciones.
- Visualización versátil: gráficos hidrológicos, mapas GIS, perfiles de cauces.
- Exportación/importación totalmente compatible con formatos estándar GIS (shapefiles, GeoTIFF, NetCDF).
- Calibración multiobjetivo basada en algoritmos genéticos NSGA-II con DEAP.
- Diseño modular y orientado a objetos para extensibilidad.

---

## Requisitos

- Python 3.10 o superior
- NumPy, Pandas, SciPy, Matplotlib
- GDAL/OGR, Rasterio, Cartopy para procesamiento y visualización geoespacial
- netCDF4 para manejo de archivos NetCDF
- DEAP para calibración genética (opcional)
- PyProj, Shapely, PyShp para trabajo con datos vectoriales

---

## Instalación

Se recomienda crear un entorno Conda a partir del archivo `environment.yml`:

conda env create -f environment.yml
conda activate wmf

Alternativamente, puede instalarse mediante Poetry con el archivo `pyproject.toml` para manejo completo de dependencias.

También es posible instalar directamente con `pip`:

```
pip install e .
```

> Nota: la compilación de archivos Fortran (`*.f90`) ya no es necesaria; todo el paquete se distribuye como código Python puro.

---

## Uso básico

1. Preparar el DEM y mapa de direcciones (DIR) de la cuenca a analizar.
2. Crear un objeto `Basin` para representar la cuenca.
3. Generar un objeto `SimuBasin` para la simulación hidrológica.
4. Ejecutar la simulación con parámetros configurables.
5. Analizar resultados con funciones de visualización y evaluación.
6. Exportar resultados a formatos GIS para integración con SIG externos.

Ejemplo básico (Python):

from wmf import Basin, SimuBasin
Crear cuenca a partir de datos raster

basin = Basin(lat=..., lon=..., DEM=mi_dem, DIR=mi_dir)
simu = SimuBasin(lat=..., lon=..., DEM=mi_dem, DIR=mi_dir)
Ejecutar simulación con parámetros

resultados = simu.run_shia(Calibracion=[...], rain_path="datos_lluvia.bin", N_intervals=1000)
Visualizar resultados

simu.plot_sim_single(resultados['Qsim'])


---

## Documentación

Para documentación completa, ejemplos y tutoriales visite el repositorio o utilice la documentación integrada en docstrings.

---

## Contribuciones

Contribuciones, reportes de bugs y solicitudes de mejora son bienvenidas. Por favor utilice issues y pull requests en el repositorio principal.

---

## Licencia

Este proyecto está licenciado bajo GNU General Public License v3.0 o superior (GPLv3+).

---

## Contacto

Autor: Julian Uran  
Email: uranzea@gmail.com  
Repositorio: https://github.com/uranzea/wmf_heavy

