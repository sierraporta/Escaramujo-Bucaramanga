# Escaramujo-Bucaramanga
Códigos para análisis de datos de Escaramujo@Bucaramanga. Escaramujo project.
El Proyecto Escaramujo es una plataforma para detección de rayos cósmicos con placas centelladoras. En este repositorio se encuentra el código utilizado para la adquisión y el procesamiento de datos.

Framework for analysis of cosmic ray data, using C++/ROOT and Python.

Sitio oficial del proyecto http://es.escaramujo.net/

Contenido:

+ El código **_qNet2root6000_GPS.exe_** se usa para el análisis de conteos de muones en el detector dado los datos que salen de la minicom. Se necesita compilar **_qNet2root6000_GPS.cc_** con **_Makefile_** y luego ejecutar **_qNet2root6000_GPS.exe_** desde la terminal de linux como:
> make

> ./qNet2root6000.exe datos_entrada.dat datos_salida.root
