#!/usr/bin/env bash

# Detener la ejecución si ocurre algún error.
set -e

# Situarnos en la carpeta principal del proyecto.
cd "$(dirname "$0")"

# Crear las carpetas necesarias.
mkdir -p results logs

# Ejecutar el análisis completo.
echo "Iniciando análisis de expresión diferencial TCGA-BRCA..."

python src/01_run_deseq2.py

echo "Análisis finalizado correctamente."