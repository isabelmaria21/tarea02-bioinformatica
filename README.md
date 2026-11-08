# Tarea 02 - Análisis de expresión diferencial con PyDESeq2

## 1. Descripción del proyecto

Este proyecto realiza un análisis de expresión diferencial de genes utilizando datos de RNA-seq del proyecto TCGA-BRCA (cáncer de mama).

El objetivo es identificar genes cuya expresión cambia significativamente entre muestras tumorales y muestras de tejido normal adyacente.

El análisis se ha realizado en Python utilizando PyDESeq2, incluyendo la variable experimental `plate` para controlar su posible efecto sobre los resultados.

## 2. Datos utilizados

Los archivos originales se encuentran en la carpeta `data/raw/`:

- `counts.csv.gz`: matriz de recuentos de expresión génica.
- `metadata.csv.gz`: información de las muestras y variables experimentales.

Se analizaron 75 muestras:
- 55 muestras tumorales.
- 20 muestras normales.

La matriz original contenía 20.502 genes.

Tras aplicar un filtrado de expresión mínima, se conservaron 18.705 genes.

## 3. Análisis de expresión diferencial

Se utilizó PyDESeq2 para comparar las muestras tumorales frente a las normales, ajustando por la variable `plate`.

El diseño estadístico utilizado fue:

`~ plate + site`

Se consideraron genes diferencialmente expresados aquellos que cumplían:

- Valor p ajustado < 0,05.
- Valor absoluto de log2FoldChange > 1.

Resultados obtenidos:

- Genes con valor p ajustado < 0,05: 10.630.
- Genes sobreexpresados en tumor (log2FoldChange > 1): 2.320.
- Genes infraexpresados en tumor (log2FoldChange < -1): 2.113.
- Total de genes diferencialmente expresados con ambos criterios: 4.433.

## 4. Visualización de resultados

Se generaron las siguientes gráficas:

- Volcano plot: `results/volcano_plot.png`
- MA plot: `results/ma_plot.png`
- Top 10 genes diferencialmente expresados: `results/top10_genes.png`

## 5. Análisis de enriquecimiento funcional

Se utilizó GSEApy para realizar análisis de enriquecimiento funcional de los 4.433 genes diferencialmente expresados.

Se analizaron dos bases de datos:

### Gene Ontology - Biological Process

Permite identificar procesos biológicos enriquecidos entre los genes diferencialmente expresados.

Resultados:
- `results/enriquecimiento_GO_BP.tsv`
- `results/go_bp_top10.png`
- `results/go_bp_top10.tsv`

### Reactome

Permite identificar rutas biológicas enriquecidas.

Resultados:
- `results/enriquecimiento_Reactome.tsv`
- `results/reactome_top10.png`
- `results/reactome_top10.tsv`

Entre las rutas identificadas destacan las relacionadas con la organización de la matriz extracelular y la señalización celular.

## 6. Organización del proyecto

- `data/raw/`: datos originales.
- `data/processed/`: datos procesados.
- `src/`: código Python del análisis.
- `results/`: tablas y gráficas obtenidas.
- `logs/`: registros e información de ejecución.
- `environment.yml`: dependencias del proyecto.
- `main.sh`: script para ejecutar el análisis.
- `README.md`: documentación del proyecto.

## 7. Ejecución del análisis

Para ejecutar el análisis desde la carpeta principal del proyecto:

```bash
bash main.sh
```

El script ejecuta el análisis de expresión diferencial y los análisis de enriquecimiento funcional, generando los archivos correspondientes en la carpeta `results/`.