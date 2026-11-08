from pathlib import Path
import pandas as pd

# Localizamos los archivos originales del experimento.
carpeta_proyecto = Path(__file__).resolve().parent.parent
carpeta_datos = carpeta_proyecto / "data" / "raw"

# Cargamos la matriz de recuentos y los metadatos.
counts = pd.read_csv(carpeta_datos / "counts.csv.gz", index_col=0)
metadata = pd.read_csv(carpeta_datos / "metadata.csv.gz", index_col=0)

# Comprobamos las dimensiones de los datos.
print("Número de genes:", counts.shape[0])
print("Número de muestras:", counts.shape[1])
print("Dimensiones de los metadatos:", metadata.shape)

# Comprobamos si las muestras coinciden entre ambos archivos.
muestras_counts = set(counts.columns)
muestras_metadata = set(metadata.index)

print("Muestras presentes solo en counts:", muestras_counts - muestras_metadata)
print("Muestras presentes solo en metadata:", muestras_metadata - muestras_counts)

# Seleccionamos las muestras que tienen datos de expresión génica.
# También colocamos los metadatos en el mismo orden que counts.
metadata = metadata.loc[counts.columns].copy()

print("Dimensiones de los metadatos alineados:", metadata.shape)
print("¿Coincide el orden de las muestras?:", metadata.index.equals(counts.columns))

# Exploramos las variables experimentales.
print("\nColumnas de los metadatos:", metadata.columns.tolist())
print("\nTipos de tejido:")
print(metadata["site"].value_counts())
print("\nPlacas experimentales:")
print(metadata["plate"].value_counts())

# Comprobamos los grupos experimentales tras alinear los datos.
print("\nGrupos experimentales que analizaremos:")
print(metadata["site"].value_counts())

# Transponemos la matriz para adaptarla a PyDESeq2.
# Ahora las filas serán muestras y las columnas serán genes.
counts_transpuestos = counts.T

print("\nDimensiones de la matriz transpuesta:")
print(counts_transpuestos.shape)

# Comprobamos que las muestras coinciden con los metadatos.
print("¿Coinciden las muestras?:", counts_transpuestos.index.equals(metadata.index))

# Filtramos genes con muy baja expresión.
# Conservamos los genes con al menos 10 recuentos en 3 muestras.

min_count = 10
min_samples = 3

genes_validos = (counts_transpuestos >= min_count).sum(axis=0) >= min_samples

counts_filtrados = counts_transpuestos.loc[:, genes_validos].copy()

print("\nGenes antes del filtrado:", counts_transpuestos.shape[1])
print("Genes después del filtrado:", counts_filtrados.shape[1])
print("Genes eliminados:", counts_transpuestos.shape[1] - counts_filtrados.shape[1])

# Preparamos las variables del diseño experimental.
metadata["site"] = pd.Categorical(
    metadata["site"],
    categories=["Solid_Tissue_Normal", "Primary_Solid_Tumor"]
)

metadata["plate"] = metadata["plate"].astype("category")

# Comprobamos que Python reconoce las variables categóricas.
print("\nTipos de variables del diseño experimental:")
print(metadata.dtypes)

print("\nCategoría de referencia para el tejido:")
print(metadata["site"].cat.categories[0])

# Importamos la herramienta principal de PyDESeq2.
from pydeseq2.dds import DeseqDataSet

# Creamos el conjunto de datos para el análisis.
dds = DeseqDataSet(
    counts=counts_filtrados,
    metadata=metadata,
    design="~ plate + site",
    refit_cooks=True,
    n_cpus=2
)

print("\nModelo estadístico creado correctamente.")
print("Número de muestras:", dds.n_obs)
print("Número de genes:", dds.n_vars)

# Ejecutamos el análisis de expresión diferencial con PyDESeq2.
print("\nIniciando el análisis DESeq2...")
print("Este proceso puede tardar varios minutos.")

dds.deseq2()

print("\n¡Análisis DESeq2 completado correctamente!")

# Importamos la herramienta para obtener los resultados estadísticos.
from pydeseq2.ds import DeseqStats

# Comparamos tejido tumoral frente a tejido normal.
stats = DeseqStats(
    dds,
    contrast=["site", "Primary_Solid_Tumor", "Solid_Tissue_Normal"],
    n_cpus=2
)

# Calculamos los resultados de expresión diferencial.
stats.summary()

print("\n¡Resultados estadísticos calculados correctamente!")

# Guardamos la tabla completa de resultados estadísticos.
carpeta_resultados = carpeta_proyecto / "results"
carpeta_resultados.mkdir(parents=True, exist_ok=True)

resultados = stats.results_df

resultados.to_csv(
    carpeta_resultados / "deseq2_resultados_completos.tsv",
    sep="\t"
)

print("\nTabla de resultados guardada correctamente.")

# Identificamos los genes con expresión diferencial significativa.
# Utilizamos un valor p ajustado inferior a 0,05.

genes_significativos = resultados[resultados["padj"] < 0.05]

# Separamos los genes según aumente o disminuya su expresión.
genes_sobreexpresados = genes_significativos[
    genes_significativos["log2FoldChange"] > 1
]

genes_infraexpresados = genes_significativos[
    genes_significativos["log2FoldChange"] < -1
]

print("\nRESUMEN DEL ANÁLISIS DE EXPRESIÓN DIFERENCIAL")
print("Genes analizados:", len(resultados))
print("Genes estadísticamente significativos:", len(genes_significativos))
print("Genes sobreexpresados en el tumor:", len(genes_sobreexpresados))
print("Genes infraexpresados en el tumor:", len(genes_infraexpresados))

# Guardamos las tablas de genes con expresión diferencial.

genes_significativos.to_csv(
    carpeta_resultados / "genes_significativos.tsv",
    sep="\t"
)

genes_sobreexpresados.to_csv(
    carpeta_resultados / "genes_sobreexpresados.tsv",
    sep="\t"
)

genes_infraexpresados.to_csv(
    carpeta_resultados / "genes_infraexpresados.tsv",
    sep="\t"
)

print("\n¡Tablas de genes significativos guardadas correctamente!")

# --------------------------------------------------
# GRÁFICA 1: VOLCANO PLOT
# --------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Preparamos los resultados para representar los genes.
datos_volcano = resultados.copy()

# Eliminamos los genes sin valores estadísticos válidos.
datos_volcano = datos_volcano.dropna(
    subset=["log2FoldChange", "padj"]
)

# Evitamos calcular el logaritmo de cero.
datos_volcano["padj_grafica"] = datos_volcano["padj"].clip(lower=1e-300)

# Calculamos el valor que aparecerá en el eje vertical.
datos_volcano["menos_log10_padj"] = -np.log10(
    datos_volcano["padj_grafica"]
)

# Clasificamos los genes según su expresión.
datos_volcano["grupo"] = "No significativo"

datos_volcano.loc[
    (datos_volcano["padj"] < 0.05) &
    (datos_volcano["log2FoldChange"] > 1),
    "grupo"
] = "Sobreexpresado"

datos_volcano.loc[
    (datos_volcano["padj"] < 0.05) &
    (datos_volcano["log2FoldChange"] < -1),
    "grupo"
] = "Infraexpresado"

# Creamos la figura.
plt.figure(figsize=(10, 7))

colores = {
    "No significativo": "lightgray",
    "Sobreexpresado": "red",
    "Infraexpresado": "royalblue"
}

for grupo, color in colores.items():
    datos_grupo = datos_volcano[
        datos_volcano["grupo"] == grupo
    ]

    plt.scatter(
        datos_grupo["log2FoldChange"],
        datos_grupo["menos_log10_padj"],
        c=color,
        s=8,
        alpha=0.6,
        label=grupo
    )

# Añadimos los límites de significación.
plt.axvline(x=1, color="black", linestyle="--")
plt.axvline(x=-1, color="black", linestyle="--")
plt.axhline(y=-np.log10(0.05), color="black", linestyle="--")

plt.xlabel("Log2 Fold Change (tumor vs normal)")
plt.ylabel("-Log10(p ajustado)")
plt.title("Volcano plot - Expresión diferencial TCGA-BRCA")

plt.legend()
plt.tight_layout()

# Guardamos la gráfica.
plt.savefig(
    carpeta_resultados / "volcano_plot.png",
    dpi=300
)

plt.close()

print("\n¡Volcano plot guardado correctamente!")

# --------------------------------------------------
# GRÁFICA 2: MA PLOT
# --------------------------------------------------

# Preparamos los datos para representar el MA plot.
datos_ma = resultados.copy()

# Eliminamos los genes sin valores válidos.
datos_ma = datos_ma.dropna(
    subset=["baseMean", "log2FoldChange", "padj"]
)

# Clasificamos los genes según su significación estadística.
genes_sig_ma = datos_ma["padj"] < 0.05
genes_no_sig_ma = datos_ma["padj"] >= 0.05

# Creamos la figura.
plt.figure(figsize=(10, 7))

# Genes no significativos en gris.
plt.scatter(
    datos_ma.loc[genes_no_sig_ma, "baseMean"],
    datos_ma.loc[genes_no_sig_ma, "log2FoldChange"],
    color="lightgray",
    s=8,
    alpha=0.5,
    label="No significativo"
)

# Genes significativos en azul.
plt.scatter(
    datos_ma.loc[genes_sig_ma, "baseMean"],
    datos_ma.loc[genes_sig_ma, "log2FoldChange"],
    color="royalblue",
    s=8,
    alpha=0.6,
    label="Significativo (padj < 0,05)"
)

# Utilizamos una escala logarítmica en el eje X.
plt.xscale("log")

# Añadimos una línea horizontal en cero.
plt.axhline(y=0, color="black", linestyle="--")

# Etiquetas y título.
plt.xlabel("Expresión media normalizada (baseMean)")
plt.ylabel("Log2 Fold Change (tumor vs normal)")
plt.title("MA plot - Expresión diferencial TCGA-BRCA")

plt.legend()
plt.tight_layout()

# Guardamos la gráfica.
plt.savefig(
    carpeta_resultados / "ma_plot.png",
    dpi=300
)

plt.close()

print("\n¡MA plot guardado correctamente!")

# ============================================================
# GRÁFICA DE LOS 10 GENES MÁS SIGNIFICATIVOS
# ============================================================

# Seleccionamos los 10 genes con menor p-valor ajustado.
top10 = resultados.dropna(subset=["padj", "log2FoldChange"]).copy()
top10 = top10.sort_values("padj").head(10)

# Ordenamos los genes para representarlos.
top10 = top10.sort_values("log2FoldChange")

# Preparamos los colores según el cambio de expresión.
colores = [
    "tomato" if valor > 0 else "royalblue"
    for valor in top10["log2FoldChange"]
]

# Creamos la gráfica.
plt.figure(figsize=(10, 6))

plt.barh(
    top10.index,
    top10["log2FoldChange"],
    color=colores
)

# Añadimos títulos y etiquetas.
plt.axvline(x=0, color="black", linestyle="--", linewidth=1)
plt.xlabel("Log2 Fold Change (tumor vs normal)")
plt.ylabel("Genes")
plt.title("Top 10 genes más significativos - TCGA-BRCA")

plt.tight_layout()

# Guardamos la gráfica.
plt.savefig(
    carpeta_resultados / "top10_genes.png",
    dpi=300
)

plt.close()

# Guardamos también la tabla con los 10 genes.
top10.to_csv(
    carpeta_resultados / "top10_genes.tsv",
    sep="\t"
)

print("\n¡Gráfica de los 10 genes principales guardada correctamente!")

# ============================================================
# ENRIQUECIMIENTO FUNCIONAL GO BIOLOGICAL PROCESS
# ============================================================

import gseapy as gp

print("\nIniciando enriquecimiento GO Biological Process...")

# Seleccionamos genes con cambios importantes de expresión.
genes_go = resultados[
    (resultados["padj"] < 0.05) &
    (resultados["log2FoldChange"].abs() > 1)
].index.tolist()

# Utilizamos como referencia los genes analizados.
universo_go = resultados.index.tolist()

print(f"Genes seleccionados para GO BP: {len(genes_go)}")

# Realizamos el enriquecimiento funcional.
go_bp = gp.enrichr(
    gene_list=genes_go,
    gene_sets="GO_Biological_Process_2023",
    organism="human",
    background=universo_go,
    outdir=None
)

# Guardamos la tabla de resultados.
resultados_go = go_bp.results

resultados_go.to_csv(
    carpeta_resultados / "enriquecimiento_GO_BP.tsv",
    sep="\t",
    index=False
)

print("\n¡Enriquecimiento GO BP guardado correctamente!")

# ============================================================
# GRÁFICA DE LOS 10 PROCESOS GO BP MÁS SIGNIFICATIVOS
# ============================================================

# Seleccionamos los procesos con significación estadística.
go_significativos = resultados_go[
    resultados_go["Adjusted P-value"] < 0.05
].copy()

# Elegimos los 10 procesos con menor p-valor ajustado.
top10_go = go_significativos.sort_values(
    "Adjusted P-value"
).head(10).copy()

# Calculamos -log10 del p-valor ajustado.
top10_go["menos_log10_padj"] = -np.log10(
    top10_go["Adjusted P-value"].clip(lower=1e-300)
)

# Ordenamos los procesos para la gráfica horizontal.
top10_go = top10_go.sort_values("menos_log10_padj")

# Acortamos los nombres para que sean más legibles.
top10_go["Nombre_corto"] = top10_go["Term"].str.replace(
    r" \(GO:\d+\)$", "", regex=True
)

# Creamos la gráfica.
plt.figure(figsize=(12, 7))

plt.barh(
    top10_go["Nombre_corto"],
    top10_go["menos_log10_padj"],
    color="mediumseagreen"
)

plt.xlabel("-log10 (p-valor ajustado)")
plt.ylabel("Procesos biológicos")
plt.title("Top 10 procesos GO Biological Process - TCGA-BRCA")

plt.tight_layout()

# Guardamos la figura.
plt.savefig(
    carpeta_resultados / "go_bp_top10.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# Guardamos también la tabla de los 10 procesos.
top10_go.to_csv(
    carpeta_resultados / "go_bp_top10.tsv",
    sep="\t",
    index=False
)

print("\n¡Gráfica GO BP guardada correctamente!")

# ============================================================
# ENRIQUECIMIENTO FUNCIONAL REACTOME
# ============================================================

print("\nIniciando enriquecimiento Reactome...")

# Utilizamos los mismos genes seleccionados para GO BP.
genes_reactome = genes_go

# Utilizamos como referencia los genes analizados.
universo_reactome = universo_go

print(f"Genes seleccionados para Reactome: {len(genes_reactome)}")

# Realizamos el enriquecimiento funcional.
reactome = gp.enrichr(
    gene_list=genes_reactome,
    gene_sets="Reactome_2022",
    organism="human",
    background=universo_reactome,
    outdir=None
)

# Guardamos los resultados.
resultados_reactome = reactome.results

resultados_reactome.to_csv(
    carpeta_resultados / "enriquecimiento_Reactome.tsv",
    sep="\t",
    index=False
)

print("\n¡Enriquecimiento Reactome guardado correctamente!")

# ============================================================
# GRÁFICA DE LAS 10 RUTAS REACTOME MÁS SIGNIFICATIVAS
# ============================================================

print("\nGenerando gráfica de Reactome...")

# Seleccionamos las rutas estadísticamente significativas.
reactome_significativos = resultados_reactome[
    resultados_reactome["Adjusted P-value"] < 0.05
].copy()

# Elegimos las 10 rutas con menor p-valor ajustado.
top10_reactome = reactome_significativos.sort_values(
    "Adjusted P-value"
).head(10).copy()

# Calculamos -log10 del p-valor ajustado.
top10_reactome["menos_log10_padj"] = -np.log10(
    top10_reactome["Adjusted P-value"].clip(lower=1e-300)
)

# Ordenamos las rutas para la gráfica horizontal.
top10_reactome = top10_reactome.sort_values("menos_log10_padj")

# Eliminamos los identificadores Reactome de los nombres.
top10_reactome["Nombre_corto"] = top10_reactome["Term"].str.replace(
    r"\s+R-HSA-\d+$", "", regex=True
)

# Creamos la gráfica.
plt.figure(figsize=(12, 7))

plt.barh(
    top10_reactome["Nombre_corto"],
    top10_reactome["menos_log10_padj"],
    color="mediumpurple"
)

plt.xlabel("-log10 (p-valor ajustado)")
plt.ylabel("Rutas biológicas")
plt.title("Top 10 rutas Reactome - TCGA-BRCA")

plt.tight_layout()

# Guardamos la gráfica.
plt.savefig(
    carpeta_resultados / "reactome_top10.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# Guardamos también la tabla de las 10 rutas.
top10_reactome.to_csv(
    carpeta_resultados / "reactome_top10.tsv",
    sep="\t",
    index=False
)

print("\n¡Gráfica Reactome guardada correctamente!")

# ============================================================
# ENRIQUECIMIENTO FUNCIONAL SEPARADO POR DIRECCIÓN DEL CAMBIO
# ============================================================

print("\nIniciando enriquecimientos separados...")

# Seleccionamos genes significativamente sobreexpresados e infraexpresados.
genes_up = resultados[
    (resultados["padj"] < 0.05) &
    (resultados["log2FoldChange"] > 1)
].index.astype(str).tolist()

genes_down = resultados[
    (resultados["padj"] < 0.05) &
    (resultados["log2FoldChange"] < -1)
].index.astype(str).tolist()

print(f"Genes sobreexpresados: {len(genes_up)}")
print(f"Genes infraexpresados: {len(genes_down)}")

# ============================================================
# GO BIOLOGICAL PROCESS - GENES SOBREEXPRESADOS
# ============================================================

print("\nAnalizando GO BP de genes sobreexpresados...")

go_up = gp.enrichr(
    gene_list=genes_up,
    gene_sets="GO_Biological_Process_2023",
    organism="human",
    background=universo_go,
    outdir=None
)

resultados_go_up = go_up.results

resultados_go_up.to_csv(
    carpeta_resultados / "GO_BP_sobreexpresados.tsv",
    sep="\t",
    index=False
)

print("¡GO BP de genes sobreexpresados guardado correctamente!")

# ============================================================
# GO BIOLOGICAL PROCESS - GENES INFRAEXPRESADOS
# ============================================================

print("\nAnalizando GO BP de genes infraexpresados...")

go_down = gp.enrichr(
    gene_list=genes_down,
    gene_sets="GO_Biological_Process_2023",
    organism="human",
    background=universo_go,
    outdir=None
)

resultados_go_down = go_down.results

resultados_go_down.to_csv(
    carpeta_resultados / "GO_BP_infraexpresados.tsv",
    sep="\t",
    index=False
)

print("¡GO BP de genes infraexpresados guardado correctamente!")

# ============================================================
# REACTOME - GENES SOBREEXPRESADOS
# ============================================================

print("\nAnalizando Reactome de genes sobreexpresados...")

reactome_up = gp.enrichr(
    gene_list=genes_up,
    gene_sets="Reactome_2022",
    organism="human",
    background=universo_reactome,
    outdir=None
)

resultados_reactome_up = reactome_up.results

resultados_reactome_up.to_csv(
    carpeta_resultados / "Reactome_sobreexpresados.tsv",
    sep="\t",
    index=False
)

print("¡Reactome de genes sobreexpresados guardado correctamente!")

# ============================================================
# REACTOME - GENES INFRAEXPRESADOS
# ============================================================

print("\nAnalizando Reactome de genes infraexpresados...")

reactome_down = gp.enrichr(
    gene_list=genes_down,
    gene_sets="Reactome_2022",
    organism="human",
    background=universo_reactome,
    outdir=None
)

resultados_reactome_down = reactome_down.results

resultados_reactome_down.to_csv(
    carpeta_resultados / "Reactome_infraexpresados.tsv",
    sep="\t",
    index=False
)

print("¡Reactome de genes infraexpresados guardado correctamente!")