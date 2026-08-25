# Harness de validación

Reproduce lo que hace SARA **sin abrir la interfaz**, para poder automatizar los
resultados y llevar a la aplicación solo lo que ya se sabe que sale bien.

Llama a las mismas funciones que el servidor de Shiny —`run_deg()`,
`run_enrichment_*()`, `run_gsea()`— y al mismo `workflow.sh`. Lo que se valida
es la aplicación, no un script paralelo que podría diferir en un detalle.

Que esto sea posible no es casualidad: SARA separa las funciones puras del
código que depende de `input`/`session`. El harness consume esa capa directamente.

## Uso

```bash
# Un dataset, todas las etapas
Rscript scripts/harness/correr.R GSE273773

# Todos los configurados
Rscript scripts/harness/correr.R --todos

# Solo algunas etapas
Rscript scripts/harness/correr.R GSE273773 --etapas deg,metricas

# Rehacer lo ya hecho
Rscript scripts/harness/correr.R GSE273773 --forzar

# Solo regenerar las tablas de defensa
Rscript scripts/harness/correr.R --tablas
```

Arráncalo con el entorno de herramientas activo, igual que la aplicación:

```bash
export PATH="$HOME/miniforge3/envs/rnaseq_ecoli/bin:$PATH"
```

## Etapas

| Etapa | Qué hace |
|---|---|
| `pipeline` | Lanza `workflow.sh` por cada ruta configurada |
| `deg` | Expresión diferencial con los tres motores, por cada ruta |
| `enriquecimiento` | ORA sobre las colecciones configuradas, más GSEA |
| `metricas` | Concordancia con lo publicado, acuerdo interno, permutación y descomposición del error |

Cada etapa **se salta si ya está hecha**, así que una ejecución interrumpida se
reanuda sin repetir lo caro. `--forzar` las rehace.

## Salida

```
validacion/
├── <dataset>/
│   ├── pipeline/<ruta>/      salida de workflow.sh
│   ├── deg/<ruta>__<motor>.rds
│   ├── enriquecimiento/<coleccion>.rds
│   └── metricas/*.rds
└── tablas/
    ├── T1..T8.md             para pegar en la memoria
    ├── T1..T8.tsv            para reprocesar
    └── TODAS.md              las ocho seguidas
```

## Las ocho tablas

Cada una responde a una pregunta que un tribunal puede hacer:

| | Pregunta |
|---|---|
| T1 | ¿Con qué datos has trabajado? |
| T2 | ¿Qué coste tiene ejecutar esto? |
| T3 | ¿Coincides con lo publicado? |
| T4 | ¿Coinciden los motores entre sí? |
| T5 | ¿Coinciden las rutas del pipeline entre sí? |
| T6 | ¿Cómo sabes que no estás inventando señal? |
| T7 | ¿De dónde viene el desacuerdo que queda? |
| T8 | ¿Con qué versiones exactas? |

T4 es la única que valida la aplicación y no el trabajo ajeno: compara los tres
motores sobre la misma matriz, de modo que la única variable es el motor.

T6 es la única que mide **especificidad**. Las demás comparan con lo publicado y
por tanto no distinguen un pipeline correcto de uno que invente señal de forma
reproducible.

## Añadir un dataset

Se añade una entrada a `datasets.R` y nada más: el harness sabe ejecutarla. Los
campos están documentados en la cabecera de ese fichero.

Dos cosas que conviene no equivocar:

- **La referencia depende de la ruta.** Las de alineamiento (`bowtie2`,
  `subjunc`) van contra el **genoma**; las de pseudoalineamiento (`salmon`,
  `kallisto`) contra el **transcriptoma**. El parámetro de `workflow.sh` se
  llama igual en los dos casos, así que es el error más fácil de cometer.
- **`bowtie2` solo es legítimo sin intrones.** En un eucariota con intrones usa
  `subjunc`, que sí reconoce las uniones exón-exón, o pseudoalineamiento.
