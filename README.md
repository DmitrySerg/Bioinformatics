# Bioinformatics [summer school](http://contest.bioinf.me/summer2017) MIPT 2017

## Single Cell Project

В данном проекте участникам предлагалось провести анализ данных single-cell RNA-Seq от 10x Genomics - 10xgenomics.com/single-cell-gene-expression/datasets/. Необходимо было провести кластеризацию данных экспрессии генов, определить присутствующие клеточные популяции и построить классифицирующие модели для предсказания по значениям экспресии генов в клетке ее принадлежность к определенной популяции.

## Data Preparation

Первичный анализ заключался в разделении ридов по клеткам и получении матрицы экспрессии генов при помощи CellRanger. Данные были отфильтрованы по содержанию митохондриальных последовательностей, а также были убраны выбросы по экспресии. При помощи PCA размерность данных была существенно снижена, статистически значимой дисперсия была до 15-й компоненты включительно, поэтому в дальнейшем использовались только они. 

![](https://github.com/DmitrySerg/Bioinformatics/blob/master/Single_cell/Pictures/Heatmap6cluster.png)

С кластеризацией данных успешно справился Shared Nearest Neighbor (SNN), а с красивой визуализацией помог tSNE. 

![](https://github.com/DmitrySerg/Bioinformatics/blob/master/Single_cell/Pictures/tSNE%20annotated.png)

## Machine Learning

Наконец, получив метки классов после интерпретации кластеров, были построены классификаторы, способные по сырым данным экспресии (без фильтрации, предобработки и PCA) предсказывать метку класса.

![](https://github.com/DmitrySerg/Bioinformatics/blob/master/Single_cell/Pictures/quality.png)

![](https://github.com/DmitrySerg/Bioinformatics/blob/master/Single_cell/Pictures/feature_importance.png)

## Feature importance

Немного хороших ссылок, которые мы использовали
https://bioconductor.org/packages/devel/bioc/vignettes/SC3/inst/doc/my-vignette.html
http://satijalab.org/seurat/pbmc3k_tutorial.html
https://hemberg-lab.github.io/scRNA.seq.course/data-visualization-reads.html