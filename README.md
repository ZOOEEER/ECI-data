# ECI-data
ECI stands for Enzyme-Chemical Interaction. The problem formulation is as below:

![problem_formulation](media\problem_formulation.png)

Datasets suitable for ECI tasks come from the literature performing high-throughput experiments or large-scale characterization.

**The high goal of ECI is to develop models which could use data of the activity of multiple compounds on related enzymes to enhance predictions of the activity of these compounds on the focused enzyme.** The problem formulation of ECI is similar to drug-target interaction(DTI), but there are some significant differences between the two:  DTI is a matrix completion problem where the problem is information insufficiency. ECI is a matrix factorization problem where the problem is information overload.



## Usage

Show all the datasets' information.

```python

```

Show the statistics for one dataset.

```python

```

Show the meta data for one dataset.

```python

```

All datasets are provided as `csv` files, so they can be easily used for your own projects.



## Data Format

Four files are provided for each dataset.

| file            | description                                                  |
| --------------- | ------------------------------------------------------------ |
| `enzymes.csv`   | Name, Sequence and other properties                          |
| `chemicals.csv` | Name, SMILES, sdf file path and other properties             |
| `activity.csv`  | raw values from literature                                   |
| `metadata.json` | the original data collection strategy, modeling strategy, and enzyme catalysis issues to be explored |

## Metadata description

This repository not only focuses on the data itself, but also focuses on the original data collection strategy, modeling strategy, and enzyme catalysis issues to be explored. The definition of `metadata.json` is listed in the table below. It is mainly divided into four categories: basic information, experiments, data, and models. One could call the function to show the help.

```python

```



| #    |      |      |
| ---- | ---- | ---- |
| 1    |      |      |
| 2    |      |      |
| 3    |      |      |
| 4    |      |      |
| 5    |      |      |
| 6    |      |      |
| 7    |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |
|      |      |      |



## Datasets



## Data Curation





## Related Sources

A very good beginning for ECI problem. A similar work based on different starting point:

[1] [samgoldman97/enzyme-datasets: Enzyme datasets used to benchmark enzyme-substrate promiscuity models (github.com)](https://github.com/samgoldman97/enzyme-datasets)

A very good introduction to DTI problem:

[2] [alizat/Chemogenomic-DTI-Prediction-Methods: Algorithms for prediction of drug-target interactions via computational (chemogenomic) methods (github.com)](https://github.com/alizat/Chemogenomic-DTI-Prediction-Methods)

