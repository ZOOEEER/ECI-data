# Workflow of data collection

## Files

| python file  | description |
| ------------ | ----------- |
| util_file.py |             |
| util_log.py  |             |
| util_func.py |             |
| util_prot.py |             |
| util_chem.py |             |



## Activity

The collation of activity data needs to solve the following problems: 

1. The number of enzymes and the determined order.
2. The number of compounds and the determined order.
3. The activity table of (enzyme, compound).

Active data may appear in the form of **long tables** or **wide tables**, and the two can be easily converted with `pandas` built-in functions.

```python

```

The processes mainly rely on the `pandas` package for data cleaning. For data provided in image format: activity values and compound 2D topology, the image processing package `cv2` was used for data cleaning. The variety of data forms makes it unpractical to define a general data processing pipeline. It is hoped that the data processing methods employed here will inspire related work in collecting data from the literature.

The framework provides the following support for data processing:

1. generate the `[raw]` folder.  (in `main.py`)
2. generate the `dev.ipynb` file in `[raw]` folder as a start. (in `main.py`)

Use automation as much as possible in each case, and common functions are defined in `util_func.py`.



## Enzyme

The acquisition of protein sequence information is a key step in the acquisition of enzyme information, followed by the generation of structure from the sequence. Before generating sequences for modeling, make sure the enzyme names and sequences are in place.

The development process is carried out in the raw folder, and the framework will provide the following support:

1. generate the  `[raw]/test/pdb` folder.  (in `main.py`)
2. provide functions to access NCBI database by `biopython` module. (`util_prot` module)
3. provide functions to access SwissModel web server by `requests` module. (`util_prot` module)
4. provide functions to generate and use the csv files for human curation.  (`util_prot` module)



## Chemical

If the raw data provided the pubchem_id(`cid`) of each compound, then collecting the data would simply call the API provided by pubchempy to download the 3D sdf file from the online database.

```python
import util_chem
util_chem.download_sdf(identifier, sdfpath)
```

However, it's much more complicated. It was an iterative process, and had to include human intervention, and here we finally present the data curation as a one-pass by documenting the human curation results as a csv file as a local database(local_db).

Before querying databases and collecting data manually, check compound names and classes, and do some filtering to build a clean compound list.

During data curation, `chemicals(pandas.DataFrames)` should be provided to `query_chemicals` function to try to fetch necessary information and download structural data from the pubchem database. The information will be fed into `chemicals`. It will return the missing compound index as a list.

```python
no_hits = util_chem.query_chemicals(chemicals) # A simple example
```

Then generate the `chemicals_curated.csv` file in the current folder.

```python
util_chem.make_local_db(chemicals)
```

Edit the file to complement the missing information, especially the `cid`  and `sdf`.

In fact, there is no guarantee that the 3D structure of every compound could be fetched from pubchem, especially for polymer mixtures and macromolecular substances. Here, you need to implement an index filter to exclude these molecules.

For molecules whose 3D structure cannot be obtained from the pubchem database, if necessary, its 3D structure must be obtained through modeling software, such as GaussView. A recommended place to save these manually modeled molecule files are the `[clean]/sdf` folder, as the final dataset folder.

After you make the curation of how to replace the info(`name`,`cid`,`SMILES`,`sdf`), call the function.

```python
util_chem.query_local(chemicals)
```

The development process (**dev**) for compound information compilation described above will determine that the required information is available. Then, a forward production process (**prod**) should be implemented. The process reads the aforementioned literature-sourced or manually curated files to generate the final compound information data in one go, which is finally available in the [clean] folder.

The development process is done in the `[raw]` folder. The  framework will provide support to:

1. generate the `[raw]`, `[raw]/test/sdf` folder.  (in `main.py`)
2. provide functions to query Pubchem by `pubchempy` module. (`util_chem` module)
3. provide functions to generate and use the csv files for human curation.  (`util_chem` module)
4. provide functions to do OSAR by  API and `request` module. (`util_chem` module)