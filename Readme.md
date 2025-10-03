# Biodeg Prediction Docker Image

This Docker image allows you to run SMILES or MS/MS spectra through a Ready biodegradation prediction model.

It is important to note that the model specifically predicts if a substance is persistent. If it is predicted as persistent (1) it means it is not predicted as being readily biodegradable and vice-versa.

## Pulling the image

First, make sure you have [Docker installed](https://docs.docker.com/get-docker/).

Then pull the image from Docker Hub:

```bash
docker pull tobiashulleman/biodeg_pred_repo:latest
```

## Running the model 

### Run with a SMILES CSV file

Mount a local folder with your input file into the container and run the model:

```bash
docker run --rm -v "/path/to/your/folder":/app/out \
    tobiashulleman/biodeg_pred_repo:latest SMILES "your_smiles.csv"
```
or using a single SMILES string

```bash
docker run --rm -v "/path/to/your/folder":/app/out \
    tobiashulleman/biodeg_pred_repo:latest SMILES "CCO"
```

or multiple SMILES

```bash
docker run --rm -v "/path/to/your/folder":/app/out \
    tobiashulleman/biodeg_pred_repo:latest SMILES "CCO,CCCl"
```
* /path/to/your/folder → replace with the folder on your computer that has your SMILES file.
* your_smiles.csv → the input CSV file containing SMILES.
* Output will be saved back into the same folder.

## Run MS2 spectra in a .msp file

```bash
docker run --rm -v "/path/to/your/folder":/app/out \
    tobiashulleman/biodeg_pred_repo:latest CNL "your_spectra.msp"
```
* /path/to/your/folder → replace with the folder on your computer that has your .msp file.
* your_spectra.ms → the input .msp MS2 data file
* Output will be saved back into the same folder.

## Run MS2 spectra in a .csv file

```bash
docker run --rm -v "C:\Users\YourName\Downloads":/app/out \
    tobiashulleman/biodeg_pred_repo:latest SMILES "molecules.csv"
```

## Notes

* Use quotes around file names if they contain spaces.

* Always mount the folder containing the input file, not just the file itself.

* The output file will be saved in the same folder you mounted.


