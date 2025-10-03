# ---------------------------
# Stage 1: Base image with Python 3.7
# ---------------------------
FROM python:3.10-slim

# Install system dependencies and Julia
RUN apt-get update && apt-get install -y wget tar xz-utils build-essential \
    && wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.9-linux-x86_64.tar.gz \
    && tar -xvzf julia-1.10.9-linux-x86_64.tar.gz -C /opt/ \
    && ln -s /opt/julia-1.10.9/bin/julia /usr/local/bin/julia \
    && rm julia-1.10.9-linux-x86_64.tar.gz \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------
# Install Python dependencies
RUN python -m pip install --upgrade pip && \
    python -m pip install \
        pandas numpy tqdm rdkit-pypi tensorflow==2.12 xgboost==1.5.2 joblib scikit-learn==1.6.1 

# ---------------------------
# Stage 2: Copy project files
# ---------------------------
WORKDIR /app

COPY Fingerprint_generation.jl .
COPY SA_scores.py .
COPY RA_score.py .
COPY main.jl .
COPY Ready_Biodegradability_model_filtered.joblib .
COPY RAscore ./RAscore
COPY fpscores.pkl.gz .
COPY sascorer.py .
COPY standalone_model_numpy.py .
COPY MACCS_keys.txt .
COPY model.ckpt-10654.as_numpy.json.gz .
COPY CNL_Headers.csv .
COPY CNL_pred.jl .
COPY Ready_Biodegradability_CNL_model_precurso.joblib .


# ---------------------------
# Stage 3: Julia dependencies
# ---------------------------
RUN julia -e 'using Pkg; \
    Pkg.add("CSV"); \
    Pkg.add("DataFrames"); \
    Pkg.add("ScikitLearn"); \
    Pkg.add("PyCall"); \
    Pkg.add("ProgressMeter"); \
    Pkg.add("Conda"); \
    Pkg.precompile()'

# Fix libmamba solver for Conda


# ---------------------------
# Stage 4: Entrypoint
# ---------------------------
ENTRYPOINT ["julia", "main.jl"]
