using CSV, DataFrames
using ScikitLearn
import ScikitLearn: predict_proba
using PyCall


# Python deps
sklearn = pyimport("sklearn")
joblib = pyimport("joblib")
np = pyimport("numpy")
pybuiltins = pyimport("builtins")
# --------------------------
#Python path
path_python = (PyCall.python)
path_python_py37 = "/usr/local/bin/python"

# --------------------------
# Main pipeline
# --------------------------



function main_SMILES(input_arg::String) 
     # Determine if input is a CSV file
    if endswith(lowercase(input_arg), ".csv")
        # Treat as CSV
        smiles_file = input_arg
          # Use full path for child Julia
    else
        # Single SMILES or comma-separated
        smiles_list = split(input_arg, ",")
        smiles_file = "temp_smiles.csv"
        df = DataFrame(SMILES = smiles_list)
        CSV.write(smiles_file, df)
        smiles_file = abspath(smiles_file)  # ensure child process sees it
    end

    

    # Step 1: Fingerprint generation
    println(">>> Running fingerprint generation...")
    run(`julia "Fingerprint_generation.jl" $smiles_file MACCS.csv`)

    # Step 2: SA scores
    println(">>> Running SA score calculation...")
    run(`$path_python SA_scores.py MACCS.csv combined_scores.csv`)

    # Step 3: RA scores
    println(">>> Running RA score calculation...")
    run(`$path_python_py37 RA_score.py combined_scores.csv final_results.csv`)

    # Step 4: Final model prediction
    println(">>> Running biodegradability prediction...")
    results = CSV.read("final_results.csv", DataFrame)

    desired_order = [
    "SA_score",
    "SCS_score",
    "RA_NN_score",
    "RA_XGB_score",
    "MC1_score",
    "MC2_score"]

    Model_data_correct_order = hcat(DataFrame(SMILES = results.SMILES),results[:,end-5:end][:,desired_order],results[:,2:end-6])

    Model_matrix = Matrix(Model_data_correct_order[:,2:end])

    Model_matrix_clean = coalesce.(Model_matrix, 0.0)  # replace missing with 0.0 (or use imputation)

    Model_matrix_np = np.array(Model_matrix_clean, dtype=pybuiltins.float)

    rf = joblib.load("Ready_Biodegradability_model_filtered.joblib")

    RB_proba = predict_proba(rf, Model_matrix_np)[:, 2]

    optimal_threshold = 0.46329756296676594

    RB_pred = map(p -> p >= optimal_threshold ? 1 : 0, RB_proba)
    
    max_probas_RB = [maximum(row) for row in eachrow(ScikitLearn.predict_proba(rf,(Model_matrix_np)))]

    RB = DataFrame(SMILES = Model_data_correct_order.SMILES,
                   Persistence = RB_pred,
                   Probability = max_probas_RB)

        # ------------------ SAVE TO MOUNTED FOLDER ------------------
    output_dir = "/app/out"  # must match docker -v mount
    CSV.write(joinpath(output_dir, "RB_predicted.csv"), RB)
    println(">>> Prediction complete! Files saved to folder: $output_dir")
end

function main_CNL(input_arg::String)
    output_dir = "/app/out"
    output_file = joinpath(output_dir, "CNL_predicted.csv")
    println(">>> Running CNL model...")
    run(`julia "CNL_pred.jl" $(joinpath("/app/out", basename(input_arg))) $output_file`)
    println(">>> Prediction complete! Results saved to $output_file")
end

# --------------------------
# CLI entrypoint
# --------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage:")
        println("  julia main.jl SMILES input.csv")
        println("  julia main.jl SMILES \"CCO\"")
        println("  julia main.jl CNL spectrum.msp")
        println("  julia main.jl CNL spectrum.csv")
        exit(1)
    end

    method = uppercase(ARGS[1])
    input_arg = ARGS[2]

    if method == "SMILES"
        main_SMILES(input_arg)
    elseif method == "CNL"
        main_CNL(input_arg)
    else
        println("Unknown method: $method. Use SMILES or CNL.")
        exit(1)
    end
end
