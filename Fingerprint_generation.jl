# maccs_scores.jl (Fingerprint_generation.jl)
using PyCall, CSV, DataFrames, ProgressMeter

# Import RDKit from Python
Chem = pyimport("rdkit.Chem")
MACCSkeys = pyimport("rdkit.Chem.MACCSkeys")

# --- Command line args ---
if length(ARGS) < 2
    println("Usage:")
    println("  julia Fingerprint_generation.jl input.csv output.csv")
    println("  julia Fingerprint_generation.jl \"C1=CC=CC=C1\" output.csv")
    println("  julia Fingerprint_generation.jl \"C1=CC=CC=C1,C1CCCCC1\" output.csv")
    exit(1)
end


input_arg = ARGS[1]
output_file = ARGS[2]

# --- Load SMILES ---
if contains(input_arg, "temp")
    df = CSV.read(input_arg, DataFrame)   # <--- use input_arg here
    SMILES = collect(df.SMILES)
else
    smiles_file = joinpath("/app/out", basename(input_arg))
    println(">>> Reading SMILES from file: $smiles_file")
    df = CSV.read(smiles_file, DataFrame)   # <--- use input_arg here
    SMILES = collect(df.SMILES)
end


# --- Load MACCS keys names ---
MACCS_keys = readlines("MACCS_keys.txt")[2:end]
rows = DataFrame(Int.(zeros(0,166)), MACCS_keys)
smiles_list = String[]
errored_smiles = String[]

# --- Batch processing ---
@showprogress for (i, smile) in enumerate(SMILES)
    try
        mol = Chem.MolFromSmiles(smile)
        if mol == PyCall.PyNULL()
            push!(errored_smiles, smile)
            continue
        end

        # --- MACCS Keys ---
        maccs = MACCSkeys.GenMACCSKeys(mol)
        maccs_array = [Int(maccs.GetBit(i)) for i in 1:166]

        push!(rows, maccs_array)
        push!(smiles_list, smile)

    catch e_inner
        @warn "Error processing SMILES: $smile" exception = e_inner
        push!(errored_smiles, smile)
    end
end

# Insert SMILES as first column
rows = insertcols!(rows, 1, :SMILES => smiles_list)

# Save output
CSV.write(output_file, rows)

println("Done! Errors for these SMILES: ", errored_smiles)
