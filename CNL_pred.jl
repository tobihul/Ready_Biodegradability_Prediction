using CSV, DataFrames
using ScikitLearn
import ScikitLearn: predict_proba
using PyCall

joblib = pyimport("joblib")
rf_CNL = joblib.load("Ready_Biodegradability_CNL_model_precurso.joblib")
youden = 0.5423759954672558

cnl_fp_headers = CSV.read("CNL_Headers.csv", DataFrame)[:,1]

# --- parsing function ---
function parse_msp_simple(file::String)
    spectra = []
    entry_id = 0
    name = ""
    peaks = Float64[]
    intensities = Float64[]
    precursor_mz = missing

    open(file, "r") do io
        for line in eachline(io)
            s = strip(line)
            if isempty(s)
                if !isempty(peaks)
                    entry_id += 1
                    if isempty(name)
                        name = "ID_$entry_id"
                    end
                    push!(spectra, (name=name, precursor_mz=precursor_mz,
                                    mz=peaks, intensity=intensities))
                    peaks = Float64[]
                    intensities = Float64[]
                    name = ""
                    precursor_mz = missing
                end
                continue
            end

            if startswith(s, "Name:")
                name = strip(split(s, ":", limit=2)[2])
            elseif startswith(s, "PrecursorMZ:")
                precursor_mz = parse(Float64, strip(split(s, ":", limit=2)[2]))
            elseif startswith(s, "Num Peaks:")
                # ignore
            elseif occursin(r"^\d", s)
                parts = split(s)
                if length(parts) == 2
                    push!(peaks, parse(Float64, parts[1]))
                    push!(intensities, parse(Float64, parts[2]))
                end
            end
        end
        if !isempty(peaks)
            entry_id += 1
            if isempty(name)
                name = "ID_$entry_id"
            end
            push!(spectra, (name=name, precursor_mz=precursor_mz,
                            mz=peaks, intensity=intensities))
        end
    end
    return spectra
end

# --- fragment → CNL vector ---
function frags_to_cnls_vec(MZ_prec_ion, MZ_frags, cnl_fp_headers)
    tolerance = 0.01
    cnls = MZ_prec_ion .- MZ_frags
    cnls_fp = zeros(length(cnl_fp_headers))
    [cnls_fp[abs.(cnl_fp_headers .- cnls[i]) .<= tolerance+0.0001] .= 1 for i in eachindex(cnls)]
    return cnls_fp
end


# --- main logic ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia CNL_pred.jl spectrum.msp [output.csv]")
        exit(1)
    end

    input_arg = ARGS[1]
    output_file = length(ARGS) > 1 ? ARGS[2] : nothing

    if endswith(lowercase(input_arg), ".msp")
        println(">>> Processing spectrum/spectra from .msp file: $input_arg")
        parsed_file = parse_msp_simple(joinpath("/app/out", basename(input_arg)))

        CNLs = zeros((length(parsed_file), 2441))
        for i in eachindex(parsed_file)
            mz_prec = parsed_file[i][:precursor_mz]
            mz_frags = parsed_file[i][:mz]
            CNL = frags_to_cnls_vec(mz_prec, mz_frags, cnl_fp_headers)
            CNL = vcat(CNL, mz_prec)
            CNLs[i,:] = CNL
        end

        proba_pred = ScikitLearn.predict_proba(rf_CNL,(CNLs))[:,2]
        max_probas = [maximum(row) for row in eachrow(ScikitLearn.predict_proba(rf_CNL,(CNLs)))]
        IDs = [entry[:name] for entry in parsed_file]
        predictions = map(p -> p >= youden ? 1 : 0, proba_pred)

        CNL_pred = DataFrame(NAME = IDs, Persistence = predictions, Probability = max_probas)

        if output_file !== nothing
            CSV.write(output_file, CNL_pred)
            println(">>> Results saved in $output_file")
        else
            println(CNL_pred)
        end

    else
        println(">>> Processing spectrum from .csv file: $input_arg")
        data = CSV.read(joinpath("/app/out", basename(input_arg)), DataFrame)
        CNLs = frags_to_cnls_vec(data[1,:Precursor_mz], data.mz, cnl_fp_headers)
        CNLs = vcat(CNLs, data[1,:Precursor_mz])

        proba_pred = ScikitLearn.predict_proba(rf_CNL,([CNLs]))[:,2]
        max_probas = [maximum(row) for row in eachrow(ScikitLearn.predict_proba(rf_CNL,([CNLs])))]

        prediction = map(p -> p >= youden ? 1 : 0, proba_pred)

        if prediction[1] == 1
            Persistence = "persistent"
        else
            Persistence = "non-persistent"
        end

        println(">>> The query chemical is predicted as $Persistence (prob = $(max_probas[1]))")
    end
end
