using ScikitLearn, CSV, DataFrames
import ScikitLearn: RandomForestClassifier

@sk_import ensemble: RandomForestClassifier

sklearn = pyimport("sklearn")
joblib = pyimport("joblib")

##Prepare Data

calculated_data = CSV.read("final_results.csv", DataFrame)

desired_order = [
    "SA_score",
    "SCS_score",
    "RA_NN_score",
    "RA_XGB_score",
    "MC1_score",
    "MC2_score"
]


Model_data_correct_order = hcat(DataFrame(SMILES = calculated_data.SMILES),calculated_data[:,end-5:end][:,desired_order],calculated_data[:,2:end-6])

Model_matrix = Matrix(Model_data_correct_order[:,2:end])

rf = joblib.load("Ready_Biodegradability_model_filtered.joblib")

RB_proba = ScikitLearn.predict_proba(rf, Matrix(Model_matrix))[:,2]
optimal_threshold = 0.46329756296676594

RB_pred = map(p -> p >= optimal_threshold ? 1 : 0, RB_proba)

RB_pred = 1 .- RB_pred

RB = DataFrame(SMILES = MACCS.SMILES, Ready_Biodegradability = RB_pred, Probability = RB_proba)

CSV.write("RB_predicted.csv", RB)