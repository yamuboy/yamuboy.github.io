# = an============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import pymc as pm
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

print(f"Running on PyMC v{pm.__version__}")
print("Setup and Imports Complete.")

# =============================================================================
# STEP 2: LOAD AND PREPARE THE DATA
# =============================================================================
print("\n--- Loading data and creating advanced features ---")
try:
    df = pd.read_csv('qb_comprehensive_advanced_stats_2020_2024.csv')
except FileNotFoundError:
    print("Error: 'qb_comprehensive_advanced_stats_2020_2024.csv' not found.")
    exit()

df = df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)

# Use the same 3-week rolling average features for a direct comparison
rolling_features = [
    'dk_fantasy_points', 'total_yards', 'epa_per_play', 'total_touchdowns',
    'fantasy_points_per_dropback', 'passer_rating'
]
for feature in rolling_features:
    df[f'{feature}_rolling_3wk_avg'] = df.groupby('player_display_name')[feature].transform(
        lambda x: x.rolling(window=3, min_periods=1).mean().shift(1)
    )

df.dropna(inplace=True)

# Define the feature set
feature_columns = [
    'is_home', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    'rolling_4wk_avg_qb_fantasy_points_allowed', 'rolling_4wk_avg_qb_sacks_forced', 'rolling_4wk_avg_qb_interceptions_forced',
    'dk_fantasy_points_rolling_3wk_avg', 'total_yards_rolling_3wk_avg', 'epa_per_play_rolling_3wk_avg',
    'total_touchdowns_rolling_3wk_avg', 'fantasy_points_per_dropback_rolling_3wk_avg', 'passer_rating_rolling_3wk_avg'
]

# Standardize features for better model stability (important for Bayesian models)
for col in feature_columns:
    mean = df[col].mean()
    std = df[col].std()
    df[col] = (df[col] - mean) / std

X = df[feature_columns].values
y_dk = df['dk_fantasy_points'].values

print(f"Data prepared with {len(feature_columns)} features.")


# =============================================================================
# STEP 3: DEFINE AND TRAIN THE BAYESIAN MODEL
# =============================================================================
print("\n--- Defining and 'training' (sampling from) the Bayesian model ---")

with pm.Model() as bayesian_model:
    # --- Priors (Our initial beliefs) ---
    # Intercept: A reasonable prior for an average QB's score
    intercept = pm.Normal('intercept', mu=y_dk.mean(), sigma=10)
    
    # Coefficients for each feature. We believe they are centered around 0.
    coeffs = pm.Normal('coeffs', mu=0, sigma=2, shape=X.shape[1])
    
    # Model error/noise
    sigma = pm.HalfNormal('sigma', sigma=5)

    # --- Data Container for Prediction ---
    # This allows us to swap in Mahomes' data later
    X_data = pm.Data('features', X, mutable=True)
    
    # --- Linear Model (The core equation) ---
    mu = intercept + pm.math.dot(X_data, coeffs)
    
    # --- Likelihood (Connecting the model to the data) ---
    y_obs = pm.Normal('y_obs', mu=mu, sigma=sigma, observed=y_dk)
    
    # --- MCMC Sampling (This is the 'training' process) ---
    # This will take longer than the LightGBM model
    idata = pm.sample(2000, tune=1000, cores=1, target_accept=0.9)

print("Bayesian model sampling complete.")


# =============================================================================
# STEP 4: MAKE THE BAYESIAN PREDICTION
# =============================================================================
print("\n--- Preparing prediction data for Patrick Mahomes, Week 16 2024 ---")
PREDICTION_YEAR = 2024
PREDICTION_WEEK = 16
PLAYER_NAME = 'Patrick Mahomes'

# Get the original, unscaled data for Week 15
prediction_input_df_unscaled = df[
    (df['player_display_name'] == PLAYER_NAME) &
    (df['season'] == PREDICTION_YEAR) &
    (df['week'] == PREDICTION_WEEK - 1)
]

if prediction_input_df_unscaled.empty:
    print(f"\nError: Could not find data for {PLAYER_NAME} in Week {PREDICTION_WEEK - 1}, {PREDICTION_YEAR}.")
else:
    # Get the scaled features for the prediction
    X_predict_scaled = prediction_input_df_unscaled[feature_columns].values

    # --- Generate the Posterior Predictive Distribution ---
    print("\n--- Generating full distribution of outcomes for Mahomes ---")
    with bayesian_model:
        # Swap in Mahomes' scaled Week 15 data
        pm.set_data({'features': X_predict_scaled})
        
        # Generate 4000 possible scores for next week based on our model's learned uncertainty
        posterior_pred = pm.sample_posterior_predictive(idata)

    # Extract the simulated scores
    simulated_scores = posterior_pred.posterior_predictive['y_obs'].values.flatten()

    # Calculate the mean, which is our new point estimate
    mahomes_dk_pred_bayesian = simulated_scores.mean()
    mahomes_fd_pred_bayesian = mahomes_dk_pred_bayesian * (df['fd_fantasy_points'].mean() / df['dk_fantasy_points'].mean()) # Simple ratio for FD
    mahomes_captain_pred_bayesian = mahomes_dk_pred_bayesian * 1.5
    
    # Get the credible interval (our floor/ceiling)
    hdi = az.hdi(simulated_scores, hdi_prob=0.8) # 80% credible interval

    print("\n" + "="*50)
    print(f"BAYESIAN FANTASY PREDICTIONS FOR: {PLAYER_NAME}")
    print(f"WEEK: {PREDICTION_WEEK}, SEASON: {PREDICTION_YEAR}")
    print("="*50)
    print(f"Actual Score:                     23.70")
    print("-" * 50)
    print(f"Bayesian Mean DK Prediction:      {mahomes_dk_pred_bayesian:.2f}")
    print(f"Bayesian Mean FD Prediction:      {mahomes_fd_pred_bayesian:.2f}")
    print(f"Bayesian Mean DK Captain:         {mahomes_captain_pred_bayesian:.2f}")
    print("-" * 50)
    print(f"80% Credible Interval (Floor/Ceiling): [{hdi[0]:.2f}, {hdi[1]:.2f}]")
    print(f"This means the model is 80% confident the true score will fall in this range.")
    print("="*50)

    # --- Visualize the Prediction ---
    plt.figure(figsize=(12, 7))
    az.plot_dist(simulated_scores, kind='hist', hist_kwargs={'bins': 50}, label="Predicted Score Distribution")
    plt.axvline(mahomes_dk_pred_bayesian, color='red', linestyle='--', label=f'Mean Prediction: {mahomes_dk_pred_bayesian:.2f}')
    plt.axvline(23.7, color='green', linestyle='-', linewidth=2, label=f'Actual Score: 23.70')
    plt.title(f"Bayesian Prediction Distribution for Patrick Mahomes", fontsize=16)
    plt.xlabel("Predicted DraftKings Points", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.legend()
    plt.show()