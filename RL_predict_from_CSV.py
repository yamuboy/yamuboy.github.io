# =============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import lightgbm as lgb
import warnings

# Suppress common warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

print("Setup and Imports Complete.")

# =============================================================================
# STEP 2: LOAD AND PREPARE THE DATA
# =============================================================================
print("\n--- Loading the comprehensive QB dataset ---")
try:
    # Load the master QB stats file you created
    df = pd.read_csv('qb_comprehensive_advanced_stats_2020_2024.csv')
except FileNotFoundError:
    print("Error: 'qb_comprehensive_advanced_stats_2020_2024.csv' not found.")
    print("Please make sure you have run the previous script to generate this file.")
    exit()

# --- Feature Engineering ---
# Create lag features. A player's performance last week is often a good predictor for this week.
df = df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)
df['dk_points_lag_1'] = df.groupby('player_display_name')['dk_fantasy_points'].shift(1)
df['fd_points_lag_1'] = df.groupby('player_display_name')['fd_fantasy_points'].shift(1)
df['total_yards_lag_1'] = df.groupby('player_display_name')['total_yards'].shift(1)
df['epa_per_play_lag_1'] = df.groupby('player_display_name')['epa_per_play'].shift(1)

# Drop rows with NaN values created by the lag (i.e., the first game for each player)
df.dropna(inplace=True)

# --- Define Features (X) and Targets (y) ---
# These are the variables our model will use to make predictions
feature_columns = [
    # Betting & Context
    'is_home', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    # Opponent Rolling Stats
    'rolling_4wk_avg_qb_fantasy_points_allowed', 'rolling_4wk_avg_qb_sacks_forced', 'rolling_4wk_avg_qb_interceptions_forced',
    # QB Lagged Performance
    'dk_points_lag_1', 'fd_points_lag_1', 'total_yards_lag_1', 'epa_per_play_lag_1'
]

# Create our training data
X = df[feature_columns]
y_dk = df['dk_fantasy_points'] # Target for DraftKings
y_fd = df['fd_fantasy_points'] # Target for FanDuel

print(f"Data prepared. Training model on {len(X)} samples.")
print(f"Using the following features: {feature_columns}")

# =============================================================================
# STEP 3: TRAIN THE PREDICTION MODELS
# =============================================================================
print("\n--- Training DraftKings and FanDuel prediction models ---")

# --- DraftKings Model ---
lgbm_dk = lgb.LGBMRegressor(
    objective='regression_l1', # L1 loss is often more robust to outliers
    n_estimators=1000,         # Number of trees
    learning_rate=0.05,
    random_state=42,
    n_jobs=-1,
    verbose=-1 # Suppress verbose output
)
lgbm_dk.fit(X, y_dk)
print("DraftKings model trained successfully.")

# --- FanDuel Model ---
lgbm_fd = lgb.LGBMRegressor(
    objective='regression_l1',
    n_estimators=1000,
    learning_rate=0.05,
    random_state=42,
    n_jobs=-1,
    verbose=-1
)
lgbm_fd.fit(X, y_fd)
print("FanDuel model trained successfully.")

# =============================================================================
# STEP 4: PREPARE THE PREDICTION DATA FOR MAHOMES
# =============================================================================
print("\n--- Preparing prediction data for Patrick Mahomes, Week 16 2024 ---")

# We need the data from Week 15 to predict Week 16
# The "rolling" stats for week 16 are based on the games up to and including week 15.
# The "lag" stats for week 16 are the actual stats from the week 15 game.

PREDICTION_YEAR = 2024
PREDICTION_WEEK = 16
PLAYER_NAME = 'Patrick Mahomes'

# Find the data for the week *before* our target prediction week
prediction_input_df = df[
    (df['player_display_name'] == PLAYER_NAME) &
    (df['season'] == PREDICTION_YEAR) &
    (df['week'] == PREDICTION_WEEK - 1)
]

if prediction_input_df.empty:
    print(f"\nError: Could not find data for {PLAYER_NAME} in Week {PREDICTION_WEEK - 1}, {PREDICTION_YEAR}.")
    print("Cannot make a prediction. Please ensure the CSV contains this data.")
else:
    # Get the features for the prediction.
    # Note: We are using the rolling stats from week 15 as inputs for week 16.
    # We must also get the actual performance from week 15 to serve as the lag features.
    X_predict = prediction_input_df[feature_columns].copy()

    # Manually update the lag features with the actual Week 15 performance
    X_predict['dk_points_lag_1'] = prediction_input_df['dk_fantasy_points'].values[0]
    X_predict['fd_points_lag_1'] = prediction_input_df['fd_fantasy_points'].values[0]
    X_predict['total_yards_lag_1'] = prediction_input_df['total_yards'].values[0]
    X_predict['epa_per_play_lag_1'] = prediction_input_df['epa_per_play'].values[0]

    # =============================================================================
    # STEP 5: MAKE AND DISPLAY THE PREDICTIONS
    # =============================================================================
    print("\n--- Generating Predictions ---")

    # Predict using the trained models
    mahomes_dk_pred = lgbm_dk.predict(X_predict)[0]
    mahomes_fd_pred = lgbm_fd.predict(X_predict)[0]
    mahomes_captain_pred = mahomes_dk_pred * 1.5

    print("\n" + "="*50)
    print(f"FANTASY PREDICTIONS FOR: {PLAYER_NAME}")
    print(f"WEEK: {PREDICTION_WEEK}, SEASON: {PREDICTION_YEAR}")
    print("="*50)
    print(f"Predicted DraftKings Points:      {mahomes_dk_pred:.2f}")
    print(f"Predicted FanDuel Points:         {mahomes_fd_pred:.2f}")
    print(f"Predicted DK Showdown Captain:    {mahomes_captain_pred:.2f}")
    print("="*50)