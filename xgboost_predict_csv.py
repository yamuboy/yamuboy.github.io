# =============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

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

# Use the same advanced 3-week rolling average features
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
all_feature_columns = [
    'is_home', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    'rolling_4wk_avg_qb_fantasy_points_allowed', 'rolling_4wk_avg_qb_sacks_forced', 'rolling_4wk_avg_qb_interceptions_forced',
    'dk_fantasy_points_rolling_3wk_avg', 'total_yards_rolling_3wk_avg', 'epa_per_play_rolling_3wk_avg',
    'total_touchdowns_rolling_3wk_avg', 'fantasy_points_per_dropback_rolling_3wk_avg', 'passer_rating_rolling_3wk_avg'
]

X = df[all_feature_columns]
y_dk = df['dk_fantasy_points']
y_fd = df['fd_fantasy_points']

print(f"Data prepared with {len(all_feature_columns)} features.")

# =============================================================================
# STEP 3: TRAIN THE XGBOOST MODELS
# =============================================================================
print("\n--- Training XGBoost models ---")

# --- DraftKings Model ---
# NOTE: We train on the full dataset now since we removed the validation set for early stopping
xgb_dk = xgb.XGBRegressor(
    objective='reg:squarederror',
    n_estimators=1000, # Using 1000 instead of 2000 as a reasonable number without early stopping
    learning_rate=0.03,
    max_depth=5,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=42,
    n_jobs=-1
)

# --- THE FIX IS HERE: The .fit() call is simplified ---
xgb_dk.fit(X, y_dk, verbose=False)
print("DraftKings XGBoost model trained successfully.")


# --- FanDuel Model ---
xgb_fd = xgb.XGBRegressor(
    objective='reg:squarederror',
    n_estimators=1000,
    learning_rate=0.03,
    max_depth=5,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=42,
    n_jobs=-1
)

# --- THE FIX IS HERE: The .fit() call is simplified ---
xgb_fd.fit(X, y_fd, verbose=False)
print("FanDuel XGBoost model trained successfully.")


# =============================================================================
# STEP 4: MAKE THE PREDICTION
# =============================================================================
print("\n--- Preparing prediction data for Patrick Mahomes, Week 16 2024 ---")
PREDICTION_YEAR = 2024
PREDICTION_WEEK = 16
PLAYER_NAME = 'Patrick Mahomes'

# Get the data from Week 15 to use as input for the Week 16 prediction
prediction_input_df = df[
    (df['player_display_name'] == PLAYER_NAME) &
    (df['season'] == PREDICTION_YEAR) &
    (df['week'] == PREDICTION_WEEK - 1)
]

if prediction_input_df.empty:
    print(f"\nError: Could not find data for {PLAYER_NAME} in Week {PREDICTION_WEEK - 1}, {PREDICTION_YEAR}.")
else:
    X_predict = prediction_input_df[all_feature_columns]

    print("\n--- Generating Predictions using the XGBoost model ---")
    mahomes_dk_pred = xgb_dk.predict(X_predict)[0]
    mahomes_fd_pred = xgb_fd.predict(X_predict)[0]
    mahomes_captain_pred = mahomes_dk_pred * 1.5

    print("\n" + "="*50)
    print(f"FANTASY PREDICTIONS FOR: {PLAYER_NAME}")
    print(f"WEEK: {PREDICTION_WEEK}, SEASON: {PREDICTION_YEAR}")
    print("="*50)
    print(f"Actual Score:                     23.70")
    print("-" * 50)
    print(f"XGBoost DK Points Prediction:     {mahomes_dk_pred:.2f}")
    print(f"XGBoost FD Points Prediction:     {mahomes_fd_pred:.2f}")
    print(f"XGBoost DK Captain Prediction:    {mahomes_captain_pred:.2f}")
    print("="*50)