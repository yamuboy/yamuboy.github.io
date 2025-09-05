# =============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import lightgbm as lgb
import optuna
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=optuna.exceptions.ExperimentalWarning)


print("Setup and Imports Complete.")

# =============================================================================
# STEP 2: ADVANCED FEATURE ENGINEERING
# =============================================================================
print("\n--- Loading data and creating advanced features ---")
try:
    df = pd.read_csv('qb_comprehensive_advanced_stats_2020_2024.csv')
except FileNotFoundError:
    print("Error: 'qb_comprehensive_advanced_stats_2020_2024.csv' not found.")
    exit()

df = df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)

# --- NEW: Create 3-week rolling averages for player's recent form ---
rolling_features = [
    'dk_fantasy_points', 'total_yards', 'epa_per_play', 'total_touchdowns',
    'fantasy_points_per_dropback', 'passer_rating'
]
for feature in rolling_features:
    # We use .shift(1) to ensure we're only using data from *before* the current week
    df[f'{feature}_rolling_3wk_avg'] = df.groupby('player_display_name')[feature].transform(
        lambda x: x.rolling(window=3, min_periods=1).mean().shift(1)
    )

df.dropna(inplace=True)

# Define all potential features for the model
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
# STEP 3: HYPERPARAMETER TUNING WITH OPTUNA
# =============================================================================
print("\n--- Finding optimal hyperparameters with Optuna ---")

def objective(trial, X, y):
    """The function Optuna tries to minimize."""
    # Define the search space for hyperparameters
    params = {
        'objective': 'regression_l1',
        'metric': 'mae',
        'n_estimators': trial.suggest_int('n_estimators', 200, 2000),
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3),
        'num_leaves': trial.suggest_int('num_leaves', 20, 300),
        'max_depth': trial.suggest_int('max_depth', 3, 12),
        'min_child_samples': trial.suggest_int('min_child_samples', 5, 100),
        'subsample': trial.suggest_float('subsample', 0.6, 1.0),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
        'random_state': 42,
        'n_jobs': -1,
        'verbose': -1
    }
    
    # Use K-Fold cross-validation to get a reliable error estimate
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    errors = []
    for train_index, val_index in kf.split(X):
        X_train, X_val = X.iloc[train_index], X.iloc[val_index]
        y_train, y_val = y.iloc[train_index], y.iloc[val_index]
        
        model = lgb.LGBMRegressor(**params)
        model.fit(X_train, y_train)
        preds = model.predict(X_val)
        errors.append(mean_absolute_error(y_val, preds))
        
    return np.mean(errors)

# Create a study object and optimize
# We do this separately for DK and FD as their optimal parameters might differ
print("Tuning DraftKings model...")
study_dk = optuna.create_study(direction='minimize')
study_dk.optimize(lambda trial: objective(trial, X, y_dk), n_trials=50) # 50 trials is a good balance of search and speed
best_params_dk = study_dk.best_params

print("Tuning FanDuel model...")
study_fd = optuna.create_study(direction='minimize')
study_fd.optimize(lambda trial: objective(trial, X, y_fd), n_trials=50)
best_params_fd = study_fd.best_params

print("\nBest hyperparameters found for DraftKings model:", best_params_dk)
print("Best hyperparameters found for FanDuel model:", best_params_fd)

# =============================================================================
# STEP 4: TRAIN FINAL MODELS WITH OPTIMIZED PARAMETERS
# =============================================================================
print("\n--- Training final models with optimized parameters ---")

# Train final DK model on ALL data using the best params
final_model_dk = lgb.LGBMRegressor(objective='regression_l1', random_state=42, n_jobs=-1, **best_params_dk)
final_model_dk.fit(X, y_dk)
print("Final DraftKings model trained.")

# Train final FD model on ALL data using the best params
final_model_fd = lgb.LGBMRegressor(objective='regression_l1', random_state=42, n_jobs=-1, **best_params_fd)
final_model_fd.fit(X, y_fd)
print("Final FanDuel model trained.")

# =============================================================================
# STEP 5: MAKE THE PREDICTION
# =============================================================================
print("\n--- Preparing prediction data for Patrick Mahomes, Week 16 2024 ---")
PREDICTION_YEAR = 2024
PREDICTION_WEEK = 16
PLAYER_NAME = 'Patrick Mahomes'

# We still need the data from Week 15 to get the rolling averages for our Week 16 prediction
prediction_input_df = df[
    (df['player_display_name'] == PLAYER_NAME) &
    (df['season'] == PREDICTION_YEAR) &
    (df['week'] == PREDICTION_WEEK - 1)
]

if prediction_input_df.empty:
    print(f"\nError: Could not find data for {PLAYER_NAME} in Week {PREDICTION_WEEK - 1}, {PREDICTION_YEAR}.")
else:
    # Select only the feature columns needed for the model
    X_predict = prediction_input_df[all_feature_columns]

    print("\n--- Generating Predictions using the OPTIMIZED model ---")
    mahomes_dk_pred = final_model_dk.predict(X_predict)[0]
    mahomes_fd_pred = final_model_fd.predict(X_predict)[0]
    mahomes_captain_pred = mahomes_dk_pred * 1.5

    print("\n" + "="*50)
    print(f"FANTASY PREDICTIONS FOR: {PLAYER_NAME}")
    print(f"WEEK: {PREDICTION_WEEK}, SEASON: {PREDICTION_YEAR}")
    print("="*50)
    print(f"Previous Model Prediction:        19.00")
    print(f"Actual Score:                     23.70")
    print("-" * 50)
    print(f"Optimized DK Points Prediction:   {mahomes_dk_pred:.2f}")
    print(f"Optimized FD Points Prediction:   {mahomes_fd_pred:.2f}")
    print(f"Optimized DK Captain Prediction:  {mahomes_captain_pred:.2f}")
    print("="*50)