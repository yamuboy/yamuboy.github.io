# =============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import lightgbm as lgb
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress common warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

print("Setup and Imports Complete.")

# =============================================================================
# STEP 2: LOAD AND PREPARE THE DATA
# =============================================================================
print("\n--- Loading the comprehensive QB dataset ---")
try:
    df = pd.read_csv('qb_comprehensive_advanced_stats_2020_2024.csv')
except FileNotFoundError:
    print("Error: 'qb_comprehensive_advanced_stats_2020_2024.csv' not found.")
    exit()

# --- Feature Engineering (Lag Features) ---
df = df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)
# Add a few more lag features for the model to consider
lag_features_to_create = [
    'dk_fantasy_points', 'fd_fantasy_points', 'total_yards', 'epa_per_play',
    'total_touchdowns', 'fantasy_points_per_dropback', 'passer_rating'
]
for feature in lag_features_to_create:
    df[f'{feature}_lag_1'] = df.groupby('player_display_name')[feature].shift(1)

df.dropna(inplace=True)

# --- Define ALL potential features ---
# We give the model everything and let it decide what's important
all_feature_columns = [
    # Betting & Context
    'is_home', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    # Opponent Rolling Stats
    'rolling_4wk_avg_qb_fantasy_points_allowed', 'rolling_4wk_avg_qb_sacks_forced',
    'rolling_4wk_avg_qb_interceptions_forced',
    # Lagged Performance Stats
    'dk_fantasy_points_lag_1', 'fd_fantasy_points_lag_1', 'total_yards_lag_1',
    'epa_per_play_lag_1', 'total_touchdowns_lag_1', 'fantasy_points_per_dropback_lag_1',
    'passer_rating_lag_1'
]

X_all = df[all_feature_columns]
y_dk = df['dk_fantasy_points']

print(f"Data prepared. Using {len(all_feature_columns)} potential features.")

# =============================================================================
# STEP 3: TRAIN INITIAL MODEL AND FIND MOST IMPORTANT FEATURES
# =============================================================================
print("\n--- Training initial model to find feature importances ---")
lgbm_initial = lgb.LGBMRegressor(random_state=42, n_jobs=-1, verbose=-1)
lgbm_initial.fit(X_all, y_dk)

# Create a DataFrame of feature importances
feature_importance_df = pd.DataFrame({
    'feature': all_feature_columns,
    'importance': lgbm_initial.feature_importances_
}).sort_values(by='importance', ascending=False)

print("\n--- Most Predictive Features for DraftKings Points ---")
print(feature_importance_df)

# Visualize the feature importances
plt.figure(figsize=(12, 8))
sns.barplot(x='importance', y='feature', data=feature_importance_df)
plt.title('Feature Importance for Predicting QB DraftKings Points', fontsize=16)
plt.xlabel('Importance Score', fontsize=12)
plt.ylabel('Feature', fontsize=12)
plt.tight_layout()
plt.show()

# =============================================================================
# STEP 4: RE-TRAIN A FOCUSED MODEL WITH ONLY THE TOP FEATURES
# =============================================================================
# Let's select the top 8 features based on the results
N_TOP_FEATURES = 8
top_features = feature_importance_df.head(N_TOP_FEATURES)['feature'].tolist()

print(f"\n--- Selecting the top {N_TOP_FEATURES} features for the final model ---")
print(top_features)

# Create new training data with only the most important features
X_top = df[top_features]

# Train the final, focused models
print("\n--- Training final models on selected features ---")
lgbm_dk_final = lgb.LGBMRegressor(objective='regression_l1', n_estimators=1000, learning_rate=0.05, random_state=42, n_jobs=-1, verbose=-1)
lgbm_dk_final.fit(X_top, df['dk_fantasy_points'])
print("Final DraftKings model trained.")

lgbm_fd_final = lgb.LGBMRegressor(objective='regression_l1', n_estimators=1000, learning_rate=0.05, random_state=42, n_jobs=-1, verbose=-1)
lgbm_fd_final.fit(X_top, df['fd_fantasy_points'])
print("Final FanDuel model trained.")

# =============================================================================
# STEP 5: MAKE THE PREDICTION USING THE FOCUSED MODEL
# =============================================================================
print("\n--- Preparing prediction data for Patrick Mahomes, Week 16 2024 ---")
PREDICTION_YEAR = 2024
PREDICTION_WEEK = 16
PLAYER_NAME = 'Patrick Mahomes'

prediction_input_df = df[
    (df['player_display_name'] == PLAYER_NAME) &
    (df['season'] == PREDICTION_YEAR) &
    (df['week'] == PREDICTION_WEEK - 1)
]

if prediction_input_df.empty:
    print(f"\nError: Could not find data for {PLAYER_NAME} in Week {PREDICTION_WEEK - 1}, {PREDICTION_YEAR}.")
else:
    # Create the prediction input using ALL features first to ensure lag columns are present
    X_predict_full = prediction_input_df[all_feature_columns].copy()

    # Manually update all lag features with Week 15's actual performance
    for feature in lag_features_to_create:
        X_predict_full[f'{feature}_lag_1'] = prediction_input_df[feature].values[0]

    # Now, select only the top features for the final prediction
    X_predict_top = X_predict_full[top_features]

    print("\n--- Generating Predictions using the focused model ---")
    mahomes_dk_pred = lgbm_dk_final.predict(X_predict_top)[0]
    mahomes_fd_pred = lgbm_fd_final.predict(X_predict_top)[0]
    mahomes_captain_pred = mahomes_dk_pred * 1.5

    print("\n" + "="*50)
    print(f"FANTASY PREDICTIONS FOR: {PLAYER_NAME}")
    print(f"WEEK: {PREDICTION_WEEK}, SEASON: {PREDICTION_YEAR}")
    print("="*50)
    print(f"Predicted DraftKings Points:      {mahomes_dk_pred:.2f}")
    print(f"Predicted FanDuel Points:         {mahomes_fd_pred:.2f}")
    print(f"Predicted DK Showdown Captain:    {mahomes_captain_pred:.2f}")
    print("="*50)