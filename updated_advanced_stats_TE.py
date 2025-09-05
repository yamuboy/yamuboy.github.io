# =============================================================================
# STEP 1: SETUP AND IMPORTS
# =============================================================================
import pandas as pd
import numpy as np
import nfl_data_py as nfl
import warnings

# Suppress common warnings for cleaner output
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

print("Setup and Imports Complete.")

# =============================================================================
# STEP 2: LOAD RAW DATA (2020-2024)
# =============================================================================
print("\n--- Loading Play-by-Play, Weekly, Schedule, and Snap Count Data (2020-2024) ---")
years = range(2020, 2025)
pbp_df = nfl.import_pbp_data(years)
weekly_df = nfl.import_weekly_data(years)
schedule_df = nfl.import_schedules(years)
snap_counts_df = nfl.import_snap_counts(years)

# Filter for regular season plays only
pbp_df_reg = pbp_df[pbp_df['season_type'] == 'REG'].copy()

# Identify all official TEs
official_te_ids = set(weekly_df[weekly_df['position'] == 'TE']['player_id'].unique())
print(f"Identified {len(official_te_ids)} unique Tight Ends. Raw data loaded.")


# =============================================================================
# STEP 3: CALCULATE ROLLING 4-WEEK DEFENSIVE STATS VS TEs
# =============================================================================
print("\n--- Calculating 4-week rolling defensive stats vs. TEs ---")
# This is a key step: filter for TE position BEFORE grouping by opponent
weekly_te_opp_stats = weekly_df[weekly_df['position'] == 'TE']
defense_stats_allowed_weekly = weekly_te_opp_stats.groupby(['opponent_team', 'season', 'week']).agg(
    te_fantasy_points_allowed=('fantasy_points', 'sum'),
    te_receiving_yards_allowed=('receiving_yards', 'sum'),
    te_receptions_allowed=('receptions', 'sum')
).reset_index()

defense_stats_allowed_weekly = defense_stats_allowed_weekly.sort_values(by=['opponent_team', 'season', 'week'])
for stat in ['te_fantasy_points_allowed', 'te_receiving_yards_allowed', 'te_receptions_allowed']:
    defense_stats_allowed_weekly[f'rolling_4wk_avg_{stat}'] = defense_stats_allowed_weekly.groupby('opponent_team')[stat].transform(
        lambda x: x.rolling(window=4, min_periods=1).mean().shift(1)
    )
rolling_defense_stats = defense_stats_allowed_weekly[['opponent_team', 'season', 'week', 'rolling_4wk_avg_te_fantasy_points_allowed', 'rolling_4wk_avg_te_receiving_yards_allowed', 'rolling_4wk_avg_te_receptions_allowed']]


# =============================================================================
# STEP 4: PREPARE BASE TE DATAFRAME AND MERGE ALL DATA
# =============================================================================
print("\n--- Preparing base TE dataframe and merging all data sources ---")

# Start with the weekly data, filtered for our official TEs
base_te_df = weekly_df[weekly_df['player_id'].isin(official_te_ids)].copy()

# --- Manually Calculate DraftKings and FanDuel Points ---
# Scoring for TEs is the same as for WRs
base_te_df['dk_fantasy_points'] = (
    base_te_df['rushing_yards'] * 0.1 + base_te_df['rushing_tds'] * 6 +
    base_te_df['receiving_yards'] * 0.1 + base_te_df['receiving_tds'] * 6 +
    base_te_df['receptions'] * 1 +
    (base_te_df['rushing_fumbles_lost'] + base_te_df['receiving_fumbles_lost']) * -1
)
base_te_df['dk_bonus'] = base_te_df['receiving_yards'].apply(lambda x: 3 if x >= 100 else 0)
base_te_df['dk_fantasy_points'] += base_te_df['dk_bonus']
base_te_df['fd_fantasy_points'] = (
    base_te_df['rushing_yards'] * 0.1 + base_te_df['rushing_tds'] * 6 +
    base_te_df['receiving_yards'] * 0.1 + base_te_df['receiving_tds'] * 6 +
    base_te_df['receptions'] * 0.5 +
    (base_te_df['rushing_fumbles_lost'] + base_te_df['receiving_fumbles_lost']) * -2
)

# Merge snap counts using the robust composite key method
snap_counts_to_merge = snap_counts_df[['player', 'team', 'season', 'week', 'offense_pct']].copy()
snap_counts_to_merge.rename(columns={'offense_pct': 'snap_percentage'}, inplace=True)
final_df = pd.merge(base_te_df, snap_counts_to_merge, left_on=['player_display_name', 'recent_team', 'season', 'week'], right_on=['player', 'team', 'season', 'week'], how='left')

# Merge rolling defensive stats
final_df = pd.merge(final_df, rolling_defense_stats, on=['opponent_team', 'season', 'week'], how='left')

# Merge schedule context
schedule_context = schedule_df[['season', 'week', 'home_team', 'away_team', 'total_line', 'spread_line', 'div_game']].copy()
schedule_context['home_spread'] = schedule_context['spread_line']
schedule_context['away_spread'] = -schedule_context['spread_line']
final_df = pd.merge(final_df, schedule_context, left_on=['season', 'week', 'recent_team', 'opponent_team'], right_on=['season', 'week', 'home_team', 'away_team'], how='left')


# =============================================================================
# STEP 5: CALCULATE FINAL METRICS AND CLEAN DATAFRAME
# =============================================================================
print("\n--- Calculating final derived metrics and cleaning data ---")

# Add context columns
final_df['is_home'] = np.where(final_df['home_team'].notna(), 1, 0)
final_df['team_spread'] = np.where(final_df['is_home'] == 1, final_df['home_spread'], final_df['away_spread'])
final_df['is_division_rival'] = final_df['div_game'].fillna(0).astype(int)
final_df['team_implied_total'] = (final_df['total_line'] / 2) - (final_df['team_spread'] / 2)

# --- NEW METRIC CALCULATIONS from PBP data ---
# Calculate Routes Run
routes_run = pbp_df_reg.groupby(['receiver_player_id', 'season', 'week'], as_index=False).agg(routes_run=('play_id', 'count'))
routes_run.rename(columns={'receiver_player_id': 'player_id'}, inplace=True)
final_df = pd.merge(final_df, routes_run, on=['player_id', 'season', 'week'], how='left')

# Calculate Red Zone Targets
rz_df = pbp_df_reg[pbp_df_reg['yardline_100'] <= 20]
rz_targets = rz_df.groupby(['receiver_player_id', 'season', 'week'], as_index=False).agg(red_zone_targets=('play_id', 'count'))
rz_targets.rename(columns={'receiver_player_id': 'player_id'}, inplace=True)
final_df = pd.merge(final_df, rz_targets, on=['player_id', 'season', 'week'], how='left')

# Calculate final derived metrics
final_df['yards_per_route_run'] = (final_df['receiving_yards'] / final_df['routes_run'])
final_df['average_depth_of_target'] = (final_df['receiving_air_yards'] / final_df['targets'])
final_df['yac_per_reception'] = (final_df['receiving_yards_after_catch'] / final_df['receptions'])
final_df['dk_showdown_captain_points'] = final_df['dk_fantasy_points'] * 1.5

# Fill NaNs for key columns before final output
fill_cols = ['snap_percentage', 'team_implied_total', 'total_line', 'team_spread', 'routes_run', 'red_zone_targets']
for col in fill_cols:
    final_df[col] = final_df[col].fillna(0)
final_df.replace([np.inf, -np.inf], 0, inplace=True)

# Drop redundant columns
final_df.drop(columns=['player', 'team'], inplace=True, errors='ignore')


# =============================================================================
# STEP 6: ORGANIZE AND EXPORT TO CSV
# =============================================================================
print("\n--- Exporting final comprehensive dataset to CSV ---")

final_columns_order = [
    # Player & Game Info
    'player_id', 'player_display_name', 'recent_team', 'season', 'week', 'opponent_team',
    # Betting & Context
    'is_home', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    # Opponent Rolling Stats vs TEs
    'rolling_4wk_avg_te_fantasy_points_allowed', 'rolling_4wk_avg_te_receiving_yards_allowed', 'rolling_4wk_avg_te_receptions_allowed',
    # Core Volume & Usage
    'snap_percentage', 'routes_run', 'targets', 'receptions', 'target_share', 'air_yards_share',
    # Red Zone Usage
    'red_zone_targets',
    # Production & Efficiency
    'receiving_yards', 'receiving_tds', 'receiving_yards_after_catch', 'yac_per_reception',
    'receiving_air_yards', 'average_depth_of_target', 'yards_per_route_run',
    # Fantasy Production
    'fantasy_points', 'fantasy_points_ppr', 'dk_fantasy_points', 'fd_fantasy_points', 'dk_showdown_captain_points'
]

final_columns_order = [col for col in final_columns_order if col in final_df.columns]
final_te_stats_df = final_df[final_columns_order].copy()
final_te_stats_df = final_te_stats_df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)

output_filename = 'te_comprehensive_advanced_stats_2020_2024.csv'
final_te_stats_df.to_csv(output_filename, index=False)

print(f"\n✅ Success! Comprehensive advanced stats for TEs have been exported to '{output_filename}'")
print("\n--- Data Preview ---")
print(final_te_stats_df.head())