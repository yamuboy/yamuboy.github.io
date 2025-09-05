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

# Identify all official RBs
official_rb_ids = set(weekly_df[weekly_df['position'] == 'RB']['player_id'].unique())
print(f"Identified {len(official_rb_ids)} unique Running Backs. Raw data loaded.")


# =============================================================================
# STEP 3: CALCULATE RED ZONE AND GOAL-LINE USAGE
# =============================================================================
print("\n--- Calculating Red Zone and Goal-Line Usage from PBP data ---")
rz_df = pbp_df_reg[pbp_df_reg['yardline_100'] <= 20].copy()
gl_df = pbp_df_reg[pbp_df_reg['yardline_100'] <= 5].copy()
def get_usage_stats(df):
    rush_usage = df.groupby(['rusher_player_id', 'season', 'week'], as_index=False).agg(touches=('rush_attempt', 'sum'))
    rush_usage.rename(columns={'rusher_player_id': 'player_id'}, inplace=True)
    rec_usage = df.groupby(['receiver_player_id', 'season', 'week'], as_index=False).agg(touches=('pass_attempt', 'sum'))
    rec_usage.rename(columns={'receiver_player_id': 'player_id'}, inplace=True)
    total_usage = pd.concat([rush_usage, rec_usage]).groupby(['player_id', 'season', 'week'], as_index=False).sum()
    return total_usage
rz_player_usage = get_usage_stats(rz_df).rename(columns={'touches': 'red_zone_touches'})
gl_player_usage = get_usage_stats(gl_df).rename(columns={'touches': 'goal_line_touches'})

# =============================================================================
# STEP 4: CALCULATE ROLLING 4-WEEK DEFENSIVE STATS VS RBs
# =============================================================================
print("\n--- Calculating 4-week rolling defensive stats vs. RBs ---")
weekly_rb_opp_stats = weekly_df[weekly_df['position'] == 'RB']
defense_stats_allowed_weekly = weekly_rb_opp_stats.groupby(['opponent_team', 'season', 'week']).agg(
    rb_fantasy_points_allowed=('fantasy_points', 'sum'), rb_rushing_yards_allowed=('rushing_yards', 'sum'),
    rb_receiving_yards_allowed=('receiving_yards', 'sum'), rb_receptions_allowed=('receptions', 'sum')
).reset_index()
defense_stats_allowed_weekly = defense_stats_allowed_weekly.sort_values(by=['opponent_team', 'season', 'week'])
for stat in ['rb_fantasy_points_allowed', 'rb_rushing_yards_allowed', 'rb_receiving_yards_allowed', 'rb_receptions_allowed']:
    defense_stats_allowed_weekly[f'rolling_4wk_avg_{stat}'] = defense_stats_allowed_weekly.groupby('opponent_team')[stat].transform(
        lambda x: x.rolling(window=4, min_periods=1).mean().shift(1)
    )
rolling_defense_stats = defense_stats_allowed_weekly[['opponent_team', 'season', 'week', 'rolling_4wk_avg_rb_fantasy_points_allowed', 'rolling_4wk_avg_rb_rushing_yards_allowed', 'rolling_4wk_avg_rb_receiving_yards_allowed', 'rolling_4wk_avg_rb_receptions_allowed']]

# =============================================================================
# STEP 5: PREPARE BASE RB DATAFRAME AND MERGE ALL DATA
# =============================================================================
print("\n--- Preparing base RB dataframe and merging all data sources ---")
base_rb_df = weekly_df[weekly_df['player_id'].isin(official_rb_ids)].copy()
base_rb_df['dk_fantasy_points'] = (
    base_rb_df['rushing_yards'] * 0.1 + base_rb_df['rushing_tds'] * 6 +
    base_rb_df['receiving_yards'] * 0.1 + base_rb_df['receiving_tds'] * 6 +
    base_rb_df['receptions'] * 1 +
    (base_rb_df['rushing_fumbles_lost'] + base_rb_df['receiving_fumbles_lost']) * -1
)
base_rb_df['dk_bonus'] = base_rb_df['rushing_yards'].apply(lambda x: 3 if x >= 100 else 0) + base_rb_df['receiving_yards'].apply(lambda x: 3 if x >= 100 else 0)
base_rb_df['dk_fantasy_points'] += base_rb_df['dk_bonus']
base_rb_df['fd_fantasy_points'] = (
    base_rb_df['rushing_yards'] * 0.1 + base_rb_df['rushing_tds'] * 6 +
    base_rb_df['receiving_yards'] * 0.1 + base_rb_df['receiving_tds'] * 6 +
    base_rb_df['receptions'] * 0.5 +
    (base_rb_df['rushing_fumbles_lost'] + base_rb_df['receiving_fumbles_lost']) * -2
)

# --- THE FIX IS HERE ---
# Use the correct column name `offense_pct` for snap percentage.
snap_counts_to_merge = snap_counts_df[['player', 'team', 'season', 'week', 'offense_pct']].copy()
snap_counts_to_merge.rename(columns={'offense_pct': 'snap_percentage'}, inplace=True)

final_df = pd.merge(
    base_rb_df,
    snap_counts_to_merge,
    left_on=['player_display_name', 'recent_team', 'season', 'week'],
    right_on=['player', 'team', 'season', 'week'],
    how='left'
)
# --- END FIX ---


# The rest of the merges will now work correctly
final_df = pd.merge(final_df, rz_player_usage, on=['player_id', 'season', 'week'], how='left')
final_df = pd.merge(final_df, gl_player_usage, on=['player_id', 'season', 'week'], how='left')
final_df = pd.merge(final_df, rolling_defense_stats, on=['opponent_team', 'season', 'week'], how='left')
schedule_context = schedule_df[['season', 'week', 'home_team', 'away_team', 'total_line', 'spread_line', 'div_game']].copy()
schedule_context['home_spread'] = schedule_context['spread_line']
schedule_context['away_spread'] = -schedule_context['spread_line']
final_df = pd.merge(final_df, schedule_context, left_on=['season', 'week', 'recent_team', 'opponent_team'], right_on=['season', 'week', 'home_team', 'away_team'], how='left')

# =============================================================================
# STEP 6: CALCULATE FINAL METRICS AND CLEAN DATAFRAME
# =============================================================================
print("\n--- Calculating final derived metrics and cleaning data ---")
final_df['is_home'] = np.where(final_df['home_team'].notna(), 1, 0)
final_df['team_spread'] = np.where(final_df['is_home'] == 1, final_df['home_spread'], final_df['away_spread'])
final_df['is_division_rival'] = final_df['div_game'].fillna(0).astype(int)
final_df['team_implied_total'] = (final_df['total_line'] / 2) - (final_df['team_spread'] / 2)
final_df['total_touches'] = final_df['carries'] + final_df['targets']
final_df['yards_per_touch'] = (final_df['rushing_yards'] + final_df['receiving_yards']) / final_df['total_touches']
final_df['dk_showdown_captain_points'] = final_df['dk_fantasy_points'] * 1.5
fill_cols = ['snap_percentage', 'red_zone_touches', 'goal_line_touches', 'team_implied_total', 'total_line', 'team_spread']
for col in fill_cols:
    final_df[col] = final_df[col].fillna(0)
final_df.replace([np.inf, -np.inf], 0, inplace=True)

# Drop redundant columns from the merge
final_df.drop(columns=['player', 'team'], inplace=True, errors='ignore')

# =============================================================================
# STEP 7: ORGANIZE AND EXPORT TO CSV
# =============================================================================
print("\n--- Exporting final comprehensive dataset to CSV ---")
final_columns_order = [
    'player_id', 'player_display_name', 'recent_team', 'season', 'week', 'opponent_team',
    'is_home', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    'rolling_4wk_avg_rb_fantasy_points_allowed', 'rolling_4wk_avg_rb_rushing_yards_allowed', 'rolling_4wk_avg_rb_receiving_yards_allowed', 'rolling_4wk_avg_rb_receptions_allowed',
    'snap_percentage', 'carries', 'targets', 'total_touches', 'target_share',
    'red_zone_touches', 'goal_line_touches',
    'rushing_yards', 'rushing_tds', 'rushing_epa',
    'receptions', 'receiving_yards', 'receiving_tds', 'receiving_epa',
    'yards_per_touch',
    'fantasy_points', 'fantasy_points_ppr', 'dk_fantasy_points', 'fd_fantasy_points', 'dk_showdown_captain_points'
]
final_columns_order = [col for col in final_columns_order if col in final_df.columns]
final_rb_stats_df = final_df[final_columns_order].copy()
final_rb_stats_df = final_rb_stats_df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)
output_filename = 'rb_comprehensive_advanced_stats_2020_2024.csv'
final_rb_stats_df.to_csv(output_filename, index=False)

print(f"\n✅ Success! Comprehensive advanced stats for RBs have been exported to '{output_filename}'")
print("\n--- Data Preview ---")
print(final_rb_stats_df.head())