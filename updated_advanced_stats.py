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
print("\n--- Loading Play-by-Play, Weekly, and Schedule Data (2020-2024) ---")
years = range(2020, 2025)
pbp_df = nfl.import_pbp_data(years)
weekly_df = nfl.import_weekly_data(years)
schedule_df = nfl.import_schedules(years)

# Filter for regular season plays only
pbp_df_reg = pbp_df[pbp_df['season_type'] == 'REG'].copy()
official_qb_ids = set(weekly_df[weekly_df['position'] == 'QB']['player_id'].unique())
print(f"Identified {len(official_qb_ids)} unique QBs. Raw data loaded.")


# =============================================================================
# STEP 3: MANUALLY CALCULATE DRAFTKINGS AND FANDUEL FANTASY POINTS
# =============================================================================
print("\n--- Manually calculating DraftKings and FanDuel fantasy points ---")

# Use a copy of the weekly_df to avoid SettingWithCopyWarning
weekly_df_enriched = weekly_df.copy()

# --- DraftKings Scoring ---
weekly_df_enriched['dk_passing_pts'] = (
    weekly_df_enriched['passing_yards'] * 0.04 +
    weekly_df_enriched['passing_tds'] * 4 -
    weekly_df_enriched['interceptions'] * 1
)
weekly_df_enriched['dk_passing_bonus'] = weekly_df_enriched['passing_yards'].apply(lambda x: 3 if x >= 300 else 0)

weekly_df_enriched['dk_rushing_pts'] = (
    weekly_df_enriched['rushing_yards'] * 0.1 +
    weekly_df_enriched['rushing_tds'] * 6
)
weekly_df_enriched['dk_rushing_bonus'] = weekly_df_enriched['rushing_yards'].apply(lambda x: 3 if x >= 100 else 0)

# Combine fumbles from sacks and rushes for a total fumble lost count
weekly_df_enriched['fumbles_lost'] = weekly_df_enriched['sack_fumbles_lost'] + weekly_df_enriched['rushing_fumbles_lost']
weekly_df_enriched['dk_fumble_pts'] = weekly_df_enriched['fumbles_lost'] * -1

weekly_df_enriched['dk_fantasy_points'] = (
    weekly_df_enriched['dk_passing_pts'] + weekly_df_enriched['dk_passing_bonus'] +
    weekly_df_enriched['dk_rushing_pts'] + weekly_df_enriched['dk_rushing_bonus'] +
    weekly_df_enriched['dk_fumble_pts']
)

# --- FanDuel Scoring ---
weekly_df_enriched['fd_passing_pts'] = (
    weekly_df_enriched['passing_yards'] * 0.04 +
    weekly_df_enriched['passing_tds'] * 4 -
    weekly_df_enriched['interceptions'] * 1
)
weekly_df_enriched['fd_rushing_pts'] = (
    weekly_df_enriched['rushing_yards'] * 0.1 +
    weekly_df_enriched['rushing_tds'] * 6
)
weekly_df_enriched['fd_fumble_pts'] = weekly_df_enriched['fumbles_lost'] * -1

weekly_df_enriched['fd_fantasy_points'] = (
    weekly_df_enriched['fd_passing_pts'] +
    weekly_df_enriched['fd_rushing_pts'] +
    weekly_df_enriched['fd_fumble_pts']
)


# =============================================================================
# STEP 4: CALCULATE QB'S OWN WEEKLY ADVANCED STATS
# =============================================================================
print("\n--- Aggregating weekly stats for each QB ---")
pbp_df_reg['a_yards'] = pbp_df_reg['air_yards'].fillna(0)
passer_stats = pbp_df_reg.groupby(['passer_player_id', 'passer_player_name', 'season', 'week'], as_index=False).agg(
    passing_attempts=('pass_attempt', 'sum'), completions=('complete_pass', 'sum'),
    passing_yards_pbp=('passing_yards', 'sum'), passing_tds_pbp=('pass_touchdown', 'sum'),
    interceptions=('interception', 'sum'), sacks=('sack', 'sum'),
    scrambles=('qb_scramble', 'sum'),
    passing_epa=('epa', lambda x: x[pbp_df_reg.loc[x.index, 'pass_attempt'] == 1].sum()),
    total_air_yards=('a_yards', lambda x: x[pbp_df_reg.loc[x.index, 'pass_attempt'] == 1].sum()),
    completed_air_yards=('a_yards', lambda x: x[pbp_df_reg.loc[x.index, 'complete_pass'] == 1].sum())
)
rusher_stats = pbp_df_reg.groupby(['rusher_player_id', 'rusher_player_name', 'season', 'week'], as_index=False).agg(
    rushing_attempts=('rush_attempt', 'sum'), rushing_yards_pbp=('rushing_yards', 'sum'),
    rushing_tds_pbp=('rush_touchdown', 'sum'),
    rushing_epa=('epa', lambda x: x[pbp_df_reg.loc[x.index, 'rush_attempt'] == 1].sum())
)
passer_stats.rename(columns={'passer_player_id': 'player_id', 'passer_player_name': 'player_name'}, inplace=True)
rusher_stats.rename(columns={'rusher_player_id': 'player_id', 'rusher_player_name': 'player_name'}, inplace=True)
weekly_stats_unfiltered = pd.merge(passer_stats, rusher_stats, on=['player_id', 'player_name', 'season', 'week'], how='outer')
qb_weekly_stats = weekly_stats_unfiltered[weekly_stats_unfiltered['player_id'].isin(official_qb_ids)].copy()

# =============================================================================
# STEP 5: CALCULATE ROLLING 4-WEEK DEFENSIVE STATS
# =============================================================================
print("\n--- Calculating 4-week rolling defensive stats for all teams ---")
# Use the enriched weekly data for defense calculations too
weekly_qb_opp_stats = weekly_df_enriched[weekly_df_enriched['position'] == 'QB']
defense_stats_allowed_weekly = weekly_qb_opp_stats.groupby(['opponent_team', 'season', 'week']).agg(
    qb_fantasy_points_allowed=('fantasy_points', 'sum'), qb_sacks_forced=('sacks', 'sum'),
    qb_interceptions_forced=('interceptions', 'sum')
).reset_index()
defense_stats_allowed_weekly = defense_stats_allowed_weekly.sort_values(by=['opponent_team', 'season', 'week'])
for stat in ['qb_fantasy_points_allowed', 'qb_sacks_forced', 'qb_interceptions_forced']:
    defense_stats_allowed_weekly[f'rolling_4wk_avg_{stat}'] = defense_stats_allowed_weekly.groupby('opponent_team')[stat].transform(
        lambda x: x.rolling(window=4, min_periods=1).mean().shift(1)
    )
rolling_defense_stats = defense_stats_allowed_weekly[['opponent_team', 'season', 'week', 'rolling_4wk_avg_qb_fantasy_points_allowed', 'rolling_4wk_avg_qb_sacks_forced', 'rolling_4wk_avg_qb_interceptions_forced']]

# =============================================================================
# STEP 6: MERGE ALL DATASETS AND ADD GAME CONTEXT
# =============================================================================
print("\n--- Merging all datasets and adding game context ---")

# Now we use our weekly_df_enriched which contains the calculated DK/FD points
weekly_info_to_merge = weekly_df_enriched[[
    'player_id', 'season', 'week', 'player_display_name', 'recent_team', 'opponent_team',
    'fantasy_points', 'fantasy_points_ppr', 'dk_fantasy_points', 'fd_fantasy_points',
    'passing_yards', 'passing_tds', 'rushing_yards', 'rushing_tds' # Official weekly stats
]]
final_df = pd.merge(qb_weekly_stats, weekly_info_to_merge, on=['player_id', 'season', 'week'], how='left')
final_df = pd.merge(final_df, rolling_defense_stats, on=['opponent_team', 'season', 'week'], how='left')

schedule_context = schedule_df[['season', 'week', 'home_team', 'away_team', 'roof', 'total_line', 'spread_line', 'div_game']].copy()
schedule_context['home_spread'] = schedule_context['spread_line']
schedule_context['away_spread'] = -schedule_context['spread_line']
final_df = pd.merge(final_df, schedule_context, left_on=['season', 'week', 'recent_team', 'opponent_team'], right_on=['season', 'week', 'home_team', 'away_team'], how='left')

final_df['is_home'] = np.where(final_df['home_team'].notna(), 1, 0)
final_df['stadium_type'] = final_df['roof'].apply(lambda x: 'indoor' if x in ['dome', 'closed', 'retractable'] else 'outdoor')
final_df['team_spread'] = np.where(final_df['is_home'] == 1, final_df['home_spread'], final_df['away_spread'])
final_df['is_division_rival'] = final_df['div_game'].fillna(0).astype(int)
final_df['team_implied_total'] = (final_df['total_line'] / 2) - (final_df['team_spread'] / 2)

final_df.drop(columns=['home_team', 'away_team', 'roof', 'home_spread', 'away_spread', 'div_game', 'spread_line'], inplace=True)

# =============================================================================
# STEP 7: CALCULATE NEW DERIVED METRICS AND FINALIZE DATAFRAME
# =============================================================================
print("\n--- Calculating new derived metrics and finalizing data ---")

stat_cols = [
    'passing_attempts', 'completions', 'sacks', 'scrambles', 'passing_epa', 'total_air_yards',
    'completed_air_yards', 'rushing_attempts', 'rushing_epa', 'fantasy_points', 'fantasy_points_ppr',
    'dk_fantasy_points', 'fd_fantasy_points'
]
for col in stat_cols:
    final_df[col] = final_df[col].fillna(0)

# Volume & Efficiency
final_df['total_touchdowns'] = final_df['passing_tds'] + final_df['rushing_tds']
final_df['total_yards'] = final_df['passing_yards'] + final_df['rushing_yards']
final_df['passing_td_percentage'] = (final_df['passing_tds'] / final_df['passing_attempts'])
final_df['completed_air_yards_per_completion'] = (final_df['completed_air_yards'] / final_df['completions'])
final_df['air_yards_per_attempt'] = (final_df['total_air_yards'] / final_df['passing_attempts'])
# Rushing Ability
final_df['fantasy_points_from_rushing'] = (final_df['rushing_yards'] * 0.1) + (final_df['rushing_tds'] * 6)
# Advanced Metrics
final_df['dropbacks'] = final_df['passing_attempts'] + final_df['sacks']
final_df['fantasy_points_per_dropback'] = (final_df['fantasy_points'] / final_df['dropbacks'])
total_plays = final_df['passing_attempts'] + final_df['rushing_attempts']
total_epa = final_df['passing_epa'] + final_df['rushing_epa']
final_df['epa_per_play'] = (total_epa / total_plays)
def calculate_passer_rating(row):
    att, comp, yds, td, int_ = row['passing_attempts'], row['completions'], row['passing_yards'], row['passing_tds'], row['interceptions']
    if att == 0: return 0.0
    a = min(2.375, max(0, ((comp / att) - 0.3) * 5))
    b = min(2.375, max(0, ((yds / att) - 3) * 0.25))
    c = min(2.375, max(0, (td / att) * 20))
    d = min(2.375, max(0, 2.375 - ((int_ / att) * 25)))
    return ((a + b + c + d) / 6) * 100
final_df['passer_rating'] = final_df.apply(calculate_passer_rating, axis=1)

# Showdown Captain Points Calculation
final_df['dk_showdown_captain_points'] = final_df['dk_fantasy_points'] * 1.5

# Clean up any potential infinite values or NaNs from division by zero
final_df.replace([np.inf, -np.inf], 0, inplace=True)
final_df = final_df.fillna(0)

# =============================================================================
# STEP 8: ORGANIZE AND EXPORT TO CSV
# =============================================================================
print("\n--- Exporting final comprehensive dataset to CSV ---")

final_columns_order = [
    'player_id', 'player_display_name', 'recent_team', 'season', 'week', 'opponent_team',
    'is_home', 'stadium_type', 'is_division_rival', 'total_line', 'team_spread', 'team_implied_total',
    'rolling_4wk_avg_qb_fantasy_points_allowed', 'rolling_4wk_avg_qb_sacks_forced', 'rolling_4wk_avg_qb_interceptions_forced',
    'passing_attempts', 'completions', 'passing_yards', 'interceptions', 'rushing_attempts', 'rushing_yards', 'total_yards',
    'passing_tds', 'rushing_tds', 'total_touchdowns', 'passing_td_percentage',
    'scrambles', 'fantasy_points_from_rushing',
    'total_air_yards', 'completed_air_yards', 'air_yards_per_attempt', 'completed_air_yards_per_completion',
    'passer_rating', 'dropbacks', 'passing_epa', 'rushing_epa', 'epa_per_play', 'fantasy_points_per_dropback',
    'fantasy_points', 'fantasy_points_ppr', 'dk_fantasy_points', 'fd_fantasy_points', 'dk_showdown_captain_points'
]
final_columns_order = [col for col in final_columns_order if col in final_df.columns]
final_qb_stats_df = final_df[final_columns_order].copy()
final_qb_stats_df = final_qb_stats_df.sort_values(by=['player_display_name', 'season', 'week']).reset_index(drop=True)
output_filename = 'qb_comprehensive_advanced_stats_2020_2024.csv'
final_qb_stats_df.to_csv(output_filename, index=False)

print(f"\n✅ Success! Comprehensive advanced stats for QBs have been exported to '{output_filename}'")
print("\n--- Data Preview ---")
print(final_qb_stats_df.head())