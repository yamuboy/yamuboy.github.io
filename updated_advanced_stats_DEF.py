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
schedule_df = nfl.import_schedules(years)

# Filter for regular season plays only
pbp_df_reg = pbp_df[pbp_df['season_type'] == 'REG'].copy()
print("Raw data loaded successfully.")

# =============================================================================
# STEP 3: CALCULATE BASE WEEKLY DEFENSIVE STATS FROM PBP
# =============================================================================
print("\n--- Aggregating weekly defensive performance stats from PBP data ---")

# Aggregate core defensive stats for each defense, each week
defensive_stats_weekly = pbp_df_reg.groupby(['defteam', 'season', 'week'], as_index=False).agg(
    sacks=('sack', 'sum'),
    interceptions=('interception', 'sum'),
    fumbles_forced=('fumble_forced', 'sum'),
    fumbles_recovered=('fumble_lost', 'sum') # A fumble lost by the offense is recovered by the defense
)
# Calculate defensive touchdowns
defensive_tds = pbp_df_reg[pbp_df_reg['td_team'] == pbp_df_reg['defteam']]
defensive_td_counts = defensive_tds.groupby(['defteam', 'season', 'week'], as_index=False).agg(
    defensive_tds=('play_id', 'count')
)
# Merge TDs back into the main defensive stats
defensive_stats_weekly = pd.merge(defensive_stats_weekly, defensive_td_counts, on=['defteam', 'season', 'week'], how='left')
defensive_stats_weekly['defensive_tds'] = defensive_stats_weekly['defensive_tds'].fillna(0)

# =============================================================================
# STEP 4: CREATE THE MAIN DST DATAFRAME FROM SCHEDULE DATA
# =============================================================================
print("\n--- Creating the main DST-week dataframe from schedule data ---")

# We need a row for each defense in each game, so we reshape the schedule
home_teams = schedule_df[['season', 'week', 'home_team', 'away_team', 'home_score', 'away_score', 'total_line', 'spread_line']].copy()
home_teams.rename(columns={'home_team': 'team', 'away_team': 'opponent', 'home_score': 'points_for', 'away_score': 'points_allowed'}, inplace=True)
home_teams['is_home'] = 1

away_teams = schedule_df[['season', 'week', 'away_team', 'home_team', 'away_score', 'home_score', 'total_line', 'spread_line']].copy()
away_teams.rename(columns={'away_team': 'team', 'home_team': 'opponent', 'away_score': 'points_for', 'home_score': 'points_allowed'}, inplace=True)
away_teams['is_home'] = 0
# For away teams, the spread must be inverted
away_teams['spread_line'] = -away_teams['spread_line']

# Combine into one DataFrame where each row is a team's weekly game
dst_weekly_df = pd.concat([home_teams, away_teams])
dst_weekly_df = dst_weekly_df.sort_values(by=['team', 'season', 'week']).reset_index(drop=True)

# =============================================================================
# STEP 5: CALCULATE ROLLING OFFENSIVE STATS FOR MATCHUP ANALYSIS
# =============================================================================
print("\n--- Calculating rolling offensive stats for matchup analysis ---")
# These are the stats of the offense a DST will face
offensive_stats_weekly = pbp_df_reg.groupby(['posteam', 'season', 'week'], as_index=False).agg(
    sacks_allowed=('sack', 'sum'),
    interceptions_thrown=('interception', 'sum'),
    fumbles_lost=('fumble_lost', 'sum')
)
offensive_stats_weekly['turnovers_committed'] = offensive_stats_weekly['interceptions_thrown'] + offensive_stats_weekly['fumbles_lost']
offensive_stats_weekly = offensive_stats_weekly.sort_values(by=['posteam', 'season', 'week'])

for stat in ['sacks_allowed', 'turnovers_committed']:
    offensive_stats_weekly[f'rolling_4wk_avg_{stat}'] = offensive_stats_weekly.groupby('posteam')[stat].transform(
        lambda x: x.rolling(window=4, min_periods=1).mean().shift(1)
    )
rolling_offensive_stats = offensive_stats_weekly[['posteam', 'season', 'week', 'rolling_4wk_avg_sacks_allowed', 'rolling_4wk_avg_turnovers_committed']]

# =============================================================================
# STEP 6: MERGE ALL DATASETS AND CALCULATE FANTASY POINTS
# =============================================================================
print("\n--- Merging all datasets and calculating fantasy points ---")

# Merge the aggregated defensive stats
final_df = pd.merge(dst_weekly_df, defensive_stats_weekly, left_on=['team', 'season', 'week'], right_on=['defteam', 'season', 'week'], how='left')

# Merge the rolling offensive stats of the opponent
final_df = pd.merge(final_df, rolling_offensive_stats, left_on=['opponent', 'season', 'week'], right_on=['posteam', 'season', 'week'], how='left')

# --- Manually Calculate DraftKings DST Fantasy Points ---
def get_dk_dst_points(row):
    score = 0
    pa = row['points_allowed']
    if pa == 0: score += 10
    elif 1 <= pa <= 6: score += 7
    elif 7 <= pa <= 13: score += 4
    elif 14 <= pa <= 20: score += 1
    elif 28 <= pa <= 34: score -= 1
    elif pa >= 35: score -= 4
    
    score += row['sacks'] * 1
    score += row['interceptions'] * 2
    score += row['fumbles_recovered'] * 2
    score += row['defensive_tds'] * 6
    # Note: Safeties and Blocked Kicks are rare and harder to parse reliably, omitted for consistency
    return score

# Fill NaNs before calculation
fill_cols = ['sacks', 'interceptions', 'fumbles_recovered', 'defensive_tds', 'points_allowed']
for col in fill_cols:
    final_df[col] = final_df[col].fillna(0)
final_df['dk_fantasy_points'] = final_df.apply(get_dk_dst_points, axis=1)

# =============================================================================
# STEP 7: FINALIZE DATAFRAME AND EXPORT
# =============================================================================
print("\n--- Finalizing dataframe and exporting to CSV ---")

# Final context calculations
final_df['team_implied_total'] = (final_df['total_line'] / 2) - (final_df['spread_line'] / 2)
final_df['opponent_implied_total'] = (final_df['total_line'] / 2) + (final_df['spread_line'] / 2)
final_df['dk_showdown_captain_points'] = final_df['dk_fantasy_points'] * 1.5

# Fill any remaining NaNs
final_df.fillna(0, inplace=True)

# Define final column order
final_columns_order = [
    # Game Info
    'team', 'season', 'week', 'opponent',
    # Betting & Context
    'is_home', 'total_line', 'spread_line', 'team_implied_total', 'opponent_implied_total',
    # Opponent Rolling Stats (Matchup)
    'rolling_4wk_avg_sacks_allowed', 'rolling_4wk_avg_turnovers_committed',
    # DST Performance Stats
    'sacks', 'interceptions', 'fumbles_recovered', 'defensive_tds', 'points_allowed',
    # Fantasy Production
    'dk_fantasy_points', 'dk_showdown_captain_points'
]

final_columns_order = [col for col in final_columns_order if col in final_df.columns]
final_dst_stats_df = final_df[final_columns_order].copy()
final_dst_stats_df = final_dst_stats_df.sort_values(by=['team', 'season', 'week']).reset_index(drop=True)

output_filename = 'dst_comprehensive_advanced_stats_2020_2024.csv'
final_dst_stats_df.to_csv(output_filename, index=False)

print(f"\n✅ Success! Comprehensive advanced stats for DSTs have been exported to '{output_filename}'")
print("\n--- Data Preview ---")
print(final_dst_stats_df.head())