# Movie Recommendation System

## Overview

This project focuses on building a **movie recommendation system** using Python. It involves the creation of two custom classes: `MovieRecommender` and `GenreRecommender`. The goal is to develop a system that suggests movies or genres to users based on the preferences of other users, using rating data and correlation analysis.

## Datasets

The recommendation engine is based on the **MovieLens dataset**, which contains:
- Over 33 million ratings
- About 86,000 movies
- From more than 330,000 users

Two datasets are used:
- `movies.csv`: contains movie IDs, titles, and genres.
- `ratings.csv`: contains user IDs, movie IDs, and ratings (converted to a 1–10 scale).

## Methodology

### Data Preprocessing
The system performs basic data cleaning, such as:
- Removing unneeded columns like `timestamp`
- Validating rating values
- Reducing matrix size using memory-efficient formats (e.g., `int8` and `float16`)

### MovieRecommender Class
This class constructs a user-movie rating matrix and uses Pearson c
