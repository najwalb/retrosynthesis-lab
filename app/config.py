"""Application configuration, read from environment."""
import os


class Config:
    DATABASE_URL = (os.environ.get('DATABASE_URL') or '').strip() or None
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'dev-insecure-change-in-production'
    SESSION_COOKIE_NAME = os.environ.get('RXNLAB_SESSION_COOKIE', 'rxnlab_session')
    SESSION_COOKIE_MAX_AGE = 60 * 60 * 24 * 365  # 1 year


DIFFALIGN_MODEL_ID = 'diffalign-align-absorbing-v1'
