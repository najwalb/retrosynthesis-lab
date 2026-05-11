"""SQLAlchemy engine, session helpers, and on-startup table creation."""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Optional

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

_engine = None
_SessionLocal: Optional[sessionmaker] = None


def init_db(database_url: str) -> None:
    global _engine, _SessionLocal
    _engine = create_engine(database_url, future=True, pool_pre_ping=True)
    _SessionLocal = sessionmaker(bind=_engine, autoflush=False, autocommit=False, future=True)

    from app.models_db import Base
    Base.metadata.create_all(_engine)
    _seed_models()


def _seed_models() -> None:
    from app.config import DIFFALIGN_MODEL_ID
    from app.models_db import Model

    assert _SessionLocal is not None
    with _SessionLocal() as session:
        if session.get(Model, DIFFALIGN_MODEL_ID) is None:
            session.add(Model(
                model_id=DIFFALIGN_MODEL_ID,
                display_name='DiffAlign (align-absorbing)',
                version='epoch760',
                description='Graph diffusion model for single-step retrosynthesis.',
                supports_inpainting=True,
                metadata_={'arch': 'graph-diffusion', 'training': 'USPTO-50k'},
                added_at=datetime.now(timezone.utc),
            ))
            session.commit()


def get_session() -> Optional[Session]:
    if _SessionLocal is None:
        return None
    return _SessionLocal()


def db_enabled() -> bool:
    return _SessionLocal is not None
