"""Chemist feedback collection endpoints."""
import uuid

from flask import Blueprint, g, jsonify, request

from app.db import db_enabled
from app.models_db import Event, Feedback, PredictionRun

bp = Blueprint('feedback', __name__)

VALID_FEEDBACK_TYPES = {'plausibility', 'synthesizable', 'issue_tag', 'free_text'}

VALID_ISSUE_TAGS = {
    'looks-correct',
    'wrong-functional-group', 'implausible-disconnection', 'wrong-stereochemistry',
    'regio-selectivity', 'chemo-selectivity', 'protecting-group-missing',
    'byproduct-ignored', 'multi-step-disguised', 'reagent-uncommon', 'other',
}

VALID_SYNTHESIZABLE = {'yes', 'no', 'maybe'}

MAX_TEXT_LEN = 5000


def _validate_payload(feedback_type: str, payload: dict):
    if feedback_type == 'plausibility':
        score = payload.get('score')
        if not isinstance(score, int) or not (1 <= score <= 5):
            return 'payload.score must be int in [1, 5]'
    elif feedback_type == 'synthesizable':
        value = payload.get('value')
        if value not in VALID_SYNTHESIZABLE:
            return f'payload.value must be one of {sorted(VALID_SYNTHESIZABLE)}'
    elif feedback_type == 'issue_tag':
        tags = payload.get('tags')
        if not isinstance(tags, list) or not tags:
            return 'payload.tags must be a non-empty list'
        if any(t not in VALID_ISSUE_TAGS for t in tags):
            return f'payload.tags entries must be from {sorted(VALID_ISSUE_TAGS)}'
    elif feedback_type == 'free_text':
        text = payload.get('text')
        if not isinstance(text, str) or not text.strip():
            return 'payload.text must be a non-empty string'
        if len(text) > MAX_TEXT_LEN:
            return f'payload.text exceeds {MAX_TEXT_LEN} chars'
    return None


@bp.route('/api/feedback', methods=['POST'])
def post_feedback():
    if not db_enabled() or g.get('db') is None:
        return jsonify({'error': 'Feedback persistence is disabled in this deployment.'}), 503
    if g.get('session_id') is None:
        return jsonify({'error': 'No active session.'}), 400

    data = request.get_json(silent=True) or {}

    try:
        run_id = uuid.UUID(data.get('run_id', ''))
    except (TypeError, ValueError, AttributeError):
        return jsonify({'error': 'run_id must be a UUID'}), 400

    prediction_index = data.get('prediction_index')
    if not isinstance(prediction_index, int) or prediction_index < 0:
        return jsonify({'error': 'prediction_index must be a non-negative int'}), 400

    feedback_type = data.get('feedback_type')
    if feedback_type not in VALID_FEEDBACK_TYPES:
        return jsonify({'error': f'feedback_type must be one of {sorted(VALID_FEEDBACK_TYPES)}'}), 400

    payload = data.get('payload')
    if not isinstance(payload, dict):
        return jsonify({'error': 'payload must be an object'}), 400

    err = _validate_payload(feedback_type, payload)
    if err:
        return jsonify({'error': err}), 400

    if g.db.get(PredictionRun, run_id) is None:
        return jsonify({'error': 'unknown run_id'}), 404

    fb = Feedback(
        run_id=run_id,
        session_id=g.session_id,
        prediction_index=prediction_index,
        feedback_type=feedback_type,
        payload=payload,
    )
    g.db.add(fb)
    g.db.commit()
    return jsonify({'feedback_id': str(fb.feedback_id)}), 201


@bp.route('/api/event', methods=['POST'])
def post_event():
    if not db_enabled() or g.get('db') is None or g.get('session_id') is None:
        return jsonify({'ok': True}), 200  # fire-and-forget when DB is off

    data = request.get_json(silent=True) or {}
    event_type = data.get('event_type')
    if not isinstance(event_type, str) or not event_type:
        return jsonify({'error': 'event_type required'}), 400

    run_id = None
    raw_run_id = data.get('run_id')
    if raw_run_id:
        try:
            run_id = uuid.UUID(raw_run_id)
        except (TypeError, ValueError):
            run_id = None

    payload = data.get('payload')
    if not isinstance(payload, dict):
        payload = None

    g.db.add(Event(
        session_id=g.session_id,
        run_id=run_id,
        event_type=event_type[:100],
        payload=payload,
    ))
    g.db.commit()
    return jsonify({'ok': True}), 201
