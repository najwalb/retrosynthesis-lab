/* Chemist feedback panel interactions.
   Event-delegated so it works on both server-rendered cards and AJAX-injected
   cards (results-container.innerHTML = html). */
(function () {
    'use strict';

    var SAVED_VISIBLE_MS = 1800;     // how long the inline "✓ Saved" stays
    var ERROR_VISIBLE_MS = 4000;
    var DROPDOWN_CLOSE_MS = 1000;    // delay before auto-closing the dropdown

    function panelOf(el) { return el.closest('.feedback-panel'); }

    function showStatus(statusEl, msg, isError) {
        if (!statusEl) return;
        // Clear any previous fade timer on this element so rapid actions don't
        // overlap (e.g., user wiggles the slider).
        if (statusEl._fadeTimer) {
            clearTimeout(statusEl._fadeTimer);
            statusEl._fadeTimer = null;
        }
        statusEl.textContent = msg;
        statusEl.classList.toggle('fb-error', !!isError);
        // Force reflow so the transition replays even if the same class is re-added.
        // eslint-disable-next-line no-unused-expressions
        statusEl.offsetHeight;
        statusEl.classList.add('fb-visible');
        statusEl._fadeTimer = setTimeout(function () {
            statusEl.classList.remove('fb-visible');
            statusEl._fadeTimer = null;
        }, isError ? ERROR_VISIBLE_MS : SAVED_VISIBLE_MS);
    }

    function markSummarySaved(detailsEl) {
        if (!detailsEl) return;
        var tick = detailsEl.querySelector('[data-summary-tick]');
        if (!tick) return;
        tick.textContent = '✓';
        tick.classList.add('fb-saved');
    }

    function postFeedback(panel, feedbackType, payload, options) {
        options = options || {};
        var runId = panel.getAttribute('data-run-id');
        var predIdx = parseInt(panel.getAttribute('data-prediction-index'), 10);
        if (!runId) return;  // disabled-mode panel; banner already explains

        var statusEl = options.statusEl || null;

        fetch('/api/feedback', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                run_id: runId,
                prediction_index: predIdx,
                feedback_type: feedbackType,
                payload: payload,
            }),
        }).then(function (r) {
            if (r.status === 201) {
                showStatus(statusEl, '✓ Saved', false);
                if (options.onSuccess) options.onSuccess();
            } else {
                return r.json().then(function (j) {
                    showStatus(statusEl, 'Error: ' + (j.error || r.statusText), true);
                });
            }
        }).catch(function (err) {
            showStatus(statusEl, 'Network error: ' + err.message, true);
        });
    }

    /* ── Plausibility slider ─────────────────────────────────────────────
       `change` fires once on release, not on every drag pixel.            */
    document.addEventListener('change', function (e) {
        var t = e.target;
        if (t.classList && t.classList.contains('fb-rate-slider')) {
            var panel = panelOf(t);
            if (!panel) return;
            var score = parseInt(t.value, 10);
            if (!Number.isInteger(score) || score < 1 || score > 5) return;
            postFeedback(panel, 'plausibility', { score: score }, {
                statusEl: panel.querySelector('.fb-rate-status'),
            });
        }
    });

    /* Live numeric readout next to the slider. */
    document.addEventListener('input', function (e) {
        var t = e.target;
        if (t.classList && t.classList.contains('fb-rate-slider')) {
            var panel = panelOf(t);
            var v = panel && panel.querySelector('.fb-rate-value');
            if (v) v.textContent = t.value;
        }
    });

    document.addEventListener('click', function (e) {
        /* ── Submit tags ───────────────────────────────────────────────── */
        var tagSubmit = e.target.closest('.fb-tag-submit');
        if (tagSubmit) {
            var panel = panelOf(tagSubmit);
            if (!panel) return;
            var checkboxes = panel.querySelectorAll('.fb-tag-body input[type=checkbox]:checked');
            var tags = Array.prototype.map.call(checkboxes, function (cb) { return cb.value; });
            var statusEl = panel.querySelector('.fb-tag-status');
            if (tags.length === 0) {
                showStatus(statusEl, 'Select at least one tag.', true);
                return;
            }
            postFeedback(panel, 'issue_tag', { tags: tags }, {
                statusEl: statusEl,
                onSuccess: function () {
                    var details = panel.querySelector('details.fb-tag');
                    markSummarySaved(details);
                    // Uncheck so a follow-up submission requires fresh intent.
                    checkboxes.forEach(function (cb) { cb.checked = false; });
                    setTimeout(function () {
                        if (details) details.open = false;
                    }, DROPDOWN_CLOSE_MS);
                },
            });
            return;
        }

        /* ── Submit comment ────────────────────────────────────────────── */
        var commentSubmit = e.target.closest('.fb-comment-submit');
        if (commentSubmit) {
            var panel2 = panelOf(commentSubmit);
            if (!panel2) return;
            var ta = panel2.querySelector('.fb-text');
            var text = ta ? ta.value.trim() : '';
            var statusEl2 = panel2.querySelector('.fb-comment-status');
            if (!text) {
                showStatus(statusEl2, 'Comment is empty.', true);
                return;
            }
            postFeedback(panel2, 'free_text', { text: text }, {
                statusEl: statusEl2,
                onSuccess: function () {
                    var details = panel2.querySelector('details.fb-comment');
                    markSummarySaved(details);
                    if (ta) ta.value = '';
                    setTimeout(function () {
                        if (details) details.open = false;
                    }, DROPDOWN_CLOSE_MS);
                },
            });
        }
    });
})();
