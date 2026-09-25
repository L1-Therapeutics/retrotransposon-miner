---
name: github-pr-edit
description: >-
  Edit pull request titles and bodies on L1-Therapeutics/retrotransposon-miner.
  Use when updating a GitHub PR description, gh pr edit, gh pr comment, or
  when gh fails with Projects (classic) or pullRequest.projectCards.
---

# GitHub pull request edits

`gh pr edit` fails on this repository. The CLI loads classic Projects via
`repository.pullRequest.projectCards`, and that GraphQL field errors:

```
GraphQL: Projects (classic) is being deprecated in favor of the new Projects experience
(repository.pullRequest.projectCards)
```

## Edit a title or body

Use the REST API. `gh` user auth can do this. The GitHub App token from
`scripts/git_app_push.sh` is for push and PR creation, not for this edit.

```bash
gh api -X PATCH repos/L1-Therapeutics/retrotransposon-miner/pulls/NUMBER \
  -f title="..." \
  -f body="$(cat <<'EOF'
## Summary
- …

## Test plan
- [ ] …
EOF
)"
```

Do not retry `gh pr edit` after the `projectCards` error.

`gh pr view` and `gh pr create` still work. Push and open a PR with
`scripts/git_app_push.sh` (deploy-key `git push` is read-only). See
`.cursor/rules/git-app-push.mdc`.
