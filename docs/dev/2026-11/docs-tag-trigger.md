# Restore tag trigger for documentation deployment

> Status: **completed**. `tags: '*'` restored under `push:` in `.github/workflows/Documenter.yml`.

Date of initial draft: 2026-11-05

Base revision: f74105bf

## Originating prompt

> Somehow stable in docs points to v0.22.1 not the latest release v0.23, which never got deployed, can you investigate why that is?

> Yes please open a PR for this and describe the problem in the description!

## Revision log

- Initial draft.

## Problem description

The `stable` documentation still points to v0.22.1. v0.23.0 was tagged by TagBot on 2026-10-08
(SSH deploy key authenticated fine), but no `Documentation` workflow run was triggered for the tag.

## Background

#1183 (50ee88ad) changed the `on:` block of `Documenter.yml`, restricting `pull_request` to `main`
and, apparently unintentionally, dropping `tags: '*'` from `push:`. Without it, tag pushes no
longer trigger the workflow, so Documenter never deploys versioned docs.

## Summary of changes

Re-add `tags: '*'` under `push:`; keep the `pull_request: branches: [main]` restriction.

## Testing and verification

Will be visible on the next release tag. v0.23.0 can be deployed retroactively via
`gh workflow run Documenter.yml --ref v0.23.0` (the workflow at that tag has `workflow_dispatch`).

## Documentation changes

None.

## Known limitations

v0.23.0 docs require the manual dispatch above; this PR only fixes future tags.

## Future work

None.
