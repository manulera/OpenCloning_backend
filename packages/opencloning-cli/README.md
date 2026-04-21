# opencloning-cli

`opencloning-cli` is the command-line companion to the OpenCloning backend. Its first iteration exists to support **Cypress-driven frontend testing**: spin up a seeded SQLite database, snapshot it, and reset it to that snapshot between tests.

The CLI is intentionally narrow. It does not talk to the API. It does not manage production databases. It only orchestrates files and DB seeding on the local filesystem.

## Install

`opencloning-cli` is a `uv` workspace member. From the repository root:

```bash
uv sync
uv run opencloning-cli --help
```
