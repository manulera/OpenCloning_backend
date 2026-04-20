# opencloning-db

`opencloning-db` is the database/API companion package for the OpenCloning backend. It provides the app and local data workflows used for OpenCloning database features.

## Run locally

From the repository root:

```bash
# Install or update workspace dependencies
uv sync

# Run the opencloning-db API
uvicorn opencloning_db.api:app --port 8001 --reload --reload-exclude='.venv'
```

If startup succeeds, the API docs should be available at [http://127.0.0.1:8001/docs](http://127.0.0.1:8001/docs).

## Related local workflows

For test and database seed commands, see the repository root [README](../../README.md#running-opencloning-db-locally).
