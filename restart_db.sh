set -e
DB_SRC_DIR="packages/opencloning-db/src"

rm -f "${DB_SRC_DIR}/sequencing_files/"*
rm -f "${DB_SRC_DIR}/sequence_files/"*

uv run --directory "${DB_SRC_DIR}" python -m opencloning_db.init_db
