from typing import Dict, Callable, Tuple, Any, Optional
from packaging import version

# Migration registry - maps version ranges to migration functions
# Format: (start_version, end_version) -> migration_function

MigrationDict = dict[tuple[str, str], Callable]


# Import migration modules to register migrations
def load_migrations() -> MigrationDict:
    """Load all migration modules to register migrations."""

    # Import migration modules
    from .transformations.v0_2_8_3_to_v_0_3_0 import migrate_0_2_8_3_to_0_3_0  # noqa: F401

    return {
        ('0.2.8.3', '0.3.0'): migrate_0_2_8_3_to_0_3_0,
    }


def migrate(data: Dict[str, Any], target_version: Optional[str] = None) -> Dict[str, Any]:
    """
    Migrate data from its current version to the target version.

    Args:
        data: Data to migrate (must contain "backend_version" field)
        target_version: Target version (defaults to latest available)

    Returns:
        Migrated data dictionary
    """
    migration_dict = load_migrations()
    current_version = data.get('backend_version')
    if not current_version:
        current_version = '0.2.8.3'

    all_end_versions = [end for _, end in migration_dict.keys()]
    target_version = max(all_end_versions, key=lambda v: version.parse(v))

    # No migration needed if already at target version
    if current_version == target_version:
        return data

    # Make a copy to avoid modifying the original
    result = data.copy()

    # Continue migrating until we reach the target version
    if 'backend_version' not in result or result['backend_version'] is None:
        current_version = '0.2.8.3'
    else:
        current_version = result['backend_version']

    while version.parse(current_version) < version.parse(target_version):
        next_migration = _find_next_migration(current_version, target_version, migration_dict)
        if not next_migration:
            break  # No more applicable migrations

        start_ver, current_version = next_migration
        migration_func = migration_dict[(start_ver, current_version)]

        # Apply the migration
        result = migration_func(result)
        result['backend_version'] = current_version  # Update the version

    return result


def _find_next_migration(
    current_version: str, target_version: str, migration_dict: MigrationDict
) -> Optional[Tuple[str, str]]:
    """Find the next applicable migration step."""
    applicable_migrations = []

    for (start_ver, end_ver), _ in migration_dict.items():
        # Migration is applicable if:
        # 1. Current version is in the range [start_ver, end_ver)
        # 2. End version is closer to or equal to the target version
        if version.parse(start_ver) <= version.parse(current_version) < version.parse(end_ver) and version.parse(
            end_ver
        ) <= version.parse(target_version):
            applicable_migrations.append((start_ver, end_ver))

    if not applicable_migrations:
        return None

    # Choose the migration that gets us closest to the target version
    return max(applicable_migrations, key=lambda x: version.parse(x[1]))
