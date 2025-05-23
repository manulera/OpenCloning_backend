from opencloning_linkml.migrations.model_archive.v0_2_9 import CloningStrategy as OldCloningStrategy


def migrate_0_2_8_3_to_0_3_0(data: dict) -> dict:
    """
    Migrate data from backed version 0.2.8.3 to 0.3.0.

    ## Error fixing
    There were two errors, fixed in https://github.com/manulera/OpenCloning_backend/pull/305

    The first error concerned gateway assemblies. `gateway_overlap` was returning the entire
    overlap, which matched regex like twtGTACAAAaaa (for attB1). That created assemblies in which
    the overlapping part may have mismatches on the w. Now, instead of returning the whole twtGTACAAAaaa
    as overlap, it returns only the common part GTACAAA.

    The second error was a bug in `SimpleSequenceLocation.from_simple_location`, where origin-spanning
    features were not being read correctly and turned into the entire sequence. This was being used in
    `generate_assemblies` and producing wrong assembly products. There is an example in the pull request.

    ## Field setting

    This migration sets the `backend_version` field to the current version of the backend.
    """
    cs = OldCloningStrategy.model_validate(data)
    cs.backend_version = '0.3.0'
    return cs.model_dump()
