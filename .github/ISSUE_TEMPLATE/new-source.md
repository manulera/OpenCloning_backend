---
name: New source
about: Adding a new type of source and corresponding functionality
title: 'New source: <name-of-source>'
labels: new-source
assignees: ''

---

## Description

<!-- Add your description here -->

## Checklist

You can do some of these tasks in parallel, specially if you need to wait for approval, but it's good to keep in mind the intended order:

* [ ] Make a branch in [OpenCloning_LinkML](https://github.com/OpenCloning/OpenCloning_LinkML) where you implement the new source (see extra docs in there) on how to do this.
* [ ] Add your branch as a dependency (via commit id) using poetry in this repository, to update `opencloning-linkml` dependency in `pyproject.toml`:
    ```bash
    poetry add git+https://github.com/OpenCloning/OpenCloning_LinkML#<commit-id>
    ```
* [ ] Implement the new source in a branch in this repository. You will need to add a new endpoint in `main.py` and a new class in `pydantic_models.py` that will handle the new source. You can use the existing sources as a template. Note that you can add extra validation or methods (see examples as well).
* [ ] If the source is a repository_id, make sure to add the appropriate redirect in `endpoints/external_import.py`, endpoint `/repository_id`.
* [ ] Write tests for the new source in `tests/test_endpoints.py`. You can use the existing tests as a template.
* [ ] Once the tests pass, merge the PR of the new source in [OpenCloning_LinkML] and make a release of the `opencloning-linkml` package.
* [ ] Update the dependency in the branch of this repository repository to the new version.
    ```bash
    # Get version from https://pypi.org/project/opencloning-linkml/
    # For instance, if the version is 0.1.6a0, then:
    poetry add opencloning-linkml@0.1.6a0
    ```
* [ ] Implement the frontend functionality in a branch, following [the frontend wiki](https://github.com/manulera/OpenCloning_frontend/wiki/Checklist-%E2%80%90--adding-a-source).
* [ ] Make tests for the frontend functionality.
* [ ] Merge the PR in this repository into master.
* [ ]  Update the backend submodule in the frontend repository to the latest master.
* [ ]  Merge the PR in the frontend repository into master.
* [ ] You are done!
