# Scripting with OpenCloning

This folder contains examples of how to use the OpenCloning API functions to create scripts to automate cloning.

## Installation

You can install the opencloning package from pypi:

```bash
pip install opencloning
```
If you want to install a particular github commit, you can do so by specifying the commit hash:

```bash
pip install git+https://github.com/OpenCloning/OpenCloning_backend.git@<commit_hash>
```
For examples, see the jupyter notebooks in this folder.

## Is OpenCloning the right tool for me?

To simulate cloning, OpenCloning uses [pydna](https://github.com/pydna-group/pydna), a python package that extends Biopython's sequence classes, providing functions to do PCR, Gibson, Restriction & Ligation, etc. You can do all the cloning of OpenCloning using pydna alone, and that would require less code (some functionality such as Gateway Cloning and Golden Gate is only available in OpenCloning now, but the code will be moved to pydna in the future).

On the other hand, OpenCloning allows you to record the cloning history in a standard format that you can also load on the web interface. If you want to know more about this, create an issue on this repository and we can help you decide what is the best tool for you.
