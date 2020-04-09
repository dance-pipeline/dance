"""Primary pipeline object for DANCE.

Provides a pipeline for generating molecule datasets from a database in various
formats.
"""
import pathlib
from typing import Union


class DancePipeline:
    """A pipeline for creating diverse molecule datasets.

    The molecules are sourced from a database. Several types of databases may be
    used. ``database_type`` indicates the database type, and ``database_info``
    contains info for accessing the database. The following table lists the
    possible values for ``database_type``, and the corresponding info that must
    be provided in ``database_info``.

    +---------------+------------------------------------------------+-----------------------------+
    | database_type | database_info                                  | type                        |
    +===============+================================================+=============================+
    | SMILES        | The name of a file listing SMILES strings.     | ``str`` or ``pathlib.Path`` |
    +---------------+------------------------------------------------+-----------------------------+
    | MOL2_DIR      | The name of a directory containing mol2 files. | ``str`` or ``pathlib.Path`` |
    +---------------+------------------------------------------------+-----------------------------+

    .. note::
        The pipeline does not load any molecules when initialized, so
        initializing a pipeline is very cheap.

    Parameters
    ----------
    database_type : str
        The type of database.
    database_info : Union[str, pathlib.Path]
        Info for the database -- refer to above table for acceptable types.

    Attributes
    ----------
    database_type : str
        The type of database.
    database_info : Union[str, pathlib.Path]
        Info for the database -- refer to the above table for acceptable types.
    filter_output_oeb : pathlib.Path
        Output OEB (Openeye Binary) filename for the :meth:`~filter` step.
    fingerprint_output_csv : pathlib.Path
        Output CSV filename for the :meth:`~assign_fingerprint` step.
    SUPPORTED_DATABASE_TYPES : FrozenSet[str]
        Set of database types supported by the pipeline.
    """

    SUPPORTED_DATABASE_TYPES = frozenset(["SMILES", "MOL2_DIR"])

    def __init__(self, database_type: str, database_info):
        if database_type not in self.SUPPORTED_DATABASE_TYPES:
            raise RuntimeError(f"`{database_type}` is not a supported database type")

        self.database_type = database_type
        self.database_info = database_info
