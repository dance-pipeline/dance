"""Primary pipeline object for DANCE.

Provides a pipeline for generating molecule datasets from a database in various
formats.
"""

_SUPPORTED_DATABASE_TYPES = {"SMILES", "MOL2_DIR"}


class DancePipeline:
    """A pipeline for creating diverse molecule datasets.

    The molecules are sourced from a database. Several types of databases may be
    used. ``database_type`` indicates the database type, and ``database_info``
    contains info for accessing the database. The following table lists the
    possible values for ``database_type``, and the corresponding info that must
    be provided in ``database_info``.

    +---------------+------------------------------------------------+
    | database_type | database_info                                  |
    +===============+================================================+
    | SMILES        | The name of a file listing SMILES strings.     |
    +---------------+------------------------------------------------+
    | MOL2_DIR      | The name of a directory containing mol2 files. |
    +---------------+------------------------------------------------+

    .. note::
        The pipeline does not load any molecules when initialized.

    Parameters
    ----------
    database_type: str
        The type of database.
    database_info
        Info for the database.

    Attributes
    ----------
    database_type : str
        The type of database.
    database_info : str
        Info for the database.
    filter_output_oeb : str
        Output OEB (Openeye Binary) filename for the :meth:`~filter` step.
    fingerprint_output_csv : str
        Output CSV filename for the :meth:`~assign_fingerprint` step.
    """
    def __init__(self, database_type: str, database_info):
        if database_type not in _SUPPORTED_DATABASE_TYPES:
            raise RuntimeError(f"`{database_type}` is not a supported database type")

        self.database_type = database_type
        self.database_info = database_info
