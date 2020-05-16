"""Provides a function which runs the entire pipeline.

This function takes in a single, all-encompassing configuration. A default
configuration for this function is also provided.
"""
from .dance_pipeline import DancePipeline

# This tag marks the code section to be included in the documentation.
# __sphinx_doc_begin__
DEFAULT_CONFIG = {
    # Parameters set to None MUST be provided in your configuration.

    # Database
    # --------

    # `database_type` is the type of database, and `database_info` is
    # information for the database. For instance, to read molecules from a
    # SMILES file, set `database_type` to "SMILES" and `database_info` to the
    # name of the SMILES file. Refer to the `dance.DancePipeline` documentation
    # for a full list of supported types.
    "database_type": None,
    "database_info": None,

    # Filter
    # ------

    # Function which takes in a single OEMol and outputs a bool telling whether
    # or not the molecule is relevant to the pipeline. For instance, `lambda
    # mol: True` would mark all molecules as relevant.
    "relevance_function": None,

    # Location for storing molecules that have been marked as relevant.
    "filter_output_oeb": "filter_output.oeb",

    # Assign Fingerprint
    # ------------------

    # Function which takes in a single OEMol and outputs a fingerprint.
    # This fingerprint must be an array-like of floats.
    "fingerprint_function": None,

    # Location for storing molecules along with their fingerprints.
    "fingerprint_output_oeb": "fingerprint_output.oeb",

    # Select
    # ------

    # How often to select molecules.
    "selection_frequency": None,

    # This is the main output of the pipeline. `dataset_type` is the type of
    # dataset to output, and `dataset_info` is the information required for that
    # output type. For instance, to output a SMILES file for the dataset, set
    # `dataset_type` to "SMILES" and set `dataset_info` to the name of the
    # SMILES file. Refer to the `dance.DancePipeline` documentation for a full
    # list of supported types.
    "dataset_type": None,
    "dataset_info": None,

    # Location for storing the molecules sorted by fingerprint.
    "sorted_by_fingerprint_oeb": "sorted_by_fingerprint.oeb",

    # Number of molecules that can be sorted in memory at once. This is a
    # low-level parameter that you will likely not need to set.
    "in_memory_sorting_threshold": 25000,
}
# __sphinx_doc_end__


def run_dance(config: {}):
    """Runs the entire DANCE pipeline.

    Parameters
    ----------
    config : dict
        All parameters for running the pipeline. Refer to the documentation for
        :data:`~dance.run.DEFAULT_CONFIG` for all options. Any parameters not
        provided in this variable will be taken from
        :data:`~dance.run.DEFAULT_CONFIG`.

    Raises
    ------
    RuntimeError
        If any necessary fields were not specified in the configuration.
    """

    # Check for errors in the configuration.
    unspecified_fields = []  # Fields the user failed to specify.
    for field in DEFAULT_CONFIG:
        if field not in config:
            if DEFAULT_CONFIG[field] is None:
                # Fields in the DEFAULT_CONFIG that are set to None must be
                # specified by the user.
                unspecified_fields.append(field)
            else:
                # Set the default value for the field.
                config[field] = DEFAULT_CONFIG[field]
    if len(unspecified_fields) > 0:
        raise RuntimeError(f"The following fields were not specified in your configuration: {unspecified_fields}")

    #
    # Run the pipeline.
    #

    pipeline = DancePipeline(
        config["database_type"],
        config["database_info"],
    )
    pipeline.filter(
        config["relevance_function"],
        config["filter_output_oeb"],
    )
    pipeline.assign_fingerprint(
        config["fingerprint_function"],
        config["fingerprint_output_oeb"],
    )
    pipeline.select(
        config["selection_frequency"],
        config["dataset_type"],
        config["dataset_info"],
        config["sorted_by_fingerprint_oeb"],
        config["in_memory_sorting_threshold"],
    )
