"""Primary pipeline object for DANCE.

Provides a pipeline for generating molecule datasets from a database in various
formats.
"""
import glob
import logging
import pathlib
import shutil
import tempfile
from typing import Tuple, Union

from openeye import oechem
from sortedcontainers import SortedList


class DancePipeline:
    """A pipeline for creating diverse molecule datasets.

    The molecules are sourced from a database. Several types of databases may be
    used. ``database_type`` indicates the database type, and ``database_info``
    contains info for accessing the database. The following table lists the
    possible values for ``database_type``, and the corresponding info that must
    be provided in ``database_info``.

    +---------------+-------------------------------------------------------------------------+
    | database_type | database_info                                                           |
    +===============+=========================================================================+
    | SMILES        | The name of a file listing SMILES strings (`str` or `pathlib.Path`)     |
    +---------------+-------------------------------------------------------------------------+
    | MOL2_DIR      | The name of a directory containing mol2 files (`str` or `pathlib.Path`) |
    +---------------+-------------------------------------------------------------------------+
    | SDF           | The name of an SDF file (`str` or `pathlib.Path`)                       |
    +---------------+-------------------------------------------------------------------------+

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
    fingerprint_output_oeb : pathlib.Path
        Output OEB filename for the :meth:`~assign_fingerprint` step.
    sorted_by_fingerprint_oeb : pathlib.Path
        Output OEB filename for the molecules sorted by fingerprint in the
        :meth:`~select` step.
    num_molecules : int
        The number of molecules in the pipeline. This value is ``None`` until
        :meth:`~filter` is called.

    Raises
    ------
    RuntimeError
        If the provided ``database_type`` is not supported.
    """

    #: Set of database types supported by the pipeline.
    SUPPORTED_DATABASE_TYPES = frozenset(["SMILES", "MOL2_DIR", "SDF"])

    #: Name of the tag used to store the fingerprint length in the molecule
    #: during the :meth:`~assign_fingerprint` step.
    FINGERPRINT_LENGTH_NAME = "dance_fingerprint_length"

    #: Prefix for the tags used to store the fingerprint values in the molecule
    #: during the :meth:`~assign_fingerprint` step.
    FINGERPRINT_VALUE_NAME = "dance_fingerprint_value"

    #: Set of types which the pipeline supports for the final dataset.
    SUPPORTED_DATASET_TYPES = frozenset(["SMILES"])

    def __init__(self, database_type: str, database_info):
        logging.info("Initializing pipeline")

        if database_type not in self.SUPPORTED_DATABASE_TYPES:
            raise RuntimeError(f"`{database_type}` is not a supported database type")

        self.database_type = database_type
        self.database_info = database_info
        self.filter_output_oeb = None
        self.fingerprint_output_oeb = None
        self.sorted_by_fingerprint_oeb = None
        self.num_molecules = None

    @staticmethod
    def get_fingerprint_from_mol(mol: oechem.OEMol) -> Tuple[float]:
        """Utility that retrieves a molecule's fingerprint and returns it as a tuple.

        Refer to :meth:`~assign_fingerprint` for how the fingerprint is stored
        in the molecule.

        Parameters
        ----------
        mol : oechem.OEMol
            The molecule from which to retrieve the fingerprint.

        Returns
        -------
        Tuple[float]
            A tuple containing the fingerprint.

        Raises
        ------
        ValueError
            If the molecule does not contain fingerprint data.
        """
        if not mol.HasData(DancePipeline.FINGERPRINT_LENGTH_NAME):
            raise ValueError("Could not retrieve fingerprint length for molecule.")
        length = mol.GetIntData(DancePipeline.FINGERPRINT_LENGTH_NAME)

        def get_fingerprint_index(i):
            name = f"{DancePipeline.FINGERPRINT_VALUE_NAME}_{i}"
            if not mol.HasData(name):
                raise ValueError(f"Unable to retrieve fingerprint value at index {i}")
            return mol.GetDoubleData(name)

        return tuple(get_fingerprint_index(i) for i in range(length))

    def filter(self, relevance_function, output_oeb: Union[str, pathlib.Path] = "filter_output.oeb"):
        """Uses the ``relevance_function`` to choose molecules from the database.

        Each molecule from the database is read into an OEMol, and the
        ``relevance_function`` is used to determine whether the molecule is
        relevant or not. If the molecule is relevant, it is saved to the
        ``output_oeb`` file.

        ``output_oeb`` is stored in the ``filter_output_oeb`` attribute as a
        ``pathlib.Path`` for future reference by the pipeline.  Also, the
        ``num_molecules`` attribute is set to the number of molecules that were
        marked as relevant.

        Parameters
        ----------
        relevance_function : function(molecule) -> bool
            A function which takes in a single OEMol and outputs a bool telling
            whether or not the molecule is relevant
        output_oeb : Union[str, pathlib.Path], optional
            Name of an OEB (Openeye Binary) file for storing the relevant
            molecules. If this file already exists, it will be overwritten!
        """
        logging.info("Filtering molecules")

        self.filter_output_oeb = pathlib.Path(output_oeb)
        self.num_molecules = 0

        filtered_molecule_stream = oechem.oemolostream(str(output_oeb))
        for idx, mol in enumerate(self._generate_molecules_from_database()):
            if relevance_function(oechem.OEMol(mol)):
                logging.debug("Molecule %d relevant!", idx + 1)
                self.num_molecules += 1
                oechem.OEWriteMolecule(filtered_molecule_stream, mol)
            else:
                logging.debug("Molecule %d not relevant", idx + 1)
        filtered_molecule_stream.close()

        logging.info("Molecules after filtering: %d", self.num_molecules)

    def _generate_molecules_from_database(self) -> oechem.OEMol:
        """Generator which yields molecules from the database.

        Information about the database is passed in via the ``database_type``
        and ``database_info`` attributes of this class.
        """
        if self.database_type == "SMILES" or self.database_type == "SDF":
            ifs = oechem.oemolistream(str(self.database_info))
            for mol in ifs.GetOEMols():
                yield mol
            ifs.close()
        elif self.database_type == "MOL2_DIR":
            for mol2file in glob.iglob(str(pathlib.Path(self.database_info) / "*.mol2")):
                ifs = oechem.oemolistream(mol2file)
                mol = oechem.OEMol()
                oechem.OEReadMolecule(ifs, mol)
                yield mol
                ifs.close()

    def assign_fingerprint(self,
                           fingerprint_function,
                           output_oeb: Union[str, pathlib.Path] = "fingerprint_output.oeb"):
        """Assigns a fingerprint to the molecules from the :meth:`~filter` step.

        Each molecule is passed to the ``fingerprint_function`` to create a
        fingerprint. This fingerprint is attached to the molecule, and the
        molecule is saved to the ``output_oeb``.

        Specifically, the fingerprints are stored as data entries in each
        molecule and may be accessed later with ``OEMol.GetIntData()`` and
        ``OEMol.GetDoubleData()``. The following data fields are stored in each
        molecule:

        - { :attr:`~FINGERPRINT_LENGTH_NAME` } (`int`) - the number of values in the
          fingerprint
        - { :attr:`~FINGERPRINT_VALUE_NAME` }_{i} (`float`) - the values of the
          fingerprint, with i = [0,1,2,...,{ :attr:`~FINGERPRINT_LENGTH_NAME` } - 1]

        Note that :attr:`~FINGERPRINT_LENGTH_NAME` and
        :attr:`~FINGERPRINT_VALUE_NAME` are attributes of this class.

        ``output_oeb`` is stored in the ``fingerprint_output_oeb`` attribute as
        a ``pathlib.Path`` for future reference by the pipeline.

        Parameters
        ----------
        fingerprint_function : function(molecule) -> array_like[float]
            A function which provides the fingerprint for a single OEMol. This
            fingerprint must be an array of floats.
        output_oeb : Union[str, pathlib.Path], optional
            Name of an OEB (Openeye Binary) file for storing the molecules along
            with their fingerprints. If this file already exists, it will be
            overwritten!
        """
        logging.info("Assigning fingerprints")

        self.fingerprint_output_oeb = pathlib.Path(output_oeb)

        filtered_molecule_stream = oechem.oemolistream(str(self.filter_output_oeb))
        fingerprinted_molecule_stream = oechem.oemolostream(str(output_oeb))
        for idx, mol in enumerate(filtered_molecule_stream.GetOEMols()):
            logging.debug("Assigning fingerprint to molecule %d / %d", idx + 1, self.num_molecules)
            fingerprint = fingerprint_function(oechem.OEMol(mol))
            self._add_single_fingerprint(mol, fingerprint)
            oechem.OEWriteMolecule(fingerprinted_molecule_stream, mol)
        filtered_molecule_stream.close()
        fingerprinted_molecule_stream.close()

        logging.info("%d fingerprints assigned", self.num_molecules)

    def _add_single_fingerprint(self, mol: oechem.OEMol, fingerprint: "array-like of float"):
        """Adds the fingerprint to the molecule in-place."""
        mol.SetIntData(self.FINGERPRINT_LENGTH_NAME, len(fingerprint))
        for i, val in enumerate(fingerprint):
            mol.SetDoubleData(f"{self.FINGERPRINT_VALUE_NAME}_{i}", val)

    def select(self,
               selection_frequency: int,
               dataset_type: str,
               dataset_info,
               sorted_oeb: Union[str, pathlib.Path] = "sorted_by_fingerprint.oeb",
               in_memory_sorting_threshold=25000):
        """Chooses molecules based on their assigned fingerprints.

        Specifically, this method sorts the molecules by their fingerprint,
        storing the sorted molecules in ``sorted_oeb``.  It then selects every
        ``selection_frequency``-th molecule, starting from 0. For example, if
        ``selection_frequency = 4``, we select molecule ``0,4,8,...``

        The selected molecules are stored in a format determined by the
        ``dataset_type`` and ``dataset_info``, similar to the ``database_type``
        and ``database_info`` passed in to ``__init__``. The following table
        shows the possible values of ``dataset_type``, and the corresponding info
        that is expected to be in ``dataset_info``

        +--------------+-------------------------------------------------------------------------+
        | dataset_type | dataset_info                                                            |
        +==============+=========================================================================+
        | SMILES       | The name of a file for writing SMILES strings (`str` or `pathlib.Path`) |
        +--------------+-------------------------------------------------------------------------+

        ``sorted_oeb`` is stored in the ``sorted_by_fingerprint_oeb`` attribute
        as a ``pathlib.Path`` for future reference by the pipeline.

        Parameters
        ----------
        selection_frequency : int
            How often to select molecules.
        dataset_type : str
            The type of dataset.
        dataset_info : Union[str, pathlib.Path]
            Info for the dataset -- refer to above table for acceptable types.
        sorted_oeb : Union[str, pathlib.Path], optional
            Name of an OEB (Openeye Binary) file for storing the molecules
            sorted by fingerprint. If this file already exists, it will be
            overwritten!
        in_memory_sorting_threshold : int
            The number of molecules that can be sorted in memory at once. This
            is a low-level parameter for the sorting algorithm, and most users
            should not have to deal with it. Note that if this threshold is
            maximized, then the maximum number of molecules that can be sorted
            by the pipeline is approximately the square of this threshold.

        Raises
        ------
        RuntimeError
            If the provided ``dataset_type`` is not supported.
        RuntimeError
            If ``selection_frequency`` is less than 1.
        """
        logging.info("Selecting molecules for final dataset")

        if dataset_type not in self.SUPPORTED_DATASET_TYPES:
            raise RuntimeError(f"`{dataset_type}` is not a supported dataset type")
        if selection_frequency < 1:
            raise RuntimeError(f"selection_frequency({selection_frequency}) must be an integer greater >= 1")

        logging.info("Sorting molecules by fingerprint")
        self.sorted_by_fingerprint_oeb = pathlib.Path(sorted_oeb)
        self._sort_molecules_by_fingerprint(self.fingerprint_output_oeb, self.sorted_by_fingerprint_oeb,
                                            in_memory_sorting_threshold)
        sorted_molstream = oechem.oemolistream(str(self.sorted_by_fingerprint_oeb))

        # Open any necessary files for outputting the dataset. This will vary by
        # dataset_type.
        if dataset_type == "SMILES":
            dataset_stream = oechem.oemolostream(str(dataset_info))
            dataset_stream.SetFormat(oechem.OEFormat_SMI)

        # Perform the selection.
        logging.info("Selecting molecules from sorted file")
        for idx, mol in enumerate(sorted_molstream.GetOEMols()):
            if idx % selection_frequency == 0:
                if dataset_type == "SMILES":
                    oechem.OEWriteMolecule(dataset_stream, mol)

        # Clean up.
        if dataset_type == "SMILES":
            dataset_stream.close()

        sorted_molstream.close()

        logging.info("Dataset selection completed")

    @staticmethod
    def _sort_molecules_by_fingerprint(input_oeb: pathlib.Path,
                                       output_oeb: pathlib.Path,
                                       in_memory_sorting_threshold: int = 25000):
        """Sorts the molecules in ``input_oeb`` and writes them to ``output_oeb``.

        This method implements an external mergesort on the molecules (see
        https://en.wikipedia.org/wiki/External_sorting). The general idea is to
        split the molecule file into several smaller files, sort each file in
        memory, and combine those smaller files into the final file.

        Parameters
        ----------
        input_oeb : pathlib.Path
            The input molecules.
        output_oeb: pathlib.Path
            The output molecules.
        in_memory_sorting_threshold : int
            The number of molecules that can be sorted in memory.
        """
        # Create a temporary directory for storing sorted chunks of molecules.
        tmpdir = pathlib.Path(tempfile.mkdtemp())

        # Break the original input file into chunks. Sort each chunk and store
        # it in a temporary file.
        input_stream = oechem.oemolistream(str(input_oeb))
        chunk_idx = 0
        molecules_in_chunk = []  # Tuples of (fingerprint, molecule)

        def sort_and_save_molecules_in_chunk(chunk_idx):
            # Save the molecules in the current chunk to a temporary file.
            chunk_file = tmpdir / f"chunk_{chunk_idx}.oeb"
            chunk_stream = oechem.oemolostream(str(chunk_file))
            molecules_in_chunk.sort()
            for _, mol in molecules_in_chunk:
                oechem.OEWriteMolecule(chunk_stream, mol)
            chunk_stream.close()

            # Prepare for the next chunk.
            molecules_in_chunk.clear()

        for idx, mol in enumerate(input_stream.GetOEMols()):
            molecules_in_chunk.append((DancePipeline.get_fingerprint_from_mol(mol), oechem.OEMol(mol)))
            if idx % in_memory_sorting_threshold == (in_memory_sorting_threshold - 1):
                sort_and_save_molecules_in_chunk(chunk_idx)
                chunk_idx += 1
        if len(molecules_in_chunk) > 0:
            # Handle any remaining molecules.
            sort_and_save_molecules_in_chunk(chunk_idx)
            chunk_idx += 1
        input_stream.close()

        # Merge the sorted chunks. This priority queue stores one entry for each
        # stream. The entry consists of the stream, the molecule most recently
        # taken from the stream, and the fingerprint of that molecule. We
        # repeatedly remove the molecule with the smallest fingerprint from the
        # queue and write it to the output stream. Refer to
        # https://en.wikipedia.org/wiki/K-way_merge_algorithm#Heap for more
        # info.
        priority_queue = SortedList()
        num_chunks = chunk_idx

        # Build the priority queue.
        for chunk_idx in range(num_chunks):
            chunk_file = tmpdir / f"chunk_{chunk_idx}.oeb"
            chunk_stream = oechem.oemolistream(str(chunk_file))
            mol = oechem.OEMol()
            oechem.OEReadMolecule(chunk_stream, mol)
            fingerprint = DancePipeline.get_fingerprint_from_mol(mol)
            priority_queue.add((fingerprint, mol, chunk_stream))

        # Write the output.
        output_stream = oechem.oemolostream(str(output_oeb))
        while len(priority_queue) > 0:
            # Retrieve the next molecule to write to the output.
            fingerprint, mol, chunk_stream = priority_queue.pop(0)
            oechem.OEWriteMolecule(output_stream, mol)

            if oechem.OEReadMolecule(chunk_stream, mol):
                # If the chunk stream still has molecules, read the next one and
                # push it into the priority queue.
                fingerprint = DancePipeline.get_fingerprint_from_mol(mol)
                priority_queue.add((fingerprint, mol, chunk_stream))
            else:
                # Otherwise, just close the stream.
                chunk_stream.close()
        output_stream.close()

        # Clean up the temporary directory.
        shutil.rmtree(str(tmpdir))
